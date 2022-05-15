use std::collections::HashMap;

use uuid::Uuid;

use crate::{utils::PanicCast, Genome, PaletteMode, Rgba, Transform};

use super::{
    rng::Flam3Rng,
    variations::{apply_xform, VariationPrecalculations},
    Flam3DeThreadHelper, Flam3ThreadHelper, RenderOps, TransformSelector,
};

const FUSE_27: u32 = 15;
const FUSE_28: u32 = 100;

#[derive(Default)]
struct TransformPrecalcs {
    precalcs: HashMap<Uuid, VariationPrecalculations>,
}

impl TransformPrecalcs {
    fn get(&mut self, transform: &Transform) -> &mut VariationPrecalculations {
        let id = transform.id;
        self.precalcs
            .entry(id)
            .or_insert_with(|| VariationPrecalculations::new(transform))
    }
}

/// Runs a set of iterations mutating samples after a set of iterations have been skipped.
fn flam3_iterate(
    cp: &Genome,
    iterations: u32,
    skip_iterations: u32,
    samples: &mut [f64],
    selector: &mut TransformSelector,
    rng: &mut Flam3Rng,
) -> u32 {
    let mut consecutive_failures = 0;
    let mut bad_iterations = 0;

    let mut p = [samples[0], samples[1], samples[2], samples[3]];

    selector.reset();

    let mut precalcs = TransformPrecalcs::default();

    let mut iteration = -(skip_iterations.i32());
    while iteration < (iterations.i32()) {
        let xform = selector.next(cp, rng);
        let precalc = precalcs.get(xform);

        let mut q = [0.0; 4];
        if !apply_xform(xform, &p, &mut q, precalc, rng) {
            consecutive_failures += 1;
            bad_iterations += 1;
            if consecutive_failures < 5 {
                p = q;
                iteration -= 1;
                continue;
            } else {
                consecutive_failures = 0;
            }
        } else {
            consecutive_failures = 0;
        }

        p = q;

        if let Some(ref xform) = cp.final_transform {
            if xform.opacity == 1.0 || rng.isaac_01() < xform.opacity {
                let precalc = precalcs.get(xform);
                apply_xform(xform, &p, &mut q, precalc, rng);
                /* Keep the opacity from the original xform */
                q[3] = p[3];
            }
        }

        /* if fuse over, store it */
        if iteration >= 0 {
            let i = (iteration * 4).usize();

            samples[i..i + 4].copy_from_slice(&q);
        }

        iteration += 1;
    }

    bad_iterations
}

pub(super) fn iter_thread<Ops: RenderOps>(
    mut fthp: Flam3ThreadHelper,
    buckets: &mut [[Ops::Bucket; 5]],
) -> Result<(), String> {
    log::trace!("Starting iteration thread");
    let ficp = &mut fthp.fic;
    let cmap_size = ficp.cmap_size.i32();
    let cmap_size_m1 = cmap_size - 1;

    //    double eta = 0.0;

    let fuse = if ficp.spec.earlyclip {
        FUSE_28
    } else {
        FUSE_27
    };
    let mut xform_distrib = TransformSelector::new(&fthp.cp)?;
    let mut iter_storage = vec![0.0; 4 * ficp.spec.sub_batch_size.usize()];

    for sub_batch in (0..ficp.batch_size).step_by(ficp.spec.sub_batch_size.usize()) {
        /* sub_batch is double so this is sketchy */
        let sub_batch_size = if sub_batch + ficp.spec.sub_batch_size > ficp.batch_size {
            ficp.batch_size - sub_batch
        } else {
            ficp.spec.sub_batch_size
        };

        log::trace!(
            "Rendering sub batches {} to {}",
            sub_batch,
            sub_batch + sub_batch_size
        );

        /* Seed iterations */
        iter_storage[0] = fthp.rng.isaac_11();
        iter_storage[1] = fthp.rng.isaac_11();
        iter_storage[2] = fthp.rng.isaac_01();
        iter_storage[3] = fthp.rng.isaac_01();

        /* Execute iterations */
        let badcount = flam3_iterate(
            &fthp.cp,
            sub_batch_size,
            fuse,
            &mut iter_storage,
            &mut xform_distrib,
            &mut fthp.rng,
        );

        /* Add the badcount to the counter */
        ficp.badvals += badcount;

        /* Put them in the bucket accumulator */
        for j in (0..(sub_batch_size * 4)).step_by(4) {
            let p = &iter_storage[j.usize()..];

            let (p0, p1) = if fthp.cp.rotate != 0.0 {
                ficp.rot
                    .transform((&[p[0] - fthp.cp.rot_center.x, p[1] - fthp.cp.rot_center.y]).into())
                    .into()
            } else {
                (p[0], p[1])
            };

            if p0 >= ficp.bounds[0]
                && p1 >= ficp.bounds[1]
                && p0 <= ficp.bounds[2]
                && p1 <= ficp.bounds[3]
            {
                let logvis = p[3];

                /* Skip if invisible */
                if logvis == 0.0 {
                    continue;
                }

                let start = (ficp.ws0 * p0 - ficp.wb0s0).usize()
                    + ficp.width.usize() * (ficp.hs1 * p1 - ficp.hb1s1).usize();
                let b = &mut buckets[start];

                let dbl_index0 = p[2] * cmap_size.f64();
                let color_index0 = dbl_index0.i32();

                let interpcolor = if PaletteMode::Linear == fthp.cp.palette_mode {
                    let (cindex, dbl_frac) = if color_index0 < 0 {
                        (0, 0.0)
                    } else if color_index0 >= cmap_size_m1 {
                        (cmap_size_m1.usize() - 1, 1.0)
                    } else {
                        /* interpolate between color_index0 and color_index0+1 */
                        (color_index0.usize(), dbl_index0 - color_index0.f64())
                    };

                    let mut interpcolor = Rgba::default();
                    for ci in 0..4 {
                        interpcolor[ci] = ficp.dmap[cindex][ci] * (1.0 - dbl_frac)
                            + ficp.dmap[cindex + 1][ci] * dbl_frac;
                    }
                    interpcolor
                } else {
                    /* Palette mode step */
                    let cindex = if color_index0 < 0 {
                        0
                    } else if color_index0 >= cmap_size_m1 {
                        cmap_size_m1.usize()
                    } else {
                        color_index0.usize()
                    };

                    ficp.dmap[cindex].clone()
                };

                Ops::bump_no_overflow(
                    b,
                    &[
                        logvis * interpcolor[0],
                        logvis * interpcolor[1],
                        logvis * interpcolor[2],
                        logvis * interpcolor[3],
                        logvis * 255.0,
                    ],
                );
            }
        }
    }

    log::trace!("Thread complete, found {} bad values", ficp.badvals);
    Ok(())
}

pub(super) fn de_thread(mut fdthp: Flam3DeThreadHelper) -> Result<(), String> {
    log::trace!("Starting density estimation thread");
    Ok(())
}
