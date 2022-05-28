use std::collections::HashMap;

use palette::{Pixel, Srgba};
use uuid::Uuid;

use crate::{render::flam3::filters::DE_THRESH, utils::PanicCast, Genome, PaletteMode, Transform};

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
        for j in (0..(sub_batch_size.usize() * 4)).step_by(4) {
            let p = &iter_storage[j..j + 4];

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

                    let mut interpcolor = [0.0; 4];
                    let first = ficp.dmap[cindex].as_raw::<[f64]>();
                    let second = ficp.dmap[cindex + 1].as_raw::<[f64]>();
                    for ci in 0..4 {
                        interpcolor[ci] = first[ci] * (1.0 - dbl_frac) + second[ci] * dbl_frac;
                    }
                    *Srgba::from_raw(&interpcolor)
                } else {
                    /* Palette mode step */
                    let cindex = if color_index0 < 0 {
                        0
                    } else if color_index0 >= cmap_size_m1 {
                        cmap_size_m1.usize()
                    } else {
                        color_index0.usize()
                    };

                    ficp.dmap[cindex]
                };

                Ops::bump_no_overflow(
                    &mut buckets[start],
                    &[
                        logvis * interpcolor.red,
                        logvis * interpcolor.green,
                        logvis * interpcolor.blue,
                        logvis * interpcolor.alpha,
                        logvis * 255.0,
                    ],
                );
            }
        }
    }

    log::trace!(
        "Iteration thread complete, found {} bad values",
        ficp.badvals
    );
    Ok(())
}

pub(super) fn de_thread<Ops: RenderOps>(
    dthp: Flam3DeThreadHelper,
    buckets: &[[Ops::Bucket; 5]],
    accumulate: &mut [[Ops::Accumulator; 4]],
) {
    let oversample = dthp.oversample;
    let ss = (oversample.f64() / 2.0).floor().i32();
    let scf = (oversample & 1) == 0;
    let scfact = (oversample.f64() / (oversample.f64() + 1.0)).powf(2.0);
    let wid = dthp.width;
    let hig = dthp.height;
    let str = (oversample - 1) + dthp.start_row;
    let enr = ((oversample - 1).i32() + dthp.end_row).u32();

    log::trace!(
        "Starting density estimation thread for rows {} to {}, width={}, height={}",
        str,
        enr,
        wid,
        hig
    );

    /* Density estimation code */
    for j in str..enr {
        for i in oversample - 1..wid - (oversample - 1) {
            let mut c = [0.0; 4];
            let mut f_select = 0.0;
            let index = (i + j * wid).usize();

            /* Don't do anything if there's no hits here */
            if buckets[index][3].f64() == 0.0 || buckets[index][4].f64() == 0.0 {
                continue;
            }

            /* Count density in ssxss area   */
            /* Scale if OS>1 for equal iters */
            for ii in -ss..=ss {
                for jj in -ss..=ss {
                    let index = ((i.i32() + ii) + (j.i32() + jj) * wid.i32()).usize();
                    f_select += buckets[index][4].f64() / 255.0;
                }
            }

            if scf {
                f_select *= scfact;
            }

            let mut f_select_int = if f_select > dthp.de.max_filtered_counts.f64() {
                dthp.de.max_filter_index
            } else if f_select <= DE_THRESH.f64() {
                f_select.ceil().u32() - 1
            } else {
                DE_THRESH + (f_select - DE_THRESH.f64()).powf(dthp.curve).floor().u32()
            };

            /* If the filter selected below the min specified clamp it to the min */
            if f_select_int > dthp.de.max_filter_index {
                f_select_int = dthp.de.max_filter_index;
            }

            /* We only have to calculate the values for ~1/8 of the square */
            let mut f_coef_idx = (f_select_int * dthp.de.kernel_size).usize();

            let arr_filt_width = (dthp.de.filter_widths[f_select_int.usize()]).ceil().u32() - 1;

            let b = &buckets[(i + j * wid).usize()..];

            for jj in 0..=arr_filt_width.i32() {
                for ii in 0..=jj {
                    /* Skip if coef is 0 */
                    if dthp.de.filter_coefs[f_coef_idx] == 0.0 {
                        f_coef_idx += 1;
                        continue;
                    }

                    c[0] = b[0][0].f64();
                    c[1] = b[0][1].f64();
                    c[2] = b[0][2].f64();
                    c[3] = b[0][3].f64();

                    let ls = dthp.de.filter_coefs[f_coef_idx]
                        * (dthp.k1 * (1.0 + c[3] * dthp.k2).ln())
                        / c[3];

                    c[0] *= ls;
                    c[1] *= ls;
                    c[2] *= ls;
                    c[3] *= ls;

                    if jj == 0 && ii == 0 {
                        Ops::add_c_to_accum(accumulate, i, ii, j, jj, wid, hig, &c);
                    } else if ii == 0 {
                        Ops::add_c_to_accum(accumulate, i, jj, j, 0, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, -jj, j, 0, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, 0, j, jj, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, 0, j, -jj, wid, hig, &c);
                    } else if jj == ii {
                        Ops::add_c_to_accum(accumulate, i, ii, j, jj, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, -ii, j, jj, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, ii, j, -jj, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, -ii, j, -jj, wid, hig, &c);
                    } else {
                        Ops::add_c_to_accum(accumulate, i, ii, j, jj, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, -ii, j, jj, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, ii, j, -jj, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, -ii, j, -jj, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, jj, j, ii, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, -jj, j, ii, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, jj, j, -ii, wid, hig, &c);
                        Ops::add_c_to_accum(accumulate, i, -jj, j, -ii, wid, hig, &c);
                    }

                    f_coef_idx += 1;
                }
            }
        }
    }

    log::trace!("Density estimation thread complete");
}
