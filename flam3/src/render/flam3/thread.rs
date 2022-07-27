use std::collections::HashMap;

use palette::{ComponentWise, Srgba};
use uuid::Uuid;

use crate::math::{ln, pow, sqr};
use crate::render::flam3::rng::Flam3Rng;
use crate::{render::flam3::filters::DE_THRESH, utils::PanicCast, Genome, PaletteMode, Transform};

use super::{
    rng::IsaacRng,
    variations::{apply_xform, VariationPrecalculations},
    Flam3DeThreadHelper, TransformSelector,
};
use super::{Accumulator, Flam3IterConstants};

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
    rng: &mut IsaacRng,
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
            if xform.opacity == 1.0 || rng.next_01() < xform.opacity {
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

pub(super) trait IterationStorage {
    fn increase_bucket(&mut self, x: usize, y: usize, color: Srgba<f64>, logvis: f64);
}

pub(super) fn iter_thread<S: IterationStorage>(
    cp: &Genome,
    ficp: &Flam3IterConstants,
    rng: &mut IsaacRng,
    storage: &mut S,
) -> Result<(), String> {
    log::trace!("Starting iteration thread");
    let cmap_size = ficp.cmap_size.i32();
    let cmap_size_m1 = cmap_size - 1;

    let fuse = if ficp.earlyclip { FUSE_28 } else { FUSE_27 };
    let mut xform_distrib = TransformSelector::new(cp)?;
    let mut iter_storage = vec![0.0; 4 * ficp.sub_batch_size.usize()];
    let mut badvals: u32 = 0;

    for sub_batch in (0..ficp.batch_size).step_by(ficp.sub_batch_size.usize()) {
        /* sub_batch is double so this is sketchy */
        let sub_batch_size = if sub_batch + ficp.sub_batch_size > ficp.batch_size {
            ficp.batch_size - sub_batch
        } else {
            ficp.sub_batch_size
        };

        /* Seed iterations */
        iter_storage[0] = rng.next_11();
        iter_storage[1] = rng.next_11();
        iter_storage[2] = rng.next_01();
        iter_storage[3] = rng.next_01();

        /* Execute iterations */
        let badcount = flam3_iterate(
            cp,
            sub_batch_size,
            fuse,
            &mut iter_storage,
            &mut xform_distrib,
            rng,
        );

        /* Add the badcount to the counter */
        badvals += badcount;

        /* Put them in the bucket accumulator */
        for j in (0..(sub_batch_size.usize() * 4)).step_by(4) {
            let p = &iter_storage[j..j + 4];

            let (p0, p1) = if cp.rotate != 0.0 {
                ficp.rot
                    .transform((&[p[0] - cp.rot_center.x, p[1] - cp.rot_center.y]).into())
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

                let dbl_index0 = p[2] * cmap_size.f64();
                let color_index0 = dbl_index0.i32();

                let interpcolor = if PaletteMode::Linear == cp.palette_mode {
                    let (cindex, dbl_frac) = if color_index0 < 0 {
                        (0, 0.0)
                    } else if color_index0 >= cmap_size_m1 {
                        (cmap_size_m1.usize() - 1, 1.0)
                    } else {
                        /* interpolate between color_index0 and color_index0+1 */
                        (color_index0.usize(), dbl_index0 - color_index0.f64())
                    };

                    ficp.dmap[cindex].component_wise(&ficp.dmap[cindex + 1], |first, second| {
                        first * (1.0 - dbl_frac) + second * dbl_frac
                    })
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

                storage.increase_bucket(
                    (ficp.ws0 * p0 - ficp.wb0s0).usize(),
                    (ficp.hs1 * p1 - ficp.hb1s1).usize(),
                    interpcolor,
                    logvis,
                );
            }
        }
    }

    log::trace!("Iteration thread complete, found {} bad values", badvals);

    Ok(())
}

pub(super) trait DensityEstimationStorage {
    fn width(&self) -> usize;
    fn height(&self) -> usize;
    fn density(&self, x: usize, y: usize) -> f64;
    fn accumulate(&self, x: usize, y: usize) -> Accumulator<f64>;
    fn increase_accumulator(&mut self, x: usize, y: usize, pixel: &Accumulator<f64>);
}

pub(super) fn empty_de_thread<S: DensityEstimationStorage>(k1: f64, k2: f64, storage: &mut S) {
    for j in 0..storage.height() {
        for i in 0..storage.width() {
            let mut pixel = storage.accumulate(i, j);

            if pixel.alpha == 0.0 {
                continue;
            }

            let ls = (k1 * ln!(1.0 + pixel.alpha * k2)) / pixel.alpha;
            pixel *= ls;

            storage.increase_accumulator(i, j, &pixel);
        }
    }
}

pub(super) fn de_thread<S: DensityEstimationStorage>(dthp: Flam3DeThreadHelper, storage: &mut S) {
    let supersample = dthp.supersample;
    let ss = (supersample.f64() / 2.0).floor().usize();
    let scf = (supersample & 1) == 0;
    let scfact = sqr!(supersample.f64() / (supersample.f64() + 1.0));
    let start_row = ((supersample - 1) + dthp.start_row).usize();
    let end_row = ((supersample - 1) + dthp.end_row).usize();
    let de_filters = dthp.de_filters;

    log::trace!(
        "Starting density estimation thread for rows {} to {}, width={}, height={}",
        start_row,
        end_row,
        storage.width(),
        storage.height()
    );

    /* Density estimation code */
    for j in start_row..end_row {
        for i in supersample.usize() - 1..storage.width() - (supersample.usize() - 1) {
            let mut f_select = 0.0;
            let current_pixel = storage.accumulate(i, j);

            /* Don't do anything if there's no hits here */
            if current_pixel.alpha == 0.0 || storage.density(i, j) == 0.0 {
                continue;
            }

            /* Count density in ssxss area   */
            /* Scale if OS>1 for equal iters */
            for ii in i - ss..=i + ss {
                for jj in j - ss..=j + ss {
                    f_select += storage.density(ii, jj);
                }
            }

            if scf {
                f_select *= scfact;
            }

            let mut f_select_int = if f_select > de_filters.max_filtered_counts {
                de_filters.filters.len() - 1
            } else if f_select <= DE_THRESH.f64() {
                f_select.ceil().usize() - 1
            } else {
                DE_THRESH + pow!(f_select - DE_THRESH.f64(), dthp.curve).floor().usize()
            };

            /* If the filter selected below the min specified clamp it to the min */
            if f_select_int > (de_filters.filters.len() - 1) {
                f_select_int = de_filters.filters.len() - 1;
            }

            /* We only have to calculate the values for ~1/8 of the square */
            let filter = &de_filters.filters[f_select_int];

            for (coef, coef_operations) in
                filter.operations(i, j, storage.width(), storage.height())
            {
                if coef_operations.is_empty() {
                    continue;
                }

                let ls = coef * (dthp.k1 * ln!(1.0 + current_pixel.alpha * dthp.k2))
                    / current_pixel.alpha;
                let pixel = &current_pixel * ls;

                for (x, y) in coef_operations {
                    storage.increase_accumulator(x, y, &pixel)
                }
            }
        }
    }

    log::trace!("Density estimation thread complete");
}
