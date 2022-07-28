use palette::{ComponentWise, Srgba};

use crate::math::{ln, pow, sqr};
use crate::render::flam3::rng::Flam3Rng;
use crate::{render::flam3::filters::DE_THRESH, utils::PanicCast, PaletteMode, Transform};

use super::variations::VariationPrecalculations;
use super::{rng::IsaacRng, variations::apply_xform, DensityEstimationContext, TransformSelector};
use super::{Accumulator, IterationContext};

/// Runs a set of iterations mutating samples after a set of iterations have been skipped.
fn flam3_iterate(
    iterations: u32,
    skip_iterations: u32,
    samples: &mut [f64],
    selector: &TransformSelector,
    final_xform: &Option<(Transform, VariationPrecalculations)>,
    rng: &mut IsaacRng,
) -> u32 {
    let mut consecutive_failures = 0;
    let mut bad_iterations = 0;

    let mut p = [samples[0], samples[1], samples[2], samples[3]];

    let mut iteration = 0;
    let mut pos = 0;
    let total_iterations = skip_iterations + iterations;
    while iteration < total_iterations {
        let (xform, precalc) = selector.next(rng);

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

        if let Some((xform, ref precalc)) = final_xform {
            if xform.opacity == 1.0 || rng.next_01() < xform.opacity {
                apply_xform(xform, &p, &mut q, precalc, rng);
                /* Keep the opacity from the original xform */
                q[3] = p[3];
            }
        }

        /* if fuse over, store it */
        if iteration >= skip_iterations {
            samples[pos..pos + 4].copy_from_slice(&q);
            pos += 4;
        }

        iteration += 1;
    }

    bad_iterations
}

pub(super) trait IterationThreadStorage {
    fn increase_bucket(&mut self, x: usize, y: usize, color: Srgba<f64>, logvis: f64);
}

pub(super) fn iter_thread<S: IterationThreadStorage>(
    context: &IterationContext,
    xform_distrib: &TransformSelector,
    final_xform: &Option<(Transform, VariationPrecalculations)>,
    rng: &mut IsaacRng,
    storage: &mut S,
) -> Result<(), String> {
    log::trace!("Starting iteration thread");
    let cmap_size = context.dmap.len().f64();
    let cmap_size_m1 = context.dmap.len().i32() - 1;

    let mut iter_storage = vec![0.0; 4 * context.sub_batch_size.usize()];
    let mut badvals: u32 = 0;

    for sub_batch in (0..context.batch_size).step_by(context.sub_batch_size.usize()) {
        /* sub_batch is double so this is sketchy */
        let sub_batch_size = if sub_batch + context.sub_batch_size > context.batch_size {
            context.batch_size - sub_batch
        } else {
            context.sub_batch_size
        };

        /* Seed iterations */
        iter_storage[0] = rng.next_11();
        iter_storage[1] = rng.next_11();
        iter_storage[2] = rng.next_01();
        iter_storage[3] = rng.next_01();

        /* Execute iterations */
        let badcount = flam3_iterate(
            sub_batch_size,
            context.skip_iterations,
            &mut iter_storage,
            xform_distrib,
            final_xform,
            rng,
        );

        /* Add the badcount to the counter */
        badvals += badcount;

        /* Put them in the bucket accumulator */
        for j in (0..(sub_batch_size.usize() * 4)).step_by(4) {
            let p = &iter_storage[j..j + 4];

            let (p0, p1) = if context.rotate != 0.0 {
                context
                    .rot
                    .transform((&[p[0] - context.rot_center.x, p[1] - context.rot_center.y]).into())
                    .into()
            } else {
                (p[0], p[1])
            };

            if p0 >= context.bounds[0].x
                && p1 >= context.bounds[0].y
                && p0 <= context.bounds[1].x
                && p1 <= context.bounds[1].y
            {
                let logvis = p[3];

                /* Skip if invisible */
                if logvis == 0.0 {
                    continue;
                }

                let dbl_index0 = p[2] * cmap_size.f64();
                let color_index0 = dbl_index0.i32();

                let interpcolor = if PaletteMode::Linear == context.palette_mode {
                    let (cindex, dbl_frac) = if color_index0 < 0 {
                        (0, 0.0)
                    } else if color_index0 >= cmap_size_m1 {
                        (cmap_size_m1.usize() - 1, 1.0)
                    } else {
                        /* interpolate between color_index0 and color_index0+1 */
                        (color_index0.usize(), dbl_index0 - color_index0.f64())
                    };

                    context.dmap[cindex]
                        .component_wise(&context.dmap[cindex + 1], |first, second| {
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

                    context.dmap[cindex]
                };

                storage.increase_bucket(
                    (context.ws0 * p0 - context.wb0s0).usize(),
                    (context.hs1 * p1 - context.hb1s1).usize(),
                    interpcolor,
                    logvis,
                );
            }
        }
    }

    log::trace!("Iteration thread complete, found {} bad values", badvals);

    Ok(())
}

pub(super) trait DensityEstimationThreadStorage {
    fn width(&self) -> usize;
    fn height(&self) -> usize;
    fn density(&self, x: usize, y: usize) -> f64;
    fn accumulate(&self, x: usize, y: usize) -> Accumulator<f64>;
    fn increase_accumulator(&mut self, x: usize, y: usize, pixel: &Accumulator<f64>);
}

pub(super) fn empty_de_thread<S: DensityEstimationThreadStorage>(
    k1: f64,
    k2: f64,
    storage: &mut S,
) {
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

pub(super) fn de_thread<S: DensityEstimationThreadStorage>(
    context: DensityEstimationContext,
    storage: &mut S,
) {
    let supersample = context.supersample;
    let ss = (supersample.f64() / 2.0).floor().usize();
    let scf = (supersample & 1) == 0;
    let scfact = sqr!(supersample.f64() / (supersample.f64() + 1.0));
    let start_row = ((supersample - 1) + context.start_row).usize();
    let end_row = ((supersample - 1) + context.end_row).usize();
    let de_filters = context.de_filters;

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
                DE_THRESH
                    + pow!(f_select - DE_THRESH.f64(), context.curve)
                        .floor()
                        .usize()
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

                let ls = coef * (context.k1 * ln!(1.0 + current_pixel.alpha * context.k2))
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
