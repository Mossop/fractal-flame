use std::f64::consts::PI;

use crate::math::{cos, exp, pow, sin, sqr, sqrt, sum_sqr};
use crate::rect::Rect;
use crate::{fastdiv, utils::PanicCast, SpatialFilter, TemporalFilterType};

use super::{Field, Flam3Frame};

const FLAM3_MITCHELL_B: f64 = 1.0 / 3.0;
const FLAM3_MITCHELL_C: f64 = 1.0 / 3.0;
pub const DE_THRESH: usize = 100;

fn filter_scale(spatial_filter: SpatialFilter) -> f64 {
    match spatial_filter {
        SpatialFilter::Gaussian => 1.5,
        SpatialFilter::Hermite => 1.0,
        SpatialFilter::Box => 0.5,
        SpatialFilter::Triangle => 1.0,
        SpatialFilter::Bell => 1.5,
        SpatialFilter::BSpline => 2.0,
        SpatialFilter::Mitchell => 2.0,
        SpatialFilter::Blackman => 1.0,
        SpatialFilter::Catrom => 2.0,
        SpatialFilter::Hanning => 1.0,
        SpatialFilter::Hamming => 1.0,
        SpatialFilter::Lanczos3 => 3.0,
        SpatialFilter::Lanczos2 => 2.0,
        SpatialFilter::Quadratic => 1.5,
    }
}

fn flam3_hermite_filter(mut t: f64) -> f64 {
    /* f(t) = 2|t|^3 - 3|t|^2 + 1, -1 <= t <= 1 */
    if t < 0.0 {
        t = -t;
    }

    if t < 1.0 {
        (2.0 * t - 3.0) * sqr!(t) + 1.0
    } else {
        0.0
    }
}

fn flam3_box_filter(t: f64) -> f64 {
    if (t > -0.5) && (t <= 0.5) {
        1.0
    } else {
        0.0
    }
}

fn flam3_triangle_filter(mut t: f64) -> f64 {
    if t < 0.0 {
        t = -t;
    }

    if t < 1.0 {
        1.0 - t
    } else {
        0.0
    }
}

fn flam3_bell_filter(mut t: f64) -> f64 {
    /* box (*) box (*) box */
    if t < 0.0 {
        t = -t;
    }

    if t < 0.5 {
        0.75 - sqr!(t)
    } else if t < 1.5 {
        t -= 1.5;
        0.5 * sqr!(t)
    } else {
        0.0
    }
}

fn flam3_b_spline_filter(mut t: f64) -> f64 {
    /* box (*) box (*) box (*) box */
    if t < 0.0 {
        t = -t;
    }

    if t < 1.0 {
        let tt = sqr!(t);
        (0.5 * tt * t) - tt + (2.0 / 3.0)
    } else if t < 2.0 {
        t = 2.0 - t;
        (1.0 / 6.0) * (t * t * t)
    } else {
        0.0
    }
}

fn flam3_sinc(mut x: f64) -> f64 {
    x *= PI;
    if x != 0.0 {
        sin!(x) / x
    } else {
        1.0
    }
}

fn flam3_blackman_filter(x: f64) -> f64 {
    0.42 + 0.5 * cos!(PI * x) + 0.08 * cos!(2.0 * PI * x)
}

fn flam3_catrom_filter(x: f64) -> f64 {
    if x < -2.0 {
        0.0
    } else if x < -1.0 {
        0.5 * (4.0 + x * (8.0 + x * (5.0 + x)))
    } else if x < 0.0 {
        0.5 * (2.0 + sqr!(x) * (-5.0 - 3.0 * x))
    } else if x < 1.0 {
        0.5 * (2.0 + sqr!(x) * (-5.0 + 3.0 * x))
    } else if x < 2.0 {
        0.5 * (4.0 + x * (-8.0 + x * (5.0 - x)))
    } else {
        0.0
    }
}

fn flam3_mitchell_filter(mut t: f64) -> f64 {
    let tt = sqr!(t);
    if t < 0.0 {
        t = -t;
    }

    if t < 1.0 {
        t = ((12.0 - 9.0 * FLAM3_MITCHELL_B - 6.0 * FLAM3_MITCHELL_C) * (t * tt))
            + ((-18.0 + 12.0 * FLAM3_MITCHELL_B + 6.0 * FLAM3_MITCHELL_C) * tt)
            + (6.0 - 2.0 * FLAM3_MITCHELL_B);
        t / 6.0
    } else if t < 2.0 {
        t = ((-1.0 * FLAM3_MITCHELL_B - 6.0 * FLAM3_MITCHELL_C) * (t * tt))
            + ((6.0 * FLAM3_MITCHELL_B + 30.0 * FLAM3_MITCHELL_C) * tt)
            + ((-12.0 * FLAM3_MITCHELL_B - 48.0 * FLAM3_MITCHELL_C) * t)
            + (8.0 * FLAM3_MITCHELL_B + 24.0 * FLAM3_MITCHELL_C);
        t / 6.0
    } else {
        0.0
    }
}

fn flam3_hanning_filter(x: f64) -> f64 {
    0.5 + 0.5 * cos!(PI * x)
}

fn flam3_hamming_filter(x: f64) -> f64 {
    0.54 + 0.46 * cos!(PI * x)
}

fn flam3_lanczos3_filter(mut t: f64) -> f64 {
    if t < 0.0 {
        t = -t;
    }

    if t < 3.0 {
        flam3_sinc(t) * flam3_sinc(t / 3.0)
    } else {
        0.0
    }
}

fn flam3_lanczos2_filter(mut t: f64) -> f64 {
    if t < 0.0 {
        t = -t;
    }

    if t < 2.0 {
        flam3_sinc(t) * flam3_sinc(t / 2.0)
    } else {
        0.0
    }
}

fn flam3_gaussian_filter(x: f64) -> f64 {
    (exp!(-2.0 * sqr!(x))) * sqrt!(2.0 / PI)
}

fn flam3_quadratic_filter(x: f64) -> f64 {
    if x < -1.5 {
        0.0
    } else if x < -0.5 {
        0.5 * (x + 1.5) * (x + 1.5)
    } else if x < 0.5 {
        0.75 - sqr!(x)
    } else if x < 1.5 {
        0.5 * (x - 1.5) * (x - 1.5)
    } else {
        0.0
    }
}

fn flam3_spatial_filter(spatial_filter: SpatialFilter, x: f64) -> f64 {
    match spatial_filter {
        SpatialFilter::Gaussian => flam3_gaussian_filter(x),
        SpatialFilter::Hermite => flam3_hermite_filter(x),
        SpatialFilter::Box => flam3_box_filter(x),
        SpatialFilter::Triangle => flam3_triangle_filter(x),
        SpatialFilter::Bell => flam3_bell_filter(x),
        SpatialFilter::BSpline => flam3_b_spline_filter(x),
        SpatialFilter::Mitchell => flam3_mitchell_filter(x),
        SpatialFilter::Blackman => flam3_sinc(x) * flam3_blackman_filter(x),
        SpatialFilter::Catrom => flam3_catrom_filter(x),
        SpatialFilter::Hanning => flam3_sinc(x) * flam3_hanning_filter(x),
        SpatialFilter::Hamming => flam3_sinc(x) * flam3_hamming_filter(x),
        SpatialFilter::Lanczos3 => flam3_lanczos3_filter(x) * flam3_sinc(x / 3.0),
        SpatialFilter::Lanczos2 => flam3_lanczos2_filter(x) * flam3_sinc(x / 2.0),
        _ => flam3_quadratic_filter(x),
    }
}

fn normalize_vector(v: &mut Rect<f64>) -> Result<(), String> {
    let mut t = 0.0;
    for val in v.iter() {
        t += val;
    }

    if t == 0.0 {
        return Err("Spatial filter value is too small".to_string());
    }

    t = 1.0 / t;
    for val in v.iter_mut() {
        *val *= t;
    }

    Ok(())
}

pub(super) fn create_spatial_filter(frame: &Flam3Frame, field: Field) -> Result<Rect<f64>, String> {
    let spatial_filter = frame.genomes[0].spatial_filter;
    let supersample = frame.genomes[0].spatial_supersample;
    let filter_radius = frame.genomes[0].spatial_filter_radius;
    let aspect_ratio = frame.pixel_aspect_ratio;
    let scale = filter_scale(spatial_filter);

    let fw = 2.0 * scale * supersample.f64() * filter_radius / aspect_ratio;
    let mut fwidth = fw.usize() + 1;

    /* Make sure the filter kernel has same parity as supersample */
    if (fwidth.u32() ^ supersample) & 1 > 0 {
        fwidth += 1;
    }

    /* Calculate the coordinate scaling factor for the kernel values */
    let adjust = if fw > 0.0 {
        scale * fwidth.f64() / fw
    } else {
        1.0
    };

    let mut filter = Rect::square(fwidth);

    /* fill in the coefs */
    for i in 0..fwidth {
        let ii = ((2 * i + 1).f64() / fwidth.f64() - 1.0) * adjust;

        for j in 0..fwidth {
            /* Calculate the function inputs for the kernel function */
            let mut jj = ((2 * j + 1).f64() / fwidth.f64() - 1.0) * adjust;

            /* Scale for scanlines */
            if field != Field::Both {
                jj *= 2.0;
            }

            /* Adjust for aspect ratio */
            jj /= aspect_ratio;

            filter[(i, j)] =
                flam3_spatial_filter(spatial_filter, ii) * flam3_spatial_filter(spatial_filter, jj);
        }
    }

    normalize_vector(&mut filter)?;

    Ok(filter)
}

pub(super) struct TemporalFilter {
    pub sumfilt: f64,
    num_temporal_samples: u32,
    filter: Vec<f64>,
    deltas: Vec<f64>,
}

impl TemporalFilter {
    pub fn new(
        filter_type: TemporalFilterType,
        filter_exp: f64,
        filter_width: f64,
        num_batches: u32,
        num_temporal_samples: u32,
    ) -> TemporalFilter {
        let num_steps = (num_batches * num_temporal_samples).usize();

        let mut maxfilt = 0.0;
        let mut sumfilt = 0.0;

        /* Allocate memory - this must be freed in the calling routine! */
        let mut deltas = vec![0.0; num_steps];
        let mut filter = vec![0.0; num_steps];

        /* Deal with only one step */
        if num_steps == 1 {
            deltas[0] = 0.0;
            filter[0] = 1.0;
            return TemporalFilter {
                sumfilt: 1.0,
                num_temporal_samples,
                filter,
                deltas,
            };
        }

        /* Define the temporal deltas */
        for (i, delta) in deltas.iter_mut().enumerate() {
            *delta = (i.f64() / (num_steps.f64() - 1.0) - 0.5) * filter_width;
        }

        /* Define the filter coefs */
        match filter_type {
            TemporalFilterType::Exp => {
                for (i, filt) in filter.iter_mut().enumerate() {
                    let slpx = if filter_exp >= 0.0 {
                        (i.f64() + 1.0) / num_steps.f64()
                    } else {
                        (num_steps - i).f64() / num_steps.f64()
                    };

                    /* Scale the color based on these values */
                    *filt = pow!(slpx, filter_exp.abs());

                    /* Keep the max */
                    if *filt > maxfilt {
                        maxfilt = *filt;
                    }
                }
            }
            TemporalFilterType::Gaussian => {
                let halfsteps = num_steps.f64() / 2.0;
                for (i, filt) in filter.iter_mut().enumerate() {
                    /* Gaussian */
                    *filt = flam3_spatial_filter(
                        SpatialFilter::Gaussian,
                        filter_scale(SpatialFilter::Gaussian) * (i.f64() - halfsteps).abs()
                            / halfsteps,
                    );
                    /* Keep the max */
                    if *filt > maxfilt {
                        maxfilt = *filt;
                    }
                }
            }
            TemporalFilterType::Box => {
                for filt in filter.iter_mut() {
                    *filt = 1.0;
                }

                maxfilt = 1.0;
            }
        }

        /* Adjust the filter so that the max is 1.0, and */
        /* calculate the K2 scaling factor  */
        for filt in filter.iter_mut() {
            *filt /= maxfilt;
            sumfilt += *filt;
        }

        sumfilt /= num_steps.f64();

        TemporalFilter {
            sumfilt,
            num_temporal_samples,
            filter,
            deltas,
        }
    }

    pub fn color_scale(&self, batch_num: u32, temporal_sample_num: u32) -> f64 {
        self.filter[(batch_num * self.num_temporal_samples + temporal_sample_num).usize()]
    }

    pub fn time_offset(&self, batch_num: u32, temporal_sample_num: u32) -> f64 {
        self.deltas[(batch_num * self.num_temporal_samples + temporal_sample_num).usize()]
    }
}

#[derive(Default, Debug, Clone)]
pub(super) struct DensityEstimatorFilter {
    pub width: usize,
    pub coefs: Vec<f64>,
}

#[derive(Default, Debug, Clone)]
pub(super) struct DensityEstimatorFilters {
    pub max_filtered_counts: f64,
    pub filters: Vec<DensityEstimatorFilter>,
}

impl DensityEstimatorFilters {
    pub fn new(
        max_radius: f64,
        min_radius: f64,
        curve: f64,
        supersample: u32,
    ) -> Result<DensityEstimatorFilters, String> {
        let mut de_filters = DensityEstimatorFilters::default();

        if curve <= 0.0 {
            return Err("Estimator curve must be > 0".to_string());
        }

        if max_radius < min_radius {
            return Err("Estimator must be larger than estimator_minimum".to_string());
        }

        // We should scale the filter width by the supersample
        // The '+1' comes from the assumed distance to the first pixel
        let comp_max_radius = max_radius * supersample.f64() + 1.0;
        let comp_min_radius = min_radius * supersample.f64() + 1.0;

        // Calculate how many filter kernels we need based on the decay function
        let filter_count = pow!(comp_max_radius / comp_min_radius, 1.0 / curve);
        if filter_count > 1e7 {
            return Err(format!(
                "Too many filters required in this configuration ({})",
                filter_count
            ));
        }
        let filter_count = filter_count.ceil().usize();

        // Condense the smaller kernels to save space
        let de_max_ind;
        if filter_count > DE_THRESH {
            de_max_ind =
                DE_THRESH + pow!((filter_count - DE_THRESH).f64(), curve).ceil().usize() + 1;
            de_filters.max_filtered_counts =
                (pow!((de_max_ind - DE_THRESH).f64(), 1.0 / curve).usize() + DE_THRESH).f64();
        } else {
            de_max_ind = filter_count;
            de_filters.max_filtered_counts = de_max_ind.f64();
        }

        /* Allocate the memory for these filters */
        /* and the hit/width lookup vector       */
        let max_radius_int = comp_max_radius.ceil().u32();
        let de_half_size = (max_radius_int - 1).i32();
        let kernel_size = (max_radius_int * (max_radius_int + 1) / 2).usize();

        /* Generate the filter coefficients */
        let mut done = false;

        for filter_index in 0..de_max_ind {
            let mut de_filt_sum = 0.0;

            /* Calculate the filter width for this number of hits in a bin */
            let mut de_filt_h = if filter_index < DE_THRESH {
                comp_max_radius / pow!((filter_index + 1).f64(), curve)
            } else {
                let adjloop = pow!((filter_index - DE_THRESH).f64(), 1.0 / curve) + DE_THRESH.f64();
                fastdiv!(comp_max_radius, pow!(adjloop + 1.0, curve))
            };

            /* Once we've reached the min radius, don't populate any more */
            if de_filt_h <= comp_min_radius {
                de_filt_h = comp_min_radius;
                done = true;
            }

            let mut filter = DensityEstimatorFilter {
                width: de_filt_h.ceil().usize() - 1,
                coefs: vec![0.0; kernel_size],
            };

            /* Calculate norm of kernel separately (easier) */
            for dej in -de_half_size..=de_half_size {
                for dek in -de_half_size..=de_half_size {
                    let de_filt_d = sqrt!(sum_sqr!(dej, dek).f64()) / de_filt_h;

                    /* Only populate the coefs within this radius */
                    if de_filt_d <= 1.0 {
                        /* Gaussian */
                        de_filt_sum += flam3_spatial_filter(
                            SpatialFilter::Gaussian,
                            filter_scale(SpatialFilter::Gaussian) * de_filt_d,
                        );

                        /* Epanichnikov */
                        // de_filt_sum += (1.0 - sqr!(de_filt_d));
                    }
                }
            }

            let mut filter_coef_idx = 0;

            /* Calculate the unique entries of the kernel */
            for dej in 0..=de_half_size {
                for dek in 0..=dej {
                    let de_filt_d = sqrt!(sum_sqr!(dej, dek).f64()) / de_filt_h;

                    /* Only populate the coefs within this radius */
                    if de_filt_d > 1.0 {
                        filter.coefs[filter_coef_idx] = 0.0;
                    } else {
                        /* Gaussian */
                        filter.coefs[filter_coef_idx] = flam3_spatial_filter(
                            SpatialFilter::Gaussian,
                            filter_scale(SpatialFilter::Gaussian) * de_filt_d,
                        ) / de_filt_sum;

                        /* Epanichnikov */
                        // de_filter_coefs[filter_coef_idx] = (1.0 - sqr!(de_filt_d))/de_filt_sum;
                    }

                    filter_coef_idx += 1;
                }
            }

            de_filters.filters.push(filter);

            if done {
                break;
            }
        }

        Ok(de_filters)
    }
}
