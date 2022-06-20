use std::f64::consts::PI;

use palette::{encoding, FromColor, Hsv, Pixel, Srgb, Srgba};

use crate::logging::state;
use crate::math::{cos, ln, pow, sin, sqr};
use crate::{
    render::flam3::{rng::Flam3Rng, Flam3DeHelper},
    utils::PanicCast,
};

use super::{
    filters::{create_spatial_filter, flam3_create_de_filters, flam3_create_temporal_filter},
    flam3_interpolate,
    thread::{de_thread, iter_thread},
    Field, Flam3DeThreadHelper, Flam3Frame, Flam3IterConstants, Flam3ThreadHelper, RenderStorage,
};

const WHITE_LEVEL: u32 = 255;
const PREFILTER_WHITE: u32 = 255;

fn flam3_calc_alpha(density: f64, gamma: f64, linrange: f64) -> f64 {
    let dnorm = density;
    let funcval = pow!(linrange, gamma);

    if dnorm > 0.0 {
        if dnorm < linrange {
            let frac = dnorm / linrange;
            (1.0 - frac) * dnorm * (funcval / linrange) + frac * pow!(dnorm, gamma)
        } else {
            pow!(dnorm, gamma)
        }
    } else {
        0.0
    }
}

fn flam3_calc_newrgb(cbuf: &[f64; 3], ls: f64, highpow: f64) -> [f64; 3] {
    let mut maxa = -1.0;
    let mut maxc = 0.0;
    let mut newrgb = [0.0; 3];

    if ls == 0.0 || (cbuf[0] == 0.0 && cbuf[1] == 0.0 && cbuf[2] == 0.0) {
        return [0.0, 0.0, 0.0];
    }

    /* Identify the most saturated channel */
    for val in cbuf {
        let a = ls * (val / PREFILTER_WHITE.f64());
        if a > maxa {
            maxa = a;
            maxc = val / PREFILTER_WHITE.f64();
        }
    }

    /* If a channel is saturated and we have a non-negative highlight power */
    /* modify the color to prevent hue shift                                */
    if maxa > 255.0 && highpow >= 0.0 {
        let newls = 255.0 / maxc;
        let lsratio = pow!(newls / ls, highpow);

        /* Calculate the max-value color (ranged 0 - 1) */
        for rgbi in 0..3 {
            newrgb[rgbi] = newls * (cbuf[rgbi] / PREFILTER_WHITE.f64()) / 255.0;
        }

        /* Reduce saturation by the lsratio */
        let mut hsv = Hsv::<encoding::Srgb, f64>::from_color(*Srgb::<f64>::from_raw(&newrgb));
        hsv.saturation *= lsratio;
        newrgb = Srgb::from_color(hsv).into_raw();

        for rgbi in newrgb.iter_mut() {
            *rgbi *= 255.0;
        }
    } else {
        let newls = 255.0 / maxc;
        let mut adjhlp = -highpow;
        if adjhlp > 1.0 {
            adjhlp = 1.0;
        }
        if maxa <= 255.0 {
            adjhlp = 1.0;
        }

        /* Calculate the max-value color (ranged 0 - 1) interpolated with the old behaviour */
        for rgbi in 0..3 {
            newrgb[rgbi] =
                ((1.0 - adjhlp) * newls + adjhlp * ls) * (cbuf[rgbi] / PREFILTER_WHITE.f64());
        }
    }
    newrgb
}

pub(super) fn render_rectangle<S: RenderStorage>(
    mut frame: Flam3Frame,
    mut buffer: &mut [u8],
    field: Field,
) -> Result<(), String> {
    let mut background = Srgba::default();
    let mut gamma = 0.0;
    let mut vib_gam_n = 0.0;
    let mut vibrancy = 0.0;

    // interpolate and get a control point
    let cp = flam3_interpolate(&frame.genomes, frame.time, 0.0)?;

    let supersample = cp.spatial_supersample;
    let highpow = cp.highlight_power;
    let num_batches = cp.passes;
    let num_temporal_samples = cp.num_temporal_samples;

    // Set up the output image dimensions, adjusted for scanline
    let image_width = cp.size.width;
    let (image_height, out_width) = if field != Field::Both {
        if field == Field::Odd {
            // Offset the buffer
            buffer =
                &mut buffer[(frame.channels * frame.bytes_per_channel * cp.size.width).usize()..];
        }

        (cp.size.height / 2, cp.size.width * 2)
    } else {
        (cp.size.height, cp.size.width)
    };

    // Spatial Filter kernel creation
    let spatial_filter = create_spatial_filter(&frame, field)?;

    //  batch filter
    //  may want to revisit this at some point
    let batch_filter = vec![1.0 / num_batches.f64(); num_batches.usize()];

    //  temporal filter - we must free temporal_filter and temporal_deltas at the end
    let (sumfilt, temporal_filter, temporal_deltas) = flam3_create_temporal_filter(
        (num_batches * num_temporal_samples).usize(),
        cp.temporal_filter,
        cp.temporal_filter_exp,
        cp.temporal_filter_width,
    )?;

    /*
       the number of additional rows of buckets we put at the edge so
       that the filter doesn't go off the edge
    */
    let mut gutter_width = (spatial_filter.width().u32() - supersample) / 2;

    /*
       Check the size of the density estimation filter.
       If the 'radius' of the density estimation filter is greater than the
       gutter width, we have to pad with more.  Otherwise, we can use the same value.
    */
    let mut max_gnm_de_fw = 0;
    for genome in frame.genomes.iter() {
        let this_width = (genome.estimator_radius * supersample.f64()).ceil().u32();
        if this_width > max_gnm_de_fw {
            max_gnm_de_fw = this_width;
        }
    }

    //  Add OS-1 for the averaging at the edges, if it's > 0 already
    if max_gnm_de_fw > 0 {
        max_gnm_de_fw += supersample - 1;
    }

    //  max_gnm_de_fw is now the number of pixels of additional gutter
    //  necessary to appropriately perform the density estimation filtering
    //  Check to see if it's greater than the gutter_width
    let de_offset = if max_gnm_de_fw > gutter_width {
        let delta = max_gnm_de_fw - gutter_width;
        gutter_width = max_gnm_de_fw;
        delta.usize()
    } else {
        0
    };

    let mut fic = Flam3IterConstants::new(frame.clone());

    //  Allocate the space required to render the image
    let storage_width = supersample * image_width + 2 * gutter_width;
    let storage_height = supersample * image_height + 2 * gutter_width;

    let mut accumulate = S::accumulator_storage(storage_width.usize(), storage_height.usize());

    //  Batch loop - outermost
    for batch_num in 0..num_batches {
        log::trace!("Rendering batch {} of {}", batch_num, num_batches);
        let mut sample_density = 0.0;
        let de_time = frame.time + temporal_deltas[(batch_num * num_temporal_samples).usize()];

        let mut buckets = S::bucket_storage(storage_width.usize(), storage_height.usize());

        //  interpolate and get a control point
        //  ONLY FOR DENSITY FILTER WIDTH PURPOSES
        //  additional interpolation will be done in the temporal_sample loop
        let cp = flam3_interpolate(&frame.genomes, de_time, 0.0)?;

        //  if instructed to by the genome, create the density estimation
        //  filter kernels.  Check boundary conditions as well.
        if cp.estimator_radius < 0.0 || cp.estimator_minimum < 0.0 {
            return Err("density estimator filter widths must be >= 0".to_string());
        }

        // Create DE filters
        let de = if cp.estimator_radius > 0.0 {
            flam3_create_de_filters(
                cp.estimator_radius,
                cp.estimator_minimum,
                cp.estimator_curve,
                supersample,
            )?
        } else {
            Flam3DeHelper::default()
        };

        let mut ppux = 0.0;
        let mut ppuy = 0.0;

        //  Temporal sample loop
        for temporal_sample_num in 0..num_temporal_samples {
            log::trace!(
                "Rendering temporal sample {} of {}",
                temporal_sample_num,
                num_temporal_samples
            );

            let run_state = state!({
                step: "iterating",
                batch: batch_num,
                temporal_sample: temporal_sample_num,
            });

            let color_scalar =
                temporal_filter[(batch_num * num_temporal_samples + temporal_sample_num).usize()];

            let temporal_sample_time = frame.time
                + temporal_deltas[(batch_num * num_temporal_samples + temporal_sample_num).usize()];

            //  Interpolate and get a control point
            let cp = flam3_interpolate(&frame.genomes, temporal_sample_time, 0.0)?;

            let dmap: Vec<Srgba<f64>> = cp
                .palette
                .iter()
                .map(|color| {
                    Srgba::<f64>::from_components((
                        color.red * WHITE_LEVEL.f64() * color_scalar,
                        color.green * WHITE_LEVEL.f64() * color_scalar,
                        color.blue * WHITE_LEVEL.f64() * color_scalar,
                        color.alpha * WHITE_LEVEL.f64() * color_scalar,
                    ))
                })
                .collect();

            //  compute camera
            if cp.sample_density <= 0.0 {
                return Err("Sample density (quality) must be greater than zero".to_string());
            }

            let scale = pow!(2.0_f64, cp.zoom);
            sample_density = cp.sample_density * sqr!(scale);

            ppux = cp.pixels_per_unit * scale;
            ppuy = if field != Field::Both {
                ppux / 2.0
            } else {
                ppux
            };
            ppux /= frame.pixel_aspect_ratio;
            let mut shift = match field {
                Field::Both => 0.0,
                // Field::Even => -0.5,
                Field::Odd => 0.5,
            };

            shift /= ppux;
            let t0 = gutter_width.f64() / (supersample.f64() * ppux).f64();
            let t1 = gutter_width.f64() / (supersample.f64() * ppuy).f64();
            let corner0 = cp.center.x - image_width.f64() / ppux / 2.0;
            let corner1 = cp.center.y - image_height.f64() / ppuy / 2.0;
            fic.bounds[0] = corner0 - t0;
            fic.bounds[1] = corner1 - t1 + shift;
            fic.bounds[2] = corner0 + image_width.f64() / ppux + t0;
            fic.bounds[3] = corner1 + image_height.f64() / ppuy + t1 + shift;
            fic.size[0] = 1.0 / (fic.bounds[2] - fic.bounds[0]);
            fic.size[1] = 1.0 / (fic.bounds[3] - fic.bounds[1]);

            let rot_x = cos!(cp.rotate * 2.0 * PI / 360.0);
            let rot_y = -sin!(cp.rotate * 2.0 * PI / 360.0);
            fic.rot = [
                [rot_x, -rot_y],
                [rot_y, rot_x],
                [cp.rot_center.x, cp.rot_center.y],
            ]
            .into();

            fic.ws0 = storage_width.f64() * fic.size[0];
            fic.wb0s0 = fic.ws0 * fic.bounds[0];
            fic.hs1 = storage_height.f64() * fic.size[1];
            fic.hb1s1 = fic.hs1 * fic.bounds[1];

            //  number of samples is based only on the output image size
            let nsamples = sample_density * image_width.f64() * image_height.f64();

            //  how many of these samples are rendered in this loop?
            let batch_size = nsamples / (num_batches * num_temporal_samples).f64();

            //  Fill in the iter constants
            fic.batch_size = (batch_size / frame.num_threads.f64()).u32();
            fic.temporal_sample_num = temporal_sample_num;
            fic.ntemporal_samples = num_temporal_samples;
            fic.batch_num = batch_num;
            fic.nbatches = num_batches;

            fic.dmap = dmap;
            fic.color_scalar = color_scalar;

            //  Initialize the thread helper structures
            let mut fth: Vec<Flam3ThreadHelper> = Vec::new();
            for _ in 0..frame.num_threads {
                fth.push(Flam3ThreadHelper {
                    rng: Flam3Rng::from_rng(&mut frame.rng),
                    cp: cp.clone(),
                    fic: fic.clone(),
                });
            }

            for (_index, helper) in fth.drain(..).enumerate() {
                // TODO actually use threads
                iter_thread::<S>(
                    helper,
                    &mut buckets,
                    state!(run_state, {
                        thread: _index,
                    }),
                )?;
            }

            if fic.aborted {
                return Err("Aborted".to_string());
            }

            vibrancy += cp.vibrancy;
            gamma += cp.gamma;
            background.red += cp.background.red;
            background.green += cp.background.green;
            background.blue += cp.background.blue;
            vib_gam_n += 1.0;
        }

        log::trace!("Temporal samples complete");

        let k1 = (cp.contrast
            * cp.brightness
            * PREFILTER_WHITE.f64()
            * 268.0
            * batch_filter[batch_num.usize()])
            / 256.0;
        let k2 = ((sqr!(supersample) * num_batches).f64() * ppux * ppuy)
            / (cp.contrast
                * image_width.f64()
                * image_height.f64()
                * WHITE_LEVEL.f64()
                * sample_density
                * sumfilt);

        if de.max_filter_index == 0 {
            for j in 0..buckets.height() {
                for i in 0..buckets.width() {
                    let mut pixel = buckets[(i, j)].accumulator();

                    if pixel.alpha == 0.0 {
                        continue;
                    }

                    let ls = (k1 * ln!(1.0 + pixel.alpha * k2)) / pixel.alpha;
                    pixel *= ls;

                    S::abump_no_overflow(&mut accumulate[(i, j)], &pixel);
                }
            }
        } else {
            let myspan = buckets.height().u32() - 2 * (supersample - 1) + 1;
            let swath = myspan / frame.num_threads.u32();

            //  Create the de helper structures
            let mut deth: Vec<Flam3DeThreadHelper> = Vec::new();

            for i in 0..frame.num_threads {
                let start_row;
                let end_row;
                if frame.num_threads > myspan.usize() {
                    //  More threads than rows
                    start_row = 0;
                    if i == frame.num_threads - 1 {
                        end_row = myspan.i32();
                    } else {
                        end_row = -1;
                    }
                } else {
                    //  Normal case
                    start_row = i.u32() * swath;
                    if i == frame.num_threads - 1 {
                        end_row = buckets.height().i32(); //myspan.i32();
                    } else {
                        end_row = ((i.u32() + 1) * swath).i32();
                    }
                }

                //  Set up the contents of the helper structure
                deth.push(Flam3DeThreadHelper {
                    supersample,
                    de: de.clone(),
                    k1,
                    k2,
                    curve: cp.estimator_curve,
                    start_row,
                    end_row,
                })
            }

            for thread_helper in deth {
                // TODO Actually use threads
                de_thread::<S>(thread_helper, &buckets, &mut accumulate);
            }
        }
    }

    log::trace!("Batches complete");

    //  filter the accumulation buffer down into the image
    let g = 1.0 / (gamma / vib_gam_n);

    let linrange = cp.gamma_threshold;

    vibrancy /= vib_gam_n.f64();
    background.red /= vib_gam_n / 256.0;
    background.green /= vib_gam_n / 256.0;
    background.blue /= vib_gam_n / 256.0;

    //  If we're in the early clip mode, perform this first step to
    //  apply the gamma correction and clipping before the spatial filter

    if frame.earlyclip {
        log::trace!("Applying gamma correction and clipping");

        for j in 0..accumulate.height() {
            for i in 0..accumulate.width() {
                let accumulator = &mut accumulate[(i, j)];

                let (alpha, ls) = if accumulator.alpha.f64() <= 0.0 {
                    (0.0, 0.0)
                } else {
                    let tmp = accumulator.alpha.f64() / PREFILTER_WHITE.f64();
                    let mut alpha = flam3_calc_alpha(tmp, g, linrange);
                    let ls = vibrancy * 256.0 * alpha / tmp;
                    if alpha < 0.0 {
                        alpha = 0.0;
                    }
                    if alpha > 1.0 {
                        alpha = 1.0;
                    }
                    (alpha, ls)
                };

                let t = [
                    accumulator.red.f64(),
                    accumulator.green.f64(),
                    accumulator.blue.f64(),
                ];
                let newrgb = flam3_calc_newrgb(&t, ls, highpow);

                let bg_color = background.as_raw::<[f64]>();
                for rgbi in 0..3 {
                    let mut a = newrgb[rgbi];
                    let t_val: f64 = t[rgbi];
                    a += (1.0 - vibrancy) * 256.0 * pow!(t_val / PREFILTER_WHITE.f64(), g);
                    if frame.channels <= 3 || !frame.transparency {
                        a += (1.0 - alpha) * bg_color[rgbi];
                    } else if alpha > 0.0 {
                        a /= alpha;
                    } else {
                        a = 0.0;
                    }

                    //  Clamp here to ensure proper filter functionality
                    if a > 255.0 {
                        a = 255.0;
                    }
                    if a < 0.0 {
                        a = 0.0;
                    }

                    //  Replace values in accumulation buffer with these new ones
                    accumulator[rgbi] = S::into_accumulator(a);
                }

                accumulator.alpha = S::into_accumulator(alpha);
            }
        }
    }

    //  Apply the spatial filter
    log::trace!("Applying spatial filter");
    let mut y = de_offset;
    for j in 0..image_height {
        let mut x = de_offset;
        for i in 0..image_width {
            let mut t = [0.0; 4];
            for ii in 0..spatial_filter.width() {
                for jj in 0..spatial_filter.height() {
                    let k = spatial_filter[(ii, jj)];
                    let ac = &mut accumulate[(x + ii, y + jj)];

                    t[0] += k * ac.red.f64();
                    t[1] += k * ac.green.f64();
                    t[2] += k * ac.blue.f64();
                    t[3] += k * ac.alpha.f64();
                }
            }

            let p = &mut buffer
                [(frame.channels * frame.bytes_per_channel * (i + j * out_width)).usize()..];

            //  The old way, spatial filter first and then clip after gamma
            if !frame.earlyclip {
                let tmp = t[3] / PREFILTER_WHITE.f64();

                let (alpha, ls) = if t[3] <= 0.0 {
                    (0.0, 0.0)
                } else {
                    let mut alpha = flam3_calc_alpha(tmp, g, linrange);
                    let ls = vibrancy * 256.0 * alpha / tmp;
                    if alpha < 0.0 {
                        alpha = 0.0;
                    }
                    if alpha > 1.0 {
                        alpha = 1.0;
                    }
                    (alpha, ls)
                };

                let newrgb = flam3_calc_newrgb(&[t[0], t[1], t[2]], ls, highpow);

                let bg_color = background.as_raw::<[f64]>();
                for rgbi in 0..3 {
                    let mut a = newrgb[rgbi];
                    a += (1.0 - vibrancy) * 256.0 * pow!(t[rgbi] / PREFILTER_WHITE.f64(), g);
                    if frame.channels <= 3 || !frame.transparency {
                        a += (1.0 - alpha) * bg_color[rgbi];
                    } else if alpha > 0.0 {
                        a /= alpha;
                    } else {
                        a = 0.0;
                    }

                    //  Clamp here to ensure proper filter functionality
                    if a > 255.0 {
                        a = 255.0;
                    }
                    if a < 0.0 {
                        a = 0.0;
                    }

                    //  Replace values in accumulation buffer with these new ones
                    t[rgbi] = a;
                }
                t[3] = alpha;
            }

            for (rgbi, val) in t.iter().enumerate().take(3) {
                let mut a = if *val > 255.0 {
                    255.0
                } else if *val < 0.0 {
                    0.0
                } else {
                    *val
                };

                let (pos, bytes) = if frame.bytes_per_channel == 2 {
                    a *= 256.0; //  Scales to [0-65280]
                    let bytes = a.u16().to_le_bytes();
                    (rgbi * 2, vec![bytes[0], bytes[1]])
                } else {
                    (rgbi, vec![a.u8()])
                };

                p[pos..pos + bytes.len()].clone_from_slice(&bytes);
            }

            if t[3] > 1.0 {
                t[3] = 1.0;
            }
            if t[3] < 0.0 {
                t[3] = 0.0;
            }

            //  alpha
            if frame.channels > 3 {
                let (pos, bytes) = if frame.transparency {
                    if frame.bytes_per_channel == 2 {
                        let bytes = (t[3] * 65535.0).u16().to_le_bytes();
                        (6, vec![bytes[0], bytes[1]])
                    } else {
                        (3, vec![(t[3] * 255.0).u8()])
                    }
                } else if frame.bytes_per_channel == 2 {
                    let bytes = 65535_u16.to_le_bytes();
                    (6, vec![bytes[0], bytes[1]])
                } else {
                    (3, vec![255])
                };

                p[pos..pos + bytes.len()].clone_from_slice(&bytes);
            }

            x += supersample.usize();
        }
        y += supersample.usize();
    }

    Ok(())
}
