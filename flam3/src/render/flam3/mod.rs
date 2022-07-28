mod filters;
mod rect;
mod rng;
mod storage;
mod thread;
mod variations;

use std::thread::available_parallelism;

use image::RgbaImage;
use rand::RngCore;

use crate::math::{ln, pow};
use crate::{utils::PanicCast, Affine, Genome, Palette, Transform};
use crate::{Coordinate, PaletteMode};
pub(crate) use storage::Accumulator;

use self::filters::DensityEstimatorFilters;
use self::storage::atomic::{
    RenderStorageAtomicDouble, RenderStorageAtomicFloat, RenderStorageAtomicInt,
};
use self::storage::sync::{RenderStorageSyncDouble, RenderStorageSyncFloat, RenderStorageSyncInt};
use self::storage::RenderStorage;
use self::variations::VariationPrecalculations;
use self::{rect::render_rectangle, rng::IsaacRng};

use super::{Buffers, RenderOptions, ThreadingMode};

pub const CHOOSE_XFORM_GRAIN: usize = 16384;
pub const CHOOSE_XFORM_GRAIN_M1: usize = 16383;

fn adjust_percentage(perc: f64) -> f64 {
    if perc == 0.0 {
        0.0
    } else {
        pow!(10.0_f64, -ln!(1.0 / perc) / ln!(2.0_f64))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LineField {
    Both,
    Odd,
}

trait ClonableRng: RngCore + Clone {}

struct Frame {
    rng: IsaacRng,
    genomes: Vec<Genome>,
    num_threads: usize,
    pixel_aspect_ratio: f64,
    time: f64,
    earlyclip: bool,
    sub_batch_size: u32,
    transparency: bool,
    channels: u32,
    bytes_per_channel: u32,
}

#[derive(Clone)]
pub(crate) struct IterationContext {
    bounds: [Coordinate<f64>; 2], //  Corner coords of viewable area
    rot: Affine,                  //  Rotation transformation
    ws0: f64,
    wb0s0: f64,
    hs1: f64,
    hb1s1: f64,    //  shortcuts for indexing
    dmap: Palette, //  palette
    batch_size: u32,
    skip_iterations: u32,
    sub_batch_size: u32,
    rotate: f64,
    rot_center: Coordinate<f64>,
    palette_mode: PaletteMode,
}

pub(crate) struct DensityEstimationContext {
    supersample: u32,
    de_filters: DensityEstimatorFilters,
    k1: f64,
    k2: f64,
    curve: f64,
    start_row: u32,
    end_row: u32,
}

pub(crate) fn render(genome: Genome, options: RenderOptions) -> Result<RgbaImage, String> {
    let rng = if let Some(ref seed) = options.isaac_seed {
        IsaacRng::from_str(seed)
    } else {
        Default::default()
    };

    let num_threads = options
        .threads
        .and_then(|v| if v == 0 { None } else { Some(v) })
        .unwrap_or_else(|| available_parallelism().map(|s| s.get()).unwrap_or(1));

    if num_threads == 1 || options.threading_mode == Some(ThreadingMode::Sync) {
        match options.buffers {
            Buffers::Int => {
                render_internal::<RenderStorageSyncInt>(genome, options, num_threads, rng)
            }
            Buffers::Float => {
                render_internal::<RenderStorageSyncFloat>(genome, options, num_threads, rng)
            }
            Buffers::Double => {
                render_internal::<RenderStorageSyncDouble>(genome, options, num_threads, rng)
            }
        }
    } else {
        match options.buffers {
            Buffers::Int => {
                render_internal::<RenderStorageAtomicInt>(genome, options, num_threads, rng)
            }
            Buffers::Float => {
                render_internal::<RenderStorageAtomicFloat>(genome, options, num_threads, rng)
            }
            Buffers::Double => {
                render_internal::<RenderStorageAtomicDouble>(genome, options, num_threads, rng)
            }
        }
    }
}

fn render_internal<S: RenderStorage>(
    mut genome: Genome,
    options: RenderOptions,
    num_threads: usize,
    rng: IsaacRng,
) -> Result<RgbaImage, String> {
    let num_strips = options.num_strips.unwrap_or(1);
    log::trace!(
        "Starting render pass, width={}, height={}, channels={}, density={}, supersample={}, threads={}, model={}, strips={}",
        genome.size.width,
        genome.size.height,
        options.channels,
        genome.sample_density,
        genome.spatial_supersample,
        num_threads,
        S::threading_model(),
        num_strips
    );

    if num_strips > genome.size.height {
        return Err(format!(
            "Cannot have more strips than rows but {}>{}",
            num_strips, genome.size.height
        ));
    }

    let img_mem =
        options.channels * genome.size.height * genome.size.width * options.bytes_per_channel;
    let mut image: Vec<u8> = vec![0; img_mem.usize()];

    let full_height = genome.size.height;
    let strip_height = (genome.size.height.f64() / num_strips.f64()).ceil().u32();
    let center_y = genome.center.y;
    let zoom_scale = pow!(2.0_f64, genome.zoom);
    let center_base = center_y
        - ((num_strips - 1).f64() * strip_height.f64())
            / (2.0 * genome.pixels_per_unit * zoom_scale);

    if let Some(index) = genome.palette_index {
        genome.palette = options.palette(index)?;
    }

    for strip in 0..num_strips {
        log::trace!("Rendering strip {} of {}", strip, num_strips);
        let mut genome = genome.clone();

        //  Force ntemporal_samples to 1 for -render
        genome.num_temporal_samples = 1;
        genome.sample_density *= num_strips.f64();
        genome.size.height = strip_height;

        let ssoff =
            strip_height * strip * genome.size.width * options.channels * options.bytes_per_channel;
        genome.center.y =
            center_base + (strip_height * strip).f64() / (genome.pixels_per_unit * zoom_scale);

        let buffer = &mut image[ssoff.usize()..];

        if strip_height * (strip + 1) > full_height {
            genome.size.height = full_height - strip_height * strip;
            genome.center.y -= (strip_height - genome.size.height).f64() * 0.5
                / (genome.pixels_per_unit * zoom_scale);
        }

        let frame = Frame {
            rng: rng.clone(),
            genomes: vec![genome.clone()],
            time: 0.0,
            pixel_aspect_ratio: options.pixel_aspect_ratio,
            num_threads,
            earlyclip: options.earlyclip,
            sub_batch_size: options.sub_batch_size,
            transparency: options.transparency,
            channels: options.channels,
            bytes_per_channel: options.bytes_per_channel,
        };

        render_rectangle::<S>(frame, buffer, LineField::Both)?;
    }

    log::trace!("Strips complete");

    Ok(RgbaImage::from_raw(genome.size.width, genome.size.height, image).unwrap())
}

fn flam3_interpolate(genomes: &[Genome], _time: f64, _stagger: f64) -> Result<Genome, String> {
    if genomes.len() == 1 {
        return Ok(genomes[0].clone());
    }

    unimplemented!();
}

fn flam3_create_chaos_distrib(
    transforms: &[Transform],
    xform_distrib: &mut [usize],
) -> Result<(), String> {
    let mut dr = 0.0;
    for xform in transforms {
        let d = xform.density;

        if d < 0.0 {
            return Err("transform weight must be non-negative".to_string());
        }

        dr += d;
    }

    if dr == 0.0 {
        return Err("cannot iterate empty flame".to_string());
    }

    dr /= CHOOSE_XFORM_GRAIN.f64();

    let mut j = 0;
    let mut t = transforms[0].density;

    let mut r = 0.0;
    for distrib_val in xform_distrib.iter_mut().take(CHOOSE_XFORM_GRAIN) {
        while r >= t {
            j += 1;

            t += transforms[j].density;
        }

        *distrib_val = j;
        r += dr;
    }

    Ok(())
}

#[derive(Clone)]
pub(crate) struct TransformSelector {
    xform_distrib: Vec<usize>,
    xforms: Vec<Transform>,
    precalcs: Vec<VariationPrecalculations>,
}

impl TransformSelector {
    fn new(cp: &Genome) -> Result<Self, String> {
        let mut xform_distrib = vec![0_usize; CHOOSE_XFORM_GRAIN];

        /* First, set up the first row of the xform_distrib (raw weights) */
        flam3_create_chaos_distrib(&cp.transforms, &mut xform_distrib)?;

        /* Check for non-unity chaos */
        let chaos_enable = false;

        if chaos_enable {
            /* Now set up a row for each of the xforms */
            for i in 0..cp.transforms.len() {
                flam3_create_chaos_distrib(
                    &cp.transforms,
                    &mut xform_distrib[CHOOSE_XFORM_GRAIN * i..],
                )?;
            }
        }

        Ok(Self {
            xform_distrib,
            xforms: cp.transforms.clone(),
            precalcs: cp
                .transforms
                .iter()
                .map(VariationPrecalculations::new)
                .collect(),
        })
    }

    fn next<'a>(&'a self, rng: &mut IsaacRng) -> (&'a Transform, &'a VariationPrecalculations) {
        let rand = rng.next_u32();
        let dist_index = rand.usize() & CHOOSE_XFORM_GRAIN_M1;
        let xform_index = self.xform_distrib[dist_index];

        (&self.xforms[xform_index], &self.precalcs[xform_index])
    }
}
