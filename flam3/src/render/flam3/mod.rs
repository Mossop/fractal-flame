mod filters;
mod rect;
mod rng;
mod thread;
mod variations;

use std::fmt::Display;

use image::RgbaImage;
use rand::RngCore;

use crate::{ln, pow, utils::PanicCast, Affine, Genome, Palette, Transform};

use self::{rect::render_rectangle, rng::Flam3Rng};

use super::{Buffers, RenderOptions};

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
enum Field {
    Both,
    Odd,
}

trait ClonableRng: RngCore + Clone {}

trait RenderOps: Clone + Default {
    type Bucket: PanicCast + Default + Display + Clone + Copy;
    type Accumulator: PanicCast + Default + Display + Clone + Copy + Into<f64>;

    fn into_accumulator(val: f64) -> Self::Accumulator;
    fn bucket_storage(nbuckets: usize) -> Vec<[Self::Bucket; 5]>;
    fn accumulator_storage(nbuckets: usize) -> Vec<[Self::Accumulator; 4]>;
    fn bump_no_overflow(dest: &mut [Self::Bucket; 5], delta: &[f64; 5]);
    fn abump_no_overflow(dest: &mut [Self::Accumulator; 4], delta: &[f64; 4]);
    fn add_c_to_accum(
        acc: &mut [[Self::Accumulator; 4]],
        i: u32,
        ii: i32,
        j: u32,
        jj: i32,
        wid: u32,
        hgt: u32,
        c: &[f64; 4],
    );
}

#[derive(Clone, Default)]
struct RenderOpsAtomicFloat {}

impl RenderOps for RenderOpsAtomicFloat {
    type Bucket = u32;
    type Accumulator = f32;

    fn bump_no_overflow(dest: &mut [Self::Bucket; 5], delta: &[f64; 5]) {
        for (index, bucket) in dest.iter_mut().enumerate() {
            let result = bucket.f64() + delta[index];
            *bucket = if result >= u32::MAX as f64 {
                u32::MAX
            } else {
                result as u32
            };
        }
    }

    fn abump_no_overflow(dest: &mut [Self::Accumulator; 4], delta: &[f64; 4]) {
        for (index, bucket) in dest.iter_mut().enumerate() {
            *bucket += delta[index].f32();
        }
    }

    fn add_c_to_accum(
        acc: &mut [[Self::Accumulator; 4]],
        i: u32,
        ii: i32,
        j: u32,
        jj: i32,
        wid: u32,
        hgt: u32,
        c: &[f64; 4],
    ) {
        let y = j.i32() + jj;
        let x = i.i32() + ii;
        let width = wid.i32();
        let height = hgt.i32();
        if y >= 0 && y < height && x >= 0 && x < width {
            Self::abump_no_overflow(&mut acc[(x + y * width).usize()], c);
        }
    }

    fn bucket_storage(nbuckets: usize) -> Vec<[Self::Bucket; 5]> {
        vec![[0; 5]; nbuckets]
    }

    fn accumulator_storage(nbuckets: usize) -> Vec<[Self::Accumulator; 4]> {
        vec![[0.0; 4]; nbuckets]
    }

    fn into_accumulator(val: f64) -> Self::Accumulator {
        val.f32()
    }
}

#[derive(Clone)]
struct Flam3Frame {
    rng: Flam3Rng,
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
struct Flam3IterConstants {
    bounds: [f64; 4], //  Corner coords of viewable area
    rot: Affine,      //  Rotation transformation
    size: [f64; 2],
    width: u32,
    height: u32,
    ws0: f64,
    wb0s0: f64,
    hs1: f64,
    hb1s1: f64, //  shortcuts for indexing
    cmap_size: usize,
    dmap: Palette,     //  palette
    color_scalar: f64, //  <1.0 if non-uniform motion blur is set
    badvals: u32,      //  accumulates all badvalue resets
    batch_size: u32,
    temporal_sample_num: u32,
    ntemporal_samples: u32,
    batch_num: u32,
    nbatches: u32,
    aborted: bool,
    spec: Flam3Frame,
}

impl Flam3IterConstants {
    fn new(frame: Flam3Frame) -> Self {
        Self {
            bounds: Default::default(),
            rot: Default::default(),
            size: Default::default(),
            width: Default::default(),
            height: Default::default(),
            ws0: Default::default(),
            wb0s0: Default::default(),
            hs1: Default::default(),
            hb1s1: Default::default(),
            cmap_size: 256,
            dmap: Default::default(),
            color_scalar: Default::default(),
            badvals: Default::default(),
            batch_size: Default::default(),
            temporal_sample_num: Default::default(),
            ntemporal_samples: Default::default(),
            batch_num: Default::default(),
            nbatches: Default::default(),
            aborted: Default::default(),
            spec: frame,
        }
    }
}

#[derive(Default, Debug, Clone)]
struct Flam3DeHelper {
    max_filtered_counts: u32,
    max_filter_index: u32,
    kernel_size: u32,
    filter_widths: Vec<f64>,
    filter_coefs: Vec<f64>,
}

struct Flam3ThreadHelper {
    rng: Flam3Rng, /* Thread-unique rng */
    cp: Genome,    /* Full copy of genome for use by the thread */
    fic: Flam3IterConstants,
}

struct Flam3DeThreadHelper {
    width: u32,
    height: u32,
    oversample: u32,
    de: Flam3DeHelper,
    k1: f64,
    k2: f64,
    curve: f64,
    start_row: u32,
    end_row: i32,
}

pub(crate) fn render(genome: Genome, options: RenderOptions) -> Result<RgbaImage, String> {
    let rng = if let Some(ref seed) = options.isaac_seed {
        Flam3Rng::from_seed(seed)
    } else {
        Default::default()
    };

    match options.buffers {
        Buffers::Int => render_internal::<RenderOpsAtomicFloat>(genome, options, rng),
        Buffers::Float => render_internal::<RenderOpsAtomicFloat>(genome, options, rng),
        Buffers::Double => render_internal::<RenderOpsAtomicFloat>(genome, options, rng),
    }
}

fn render_internal<Ops: RenderOps>(
    genome: Genome,
    options: RenderOptions,
    rng: Flam3Rng,
) -> Result<RgbaImage, String> {
    let num_strips = options.num_strips.unwrap_or(1);
    let num_threads = options.threads.unwrap_or(1);
    log::trace!(
        "Starting render pass, width={}, height={}, channels={}, density={}, oversample={}, threads={}, strips={}",
        genome.size.width,
        genome.size.height,
        options.channels,
        genome.sample_density,
        genome.spatial_oversample,
        num_threads,
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

        let frame = Flam3Frame {
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

        render_rectangle::<Ops>(frame, buffer, Field::Both)?;
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

fn flam3_create_chaos_distrib(cp: &Genome, xform_distrib: &mut [usize]) -> Result<(), String> {
    let num_std = cp.transforms.len();

    let mut dr = 0.0;
    for i in 0..num_std {
        let d = cp.transforms[i].density;

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
    let mut t = cp.transforms[0].density;

    let mut r = 0.0;
    for distrib_val in xform_distrib.iter_mut().take(CHOOSE_XFORM_GRAIN) {
        while r >= t {
            j += 1;

            t += cp.transforms[j].density;
        }

        *distrib_val = j;
        r += dr;
    }

    Ok(())
}

struct TransformSelector {
    xform_distrib: Vec<usize>,
}

impl TransformSelector {
    fn new(cp: &Genome) -> Result<Self, String> {
        let mut xform_distrib = vec![0_usize; CHOOSE_XFORM_GRAIN];

        /* First, set up the first row of the xform_distrib (raw weights) */
        flam3_create_chaos_distrib(cp, &mut xform_distrib)?;

        /* Check for non-unity chaos */
        let chaos_enable = false;

        if chaos_enable {
            /* Now set up a row for each of the xforms */
            for i in 0..cp.transforms.len() {
                flam3_create_chaos_distrib(cp, &mut xform_distrib[CHOOSE_XFORM_GRAIN * i..])?;
            }
        }

        Ok(Self { xform_distrib })
    }

    fn next<'a>(&mut self, cp: &'a Genome, rng: &mut Flam3Rng) -> &'a Transform {
        let rand = rng.irand();
        let dist_index = rand.usize() & CHOOSE_XFORM_GRAIN_M1;
        let xform_index = self.xform_distrib[dist_index];

        &cp.transforms[xform_index]
    }
}
