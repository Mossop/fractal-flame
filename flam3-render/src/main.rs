use std::{
    error::Error,
    fs::File,
    io::{stdin, BufReader},
};

use clap::{Parser, ValueEnum};
use flam3::{flam3_from_reader, render, RenderOptions};
use flexi_logger::Logger;
use image::ImageFormat;

#[derive(Default, Copy, Clone, ValueEnum)]
enum SizeMode {
    #[default]
    Contain,
    Crop,
}

/// A command line tool for rendering flam3 fractal flames.
///
/// Attempts to maintain some semblence of compatibility with the original
/// flam3-render tool from Scott Draves. Many of the same options can be
/// set using the same environment variables however command line argument
/// support is also available.
#[derive(Parser)]
#[clap(author, version, rename_all_env = "verbatim")]
struct Args {
    /// The file to render.
    #[clap(env = "in")]
    file: Option<String>,

    /// File to output to.
    #[clap(long, env)]
    out: Option<String>,

    /// The initial seed for random numbers. Generally only useful for testing.
    #[clap(long, env = "isaac_seed")]
    seed: Option<String>,

    /// How many threads to use. Uses a sane default if absent.
    #[clap(long = "nthreads", env = "nthreads")]
    threads: Option<usize>,

    /// Enables early clipping of rbg values for better antialiasing and resizing.
    #[clap(long, env)]
    earlyclip: bool,

    /// Magnifies the final image by this amount. Applies after any height/width setting.
    #[clap(long = "ss", env = "ss")]
    size_scale: Option<f64>,

    /// Magnifies the quality by this amount.
    #[clap(long = "qs", env = "qs")]
    quality_scale: Option<f64>,

    /// Scales the image to this width.
    #[clap(long)]
    width: Option<u32>,

    /// Scales the image to this height.
    #[clap(long)]
    height: Option<u32>,

    /// How to scale the image when both height and width are given.
    ///
    /// Either the final image will contain the entire bounds of the fractal
    /// plus some more if the aspect is different, or the fractal will be
    /// cropped.
    #[clap(long, value_enum, default_value = "contain")]
    size_mode: SizeMode,
}

fn make_single_name(out: &Option<String>) -> String {
    if let Some(name) = out {
        name.clone()
    } else {
        "out.png".to_string()
    }
}

fn make_name(index: usize, base: &Option<String>) -> String {
    if let Some(name) = base {
        let base = name.trim_end_matches(".png");
        format!("{}{:05}.png", base, index)
    } else {
        format!("{:05}.png", index)
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    if let Err(e) = Logger::try_with_env_or_str("info").and_then(|logger| logger.start()) {
        panic!("Failed to start logging: {}", e);
    }

    let args = Args::parse();

    let genomes = if let Some(name) = args.file {
        flam3_from_reader(BufReader::new(File::open(name)?))?
    } else {
        flam3_from_reader(BufReader::new(stdin()))?
    };
    let count = genomes.len();

    for (index, mut genome) in genomes.into_iter().enumerate() {
        match (args.width, args.height, args.size_mode) {
            (Some(width), None, _) => {
                genome.scale(width as f64 / genome.size.width as f64);
            }
            (None, Some(height), _) => {
                genome.scale(height as f64 / genome.size.height as f64);
            }
            (Some(width), Some(height), SizeMode::Crop) => {
                let wscale = width as f64 / genome.size.width as f64;
                let hscale = height as f64 / genome.size.height as f64;
                if wscale > hscale {
                    genome.scale(wscale);
                    genome.size.height = height;
                } else {
                    genome.scale(hscale);
                    genome.size.width = width;
                }
            }
            (Some(width), Some(height), SizeMode::Contain) => {
                let wscale = width as f64 / genome.size.width as f64;
                let hscale = height as f64 / genome.size.height as f64;
                if wscale < hscale {
                    genome.scale(wscale);
                    genome.size.height = height;
                } else {
                    genome.scale(hscale);
                    genome.size.width = width;
                }
            }
            (None, None, _) => {}
        }

        if let Some(scale) = args.size_scale {
            genome.scale(scale);
        }

        if let Some(scale) = args.quality_scale {
            genome.sample_density *= scale;
        }

        let name = genome
            .name
            .clone()
            .unwrap_or_else(|| format!("genome {}", index));
        eprintln!("Rendering {}...", name);

        let data = render(
            genome,
            RenderOptions {
                isaac_seed: args.seed.clone(),
                earlyclip: args.earlyclip,
                threads: args.threads,
                ..Default::default()
            },
        )?;

        let name = if count == 1 {
            make_single_name(&args.out)
        } else {
            make_name(index, &args.out)
        };
        eprintln!("Writing {}...", name);

        data.save_with_format(name, ImageFormat::Png)?;
    }

    Ok(())
}
