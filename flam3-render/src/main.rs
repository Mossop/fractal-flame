use std::{
    error::Error,
    fs::File,
    io::{stdin, BufReader},
};

use clap::Parser;
use flam3::{flam3_from_reader, render, RenderOptions};
use flexi_logger::Logger;
use image::ImageFormat;

/// A command line tool for rendering flam3 fractal flames.
///
/// Attempts to maintain some semblence of compatibility with the original
/// flam3-render tool from Scott Draves. Many of the same options can be
/// set using the same environment variables however command line argument
/// support is also available.
#[derive(Parser, Debug)]
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

    /// Magnifies the final image by this amount.
    #[clap(long = "ss", env = "ss")]
    size_scale: Option<f64>,

    /// Magnifies the quality by this amount.
    #[clap(long = "qs", env = "qs")]
    quality_scale: Option<f64>,
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
