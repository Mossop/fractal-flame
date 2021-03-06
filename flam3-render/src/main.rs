use std::{
    error::Error,
    fs::File,
    io::{stdin, BufReader},
};

use clap::Parser;
use flam3::{flam3_from_reader, render, RenderOptions};
use flexi_logger::Logger;
use image::ImageFormat;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    file: Option<String>,

    #[clap(long)]
    out: Option<String>,

    #[clap(long)]
    seed: Option<String>,

    #[clap(long)]
    earlyclip: bool,
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

    for (index, genome) in genomes.into_iter().enumerate() {
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
