use std::{borrow::Borrow, collections::HashMap, str::FromStr};

use image::RgbaImage;
use itertools::Itertools;
use palette::{IntoColor, Srgb};
use rust_embed::RustEmbed;
use serde_xml_rs::{EventReader, ParserConfig};
use xml::reader;

use crate::{utils::attr_hash, Genome, Palette};

pub mod flam3;

#[derive(RustEmbed)]
#[folder = "assets/"]
struct Assets;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Buffers {
    Int = 32,
    Double = 64,
    #[default]
    Float = 33,
}

#[derive(Debug, PartialEq, strum_macros::EnumString, strum_macros::Display)]
pub enum ThreadingMode {
    #[strum(serialize = "sync")]
    Sync,
    #[strum(serialize = "atomic")]
    Atomic,
}

#[derive(Debug)]
pub struct RenderOptions {
    pub isaac_seed: Option<String>,
    pub buffers: Buffers,
    pub bytes_per_channel: u32,
    pub channels: u32,
    pub threads: Option<usize>,
    pub threading_mode: Option<ThreadingMode>,
    pub pixel_aspect_ratio: f64,
    pub num_strips: Option<u32>,
    pub sub_batch_size: u32,
    pub earlyclip: bool,
    pub transparency: bool,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            buffers: Default::default(),
            bytes_per_channel: 1,
            channels: 4,
            threads: None,
            threading_mode: None,
            isaac_seed: None,
            pixel_aspect_ratio: 1.0,
            num_strips: None,
            sub_batch_size: 10000,
            earlyclip: false,
            transparency: false,
        }
    }
}

macro_rules! read_xml_event {
    ($parser:ident) => {
        $parser
            .next()
            .map_err(|e| format!("Xml parsing error: {}", e))?
    };
}

fn parse_palettes(path: &str) -> Result<HashMap<usize, Palette>, String> {
    let mut palettes = HashMap::new();

    let file = Assets::get(path).unwrap();
    let mut parser = EventReader::<&[u8]>::new_with_config(
        file.data.borrow(),
        ParserConfig {
            whitespace_to_characters: true,
            cdata_to_characters: true,
            ..Default::default()
        },
    );

    loop {
        match read_xml_event!(parser) {
            reader::XmlEvent::StartElement {
                name,
                attributes: element_attrs,
                ..
            } => {
                if name.local_name == "palette" {
                    let attrs = attr_hash(element_attrs);
                    let index = attrs
                        .get("number")
                        .ok_or_else(|| "Missing palette index".to_string())
                        .and_then(|v| {
                            usize::from_str(v)
                                .map_err(|e| format!("Failed to parse palette index: {}", e))
                        })?;
                    let data = attrs
                        .get("data")
                        .ok_or_else(|| "Missing palette data".to_string())?
                        .split_ascii_whitespace()
                        .join("");

                    let mut palette = Palette::new();

                    for i in (0..data.len()).step_by(8) {
                        let color = Srgb::<u8>::from_str(&data[i + 2..i + 8])
                            .map_err(|e| format!("Failed to convert palette color: {}", e))?;
                        palette.push(color.into_format().into_color());
                    }

                    palettes.insert(index, palette);
                }
            }
            reader::XmlEvent::EndDocument => {
                break;
            }
            _ => {}
        }
    }

    log::trace!("Loaded {} palettes from file", palettes.len());

    Ok(palettes)
}

impl RenderOptions {
    pub fn palette(&self, index: usize) -> Result<Palette, String> {
        let palettes = parse_palettes("flam3-palettes.xml")?;
        if let Some(palette) = palettes.get(&index) {
            Ok(palette.clone())
        } else {
            Err(format!("Unknown palette index {}", index))
        }
    }
}

pub fn render(genome: Genome, options: RenderOptions) -> Result<RgbaImage, String> {
    flam3::render(genome, options)
}
