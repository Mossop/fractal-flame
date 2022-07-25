use image::RgbaImage;

use crate::Genome;

pub mod flam3;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Buffers {
    Int = 32,
    Double = 64,
    #[default]
    Float = 33,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RenderOptions {
    pub isaac_seed: Option<String>,
    pub buffers: Buffers,
    pub bytes_per_channel: u32,
    pub channels: u32,
    pub threads: Option<usize>,
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
            isaac_seed: None,
            pixel_aspect_ratio: 1.0,
            num_strips: None,
            sub_batch_size: 10000,
            earlyclip: false,
            transparency: false,
        }
    }
}

pub fn render(genome: Genome, options: RenderOptions) -> Result<RgbaImage, String> {
    flam3::render(genome, options)
}
