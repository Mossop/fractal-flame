pub mod file;
pub(crate) mod math;
pub(crate) mod rect;
pub mod render;
pub(crate) mod utils;
pub mod variations;

use std::{cmp::Ordering, f64::consts::PI, ops::Index};

pub use file::flam3::{flam3_from_reader, flam3_to_writer};
use math::{atan2, cos, sin};
use palette::Srgba;
pub use render::{render, RenderOptions};
use utils::{parse, PanicCast};
use uuid::Uuid;
use variations::Var;

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Dimension<T> {
    pub width: T,
    pub height: T,
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Coordinate<T> {
    pub x: T,
    pub y: T,
}

impl<T> From<&[T; 2]> for Coordinate<T>
where
    T: Copy,
{
    fn from(tuple: &[T; 2]) -> Self {
        Self {
            x: tuple[0],
            y: tuple[1],
        }
    }
}

impl<T> From<&[T; 4]> for Coordinate<T>
where
    T: Copy,
{
    fn from(tuple: &[T; 4]) -> Self {
        Self {
            x: tuple[0],
            y: tuple[1],
        }
    }
}

impl<T> From<Coordinate<T>> for (T, T) {
    fn from(c: Coordinate<T>) -> Self {
        (c.x, c.y)
    }
}

#[derive(
    Debug, Clone, Copy, PartialEq, Default, strum_macros::EnumString, strum_macros::Display,
)]
pub enum Interpolation {
    #[default]
    #[strum(serialize = "linear")]
    Linear,
    #[strum(serialize = "smooth")]
    Smooth,
}

#[derive(
    Debug, Clone, Copy, PartialEq, Default, strum_macros::EnumString, strum_macros::Display,
)]
pub enum PaletteInterpolation {
    #[strum(serialize = "hsv")]
    Hsv,
    #[strum(serialize = "sweep")]
    Sweep,
    #[default]
    #[strum(serialize = "hsv_circular")]
    HsvCircular,
    #[strum(serialize = "rgb")]
    Rgb,
}

#[derive(
    Debug, Clone, Copy, PartialEq, Default, strum_macros::EnumString, strum_macros::Display,
)]
pub enum InterpolationType {
    #[strum(serialize = "linear")]
    Linear,
    #[default]
    #[strum(serialize = "log")]
    Log,
    #[strum(serialize = "old")]
    Old,
    #[strum(serialize = "older")]
    Older,
}

#[derive(
    Debug, Clone, Copy, PartialEq, Default, strum_macros::EnumString, strum_macros::Display,
)]
pub enum SpatialFilter {
    #[default]
    #[strum(serialize = "gaussian")]
    Gaussian,
    #[strum(serialize = "hermite")]
    Hermite,
    #[strum(serialize = "box")]
    Box,
    #[strum(serialize = "triangle")]
    Triangle,
    #[strum(serialize = "bell")]
    Bell,
    #[strum(serialize = "bspline")]
    BSpline,
    #[strum(serialize = "lanczos3")]
    Lanczos3,
    #[strum(serialize = "lanczos2")]
    Lanczos2,
    #[strum(serialize = "mitchell")]
    Mitchell,
    #[strum(serialize = "blackman")]
    Blackman,
    #[strum(serialize = "catrom")]
    Catrom,
    #[strum(serialize = "hamming")]
    Hamming,
    #[strum(serialize = "hanning")]
    Hanning,
    #[strum(serialize = "quadratic")]
    Quadratic,
}

#[derive(
    Debug, Clone, Copy, PartialEq, Default, strum_macros::EnumString, strum_macros::Display,
)]
pub enum TemporalFilterType {
    #[default]
    #[strum(serialize = "box")]
    Box,
    #[strum(serialize = "gaussian")]
    Gaussian,
    #[strum(serialize = "exp")]
    Exp,
}

#[derive(
    Debug, Clone, Copy, PartialEq, Default, strum_macros::EnumString, strum_macros::Display,
)]
pub enum PaletteMode {
    #[default]
    #[strum(serialize = "step")]
    Step,
    #[strum(serialize = "linear")]
    Linear,
}

#[derive(
    Debug, Clone, Copy, PartialEq, Default, strum_macros::EnumString, strum_macros::Display,
)]
pub enum MotionFunction {
    #[default]
    #[strum(serialize = "sin")]
    Sin,
    #[strum(serialize = "triangle")]
    Triangle,
    #[strum(serialize = "hill")]
    Hill,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Affine {
    coefficients: [[f64; 2]; 3],
}

impl Affine {
    pub fn transform(&self, p: Coordinate<f64>) -> Coordinate<f64> {
        Coordinate {
            x: self.coefficients[0][0] * p.x
                + self.coefficients[1][0] * p.y
                + self.coefficients[2][0],
            y: self.coefficients[0][1] * p.x
                + self.coefficients[1][1] * p.y
                + self.coefficients[2][1],
        }
    }

    fn is_identity(&self) -> bool {
        (self.coefficients[0][0] == 1.0)
            && (self.coefficients[0][1] == 0.0)
            && (self.coefficients[1][0] == 0.0)
            && (self.coefficients[1][1] == 1.0)
            && (self.coefficients[2][0] == 0.0)
            && (self.coefficients[2][1] == 0.0)
    }
}

impl From<[[f64; 2]; 3]> for Affine {
    fn from(coefficients: [[f64; 2]; 3]) -> Self {
        Self { coefficients }
    }
}

impl Default for Affine {
    fn default() -> Self {
        Self {
            coefficients: [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]],
        }
    }
}

impl Index<usize> for Affine {
    type Output = [f64; 2];

    fn index(&self, index: usize) -> &Self::Output {
        &self.coefficients[index]
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Transform {
    id: Uuid,
    /// The probability that this transform will be selected.
    ///
    /// Ignored for final transforms.
    density: f64,
    color_speed: f64,
    // Is this a bool?
    animate: f64,
    motion_frequency: f64,
    motion_function: Option<MotionFunction>,
    color: f64,
    opacity: f64,

    /// Defines the affine transform for the coordinates before being passed to
    /// the variation functions.
    coefficients: Affine,
    /// Defines the affine transform for the coordinates generated after
    /// applying all of the variation functions.
    post: Affine,

    variations: Vec<Var>,
}

impl Default for Transform {
    fn default() -> Self {
        Self {
            id: Uuid::new_v4(),
            density: 0.0,
            color_speed: 0.5,
            animate: 1.0,
            motion_frequency: 0.0,
            motion_function: None,
            color: 0.0,
            opacity: 1.0,
            coefficients: Default::default(),
            post: Default::default(),
            variations: Vec::new(),
        }
    }
}

pub type Palette = Vec<Srgba<f64>>;

#[derive(Debug, Clone, PartialEq)]
pub struct Genome {
    pub name: Option<String>,
    pub time: f64,
    pub palette: Palette,
    pub hsv_rgb_palette_blend: f64,
    pub interpolation: Interpolation,
    pub palette_interpolation: PaletteInterpolation,
    pub interpolation_type: InterpolationType,
    pub palette_index: Option<usize>,

    /// The size of the generated image in pixels.
    pub size: Dimension<u32>,
    /// The center of the image, coordinates must be in the range [-1,1].
    pub center: Coordinate<f64>,
    /// The center point of any rotation.
    pub rot_center: Coordinate<f64>,
    /// Rotation in degrees.
    pub rotate: f64,
    pub pixels_per_unit: f64,
    pub zoom: f64,
    /// The number of points to capture per pixel
    pub spatial_supersample: u32,
    pub spatial_filter: SpatialFilter,
    pub spatial_filter_radius: f64,
    pub temporal_filter: TemporalFilterType,
    pub temporal_filter_width: f64,
    pub temporal_filter_exp: f64,
    pub palette_mode: PaletteMode,
    pub sample_density: f64,
    pub passes: u32,
    pub num_temporal_samples: u32,
    pub background: Srgba<f64>,
    pub brightness: f64,
    pub gamma: f64,
    pub highlight_power: f64,
    pub contrast: f64,
    pub vibrancy: f64,
    pub hue_rotation: f64,
    pub estimator_radius: f64,
    pub estimator_minimum: f64,
    pub estimator_curve: f64,
    pub gamma_threshold: f64,

    /// A list of transforms. Each iteration a random transform is selected
    /// with its density indicating the probability of selection.
    pub transforms: Vec<Transform>,
    /// An optional final transform that is applied every iteration after the
    /// randomly selected transform.
    pub final_transform: Option<Transform>,
}

fn round6(mut x: f64) -> f64 {
    x *= 1e6;

    if x < 0.0 {
        x -= 1.0;
    }

    1e-6 * (x + 0.5).trunc()
}

fn det_matrix(c: &[[f64; 2]; 3]) -> f64 {
    c[0][0] * c[1][1] - c[0][1] * c[1][0]
}

fn compare_xforms(a: &Transform, b: &Transform) -> Ordering {
    let mut ad = det_matrix(&a.coefficients.coefficients);
    let mut bd = det_matrix(&b.coefficients.coefficients);

    if a.color_speed > b.color_speed {
        return Ordering::Greater;
    }

    if a.color_speed < b.color_speed {
        return Ordering::Less;
    }

    if a.color_speed != 0.0 {
        if ad < 0.0 {
            return Ordering::Less;
        }
        if bd < 0.0 {
            return Ordering::Greater;
        }
        ad = atan2!(a.coefficients[0][0], a.coefficients[0][1]);
        bd = atan2!(b.coefficients[0][0], b.coefficients[0][1]);
    }

    if ad < bd {
        Ordering::Less
    } else if ad > bd {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

impl Genome {
    pub fn scale(&mut self, scale: f64) {
        self.size = Dimension {
            width: (self.size.width as f64 * scale) as u32,
            height: (self.size.height as f64 * scale) as u32,
        };
        self.pixels_per_unit *= scale;
    }

    pub fn add_transform(&mut self, transform: Transform) {
        self.transforms.push(transform);
    }

    pub fn add_symmetry(&mut self, mut kind: i32) {
        let mut new_transforms = Vec::new();

        if kind < 0 {
            new_transforms.push(Transform {
                density: 1.0,
                color_speed: 0.0,
                animate: 0.0,
                variations: vec![Var::Linear(Default::default())],
                color: 1.0,
                coefficients: [[-1.0, 0.0], [0.0, 1.0], [0.0, 0.0]].into(),
                ..Default::default()
            });

            kind = -kind;
        }

        let a = 2.0 * PI / kind.f64();

        for k in 1..kind {
            let x = round6(cos!(k.f64() * a));
            let y = round6(sin!(k.f64() * a));

            new_transforms.push(Transform {
                density: 1.0,
                color_speed: 0.0,
                animate: 0.0,
                variations: vec![Var::Linear(Default::default())],
                color: if kind < 3 {
                    0.0
                } else {
                    (k - 1).f64() / (kind - 2).f64()
                },
                coefficients: [[x, y], [round6(-y), x], [0.0, 0.0]].into(),
                ..Default::default()
            });
        }

        new_transforms.sort_by(compare_xforms);

        self.transforms.extend(new_transforms);
    }
}

impl Default for Genome {
    fn default() -> Self {
        Self {
            name: None,
            palette: vec![Default::default(); 256],
            time: Default::default(),
            hsv_rgb_palette_blend: Default::default(),
            interpolation: Default::default(),
            palette_interpolation: Default::default(),
            interpolation_type: Default::default(),
            palette_index: Default::default(),
            size: Dimension {
                width: 100,
                height: 100,
            },
            center: Default::default(),
            rot_center: Default::default(),
            pixels_per_unit: 50.0,
            rotate: Default::default(),
            zoom: Default::default(),
            spatial_supersample: 1,
            spatial_filter_radius: 0.5,
            spatial_filter: Default::default(),
            temporal_filter: Default::default(),
            temporal_filter_width: 1.0,
            temporal_filter_exp: Default::default(),
            palette_mode: Default::default(),
            sample_density: 1.0,
            passes: Default::default(),
            num_temporal_samples: 1000,
            background: Default::default(),
            brightness: 4.0,
            gamma: 4.0,
            highlight_power: -1.0,
            vibrancy: 1.0,
            contrast: 1.0,
            hue_rotation: Default::default(),
            estimator_radius: 9.0,
            estimator_minimum: Default::default(),
            estimator_curve: 0.4,
            gamma_threshold: 0.01,
            transforms: Vec::new(),
            final_transform: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Dimension, Genome};

    #[test]
    fn scaling() {
        let mut genome = Genome {
            size: Dimension {
                width: 320,
                height: 240,
            },
            pixels_per_unit: 20.0,
            ..Default::default()
        };

        genome.scale(1.5);

        assert_eq!(genome.size.width, 480);
        assert_eq!(genome.size.height, 360);
        assert_eq!(genome.pixels_per_unit, 30.0);
    }
}
