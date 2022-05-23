pub mod file;
pub mod render;
mod utils;
pub mod variations;

use std::{
    cmp::Ordering,
    f64::consts::PI,
    fmt::Display,
    ops::{Index, IndexMut},
    str::FromStr,
};

use educe::Educe;
pub use file::flam3::{flam3_from_reader, flam3_to_writer};
pub use render::{render, RenderOptions};
use utils::PanicCast;
use uuid::Uuid;
use variations::Var;

#[derive(Debug, Clone, Copy)]
pub(crate) enum ColorType {
    Byte,
    Hex,
}

fn try_map<T, C, F, R>(items: C, mapper: F) -> Result<Vec<R>, String>
where
    C: IntoIterator<Item = T>,
    F: Fn(T) -> Result<R, String>,
{
    let mut result = Vec::new();

    for item in items {
        result.push(mapper(item)?);
    }

    Ok(result)
}

fn parse<T>(val: &str) -> Result<T, String>
where
    T: FromStr,
    T::Err: Display,
{
    T::from_str(val).map_err(|e| format!("Unable to parse: {}", e))
}

fn color_from_str(str: &str, color_type: ColorType) -> Result<f64, String> {
    match color_type {
        ColorType::Byte => {
            let byte = f64::from_str(str)
                .map_err(|e| format!("Could not convert value '{}' to color: {}", str, e))?;
            Ok(fastdiv!(byte, 255.0))
        }
        ColorType::Hex => {
            let byte = u8::from_str_radix(str.trim_start(), 16)
                .map_err(|e| format!("Could not convert value '{}' to color: {}", str, e))?;
            // This looks odd but it matches what the original flam3 code does when the
            // -freciprocal-math optimisation is enabled (as it is by default).
            Ok(fastdiv!(byte as f64, 255.0))
        }
    }
}

fn color_to_str(color: f64, color_type: ColorType) -> String {
    match color_type {
        ColorType::Byte => {
            format!("{:.0}", color * 255.0)
        }
        ColorType::Hex => {
            format!("{:02x}", (color * 255.0) as u8)
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Rgba {
    pub red: f64,
    pub green: f64,
    pub blue: f64,
    pub alpha: f64,
}

impl Default for Rgba {
    fn default() -> Self {
        Self {
            red: 0.0,
            green: 0.0,
            blue: 0.0,
            alpha: 1.0,
        }
    }
}

impl Index<usize> for Rgba {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.red,
            1 => &self.green,
            2 => &self.blue,
            3 => &self.alpha,
            _ => panic!("Out of index error"),
        }
    }
}

impl IndexMut<usize> for Rgba {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.red,
            1 => &mut self.green,
            2 => &mut self.blue,
            3 => &mut self.alpha,
            _ => panic!("Out of index error"),
        }
    }
}

impl Rgba {
    pub(crate) fn from_str_list(list: &str, color_type: ColorType) -> Result<Self, String> {
        let values = try_map(list.split(' '), |s| color_from_str(s, color_type))?;

        if values.len() < 3 || values.len() > 4 {
            return Err(format!("Unexpected number of color values: '{}'", list));
        }

        let alpha: f64 = *values.get(3).unwrap_or(&1.0);

        Ok(Rgba {
            red: values[0],
            green: values[1],
            blue: values[2],
            alpha,
        })
    }

    pub fn has_opacity(&self) -> bool {
        self.alpha < 1.0
    }

    pub(crate) fn to_str_list(&self, color_type: ColorType) -> String {
        if self.has_opacity() {
            format!(
                "{} {} {} {}",
                color_to_str(self.red, color_type),
                color_to_str(self.green, color_type),
                color_to_str(self.blue, color_type),
                color_to_str(self.alpha, color_type),
            )
        } else {
            format!(
                "{} {} {}",
                color_to_str(self.red, color_type),
                color_to_str(self.green, color_type),
                color_to_str(self.blue, color_type),
            )
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Dimension {
    pub width: u32,
    pub height: u32,
}

impl Dimension {
    pub(crate) fn from_str_list(list: &str) -> Result<Self, String> {
        let values = try_map(list.split(' '), parse)?;

        if values.len() != 2 {
            return Err(format!("Unexpected number of dimensions: '{}'", list));
        }

        Ok(Dimension {
            width: values[0],
            height: values[1],
        })
    }

    pub(crate) fn to_str_list(&self) -> String {
        format!("{} {}", self.width, self.height)
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Coordinate {
    pub x: f64,
    pub y: f64,
}

impl Coordinate {
    pub(crate) fn from_str_list(list: &str) -> Result<Self, String> {
        let values = try_map(list.split(' '), parse)?;

        if values.len() != 2 {
            return Err(format!("Unexpected number of dimensions: '{}'", list));
        }

        Ok(Coordinate {
            x: values[0],
            y: values[1],
        })
    }

    pub(crate) fn to_str_list(&self) -> String {
        format!("{} {}", self.x, self.y)
    }
}

impl From<&[f64; 2]> for Coordinate {
    fn from(tuple: &[f64; 2]) -> Self {
        Self {
            x: tuple[0],
            y: tuple[1],
        }
    }
}

impl From<&[f64; 4]> for Coordinate {
    fn from(tuple: &[f64; 4]) -> Self {
        Self {
            x: tuple[0],
            y: tuple[1],
        }
    }
}

impl From<Coordinate> for (f64, f64) {
    fn from(c: Coordinate) -> Self {
        (c.x, c.y)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, strum_macros::EnumString, strum_macros::Display, Educe)]
#[educe(Default)]
pub enum Interpolation {
    #[educe(Default)]
    #[strum(serialize = "linear")]
    Linear,
    #[strum(serialize = "smooth")]
    Smooth,
}

#[derive(Debug, Clone, Copy, PartialEq, strum_macros::EnumString, strum_macros::Display, Educe)]
#[educe(Default)]
pub enum PaletteInterpolation {
    #[strum(serialize = "hsv")]
    Hsv,
    #[strum(serialize = "sweep")]
    Sweep,
    #[educe(Default)]
    #[strum(serialize = "hsv_circular")]
    HsvCircular,
    #[strum(serialize = "rgb")]
    Rgb,
}

#[derive(Debug, Clone, Copy, PartialEq, strum_macros::EnumString, strum_macros::Display, Educe)]
#[educe(Default)]
pub enum InterpolationType {
    #[strum(serialize = "linear")]
    Linear,
    #[educe(Default)]
    #[strum(serialize = "log")]
    Log,
    #[strum(serialize = "old")]
    Old,
    #[strum(serialize = "older")]
    Older,
}

#[derive(Debug, Clone, Copy, PartialEq, strum_macros::EnumString, strum_macros::Display, Educe)]
#[educe(Default)]
pub enum SpatialFilter {
    #[educe(Default)]
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

#[derive(Debug, Clone, Copy, PartialEq, strum_macros::EnumString, strum_macros::Display, Educe)]
#[educe(Default)]
pub enum TemporalFilter {
    #[educe(Default)]
    #[strum(serialize = "box")]
    Box,
    #[strum(serialize = "gaussian")]
    Gaussian,
    #[strum(serialize = "exp")]
    Exp,
}

#[derive(Debug, Clone, Copy, PartialEq, strum_macros::EnumString, strum_macros::Display, Educe)]
#[educe(Default)]
pub enum PaletteMode {
    #[educe(Default)]
    #[strum(serialize = "step")]
    Step,
    #[strum(serialize = "linear")]
    Linear,
}

#[derive(Debug, Clone, Copy, PartialEq, strum_macros::EnumString, strum_macros::Display, Educe)]
#[educe(Default)]
pub enum MotionFunction {
    #[educe(Default)]
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
    pub fn transform(&self, p: Coordinate) -> Coordinate {
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

    pub(crate) fn from_str_list(list: &str) -> Result<Self, String> {
        let values: Vec<f64> = try_map(list.split(' '), parse)?;
        if values.len() != 6 {
            return Err("Unexpected number of coefficients".to_string());
        }

        let mut affine = Self::default();

        for i in 0..3 {
            for k in 0..2 {
                affine.coefficients[i][k] = values[i * 2 + k];
            }
        }

        Ok(affine)
    }

    pub(crate) fn to_str_list(&self) -> String {
        format!(
            "{} {} {} {} {} {}",
            self.coefficients[0][0],
            self.coefficients[0][1],
            self.coefficients[1][0],
            self.coefficients[1][1],
            self.coefficients[2][0],
            self.coefficients[2][1],
        )
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

pub type Palette = Vec<Rgba>;

#[derive(Debug, Clone, PartialEq)]
pub struct Genome {
    pub name: Option<String>,
    pub time: f64,
    pub palette: Palette,
    pub hsv_rgb_palette_blend: f64,
    pub interpolation: Interpolation,
    pub palette_interpolation: PaletteInterpolation,
    pub interpolation_type: InterpolationType,
    pub palette_index: u32,

    /// The size of the generated image in pixels.
    pub size: Dimension,
    /// The center of the image, coordinates must be in the range [-1,1]
    pub center: Coordinate,
    pub rot_center: Coordinate,
    pub pixels_per_unit: f64,
    pub rotate: f64,
    pub zoom: f64,
    pub spatial_oversample: u32,
    pub spatial_filter_radius: f64,
    pub spatial_filter_select: SpatialFilter,
    pub temporal_filter: TemporalFilter,
    pub temporal_filter_width: f64,
    pub temporal_filter_exp: f64,
    pub palette_mode: PaletteMode,
    pub sample_density: f64,
    pub passes: u32,
    pub num_temporal_samples: u32,
    pub background: Rgba,
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
        ad = a.coefficients[0][0].atan2(a.coefficients[0][1]);
        bd = b.coefficients[0][0].atan2(b.coefficients[0][1]);
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
            let x = round6((k.f64() * a).cos());
            let y = round6((k.f64() * a).sin());

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
            spatial_oversample: 1,
            spatial_filter_radius: 0.5,
            spatial_filter_select: Default::default(),
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
    use crate::Coordinate;

    #[test]
    fn coordinate() {
        let c = Coordinate::from_str_list("0.353453 -0.235345").unwrap();
        assert_eq!(c.x, 0.353453);
        assert_eq!(c.y, -0.235345);
    }
}
