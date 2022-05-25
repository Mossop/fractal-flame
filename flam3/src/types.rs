use std::ops::{Index, IndexMut};

use crate::utils::{color_from_str, color_to_str, try_map};

/// Holds a pixel colour in red, green, blue and alpha components.
#[derive(Debug, Clone, PartialEq)]
pub struct Rgba {
    // Red value, ranges 0-1.
    pub red: f64,
    // Green value, ranges 0-1.
    pub green: f64,
    // Blue value, ranges 0-1.
    pub blue: f64,
    // Alpha value, ranges 0-1.
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
    pub(crate) fn from_str_list(list: &str) -> Result<Self, String> {
        try_map(list.split(' '), color_from_str).and_then(|l| Rgba::try_from(&l[..]))
    }

    pub fn has_opacity(&self) -> bool {
        self.alpha < 1.0
    }

    pub(crate) fn to_str_list(&self) -> String {
        if self.has_opacity() {
            format!(
                "{} {} {} {}",
                color_to_str(self.red),
                color_to_str(self.green),
                color_to_str(self.blue),
                color_to_str(self.alpha),
            )
        } else {
            format!(
                "{} {} {}",
                color_to_str(self.red),
                color_to_str(self.green),
                color_to_str(self.blue),
            )
        }
    }
}

impl From<&Hsva> for Rgba {
    fn from(hsv: &Hsva) -> Self {
        let mut hue = hsv.hue;
        let saturation = hsv.saturation;
        let value = hsv.value;

        while hue >= 6.0 {
            hue -= 6.0;
        }
        while hue < 0.0 {
            hue += 6.0;
        }
        let j = hue.floor();
        let f = hue - j;
        let p = value * (1.0 - saturation);
        let q = value * (1.0 - (saturation * f));
        let t = value * (1.0 - (saturation * (1.0 - f)));

        let [red, green, blue] = match j as u32 {
            0 => [value, t, p],
            1 => [q, value, p],
            2 => [p, value, t],
            3 => [p, q, value],
            4 => [t, p, value],
            5 => [value, p, q],
            _ => [value, t, p],
        };

        Self {
            red,
            green,
            blue,
            alpha: hsv.alpha,
        }
    }
}

impl From<Rgba> for [f64; 3] {
    fn from(rgb: Rgba) -> [f64; 3] {
        [rgb.red, rgb.green, rgb.blue]
    }
}

impl From<&[f64; 4]> for Rgba {
    fn from(v: &[f64; 4]) -> Self {
        TryFrom::<&[f64]>::try_from(v).unwrap()
    }
}

impl From<&[f64; 3]> for Rgba {
    fn from(v: &[f64; 3]) -> Self {
        TryFrom::<&[f64]>::try_from(v).unwrap()
    }
}

impl TryFrom<&[f64]> for Rgba {
    type Error = String;

    fn try_from(values: &[f64]) -> Result<Self, String> {
        if values.len() < 3 || values.len() > 4 {
            return Err(format!(
                "Unexpected number of color values ({})",
                values.len()
            ));
        }

        Ok(Self {
            red: values[0],
            green: values[1],
            blue: values[2],
            alpha: if values.len() == 3 { 1.0 } else { values[3] },
        })
    }
}

impl From<&[u8; 4]> for Rgba {
    fn from(v: &[u8; 4]) -> Self {
        TryFrom::<&[u8]>::try_from(v).unwrap()
    }
}

impl From<&[u8; 3]> for Rgba {
    fn from(v: &[u8; 3]) -> Self {
        TryFrom::<&[u8]>::try_from(v).unwrap()
    }
}

impl TryFrom<&[u8]> for Rgba {
    type Error = String;

    fn try_from(values: &[u8]) -> Result<Self, String> {
        if values.len() < 3 || values.len() > 4 {
            return Err(format!(
                "Unexpected number of color values ({})",
                values.len()
            ));
        }

        Ok(Self {
            red: values[0] as f64 / 255.0,
            green: values[1] as f64 / 255.0,
            blue: values[2] as f64 / 255.0,
            alpha: if values.len() == 3 {
                1.0
            } else {
                values[3] as f64 / 255.0
            },
        })
    }
}

/// Holds a pixel colour in hue, saturation, value and alpha components.
#[derive(Debug, Clone, PartialEq)]
pub struct Hsva {
    // Hue value, ranges 0-6.
    pub hue: f64,
    // Saturation value, ranges 0-1.
    pub saturation: f64,
    // Value, ranges 0-1.
    pub value: f64,
    // Alpha value, ranges 0-1.
    pub alpha: f64,
}

impl Default for Hsva {
    fn default() -> Self {
        Self {
            hue: 0.0,
            saturation: 0.0,
            value: 0.0,
            alpha: 1.0,
        }
    }
}

impl Index<usize> for Hsva {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.hue,
            1 => &self.saturation,
            2 => &self.value,
            3 => &self.alpha,
            _ => panic!("Out of index error"),
        }
    }
}

impl IndexMut<usize> for Hsva {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.hue,
            1 => &mut self.saturation,
            2 => &mut self.value,
            3 => &mut self.alpha,
            _ => panic!("Out of index error"),
        }
    }
}

impl From<&Rgba> for Hsva {
    fn from(rgb: &Rgba) -> Self {
        let red = rgb.red;
        let green = rgb.green;
        let blue = rgb.blue;

        /* compute maximum of rd,gd,bd */
        let max = red.max(blue).max(green);

        /* compute minimum of rd,gd,bd */
        let min = red.min(green).min(blue);

        let del = max - min;
        let value = max;
        let saturation = if max != 0.0 { del / max } else { 0.0 };

        let mut hue = 0.0;
        if saturation != 0.0 {
            let rc = (max - red) / del;
            let gc = (max - green) / del;
            let bc = (max - blue) / del;

            hue = if red == max {
                bc - gc
            } else if green == max {
                2.0 + rc - bc
            } else {
                4.0 + gc - rc
            };

            if hue < 0.0 {
                hue += 6.0
            };
        }

        Self {
            hue,
            saturation,
            value,
            alpha: rgb.alpha,
        }
    }
}

impl From<&[f64; 4]> for Hsva {
    fn from(v: &[f64; 4]) -> Self {
        TryFrom::<&[f64]>::try_from(v).unwrap()
    }
}

impl From<&[f64; 3]> for Hsva {
    fn from(v: &[f64; 3]) -> Self {
        TryFrom::<&[f64]>::try_from(v).unwrap()
    }
}

impl TryFrom<&[f64]> for Hsva {
    type Error = String;

    fn try_from(values: &[f64]) -> Result<Self, String> {
        if values.len() < 3 || values.len() > 4 {
            return Err(format!(
                "Unexpected number of color values ({})",
                values.len()
            ));
        }

        Ok(Self {
            hue: values[0],
            saturation: values[1],
            value: values[2],
            alpha: if values.len() == 3 { 1.0 } else { values[3] },
        })
    }
}

impl From<&[u8; 4]> for Hsva {
    fn from(v: &[u8; 4]) -> Self {
        TryFrom::<&[u8]>::try_from(v).unwrap()
    }
}

impl From<&[u8; 3]> for Hsva {
    fn from(v: &[u8; 3]) -> Self {
        TryFrom::<&[u8]>::try_from(v).unwrap()
    }
}

impl TryFrom<&[u8]> for Hsva {
    type Error = String;

    fn try_from(values: &[u8]) -> Result<Self, String> {
        if values.len() < 3 || values.len() > 4 {
            return Err(format!(
                "Unexpected number of color values ({})",
                values.len()
            ));
        }

        Ok(Self {
            hue: 6.0 * (values[0] as f64) / 255.0,
            saturation: values[1] as f64 / 255.0,
            value: values[2] as f64 / 255.0,
            alpha: if values.len() == 3 {
                1.0
            } else {
                values[3] as f64 / 255.0
            },
        })
    }
}
