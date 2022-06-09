use std::{
    collections::HashMap,
    io::{Read, Write},
    str::FromStr,
};

use palette::{Pixel, Srgb, Srgba, WithAlpha};
use xml::{
    attribute::OwnedAttribute,
    reader::{self, EventReader, ParserConfig},
    writer::{self, EmitterConfig, EventWriter},
};

use crate::{
    parse,
    utils::{color_from_str, color_to_str, try_map},
    variations::{self, Var},
    Affine, Coordinate, Dimension, Genome, Interpolation, MotionFunction, Palette, TemporalFilter,
    Transform,
};

const VARIATION_COUNT: usize = 99;

trait XmlAttribute: Sized {
    fn to_attribute(&self) -> String;

    fn from_attribute(attr: &str) -> Result<Self, String>;
}

impl XmlAttribute for Affine {
    fn from_attribute(list: &str) -> Result<Self, String> {
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

    fn to_attribute(&self) -> String {
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

impl XmlAttribute for Dimension {
    fn from_attribute(list: &str) -> Result<Self, String> {
        let values = try_map(list.split(' '), parse)?;

        if values.len() != 2 {
            return Err(format!("Unexpected number of dimensions: '{}'", list));
        }

        Ok(Dimension {
            width: values[0],
            height: values[1],
        })
    }

    fn to_attribute(&self) -> String {
        format!("{} {}", self.width, self.height)
    }
}

impl XmlAttribute for Coordinate {
    fn from_attribute(list: &str) -> Result<Self, String> {
        let values = try_map(list.split(' '), parse)?;

        if values.len() != 2 {
            return Err(format!("Unexpected number of dimensions: '{}'", list));
        }

        Ok(Coordinate {
            x: values[0],
            y: values[1],
        })
    }

    fn to_attribute(&self) -> String {
        format!("{} {}", self.x, self.y)
    }
}

impl XmlAttribute for Srgba<f64> {
    fn from_attribute(list: &str) -> Result<Self, String> {
        let components = try_map(list.split(' '), color_from_str)?;

        match components.len() {
            3 => Ok(Srgb::<f64>::from_raw(&components[..]).with_alpha(1.0)),
            4 => Ok(*Srgba::<f64>::from_raw(&components[..])),
            _ => Err("Unexpected number of color components.".to_string()),
        }
    }

    fn to_attribute(&self) -> String {
        if self.alpha < 1.0 {
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

macro_rules! read_xml_event {
    ($parser:ident) => {
        $parser
            .next()
            .map_err(|e| format!("Xml parsing error: {}", e))?
    };
}

macro_rules! write_xml_event {
    ($writer:ident, $event:expr) => {
        $writer
            .write($event)
            .map_err(|e| format!("Xml writing error: {}", e))?
    };
}

macro_rules! setp {
    ($attrs:expr, $field:expr, $name:literal, $conversion:expr) => {
        if let Some(val) = $attrs.remove($name) {
            $field = $conversion(&val).map_err(|e| {
                format!("Failed to convert value '{}' for \"{}\": {}", val, $name, e)
            })?;
        }
    };
    ($attrs:expr, $field:expr, $name:literal) => {
        if let Some(val) = $attrs.remove($name) {
            $field = FromStr::from_str(&val).map_err(|e| {
                format!("Failed to convert value '{}' for \"{}\": {}", val, $name, e)
            })?;
        }
    };
}

macro_rules! writep {
    ($attrs:expr, $field:expr, $name:literal, $map:expr) => {
        $attrs.push(($name.to_string(), $map($field).to_string()));
    };
    ($attrs:expr, $field:expr, $name:literal) => {
        $attrs.push(($name.to_string(), $field.to_string()));
    };
}

macro_rules! variations {
    {$(($index:expr, $var:ident, $name:expr, {
        $(
            $field:ident: $attr:expr
        ),*$(,)?
    })),*$(,)?} => {
        fn var_name_from_index(index: usize) -> &'static str {
            match index {
                $(
                    $index => $name,
                )*
                _ => panic!("Unknown variation index.")
            }
        }

        // fn var_index_from_name(name: &str) -> usize {
        //     match name {
        //         $(
        //             $name => $index,
        //         )*
        //         _ => panic!("Unknown variation name.")
        //     }
        // }

        // fn var_index(var: Var) -> usize {
        //     $(
        //         if let Var::$var(_) = var {
        //             return $index;
        //         }
        //     )*
        //     unreachable!("Should have tested every variation");
        // }

        fn write_var(var: &Var, attrs: &mut Vec<(String, String)>) {
            match var {
                $(
                    Var::$var(v) => {
                        attrs.push(($name.to_string(), v.weight.to_string()));
                        $(
                            attrs.push(($attr.to_string(), v.$field.to_string()));
                        )*
                    }
                )*
            }
        }

        fn var_from_index(index: usize, weight: f64, attrs: Option<&mut HashMap<String, String>>) -> Result<Var, String> {
            $(
                if index == $index {
                    let mut variation = variations::$var::default();
                    variation.weight = weight;

                    variation = if let Some(_attrs) = attrs {
                        variations::$var {
                            weight,
                            $(
                                $field: if let Some(val) = _attrs.remove($attr) {
                                    parse(&val)?
                                } else {
                                    variation.$field
                                },
                            )*
                        }
                    } else {
                        variation.weight = weight;
                        variation
                    };

                    return Ok(Var::$var(variation));
                }
            )*
            unreachable!("Unknown variation index.");
        }
    }
}

variations! {
    (0, Linear, "linear", {}),
    (1, Sinusoidal, "sinusoidal", {}),
    (2, Spherical, "spherical", {}),
    (3, Swirl, "swirl", {}),
    (4, Horseshoe, "horseshoe", {}),
    (5, Polar, "polar", {}),
    (6, Handkerchief, "handkerchief", {}),
    (7, Heart, "heart", {}),
    (8, Disc, "disc", {}),
    (9, Spiral, "spiral", {}),
    (10, Hyperbolic, "hyperbolic", {}),
    (11, Diamond, "diamond", {}),
    (12, Ex, "ex", {}),
    (13, Julia, "julia", {}),
    (14, Bent, "bent", {}),
    (15, Waves, "waves", {}),
    (16, Fisheye, "fisheye", {}),
    (17, Popcorn, "popcorn", {}),
    (18, Exponential, "exponential", {}),
    (19, Power, "power", {}),
    (20, Cosine, "cosine", {}),
    (21, Rings, "rings", {}),
    (22, Fan, "fan", {}),
    (23, Blob, "blob", {
        low: "blob_low",
        high: "blob_high",
        waves: "blob_waves",
    }),
    (24, Pdj, "pdj", {
        a: "pdj_a",
        b: "pdj_b",
        c: "pdj_c",
        d: "pdj_d",
    }),
    (25, Fan2, "fan2", {
        x: "fan2_x",
        y: "fan2_y",
    }),
    (26, Rings2, "rings2", {
        val: "rings2_val",
    }),
    (27, Eyefish, "eyefish", {}),
    (28, Bubble, "bubble", {}),
    (29, Cylinder, "cylinder", {}),
    (30, Perspective, "perspective", {
        angle: "perspective_angle",
        distance: "perspective_dist",
    }),
    (31, Noise, "noise", {}),
    (32, Julian, "julian", {
        power: "julian_power",
        distance: "julian_dist",
    }),
    (33, Juliascope, "juliascope", {
        power: "juliascope_power",
        distance: "juliascope_dist",
    }),
    (34, Blur, "blur", {}),
    (35, GaussianBlur, "gaussian_blur", {}),
    (36, RadialBlur, "radial_blur", {
        angle: "radial_blur_angle",
    }),
    (37, Pie, "pie", {
        slices: "pie_slices",
        rotation: "pie_rotation",
        thickness: "pie_thickness",
    }),
    (38, Ngon, "ngon", {
        sides: "ngon_sides",
        power: "ngon_power",
        circle: "ngon_circle",
        corners: "ngon_corners",
    }),
    (39, Curl, "curl", {
        c1: "curl_c1",
        c2: "curl_c2",
    }),
    (40, Rectangles, "rectangles", {
        x: "rectangles_x",
        y: "rectangles_y",
    }),
    (41, Arch, "arch", {}),
    (42, Tangent, "tangent", {}),
    (43, Square, "square", {}),
    (44, Rays, "rays", {}),
    (45, Blade, "blade", {}),
    (46, Secant2, "secant2", {}),
    (47, Twintrian, "twintrian", {}),
    (48, Cross, "cross", {}),
    (49, Disc2, "disc2", {
        rotate: "disc2_rot",
        twist: "disc2_twist",
    }),
    (50, SuperShape, "super_shape", {
        rnd: "super_shape_rnd",
        m: "super_shape_m",
        n1: "super_shape_n1",
        n2: "super_shape_n2",
        n3: "super_shape_n3",
        holes: "super_shape_holes",
    }),
    (51, Flower, "flower", {
        petals: "flower_petals",
        holes: "flower_holes",
    }),
    (52, Conic, "conic", {
        eccentricity: "conic_eccentricity",
        holes: "conic_holes",
    }),
    (53, Parabola, "parabola", {
        height: "parabola_height",
        width: "parabola_width",
    }),
    (54, Bent2, "bent2", {
        x: "bent2_x",
        y: "bent2_y",
    }),
    (55, Bipolar, "bipolar", {
        shift: "bipolar_shift",
    }),
    (56, Boarders, "boarders", {}),
    (57, Butterfly, "butterfly", {}),
    (58, Cell, "cell", {
        size: "cell_size",
    }),
    (59, Cpow, "cpow", {
        r: "cpow_r",
        i: "cpow_i",
        power: "cpow_power",
    }),
    (60, Curve, "curve", {
        x_amp: "curve_xamp",
        y_amp: "curve_yamp",
        x_length: "curve_xlength",
        y_length: "curve_ylength",
    }),
    (61, Edisc, "edisc", {}),
    (62, Elliptic, "elliptic", {}),
    (63, Escher, "escher", {
        beta: "escher_beta",
    }),
    (64, Foci, "foci", {}),
    (65, Lazysusan, "lazysusan", {
        spin: "lazysusan_spin",
        space: "lazysusan_space",
        twist: "lazysusan_twist",
        x: "lazysusan_x",
        y: "lazysusan_y",
    }),
    (66, Loonie, "loonie", {}),
    (67, PreBlur, "pre_blur", {}),
    (68, Modulus, "modulus", {
        x: "modulus_x",
        y: "modulus_y",
    }),
    (69, Oscilloscope, "oscilloscope", {
        separation: "oscilloscope_separation",
        frequency: "oscilloscope_frequency",
        amplitude: "oscilloscope_amplitude",
        damping: "oscilloscope_damping",
    }),
    (70, Polar2, "polar2", {}),
    (71, Popcorn2, "popcorn2", {
        x: "popcorn2_x",
        y: "popcorn2_y",
        c: "popcorn2_c",
    }),
    (72, Scry, "scry", {}),
    (73, Separation, "separation", {
        x: "separation_x",
        y: "separation_y",
        x_inside: "separation_xinside",
        y_inside: "separation_yinside",
    }),
    (74, Split, "split", {
        x_size: "split_xsize",
        y_size: "split_ysize",
    }),
    (75, Splits, "splits", {
        x: "splits_x",
        y: "splits_y",
    }),
    (76, Stripes, "stripes", {
        space: "stripes_space",
        warp: "stripes_warp",
    }),
    (77, Wedge, "wedge", {
        angle: "wedge_angle",
        hole: "wedge_hole",
        count: "wedge_count",
        swirl: "wedge_swirl",
    }),
    (78, WedgeJulia, "wedge_julia", {
        angle: "wedge_julia_angle",
        count: "wedge_julia_count",
        power: "wedge_julia_power",
        dist: "wedge_julia_dist",
    }),
    (79, WedgeSph, "wedge_sph", {
        angle: "wedge_sph_angle",
        count: "wedge_sph_count",
        hole: "wedge_sph_hole",
        swirl: "wedge_sph_swirl",
    }),
    (80, Whorl, "whorl", {
        inside: "whorl_inside",
        outside: "whorl_outside",
    }),
    (81, Waves2, "waves2", {
        x_frequency: "waves2_freqx",
        y_frequency: "waves2_freqy",
        x_scale: "waves2_scalex",
        y_scale: "waves2_scaley",
    }),
    (82, Exp, "exp", {}),
    (83, Log, "log", {}),
    (84, Sin, "sin", {}),
    (85, Cos, "cos", {}),
    (86, Tan, "tan", {}),
    (87, Sec, "sec", {}),
    (88, Csc, "csc", {}),
    (89, Cot, "cot", {}),
    (90, Sinh, "sinh", {}),
    (91, Cosh, "cosh", {}),
    (92, Tanh, "tanh", {}),
    (93, Sech, "sech", {}),
    (94, Csch, "csch", {}),
    (95, Coth, "coth", {}),
    (96, Auger, "auger", {
        symmetry: "auger_sym",
        strength: "auger_weight",
        frequency: "auger_freq",
        scale: "auger_scale",
    }),
    (97, Flux, "flux", {
        spread: "flux_spread",
    }),
    (98, Mobius, "mobius", {
        re_a: "mobius_re_a",
        re_b: "mobius_re_b",
        re_c: "mobius_re_c",
        re_d: "mobius_re_d",
        im_a: "mobius_im_a",
        im_b: "mobius_im_b",
        im_c: "mobius_im_c",
        im_d: "mobius_im_d",
    }),
}

fn write_start_element<W: Write>(
    writer: &mut EventWriter<W>,
    name: &str,
    attrs: Vec<(String, String)>,
) -> Result<(), String> {
    let mut se = writer::XmlEvent::start_element(name);

    for (name, value) in attrs.iter() {
        se = se.attr(name.as_str(), value);
    }

    write_xml_event!(writer, se);
    Ok(())
}

fn finish_element<R: Read>(parser: &mut EventReader<R>) -> Result<(), String> {
    loop {
        if let reader::XmlEvent::EndElement { .. } = read_xml_event!(parser) {
            return Ok(());
        }
    }
}

fn attr_hash(attributes: Vec<OwnedAttribute>) -> HashMap<String, String> {
    HashMap::from_iter(attributes.into_iter().map(|a| (a.name.local_name, a.value)))
}

fn skip_whitespace(chars: &[u8], mut pos: usize) -> usize {
    while pos < chars.len() && chars[pos].is_ascii_whitespace() {
        pos += 1;
    }
    pos
}

fn read_hex(chars: &[u8], pos: usize) -> Result<(f64, usize), String> {
    fn as_digit(ch: u8) -> Result<u8, String> {
        if ch.is_ascii_digit() {
            Ok(ch - b'0')
        } else if ch.is_ascii_hexdigit() {
            if ch.is_ascii_lowercase() {
                Ok(ch - b'a' + 10)
            } else {
                Ok(ch - b'A' + 10)
            }
        } else {
            Err(format!("Unexpected character '{}' in palette data", ch))
        }
    }

    let start = skip_whitespace(chars, pos);
    if start >= chars.len() {
        Err("Ran out of characters parsing palette data".to_string())
    } else if start == chars.len() - 1 || chars[start + 1].is_ascii_whitespace() {
        let val = as_digit(chars[start])?;
        Ok((val as f64 / 255.0, start + 1))
    } else {
        let val1 = as_digit(chars[start])?;
        let val2 = as_digit(chars[start + 1])?;
        Ok(((val1 * 16 + val2) as f64 / 255.0, start + 2))
    }
}

fn parse_palette_data(str: &str, bytes: usize, palette: &mut Palette) -> Result<(), String> {
    let chars: Vec<u8> = str.bytes().collect();
    let mut pos = skip_whitespace(&chars, 0);

    for entry in palette.iter_mut() {
        let (r, next) = read_hex(&chars, pos)?;
        let (g, next) = read_hex(&chars, next)?;
        let (b, next) = read_hex(&chars, next)?;

        let a = if bytes == 4 {
            let (a, _) = read_hex(&chars, next)?;
            a
        } else {
            1.0
        };

        *entry = Srgba::<f64>::from_components((r, g, b, a));
        pos = skip_whitespace(&chars, pos + bytes * 2);
    }

    Ok(())
}

fn parse_palette<R: Read>(
    parser: &mut EventReader<R>,
    mut attrs: HashMap<String, String>,
    genome: &mut Genome,
) -> Result<(), String> {
    let colors = attrs
        .remove("count")
        .ok_or_else(|| "Missing number of colors in palette".to_string())
        .and_then(|s| parse::<usize>(&s))?;

    if colors != genome.palette.len() {
        return Err(format!(
            "Unexpected number of colors in palette: {}",
            colors
        ));
    }

    let bytes = match attrs.remove("format").as_deref() {
        Some("RGB") => 3,
        Some("RGBA") => 4,
        _ => return Err("Incorrect palette format".to_string()),
    };

    loop {
        match read_xml_event!(parser) {
            reader::XmlEvent::EndElement { .. } => {
                return Ok(());
            }
            reader::XmlEvent::Characters(palette) => {
                parse_palette_data(&palette, bytes, &mut genome.palette)?;
            }
            _ => {}
        }
    }
}

fn serialize_colors<W: Write>(genome: &Genome, writer: &mut EventWriter<W>) -> Result<(), String> {
    for (index, color) in genome.palette.iter().enumerate() {
        let mut attrs: Vec<(String, String)> = vec![("index".to_string(), index.to_string())];

        if color.alpha < 1.0 {
            attrs.push(("rgba".to_string(), color.to_attribute()));
        } else {
            attrs.push(("rgb".to_string(), color.to_attribute()));
        }

        write_start_element(writer, "color", attrs)?;
        write_xml_event!(writer, writer::XmlEvent::end_element());
    }

    Ok(())
}

fn parse_color<R: Read>(
    parser: &mut EventReader<R>,
    attrs: HashMap<String, String>,
    genome: &mut Genome,
) -> Result<(), String> {
    let index = attrs
        .get("index")
        .ok_or_else(|| "Missing color index".to_string())
        .and_then(|s| parse::<usize>(s))?;

    if index > 255 {
        return Err(format!("Color index {} is out of range", index));
    }

    let color = attrs
        .get("rgba")
        .or_else(|| attrs.get("rgb"))
        .ok_or_else(|| "Missing color value".to_string())
        .and_then(|s| Srgba::from_attribute(s))?;

    genome.palette[index] = color;

    finish_element(parser)
}

fn parse_symmetry<R: Read>(
    parser: &mut EventReader<R>,
    attrs: HashMap<String, String>,
    genome: &mut Genome,
) -> Result<(), String> {
    let kind = attrs
        .get("kind")
        .ok_or_else(|| "Missing symmetry kind".to_string())
        .and_then(|s| parse::<i32>(s))?;

    if matches!(kind, 0 | 1) {
        return Err(format!("Unexpected symmetry kind {}", kind));
    }

    genome.add_symmetry(kind);

    finish_element(parser)
}

fn serialize_transform<W: Write>(
    transform: &Transform,
    writer: &mut EventWriter<W>,
    is_final: bool,
) -> Result<(), String> {
    let mut attrs: Vec<(String, String)> = Vec::new();

    if !is_final {
        writep!(attrs, transform.density, "weight");
    }
    writep!(attrs, transform.color, "color");
    if transform.color_speed != 0.0 {
        writep!(attrs, transform.color_speed, "color_speed");
    }
    if !is_final {
        writep!(attrs, transform.animate, "animate");
    }

    if let Some(func) = transform.motion_function {
        attrs.push(("motion_function".to_string(), func.to_string()));
    }

    for var in &transform.variations {
        write_var(var, &mut attrs);
    }

    writep!(attrs, &transform.coefficients, "coefs", |a: &Affine| a
        .to_attribute());

    writep!(attrs, &transform.post, "post", |a: &Affine| a
        .to_attribute());

    writep!(attrs, transform.opacity, "opacity");

    let tag_name = if is_final { "finalxform" } else { "xform" };
    write_start_element(writer, tag_name, attrs)?;

    write_xml_event!(writer, writer::XmlEvent::end_element());

    Ok(())
}

fn parse_transform<R: Read>(
    parser: &mut EventReader<R>,
    mut attrs: HashMap<String, String>,
) -> Result<Transform, String> {
    let mut transform = Transform::default();

    // Deprecated form
    if let Some(sym) = attrs.remove("symmetry") {
        let sym: f64 = parse(&sym)?;
        transform.color_speed = (1.0 - sym) / 2.0;
        transform.animate = if sym > 0.0 { 0.0 } else { 1.0 };
    }

    if let Some(color) = attrs.remove("color") {
        let nums: Vec<&str> = color.split(' ').collect();
        if matches!(nums.len(), 1 | 2) {
            transform.color = parse(nums[0])?;
        } else {
            return Err(format!("Malformed color attribute: '{}'", color));
        }
    }

    if let Some(coefs) = attrs.remove("coefs") {
        transform.coefficients = Affine::from_attribute(&coefs)?;
    }

    if let Some(post) = attrs.remove("post") {
        transform.post = Affine::from_attribute(&post)?;
    }

    setp!(attrs, transform.density, "weight");
    setp!(attrs, transform.color_speed, "color_speed");
    setp!(attrs, transform.animate, "animate");
    setp!(attrs, transform.opacity, "opacity");
    setp!(attrs, transform.motion_frequency, "motion_frequency");
    setp!(attrs, transform.motion_function, "motion_function", |s| {
        MotionFunction::from_str(s).map(Some)
    });

    for index in 0..VARIATION_COUNT {
        if let Some(st) = attrs.remove(var_name_from_index(index)) {
            let variation = var_from_index(index, parse(&st)?, Some(&mut attrs))?;
            transform.variations.push(variation);
        }
    }

    for (key, _) in attrs {
        log::warn!("Unknown transform attribute: {}", key);
    }

    finish_element(parser)?;
    Ok(transform)
}

fn serialize_genome<W: Write>(genome: &Genome, writer: &mut EventWriter<W>) -> Result<(), String> {
    let mut attrs: Vec<(String, String)> = Vec::new();

    if let Some(ref name) = genome.name {
        attrs.push(("name".to_string(), name.clone()));
    }

    writep!(attrs, genome.time, "time");
    writep!(attrs, &genome.size, "size", XmlAttribute::to_attribute);
    writep!(attrs, &genome.center, "center", XmlAttribute::to_attribute);
    writep!(attrs, genome.pixels_per_unit, "scale");

    if genome.zoom > 0.0 {
        writep!(attrs, genome.zoom, "zoom");
    }
    writep!(attrs, genome.rotate, "rotate");
    writep!(attrs, genome.spatial_supersample, "supersample");
    writep!(attrs, genome.spatial_filter_radius, "filter");
    writep!(attrs, genome.spatial_filter, "filter_shape");
    writep!(attrs, genome.temporal_filter, "temporal_filter_type");
    if genome.temporal_filter == TemporalFilter::Exp {
        writep!(attrs, genome.temporal_filter_exp, "temporal_filter_exp");
    }
    writep!(attrs, genome.temporal_filter_width, "temporal_filter_width");
    writep!(attrs, genome.sample_density, "quality");
    writep!(attrs, genome.passes, "passes");
    writep!(attrs, genome.num_temporal_samples, "temporal_samples");
    writep!(
        attrs,
        &genome.background,
        "background",
        XmlAttribute::to_attribute
    );
    writep!(attrs, genome.brightness, "brightness");
    writep!(attrs, genome.gamma, "gamma");
    writep!(attrs, genome.highlight_power, "highlight_power");
    writep!(attrs, genome.vibrancy, "vibrancy");
    writep!(attrs, genome.estimator_radius, "estimator_radius");
    writep!(attrs, genome.estimator_minimum, "estimator_minimum");
    writep!(attrs, genome.estimator_curve, "estimator_curve");
    writep!(attrs, genome.gamma_threshold, "gamma_threshold");
    writep!(attrs, genome.palette_mode, "palette_mode");

    if genome.interpolation != Interpolation::Linear {
        writep!(attrs, genome.interpolation, "interpolation");
    }
    writep!(attrs, genome.interpolation_type, "interpolation_type");

    writep!(attrs, genome.palette_interpolation, "palette_interpolation");
    if genome.hsv_rgb_palette_blend > 0.0 {
        writep!(attrs, genome.hsv_rgb_palette_blend, "hsv_rgb_palette_blend");
    }

    write_start_element(writer, "flame", attrs)?;

    for transform in genome.transforms.iter() {
        serialize_transform(transform, writer, false)?;
    }

    if let Some(ref transform) = genome.final_transform {
        serialize_transform(transform, writer, true)?;
    }

    serialize_colors(genome, writer)?;

    write_xml_event!(writer, writer::XmlEvent::end_element());

    Ok(())
}

fn parse_genome<R: Read>(
    parser: &mut EventReader<R>,
    mut attrs: HashMap<String, String>,
) -> Result<Genome, String> {
    let mut genome = Genome {
        name: attrs.remove("name"),
        ..Default::default()
    };

    attrs.remove("version");

    setp!(attrs, genome.time, "time");
    setp!(attrs, genome.hsv_rgb_palette_blend, "hsv_rgb_palette_blend");
    setp!(attrs, genome.interpolation, "interpolation");
    setp!(attrs, genome.palette_interpolation, "palette_interpolation");
    setp!(attrs, genome.interpolation_type, "interpolation_space");
    setp!(attrs, genome.interpolation_type, "interpolation_type");
    setp!(attrs, genome.palette_index, "palette");
    setp!(attrs, genome.size, "size", XmlAttribute::from_attribute);
    setp!(attrs, genome.center, "center", XmlAttribute::from_attribute);
    setp!(attrs, genome.pixels_per_unit, "scale");
    setp!(attrs, genome.rotate, "rotate");
    setp!(attrs, genome.zoom, "zoom");
    setp!(attrs, genome.spatial_supersample, "supersample");
    setp!(attrs, genome.spatial_supersample, "oversample");
    setp!(attrs, genome.spatial_filter_radius, "filter");
    setp!(attrs, genome.spatial_filter, "filter_shape");
    setp!(attrs, genome.temporal_filter, "temporal_filter_type");
    setp!(attrs, genome.temporal_filter_width, "temporal_filter_width");
    setp!(attrs, genome.temporal_filter_exp, "temporal_filter_exp");
    setp!(attrs, genome.palette_mode, "palette_mode");
    setp!(attrs, genome.sample_density, "quality");
    setp!(attrs, genome.passes, "passes");
    setp!(attrs, genome.num_temporal_samples, "temporal_samples");
    setp!(
        attrs,
        genome.background,
        "background",
        XmlAttribute::from_attribute
    );
    setp!(attrs, genome.brightness, "brightness");
    setp!(attrs, genome.gamma, "gamma");
    setp!(attrs, genome.highlight_power, "highlight_power");
    setp!(attrs, genome.vibrancy, "vibrancy");
    setp!(attrs, genome.hue_rotation, "hue", |s| f64::from_str(s)
        .map(|v| v % 1.0));
    setp!(attrs, genome.estimator_radius, "estimator_radius");
    setp!(attrs, genome.estimator_minimum, "estimator_minimum");
    setp!(attrs, genome.estimator_curve, "estimator_curve");
    setp!(attrs, genome.gamma_threshold, "gamma_threshold");

    genome.rot_center = genome.center.clone();

    for (key, _) in attrs {
        log::warn!("Unknown genome attribute: {}", key);
    }

    loop {
        match read_xml_event!(parser) {
            reader::XmlEvent::StartElement {
                name,
                attributes: element_attrs,
                ..
            } => match name.local_name.as_str() {
                "color" => {
                    parse_color(parser, attr_hash(element_attrs), &mut genome)?;
                }
                "palette" => parse_palette(parser, attr_hash(element_attrs), &mut genome)?,
                "symmetry" => {
                    parse_symmetry(parser, attr_hash(element_attrs), &mut genome)?;
                }
                "xform" => {
                    let transform = parse_transform(parser, attr_hash(element_attrs))?;
                    genome.add_transform(transform);
                }
                "finalxform" => {
                    let transform = parse_transform(parser, attr_hash(element_attrs))?;
                    if transform.density != 0.0 {
                        return Err(
                            "Final transforms should not have weight specified.".to_string()
                        );
                    }
                    genome.final_transform = Some(transform);
                }
                "edit" => {}
                unknown => log::warn!("Unknown flame element {}.", unknown),
            },
            reader::XmlEvent::EndElement { .. } => {
                break;
            }
            _ => {}
        }
    }

    Ok(genome)
}

pub fn flam3_from_reader<R: Read>(reader: R) -> Result<Vec<Genome>, String> {
    let mut parser = EventReader::new_with_config(
        reader,
        ParserConfig {
            whitespace_to_characters: true,
            cdata_to_characters: true,
            ..Default::default()
        },
    );
    let mut genomes = Vec::new();

    let mut index = 0;
    loop {
        match read_xml_event!(parser) {
            reader::XmlEvent::StartElement {
                name,
                attributes: element_attrs,
                ..
            } => {
                if name.local_name == "flame" {
                    let genome = parse_genome(&mut parser, attr_hash(element_attrs))
                        .map_err(|e| format!("Failed to parse genome at index {}: {}", index, e))?;
                    genomes.push(genome);
                    index += 1;
                }
            }
            reader::XmlEvent::EndDocument => {
                break;
            }
            _ => {}
        }
    }

    Ok(genomes)
}

pub fn flam3_to_writer<W: Write>(genomes: &[Genome], sink: W) -> Result<(), String> {
    let mut writer = EventWriter::new_with_config(
        sink,
        EmitterConfig {
            perform_indent: true,
            write_document_declaration: false,
            pad_self_closing: false,
            ..Default::default()
        },
    );

    write_xml_event!(writer, writer::XmlEvent::start_element("flames"));

    for genome in genomes {
        serialize_genome(genome, &mut writer)?;
    }

    write_xml_event!(writer, writer::XmlEvent::end_element());

    Ok(())
}

#[cfg(test)]
mod tests {
    use palette::Srgba;

    use crate::file::flam3::XmlAttribute;

    use super::parse_palette_data;

    #[test]
    fn palette() {
        let palette_data = "
      baa0b6be85a6bd799ebc6e97c26897c96398ca6699cc699b
      ca85a9d094b4d6a4c0d2aec2cfb9c4c8c1c4c1cac4c3cdc7
      c6d1cbc6cecbbfc6c3b9bebcbbbdbebdbdc0c1bbc2c5bac5
      d3c5d1d9cdd8e0d6e0e0d4dfe0d2dedcd0dad9cfd7cecbcf
      c2c6c6a5a8a8969c9a87908c808b867a87807a867f7a867f
      8189897f8f8f7d969678939874909b708b996c8697607690
      54668346476d443963432c5a4522504819464b1b464e1d47
      6037566643616c506c695e75666c7f6470816374835e7483
      556f804f5c73564e695d40605d395b5e32575a2850582850
      5d385e5e3c60604062584a6050545f4b5b634662683e6971
      3a6c75386b6e3b656f3e60703f5a714055733f49753b4171
      32325c2f28522d1e492c1a462c17432c12442a f422c1243
      311a4c4b31665b38776c3f89724590784b97835b9e896ca3
      9882af9b86b09e8bb29c8ab09b89af9a80a89a7aa29a739b
      9c6a93934c899046878d40868b37858b36838b3d81873e7f
      80437f81498283508583548584588583628380657e7f657a
      7c5d736e466568405e623a5856334e4729403d1d3736132f
      38 62438 c253913273a192b3c20303f2a364d37445d3e51
      674c60756b7678737b7b7b818284888a8a8f909196949499
      999b9f9a9c9f9b9da09a9da2989ca1939aa08c98a18695a0
      81939f818d99818a97828896867e958774938a65898e557e
      96416e97416f9942719b4172963c6c8d31627a2e596c3353
      603f5253494e4a464439413b2d3e322840312c4b38395540
      596f5e6578697182758a9890a2aaa7b4b8b8bdbfc0c2c4c8
      c0c6ccb9c5ccacbec999b1bf88a3b67698ae6890a85e8ca5
      548aa35187a25188a45189a55187a54e839f4b74974b648f
      4c57884b4c85474d80444b7e4448804e4c835e528c6d6594
      7a7a9d888eae99a7bdafbacdc3cadaced8ddd1d9dccbd7d4
      c0cdc8b2b9b89da7a7848e916b767754626343464f393249
      35254b39204d4327544c3158513a635647715851805b5c8b
      5863915168914d6f914a71904d778e5276915a7b9465849c
      778ea38ba0ac9dacb5a9b3b9adb5baadb2b6aeafb6b3aab8
";
        let mut palette = vec![Srgba::default(); 256];
        parse_palette_data(palette_data, 3, &mut palette).unwrap();
        let mut i = palette.into_iter().map(|c| c.to_attribute());

        assert_eq!(i.next().unwrap(), "186 160 182");
        assert_eq!(i.next().unwrap(), "190 133 166");
        assert_eq!(i.next().unwrap(), "189 121 158");
        assert_eq!(i.next().unwrap(), "188 110 151");
        assert_eq!(i.next().unwrap(), "194 104 151");
        assert_eq!(i.next().unwrap(), "201 99 152");
        assert_eq!(i.next().unwrap(), "202 102 153");
        assert_eq!(i.next().unwrap(), "204 105 155");
        assert_eq!(i.next().unwrap(), "202 133 169");
        assert_eq!(i.next().unwrap(), "208 148 180");
        assert_eq!(i.next().unwrap(), "214 164 192");
        assert_eq!(i.next().unwrap(), "210 174 194");
        assert_eq!(i.next().unwrap(), "207 185 196");
        assert_eq!(i.next().unwrap(), "200 193 196");
        assert_eq!(i.next().unwrap(), "193 202 196");
        assert_eq!(i.next().unwrap(), "195 205 199");
        assert_eq!(i.next().unwrap(), "198 209 203");
        assert_eq!(i.next().unwrap(), "198 206 203");
        assert_eq!(i.next().unwrap(), "191 198 195");
        assert_eq!(i.next().unwrap(), "185 190 188");
        assert_eq!(i.next().unwrap(), "187 189 190");
        assert_eq!(i.next().unwrap(), "189 189 192");
        assert_eq!(i.next().unwrap(), "193 187 194");
        assert_eq!(i.next().unwrap(), "197 186 197");
        assert_eq!(i.next().unwrap(), "211 197 209");
        assert_eq!(i.next().unwrap(), "217 205 216");
        assert_eq!(i.next().unwrap(), "224 214 224");
        assert_eq!(i.next().unwrap(), "224 212 223");
        assert_eq!(i.next().unwrap(), "224 210 222");
        assert_eq!(i.next().unwrap(), "220 208 218");
        assert_eq!(i.next().unwrap(), "217 207 215");
        assert_eq!(i.next().unwrap(), "206 203 207");
        assert_eq!(i.next().unwrap(), "194 198 198");
        assert_eq!(i.next().unwrap(), "165 168 168");
        assert_eq!(i.next().unwrap(), "150 156 154");
        assert_eq!(i.next().unwrap(), "135 144 140");
        assert_eq!(i.next().unwrap(), "128 139 134");
        assert_eq!(i.next().unwrap(), "122 135 128");
        assert_eq!(i.next().unwrap(), "122 134 127");
        assert_eq!(i.next().unwrap(), "122 134 127");
        assert_eq!(i.next().unwrap(), "129 137 137");
        assert_eq!(i.next().unwrap(), "127 143 143");
        assert_eq!(i.next().unwrap(), "125 150 150");
        assert_eq!(i.next().unwrap(), "120 147 152");
        assert_eq!(i.next().unwrap(), "116 144 155");
        assert_eq!(i.next().unwrap(), "112 139 153");
        assert_eq!(i.next().unwrap(), "108 134 151");
        assert_eq!(i.next().unwrap(), "96 118 144");
        assert_eq!(i.next().unwrap(), "84 102 131");
        assert_eq!(i.next().unwrap(), "70 71 109");
        assert_eq!(i.next().unwrap(), "68 57 99");
        assert_eq!(i.next().unwrap(), "67 44 90");
        assert_eq!(i.next().unwrap(), "69 34 80");
        assert_eq!(i.next().unwrap(), "72 25 70");
        assert_eq!(i.next().unwrap(), "75 27 70");
        assert_eq!(i.next().unwrap(), "78 29 71");
        assert_eq!(i.next().unwrap(), "96 55 86");
        assert_eq!(i.next().unwrap(), "102 67 97");
        assert_eq!(i.next().unwrap(), "108 80 108");
        assert_eq!(i.next().unwrap(), "105 94 117");
        assert_eq!(i.next().unwrap(), "102 108 127");
        assert_eq!(i.next().unwrap(), "100 112 129");
        assert_eq!(i.next().unwrap(), "99 116 131");
        assert_eq!(i.next().unwrap(), "94 116 131");
        assert_eq!(i.next().unwrap(), "85 111 128");
        assert_eq!(i.next().unwrap(), "79 92 115");
        assert_eq!(i.next().unwrap(), "86 78 105");
        assert_eq!(i.next().unwrap(), "93 64 96");
        assert_eq!(i.next().unwrap(), "93 57 91");
        assert_eq!(i.next().unwrap(), "94 50 87");
        assert_eq!(i.next().unwrap(), "90 40 80");
        assert_eq!(i.next().unwrap(), "88 40 80");
        assert_eq!(i.next().unwrap(), "93 56 94");
        assert_eq!(i.next().unwrap(), "94 60 96");
        assert_eq!(i.next().unwrap(), "96 64 98");
        assert_eq!(i.next().unwrap(), "88 74 96");
        assert_eq!(i.next().unwrap(), "80 84 95");
        assert_eq!(i.next().unwrap(), "75 91 99");
        assert_eq!(i.next().unwrap(), "70 98 104");
        assert_eq!(i.next().unwrap(), "62 105 113");
        assert_eq!(i.next().unwrap(), "58 108 117");
        assert_eq!(i.next().unwrap(), "56 107 110");
        assert_eq!(i.next().unwrap(), "59 101 111");
        assert_eq!(i.next().unwrap(), "62 96 112");
        assert_eq!(i.next().unwrap(), "63 90 113");
        assert_eq!(i.next().unwrap(), "64 85 115");
        assert_eq!(i.next().unwrap(), "63 73 117");
        assert_eq!(i.next().unwrap(), "59 65 113");
        assert_eq!(i.next().unwrap(), "50 50 92");
        assert_eq!(i.next().unwrap(), "47 40 82");
        assert_eq!(i.next().unwrap(), "45 30 73");
        assert_eq!(i.next().unwrap(), "44 26 70");
        assert_eq!(i.next().unwrap(), "44 23 67");
        assert_eq!(i.next().unwrap(), "44 18 68");
        assert_eq!(i.next().unwrap(), "42 244 34");
        assert_eq!(i.next().unwrap(), "44 18 67");
        assert_eq!(i.next().unwrap(), "49 26 76");
        assert_eq!(i.next().unwrap(), "75 49 102");
        assert_eq!(i.next().unwrap(), "91 56 119");
        assert_eq!(i.next().unwrap(), "108 63 137");
        assert_eq!(i.next().unwrap(), "114 69 144");
        assert_eq!(i.next().unwrap(), "120 75 151");
        assert_eq!(i.next().unwrap(), "131 91 158");
        assert_eq!(i.next().unwrap(), "137 108 163");
        assert_eq!(i.next().unwrap(), "152 130 175");
        assert_eq!(i.next().unwrap(), "155 134 176");
        assert_eq!(i.next().unwrap(), "158 139 178");
        assert_eq!(i.next().unwrap(), "156 138 176");
        assert_eq!(i.next().unwrap(), "155 137 175");
        assert_eq!(i.next().unwrap(), "154 128 168");
        assert_eq!(i.next().unwrap(), "154 122 162");
        assert_eq!(i.next().unwrap(), "154 115 155");
        assert_eq!(i.next().unwrap(), "156 106 147");
        assert_eq!(i.next().unwrap(), "147 76 137");
        assert_eq!(i.next().unwrap(), "144 70 135");
        assert_eq!(i.next().unwrap(), "141 64 134");
        assert_eq!(i.next().unwrap(), "139 55 133");
        assert_eq!(i.next().unwrap(), "139 54 131");
        assert_eq!(i.next().unwrap(), "139 61 129");
        assert_eq!(i.next().unwrap(), "135 62 127");
        assert_eq!(i.next().unwrap(), "128 67 127");
        assert_eq!(i.next().unwrap(), "129 73 130");
        assert_eq!(i.next().unwrap(), "131 80 133");
        assert_eq!(i.next().unwrap(), "131 84 133");
        assert_eq!(i.next().unwrap(), "132 88 133");
        assert_eq!(i.next().unwrap(), "131 98 131");
        assert_eq!(i.next().unwrap(), "128 101 126");
        assert_eq!(i.next().unwrap(), "127 101 122");
        assert_eq!(i.next().unwrap(), "124 93 115");
        assert_eq!(i.next().unwrap(), "110 70 101");
        assert_eq!(i.next().unwrap(), "104 64 94");
        assert_eq!(i.next().unwrap(), "98 58 88");
        assert_eq!(i.next().unwrap(), "86 51 78");
        assert_eq!(i.next().unwrap(), "71 41 64");
        assert_eq!(i.next().unwrap(), "61 29 55");
        assert_eq!(i.next().unwrap(), "54 19 47");
        assert_eq!(i.next().unwrap(), "56 98 67");
        assert_eq!(i.next().unwrap(), "56 194 83");
        assert_eq!(i.next().unwrap(), "57 19 39");
        assert_eq!(i.next().unwrap(), "58 25 43");
        assert_eq!(i.next().unwrap(), "60 32 48");
        assert_eq!(i.next().unwrap(), "63 42 54");
        assert_eq!(i.next().unwrap(), "77 55 68");
        assert_eq!(i.next().unwrap(), "93 62 81");
        assert_eq!(i.next().unwrap(), "103 76 96");
        assert_eq!(i.next().unwrap(), "117 107 118");
        assert_eq!(i.next().unwrap(), "120 115 123");
        assert_eq!(i.next().unwrap(), "123 123 129");
        assert_eq!(i.next().unwrap(), "130 132 136");
        assert_eq!(i.next().unwrap(), "138 138 143");
        assert_eq!(i.next().unwrap(), "144 145 150");
        assert_eq!(i.next().unwrap(), "148 148 153");
        assert_eq!(i.next().unwrap(), "153 155 159");
        assert_eq!(i.next().unwrap(), "154 156 159");
        assert_eq!(i.next().unwrap(), "155 157 160");
        assert_eq!(i.next().unwrap(), "154 157 162");
        assert_eq!(i.next().unwrap(), "152 156 161");
        assert_eq!(i.next().unwrap(), "147 154 160");
        assert_eq!(i.next().unwrap(), "140 152 161");
        assert_eq!(i.next().unwrap(), "134 149 160");
        assert_eq!(i.next().unwrap(), "129 147 159");
        assert_eq!(i.next().unwrap(), "129 141 153");
        assert_eq!(i.next().unwrap(), "129 138 151");
        assert_eq!(i.next().unwrap(), "130 136 150");
        assert_eq!(i.next().unwrap(), "134 126 149");
        assert_eq!(i.next().unwrap(), "135 116 147");
        assert_eq!(i.next().unwrap(), "138 101 137");
        assert_eq!(i.next().unwrap(), "142 85 126");
        assert_eq!(i.next().unwrap(), "150 65 110");
        assert_eq!(i.next().unwrap(), "151 65 111");
        assert_eq!(i.next().unwrap(), "153 66 113");
        assert_eq!(i.next().unwrap(), "155 65 114");
        assert_eq!(i.next().unwrap(), "150 60 108");
        assert_eq!(i.next().unwrap(), "141 49 98");
        assert_eq!(i.next().unwrap(), "122 46 89");
        assert_eq!(i.next().unwrap(), "108 51 83");
        assert_eq!(i.next().unwrap(), "96 63 82");
        assert_eq!(i.next().unwrap(), "83 73 78");
        assert_eq!(i.next().unwrap(), "74 70 68");
        assert_eq!(i.next().unwrap(), "57 65 59");
        assert_eq!(i.next().unwrap(), "45 62 50");
        assert_eq!(i.next().unwrap(), "40 64 49");
        assert_eq!(i.next().unwrap(), "44 75 56");
        assert_eq!(i.next().unwrap(), "57 85 64");
        assert_eq!(i.next().unwrap(), "89 111 94");
        assert_eq!(i.next().unwrap(), "101 120 105");
        assert_eq!(i.next().unwrap(), "113 130 117");
        assert_eq!(i.next().unwrap(), "138 152 144");
        assert_eq!(i.next().unwrap(), "162 170 167");
        assert_eq!(i.next().unwrap(), "180 184 184");
        assert_eq!(i.next().unwrap(), "189 191 192");
        assert_eq!(i.next().unwrap(), "194 196 200");
        assert_eq!(i.next().unwrap(), "192 198 204");
        assert_eq!(i.next().unwrap(), "185 197 204");
        assert_eq!(i.next().unwrap(), "172 190 201");
        assert_eq!(i.next().unwrap(), "153 177 191");
        assert_eq!(i.next().unwrap(), "136 163 182");
        assert_eq!(i.next().unwrap(), "118 152 174");
        assert_eq!(i.next().unwrap(), "104 144 168");
        assert_eq!(i.next().unwrap(), "94 140 165");
        assert_eq!(i.next().unwrap(), "84 138 163");
        assert_eq!(i.next().unwrap(), "81 135 162");
        assert_eq!(i.next().unwrap(), "81 136 164");
        assert_eq!(i.next().unwrap(), "81 137 165");
        assert_eq!(i.next().unwrap(), "81 135 165");
        assert_eq!(i.next().unwrap(), "78 131 159");
        assert_eq!(i.next().unwrap(), "75 116 151");
        assert_eq!(i.next().unwrap(), "75 100 143");
        assert_eq!(i.next().unwrap(), "76 87 136");
        assert_eq!(i.next().unwrap(), "75 76 133");
        assert_eq!(i.next().unwrap(), "71 77 128");
        assert_eq!(i.next().unwrap(), "68 75 126");
        assert_eq!(i.next().unwrap(), "68 72 128");
        assert_eq!(i.next().unwrap(), "78 76 131");
        assert_eq!(i.next().unwrap(), "94 82 140");
        assert_eq!(i.next().unwrap(), "109 101 148");
        assert_eq!(i.next().unwrap(), "122 122 157");
        assert_eq!(i.next().unwrap(), "136 142 174");
        assert_eq!(i.next().unwrap(), "153 167 189");
        assert_eq!(i.next().unwrap(), "175 186 205");
        assert_eq!(i.next().unwrap(), "195 202 218");
        assert_eq!(i.next().unwrap(), "206 216 221");
        assert_eq!(i.next().unwrap(), "209 217 220");
        assert_eq!(i.next().unwrap(), "203 215 212");
        assert_eq!(i.next().unwrap(), "192 205 200");
        assert_eq!(i.next().unwrap(), "178 185 184");
        assert_eq!(i.next().unwrap(), "157 167 167");
        assert_eq!(i.next().unwrap(), "132 142 145");
        assert_eq!(i.next().unwrap(), "107 118 119");
        assert_eq!(i.next().unwrap(), "84 98 99");
        assert_eq!(i.next().unwrap(), "67 70 79");
        assert_eq!(i.next().unwrap(), "57 50 73");
        assert_eq!(i.next().unwrap(), "53 37 75");
        assert_eq!(i.next().unwrap(), "57 32 77");
        assert_eq!(i.next().unwrap(), "67 39 84");
        assert_eq!(i.next().unwrap(), "76 49 88");
        assert_eq!(i.next().unwrap(), "81 58 99");
        assert_eq!(i.next().unwrap(), "86 71 113");
        assert_eq!(i.next().unwrap(), "88 81 128");
        assert_eq!(i.next().unwrap(), "91 92 139");
        assert_eq!(i.next().unwrap(), "88 99 145");
        assert_eq!(i.next().unwrap(), "81 104 145");
        assert_eq!(i.next().unwrap(), "77 111 145");
        assert_eq!(i.next().unwrap(), "74 113 144");
        assert_eq!(i.next().unwrap(), "77 119 142");
        assert_eq!(i.next().unwrap(), "82 118 145");
        assert_eq!(i.next().unwrap(), "90 123 148");
        assert_eq!(i.next().unwrap(), "101 132 156");
        assert_eq!(i.next().unwrap(), "119 142 163");
        assert_eq!(i.next().unwrap(), "139 160 172");
        assert_eq!(i.next().unwrap(), "157 172 181");
        assert_eq!(i.next().unwrap(), "169 179 185");
        assert_eq!(i.next().unwrap(), "173 181 186");
        assert_eq!(i.next().unwrap(), "173 178 182");
        assert_eq!(i.next().unwrap(), "174 175 182");
        assert_eq!(i.next().unwrap(), "179 170 184");

        let palette_data = "
      dac1add9bba3d4b194d0a8857cc14929da e28e3 827ed 2
      37c53183a849cf8c6287c53240ff 23aff 234ff 333ff 1
      32ff 030f3 17bb826c77d4cc56733c4511add49 df741 0
      ff4f 1ff57 1ff5f 1ff7015ff812ae4803dca8051cb8154
      c980533dff 01fea18 1d5312690644c4b974f33c2521cee
      6216ff5a18fa531af5573ed95b62bd9784a3d3a689d6bba0
      dabfaadbd0bedacab6d9c4afd8bda5d8b69bd7b195d6ad8f
      d09269cd895dca8051c66734c34e18c14811c043 bf126 0
      f6 d 0db 05fbe 077a2 18f82 18962 1849f 18cae127f
      d39066d39974d4a383d6b197d9c0acdac5b4dbcbbcdcdbd6
      e8fff2dde7dfdad4c5d7c1acd4b194d1a17dcb8154c86c39
      f441 0f320 0f2 0 1f1 0 0f1 1 0a5 1 0a2 0 09e 01b
      75 14c721bfd75 dfe79 0ff7b 0fe7e 0feae 0b9c1 084
      cb8154cf8d64d39974d6a484d9b094d8bca4dac0a9ffda91
      ffb974d4ad8cd4aa89d5a887d2a1815a6fb25d4ecb5f44d3
       192ed 0a7c9 0bca6 0ba8b 1b971 0c158 0ca3a 0c841
      4ca762d7b599d8baa1d9c0aadac4afdec7b5ffe29cfed88f
      d09269d19b75d2a482d4ac8dd7b599d8c2abe9ebc3eaffe0
      f6fdd3fff1acffeeabffecaafee7a3d7b697efcea3e18657
      5e40d6642de96b1afd7515ff6317ff531af2742f70c87243
      c96d3a8d82305eaa262fd31c22e7 1 1ff 1 0ff 7 0e226
      27adac80bfb7dad1c2dfe7cfe9ffeaedf3cfdbcbbbd4b193
      509f70 1acc8 09fdb1f70ff2740ff222efc2125fa29 6dc
      5b 1947d 09aa0 0a0af 088a7 086a4385ac87644ce875b
      ce875bd18c62ff9742ff8a32fe70 eff6c cfe6e eca6c38
      c96d3ac76936c76935c5602ac5521bff4e 2fe4f 0ff58 0
      ff62 3fe7f22ffa457d5a785d4af929da8ba5c53ca561af4
      6d1afe76 1ff5611ff25 0db27 0d527 0cd39 1bc43 0af
      6c 26795 038ac1f 1ba2c 0c44b16c3541cc65f25c55e24
      c24f18c34e18c04710bf40 9bc36 1b323 0b421 0ec 0 2
      e6 024e1 03fe2 04cc8 07aa7 09d89 0ee84 0f7a1 09c
      b7 08ccb 078db 05ecb7041d0895fd2a07dd6b194d9bea9
";
        let mut palette = vec![Srgba::default(); 256];
        parse_palette_data(palette_data, 3, &mut palette).unwrap();
        let mut i = palette.into_iter().map(|c| c.to_attribute());

        assert_eq!(i.next().unwrap(), "218 193 173");
        assert_eq!(i.next().unwrap(), "217 187 163");
        assert_eq!(i.next().unwrap(), "212 177 148");
        assert_eq!(i.next().unwrap(), "208 168 133");
        assert_eq!(i.next().unwrap(), "124 193 73");
        assert_eq!(i.next().unwrap(), "41 218 226");
        assert_eq!(i.next().unwrap(), "40 227 130");
        assert_eq!(i.next().unwrap(), "39 237 2");
        assert_eq!(i.next().unwrap(), "55 197 49");
        assert_eq!(i.next().unwrap(), "131 168 73");
        assert_eq!(i.next().unwrap(), "207 140 98");
        assert_eq!(i.next().unwrap(), "135 197 50");
        assert_eq!(i.next().unwrap(), "64 255 35");
        assert_eq!(i.next().unwrap(), "58 255 35");
        assert_eq!(i.next().unwrap(), "52 255 51");
        assert_eq!(i.next().unwrap(), "51 255 1");
        assert_eq!(i.next().unwrap(), "50 255 3");
        assert_eq!(i.next().unwrap(), "48 243 23");
        assert_eq!(i.next().unwrap(), "123 184 38");
        assert_eq!(i.next().unwrap(), "199 125 76");
        assert_eq!(i.next().unwrap(), "197 103 51");
        assert_eq!(i.next().unwrap(), "196 81 26");
        assert_eq!(i.next().unwrap(), "221 73 223");
        assert_eq!(i.next().unwrap(), "247 65 0");
        assert_eq!(i.next().unwrap(), "255 79 31");
        assert_eq!(i.next().unwrap(), "255 87 31");
        assert_eq!(i.next().unwrap(), "255 95 31");
        assert_eq!(i.next().unwrap(), "255 112 21");
        assert_eq!(i.next().unwrap(), "255 129 42");
        assert_eq!(i.next().unwrap(), "228 128 61");
        assert_eq!(i.next().unwrap(), "202 128 81");
        assert_eq!(i.next().unwrap(), "203 129 84");
        assert_eq!(i.next().unwrap(), "201 128 83");
        assert_eq!(i.next().unwrap(), "61 255 1");
        assert_eq!(i.next().unwrap(), "31 234 24");
        assert_eq!(i.next().unwrap(), "29 83 18");
        assert_eq!(i.next().unwrap(), "105 6 68");
        assert_eq!(i.next().unwrap(), "196 185 116");
        assert_eq!(i.next().unwrap(), "243 60 37");
        assert_eq!(i.next().unwrap(), "33 206 14");
        assert_eq!(i.next().unwrap(), "98 22 255");
        assert_eq!(i.next().unwrap(), "90 24 250");
        assert_eq!(i.next().unwrap(), "83 26 245");
        assert_eq!(i.next().unwrap(), "87 62 217");
        assert_eq!(i.next().unwrap(), "91 98 189");
        assert_eq!(i.next().unwrap(), "151 132 163");
        assert_eq!(i.next().unwrap(), "211 166 137");
        assert_eq!(i.next().unwrap(), "214 187 160");
        assert_eq!(i.next().unwrap(), "218 191 170");
        assert_eq!(i.next().unwrap(), "219 208 190");
        assert_eq!(i.next().unwrap(), "218 202 182");
        assert_eq!(i.next().unwrap(), "217 196 175");
        assert_eq!(i.next().unwrap(), "216 189 165");
        assert_eq!(i.next().unwrap(), "216 182 155");
        assert_eq!(i.next().unwrap(), "215 177 149");
        assert_eq!(i.next().unwrap(), "214 173 143");
        assert_eq!(i.next().unwrap(), "208 146 105");
        assert_eq!(i.next().unwrap(), "205 137 93");
        assert_eq!(i.next().unwrap(), "202 128 81");
        assert_eq!(i.next().unwrap(), "198 103 52");
        assert_eq!(i.next().unwrap(), "195 78 24");
        assert_eq!(i.next().unwrap(), "193 72 17");
        assert_eq!(i.next().unwrap(), "192 67 191");
        assert_eq!(i.next().unwrap(), "241 38 0");
        assert_eq!(i.next().unwrap(), "246 13 13");
        assert_eq!(i.next().unwrap(), "219 5 251");
        assert_eq!(i.next().unwrap(), "190 7 122");
        assert_eq!(i.next().unwrap(), "162 24 248");
        assert_eq!(i.next().unwrap(), "130 24 150");
        assert_eq!(i.next().unwrap(), "98 24 73");
        assert_eq!(i.next().unwrap(), "159 24 202");
        assert_eq!(i.next().unwrap(), "174 18 127");
        assert_eq!(i.next().unwrap(), "211 144 102");
        assert_eq!(i.next().unwrap(), "211 153 116");
        assert_eq!(i.next().unwrap(), "212 163 131");
        assert_eq!(i.next().unwrap(), "214 177 151");
        assert_eq!(i.next().unwrap(), "217 192 172");
        assert_eq!(i.next().unwrap(), "218 197 180");
        assert_eq!(i.next().unwrap(), "219 203 188");
        assert_eq!(i.next().unwrap(), "220 219 214");
        assert_eq!(i.next().unwrap(), "232 255 242");
        assert_eq!(i.next().unwrap(), "221 231 223");
        assert_eq!(i.next().unwrap(), "218 212 197");
        assert_eq!(i.next().unwrap(), "215 193 172");
        assert_eq!(i.next().unwrap(), "212 177 148");
        assert_eq!(i.next().unwrap(), "209 161 125");
        assert_eq!(i.next().unwrap(), "203 129 84");
        assert_eq!(i.next().unwrap(), "200 108 57");
        assert_eq!(i.next().unwrap(), "244 65 15");
        assert_eq!(i.next().unwrap(), "243 32 15");
        assert_eq!(i.next().unwrap(), "242 0 31");
        assert_eq!(i.next().unwrap(), "241 0 15");
        assert_eq!(i.next().unwrap(), "241 1 10");
        assert_eq!(i.next().unwrap(), "165 1 10");
        assert_eq!(i.next().unwrap(), "162 0 9");
        assert_eq!(i.next().unwrap(), "158 1 11");
        assert_eq!(i.next().unwrap(), "117 20 199");
        assert_eq!(i.next().unwrap(), "114 27 253");
        assert_eq!(i.next().unwrap(), "117 223 231");
        assert_eq!(i.next().unwrap(), "121 15 247");
        assert_eq!(i.next().unwrap(), "123 15 231");
        assert_eq!(i.next().unwrap(), "126 15 234");
        assert_eq!(i.next().unwrap(), "174 11 156");
        assert_eq!(i.next().unwrap(), "193 8 4");
        assert_eq!(i.next().unwrap(), "203 129 84");
        assert_eq!(i.next().unwrap(), "207 141 100");
        assert_eq!(i.next().unwrap(), "211 153 116");
        assert_eq!(i.next().unwrap(), "214 164 132");
        assert_eq!(i.next().unwrap(), "217 176 148");
        assert_eq!(i.next().unwrap(), "216 188 164");
        assert_eq!(i.next().unwrap(), "218 192 169");
        assert_eq!(i.next().unwrap(), "255 218 145");
        assert_eq!(i.next().unwrap(), "255 185 116");
        assert_eq!(i.next().unwrap(), "212 173 140");
        assert_eq!(i.next().unwrap(), "212 170 137");
        assert_eq!(i.next().unwrap(), "213 168 135");
        assert_eq!(i.next().unwrap(), "210 161 129");
        assert_eq!(i.next().unwrap(), "90 111 178");
        assert_eq!(i.next().unwrap(), "93 78 203");
        assert_eq!(i.next().unwrap(), "95 68 211");
        assert_eq!(i.next().unwrap(), "25 46 13");
        assert_eq!(i.next().unwrap(), "10 124 9");
        assert_eq!(i.next().unwrap(), "11 202 6");
        assert_eq!(i.next().unwrap(), "11 168 11");
        assert_eq!(i.next().unwrap(), "27 151 1");
        assert_eq!(i.next().unwrap(), "12 21 8");
        assert_eq!(i.next().unwrap(), "12 163 10");
        assert_eq!(i.next().unwrap(), "12 132 1");
        assert_eq!(i.next().unwrap(), "76 167 98");
        assert_eq!(i.next().unwrap(), "215 181 153");
        assert_eq!(i.next().unwrap(), "216 186 161");
        assert_eq!(i.next().unwrap(), "217 192 170");
        assert_eq!(i.next().unwrap(), "218 196 175");
        assert_eq!(i.next().unwrap(), "222 199 181");
        assert_eq!(i.next().unwrap(), "255 226 156");
        assert_eq!(i.next().unwrap(), "254 216 143");
        assert_eq!(i.next().unwrap(), "208 146 105");
        assert_eq!(i.next().unwrap(), "209 155 117");
        assert_eq!(i.next().unwrap(), "210 164 130");
        assert_eq!(i.next().unwrap(), "212 172 141");
        assert_eq!(i.next().unwrap(), "215 181 153");
        assert_eq!(i.next().unwrap(), "216 194 171");
        assert_eq!(i.next().unwrap(), "233 235 195");
        assert_eq!(i.next().unwrap(), "234 255 224");
        assert_eq!(i.next().unwrap(), "246 253 211");
        assert_eq!(i.next().unwrap(), "255 241 172");
        assert_eq!(i.next().unwrap(), "255 238 171");
        assert_eq!(i.next().unwrap(), "255 236 170");
        assert_eq!(i.next().unwrap(), "254 231 163");
        assert_eq!(i.next().unwrap(), "215 182 151");
        assert_eq!(i.next().unwrap(), "239 206 163");
        assert_eq!(i.next().unwrap(), "225 134 87");
        assert_eq!(i.next().unwrap(), "94 64 214");
        assert_eq!(i.next().unwrap(), "100 45 233");
        assert_eq!(i.next().unwrap(), "107 26 253");
        assert_eq!(i.next().unwrap(), "117 21 255");
        assert_eq!(i.next().unwrap(), "99 23 255");
        assert_eq!(i.next().unwrap(), "83 26 242");
        assert_eq!(i.next().unwrap(), "116 47 112");
        assert_eq!(i.next().unwrap(), "200 114 67");
        assert_eq!(i.next().unwrap(), "201 109 58");
        assert_eq!(i.next().unwrap(), "141 130 48");
        assert_eq!(i.next().unwrap(), "94 170 38");
        assert_eq!(i.next().unwrap(), "47 211 28");
        assert_eq!(i.next().unwrap(), "34 231 1");
        assert_eq!(i.next().unwrap(), "31 15 1");
        assert_eq!(i.next().unwrap(), "15 15 7");
        assert_eq!(i.next().unwrap(), "14 34 6");
        assert_eq!(i.next().unwrap(), "39 173 172");
        assert_eq!(i.next().unwrap(), "128 191 183");
        assert_eq!(i.next().unwrap(), "218 209 194");
        assert_eq!(i.next().unwrap(), "223 231 207");
        assert_eq!(i.next().unwrap(), "233 255 234");
        assert_eq!(i.next().unwrap(), "237 243 207");
        assert_eq!(i.next().unwrap(), "219 203 187");
        assert_eq!(i.next().unwrap(), "212 177 147");
        assert_eq!(i.next().unwrap(), "80 159 112");
        assert_eq!(i.next().unwrap(), "26 204 8");
        assert_eq!(i.next().unwrap(), "9 253 177");
        assert_eq!(i.next().unwrap(), "247 15 242");
        assert_eq!(i.next().unwrap(), "116 15 242");
        assert_eq!(i.next().unwrap(), "34 239 194");
        assert_eq!(i.next().unwrap(), "18 95 162");
        assert_eq!(i.next().unwrap(), "9 109 12");
        assert_eq!(i.next().unwrap(), "91 25 71");
        assert_eq!(i.next().unwrap(), "125 9 170");
        assert_eq!(i.next().unwrap(), "160 10 10");
        assert_eq!(i.next().unwrap(), "175 8 138");
        assert_eq!(i.next().unwrap(), "167 8 106");
        assert_eq!(i.next().unwrap(), "164 56 90");
        assert_eq!(i.next().unwrap(), "200 118 68");
        assert_eq!(i.next().unwrap(), "206 135 91");
        assert_eq!(i.next().unwrap(), "206 135 91");
        assert_eq!(i.next().unwrap(), "209 140 98");
        assert_eq!(i.next().unwrap(), "255 151 66");
        assert_eq!(i.next().unwrap(), "255 138 50");
        assert_eq!(i.next().unwrap(), "254 112 239");
        assert_eq!(i.next().unwrap(), "255 108 207");
        assert_eq!(i.next().unwrap(), "254 110 236");
        assert_eq!(i.next().unwrap(), "202 108 56");
        assert_eq!(i.next().unwrap(), "201 109 58");
        assert_eq!(i.next().unwrap(), "199 105 54");
        assert_eq!(i.next().unwrap(), "199 105 53");
        assert_eq!(i.next().unwrap(), "197 96 42");
        assert_eq!(i.next().unwrap(), "197 82 27");
        assert_eq!(i.next().unwrap(), "255 78 47");
        assert_eq!(i.next().unwrap(), "254 79 15");
        assert_eq!(i.next().unwrap(), "255 88 0");
        assert_eq!(i.next().unwrap(), "255 98 63");
        assert_eq!(i.next().unwrap(), "254 127 34");
        assert_eq!(i.next().unwrap(), "255 164 87");
        assert_eq!(i.next().unwrap(), "213 167 133");
        assert_eq!(i.next().unwrap(), "212 175 146");
        assert_eq!(i.next().unwrap(), "157 168 186");
        assert_eq!(i.next().unwrap(), "92 83 202");
        assert_eq!(i.next().unwrap(), "86 26 244");
        assert_eq!(i.next().unwrap(), "109 26 254");
        assert_eq!(i.next().unwrap(), "118 31 245");
        assert_eq!(i.next().unwrap(), "86 17 255");
        assert_eq!(i.next().unwrap(), "37 13 178");
        assert_eq!(i.next().unwrap(), "39 13 82");
        assert_eq!(i.next().unwrap(), "39 12 211");
        assert_eq!(i.next().unwrap(), "57 27 196");
        assert_eq!(i.next().unwrap(), "67 10 15");
        assert_eq!(i.next().unwrap(), "108 38 121");
        assert_eq!(i.next().unwrap(), "149 3 138");
        assert_eq!(i.next().unwrap(), "172 31 27");
        assert_eq!(i.next().unwrap(), "186 44 12");
        assert_eq!(i.next().unwrap(), "196 75 22");
        assert_eq!(i.next().unwrap(), "195 84 28");
        assert_eq!(i.next().unwrap(), "198 95 37");
        assert_eq!(i.next().unwrap(), "197 94 36");
        assert_eq!(i.next().unwrap(), "194 79 24");
        assert_eq!(i.next().unwrap(), "195 78 24");
        assert_eq!(i.next().unwrap(), "192 71 16");
        assert_eq!(i.next().unwrap(), "191 64 155");
        assert_eq!(i.next().unwrap(), "188 54 27");
        assert_eq!(i.next().unwrap(), "179 35 11");
        assert_eq!(i.next().unwrap(), "180 33 14");
        assert_eq!(i.next().unwrap(), "236 0 2");
        assert_eq!(i.next().unwrap(), "230 2 78");
        assert_eq!(i.next().unwrap(), "225 3 254");
        assert_eq!(i.next().unwrap(), "226 4 204");
        assert_eq!(i.next().unwrap(), "200 7 170");
        assert_eq!(i.next().unwrap(), "167 9 216");
        assert_eq!(i.next().unwrap(), "137 14 232");
        assert_eq!(i.next().unwrap(), "132 15 122");
        assert_eq!(i.next().unwrap(), "161 9 12");
        assert_eq!(i.next().unwrap(), "183 8 204");
        assert_eq!(i.next().unwrap(), "203 7 141");
        assert_eq!(i.next().unwrap(), "219 5 236");
        assert_eq!(i.next().unwrap(), "203 112 65");
        assert_eq!(i.next().unwrap(), "208 137 95");
        assert_eq!(i.next().unwrap(), "210 160 125");
        assert_eq!(i.next().unwrap(), "214 177 148");
        assert_eq!(i.next().unwrap(), "217 190 169");
    }
}
