use std::{
    collections::HashMap,
    io::{Read, Write},
    str::FromStr,
};

use xml::{
    attribute::OwnedAttribute,
    reader::{self, EventReader, ParserConfig},
    writer::{self, EmitterConfig, EventWriter},
};

use crate::{
    color_from_str, parse,
    variations::{self, Var},
    Affine, ColorType, Coordinate, Dimension, Genome, Interpolation, MotionFunction, Rgba,
    TemporalFilter, Transform,
};

const VARIATION_COUNT: usize = 99;

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
                        attrs.push(($name.to_string(), v.variation_weight.to_string()));
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
                    variation.variation_weight = weight;

                    variation = if let Some(_attrs) = attrs {
                        variations::$var {
                            variation_weight: weight,
                            $(
                                $field: if let Some(val) = _attrs.remove($attr) {
                                    parse(&val)?
                                } else {
                                    variation.$field
                                },
                            )*
                        }
                    } else {
                        variation.variation_weight = weight;
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
        weight: "auger_weight",
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

fn strip_whitespace(data: String) -> String {
    // Strips out newlines and any whitespace at the start or end of the lines.
    // Allows a single space at the start of each line to allow for <16 hex
    // values.
    data.split('\n')
        .map(|s| {
            let trimmed = s.trim();
            if trimmed.len() % 2 != 0 {
                format!(" {}", trimmed)
            } else {
                trimmed.to_owned()
            }
        })
        .collect()
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

    if colors != 256 {
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
                let palette = strip_whitespace(palette);

                if palette.len() != colors * bytes * 2 {
                    return Err(format!(
                        "Unexpected length of palette data for {} colors with {} bytes per color: {}",
                        colors,
                        bytes,
                        palette.len()
                    ));
                }

                for index in 0..colors {
                    let pos = index * bytes * 2;
                    genome.palette[index] = Rgba {
                        red: color_from_str(&palette[pos..pos + 2], ColorType::Hex)?,
                        green: color_from_str(&palette[pos + 2..pos + 4], ColorType::Hex)?,
                        blue: color_from_str(&palette[pos + 4..pos + 6], ColorType::Hex)?,
                        alpha: if bytes == 4 {
                            color_from_str(&palette[pos + 4..pos + 6], ColorType::Hex)?
                        } else {
                            1.0
                        },
                    }
                }
            }
            _ => {}
        }
    }
}

fn serialize_colors<W: Write>(genome: &Genome, writer: &mut EventWriter<W>) -> Result<(), String> {
    for (index, color) in genome.palette.iter().enumerate() {
        let mut attrs: Vec<(String, String)> = vec![("index".to_string(), index.to_string())];

        if color.has_opacity() {
            attrs.push(("rgba".to_string(), color.to_str_list(ColorType::Byte)));
        } else {
            attrs.push(("rgb".to_string(), color.to_str_list(ColorType::Byte)));
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
        .and_then(|s| Rgba::from_str_list(s, crate::ColorType::Byte))?;

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
        .to_str_list());

    writep!(attrs, &transform.post, "post", |a: &Affine| a.to_str_list());

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
        transform.coefficients = Affine::from_str_list(&coefs)?;
    }

    if let Some(post) = attrs.remove("post") {
        transform.post = Affine::from_str_list(&post)?;
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
    writep!(attrs, &genome.size, "size", Dimension::to_str_list);
    writep!(attrs, &genome.center, "center", Coordinate::to_str_list);
    writep!(attrs, genome.pixels_per_unit, "scale");

    if genome.zoom > 0.0 {
        writep!(attrs, genome.zoom, "zoom");
    }
    writep!(attrs, genome.rotate, "rotate");
    writep!(attrs, genome.spatial_oversample, "supersample");
    writep!(attrs, genome.spatial_filter_radius, "filter");
    writep!(attrs, genome.spatial_filter_select, "filter_shape");
    writep!(attrs, genome.temporal_filter, "temporal_filter_type");
    if genome.temporal_filter == TemporalFilter::Exp {
        writep!(attrs, genome.temporal_filter_exp, "temporal_filter_exp");
    }
    writep!(attrs, genome.temporal_filter_width, "temporal_filter_width");
    writep!(attrs, genome.sample_density, "quality");
    writep!(attrs, genome.passes, "passes");
    writep!(attrs, genome.num_temporal_samples, "temporal_samples");
    writep!(attrs, &genome.background, "background", |c: &Rgba| c
        .to_str_list(crate::ColorType::Byte));
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
    setp!(attrs, genome.size, "size", Dimension::from_str_list);
    setp!(attrs, genome.center, "center", Coordinate::from_str_list);
    setp!(attrs, genome.pixels_per_unit, "scale");
    setp!(attrs, genome.rotate, "rotate");
    setp!(attrs, genome.zoom, "zoom");
    setp!(attrs, genome.spatial_oversample, "supersample");
    setp!(attrs, genome.spatial_oversample, "oversample");
    setp!(attrs, genome.spatial_filter_radius, "filter");
    setp!(attrs, genome.spatial_filter_select, "filter_shape");
    setp!(attrs, genome.temporal_filter, "temporal_filter_type");
    setp!(attrs, genome.temporal_filter_width, "temporal_filter_width");
    setp!(attrs, genome.temporal_filter_exp, "temporal_filter_exp");
    setp!(attrs, genome.palette_mode, "palette_mode");
    setp!(attrs, genome.sample_density, "quality");
    setp!(attrs, genome.passes, "passes");
    setp!(attrs, genome.num_temporal_samples, "temporal_samples");
    setp!(attrs, genome.background, "background", |s| {
        Rgba::from_str_list(s, crate::ColorType::Byte)
    });
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
