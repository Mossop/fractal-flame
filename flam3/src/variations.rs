use std::f64::consts::PI;

pub trait Variation {
    fn weight(&self) -> f64;
}

macro_rules! variation {
    ($typ:ident) => {
        #[derive(Debug, Default, Clone, PartialEq)]
        pub struct $typ {
            pub variation_weight: f64,
        }

        impl Variation for $typ {
            fn weight(&self) -> f64 {
                self.variation_weight
            }
        }
    };
    ($typ:ident, {
        $(
            $field_vis:vis $field_name:ident : $field_type:ty = $default:expr
        ),*$(,)+
    }) => {
        #[derive(Debug, Clone, PartialEq)]
        pub struct $typ {
            pub variation_weight: f64,
            $(
                $field_vis $field_name: $field_type,
            )*
        }

        impl Default for $typ {
            fn default() -> Self {
                Self {
                    variation_weight: 1.0,
                    $(
                        $field_name: $default,
                    )*
                }
            }
        }

        impl Variation for $typ {
            fn weight(&self) -> f64 {
                self.variation_weight
            }
        }
    };
}

variation!(Linear);
variation!(Sinusoidal);
variation!(Spherical);
variation!(Swirl);
variation!(Horseshoe);
variation!(Polar);
variation!(Handkerchief);
variation!(Heart);
variation!(Disc);
variation!(Spiral);
variation!(Hyperbolic);
variation!(Diamond);
variation!(Ex);
variation!(Julia);
variation!(Bent);
variation!(Waves);
variation!(Fisheye);
variation!(Popcorn);
variation!(Exponential);
variation!(Power);
variation!(Cosine);
variation!(Rings);
variation!(Fan);
variation!(Blob, {
    pub low: f64 = 0.0,
    pub high: f64 = 1.0,
    pub waves: f64 = 1.0,
});
variation!(Pdj, {
    pub a: f64 = 0.0,
    pub b: f64 = 0.0,
    pub c: f64 = 0.0,
    pub d: f64 = 0.0,
});
variation!(Fan2, {
    pub x: f64 = 0.0,
    pub y: f64 = 0.0,
});
variation!(Rings2, {
    pub val: f64 = 0.0,
});
variation!(Eyefish);
variation!(Bubble);
variation!(Cylinder);
variation!(Perspective, {
    pub angle: f64 = 0.0,
    pub distance: f64 = 0.0,
});
variation!(Noise);
variation!(Julian, {
    pub power: f64 = 1.0,
    pub distance: f64 = 1.0,
});
variation!(Juliascope, {
    pub power: f64 = 1.0,
    pub distance: f64 = 1.0,
});
variation!(Blur);
variation!(GaussianBlur);
variation!(RadialBlur, {
    pub angle: f64 = 0.0,
});
variation!(Pie, {
    pub slices: f64 = 6.0,
    pub rotation: f64 = 0.0,
    pub thickness: f64 = 0.5,
});
variation!(Ngon, {
    pub sides: f64 = 5.0,
    pub power: f64 = 3.0,
    pub circle: f64 = 1.0,
    pub corners: f64 = 2.0,
});
variation!(Curl, {
    pub c1: f64 = 1.0,
    pub c2: f64 = 0.0,
});
variation!(Rectangles, {
    pub x: f64 = 1.0,
    pub y: f64 = 1.0,
});
variation!(Arch);
variation!(Tangent);
variation!(Square);
variation!(Rays);
variation!(Blade);
variation!(Secant2);
variation!(Twintrian);
variation!(Cross);
variation!(Disc2, {
    pub rotate: f64 = 0.0,
    pub twist: f64 = 0.0,
});
variation!(SuperShape, {
    pub rnd: f64 = 0.0,
    pub m: f64 = 0.0,
    pub n1: f64 = 1.0,
    pub n2: f64 = 1.0,
    pub n3: f64 = 1.0,
    pub holes: f64 = 0.0,
});
variation!(Flower, {
    pub petals: f64 = 0.0,
    pub holes: f64 = 0.0,
});
variation!(Conic, {
    pub eccentricity: f64 = 1.0,
    pub holes: f64 = 0.0,
});
variation!(Parabola, {
    pub height: f64 = 0.0,
    pub width: f64 = 0.0,
});
variation!(Bent2, {
    pub x: f64 = 1.0,
    pub y: f64 = 1.0,
});
variation!(Bipolar, {
    pub shift: f64 = 0.0,
});
variation!(Boarders);
variation!(Butterfly);
variation!(Cell, {
    pub size: f64 = 1.0,
});
variation!(Cpow, {
    pub r: f64 = 1.0,
    pub i: f64 = 0.0,
    pub power: f64 = 1.0,
});
variation!(Curve, {
    pub x_amp: f64 = 0.0,
    pub y_amp: f64 = 0.0,
    pub x_length: f64 = 1.0,
    pub y_length: f64 = 1.0,
});
variation!(Edisc);
variation!(Elliptic);
variation!(Escher, {
    pub beta: f64 = 0.0,
});
variation!(Foci);
variation!(Lazysusan, {
    pub spin: f64 = 0.0,
    pub space: f64 = 0.0,
    pub twist: f64 = 0.0,
    pub x: f64 = 0.0,
    pub y: f64 = 0.0,
});
variation!(Loonie);
variation!(PreBlur);
variation!(Modulus, {
    pub x: f64 = 0.0,
    pub y: f64 = 0.0,
});
variation!(Oscilloscope, {
    pub separation: f64 = 1.0,
    pub frequency: f64 = PI,
    pub amplitude: f64 = 1.0,
    pub damping: f64 = 0.0,
});
variation!(Polar2);
variation!(Popcorn2, {
    pub x: f64 = 0.0,
    pub y: f64 = 0.0,
    pub c: f64 = 0.0,
});
variation!(Scry);
variation!(Separation, {
    pub x: f64 = 0.0,
    pub y: f64 = 0.0,
    pub x_inside: f64 = 0.0,
    pub y_inside: f64 = 0.0,
});
variation!(Split, {
    pub x_size: f64 = 0.0,
    pub y_size: f64 = 0.0,
});
variation!(Splits, {
    pub x: f64 = 0.0,
    pub y: f64 = 0.0,
});
variation!(Stripes, {
    pub space: f64 = 0.0,
    pub warp: f64 = 0.0,
});
variation!(Wedge, {
    pub angle: f64 = 0.0,
    pub hole: f64 = 0.0,
    pub count: f64 = 1.0,
    pub swirl: f64 = 0.0,
});
variation!(WedgeJulia, {
    pub angle: f64 = 0.0,
    pub count: f64 = 1.0,
    pub power: f64 = 1.0,
    pub dist: f64 = 0.0,
});
variation!(WedgeSph, {
    pub angle: f64 = 0.0,
    pub count: f64 = 1.0,
    pub hole: f64 = 0.0,
    pub swirl: f64 = 0.0,
});
variation!(Whorl, {
    pub inside: f64 = 0.0,
    pub outside: f64 = 0.0,
});
variation!(Waves2, {
    pub x_frequency: f64 = 0.0,
    pub x_scale: f64 = 0.0,
    pub y_frequency: f64 = 0.0,
    pub y_scale: f64 = 0.0,
});
variation!(Exp);
variation!(Log);
variation!(Sin);
variation!(Cos);
variation!(Tan);
variation!(Sec);
variation!(Csc);
variation!(Cot);
variation!(Sinh);
variation!(Cosh);
variation!(Tanh);
variation!(Sech);
variation!(Csch);
variation!(Coth);
variation!(Auger, {
    pub symmetry: f64 = 0.0,
    pub weight: f64 = 0.5,
    pub frequency: f64 = 1.0,
    pub scale: f64 = 1.0,
});
variation!(Flux, {
    pub spread: f64 = 0.0,
});
variation!(Mobius, {
    pub re_a: f64 = 0.0,
    pub re_b: f64 = 0.0,
    pub re_c: f64 = 0.0,
    pub re_d: f64 = 0.0,
    pub im_a: f64 = 0.0,
    pub im_b: f64 = 0.0,
    pub im_c: f64 = 0.0,
    pub im_d: f64 = 0.0,
});

macro_rules! gen_var {
    {$($var:ident),*$(,)?} => {
        #[derive(Debug, Clone, PartialEq)]
        pub enum Var {
            $(
                $var($var),
            )*
        }

        $(
            impl From<$var> for Var {
                fn from(var: $var) -> Var {
                    Var::$var(var)
                }
            }
        )*
    }
}

gen_var! {
    Linear,
    Sinusoidal,
    Spherical,
    Swirl,
    Horseshoe,
    Polar,
    Handkerchief,
    Heart,
    Disc,
    Spiral,
    Hyperbolic,
    Diamond,
    Ex,
    Julia,
    Bent,
    Waves,
    Fisheye,
    Popcorn,
    Exponential,
    Power,
    Cosine,
    Rings,
    Fan,
    Blob,
    Pdj,
    Fan2,
    Rings2,
    Eyefish,
    Bubble,
    Cylinder,
    Perspective,
    Noise,
    Julian,
    Juliascope,
    Blur,
    GaussianBlur,
    RadialBlur,
    Pie,
    Ngon,
    Curl,
    Rectangles,
    Arch,
    Tangent,
    Square,
    Rays,
    Blade,
    Secant2,
    Twintrian,
    Cross,
    Disc2,
    SuperShape,
    Flower,
    Conic,
    Parabola,
    Bent2,
    Bipolar,
    Boarders,
    Butterfly,
    Cell,
    Cpow,
    Curve,
    Edisc,
    Elliptic,
    Escher,
    Foci,
    Lazysusan,
    Loonie,
    PreBlur,
    Modulus,
    Oscilloscope,
    Polar2,
    Popcorn2,
    Scry,
    Separation,
    Split,
    Splits,
    Stripes,
    Wedge,
    WedgeJulia,
    WedgeSph,
    Whorl,
    Waves2,
    Exp,
    Log,
    Sin,
    Cos,
    Tan,
    Sec,
    Csc,
    Cot,
    Sinh,
    Cosh,
    Tanh,
    Sech,
    Csch,
    Coth,
    Auger,
    Flux,
    Mobius,
}

#[macro_export]
macro_rules! with_var {
    ($var:ident, $with:ident, $body:expr) => {
        match $var {
            Var::Linear(ref $with) => $body,
            Var::Sinusoidal(ref $with) => $body,
            Var::Spherical(ref $with) => $body,
            Var::Swirl(ref $with) => $body,
            Var::Horseshoe(ref $with) => $body,
            Var::Polar(ref $with) => $body,
            Var::Handkerchief(ref $with) => $body,
            Var::Heart(ref $with) => $body,
            Var::Disc(ref $with) => $body,
            Var::Spiral(ref $with) => $body,
            Var::Hyperbolic(ref $with) => $body,
            Var::Diamond(ref $with) => $body,
            Var::Ex(ref $with) => $body,
            Var::Julia(ref $with) => $body,
            Var::Bent(ref $with) => $body,
            Var::Waves(ref $with) => $body,
            Var::Fisheye(ref $with) => $body,
            Var::Popcorn(ref $with) => $body,
            Var::Exponential(ref $with) => $body,
            Var::Power(ref $with) => $body,
            Var::Cosine(ref $with) => $body,
            Var::Rings(ref $with) => $body,
            Var::Fan(ref $with) => $body,
            Var::Blob(ref $with) => $body,
            Var::Pdj(ref $with) => $body,
            Var::Fan2(ref $with) => $body,
            Var::Rings2(ref $with) => $body,
            Var::Eyefish(ref $with) => $body,
            Var::Bubble(ref $with) => $body,
            Var::Cylinder(ref $with) => $body,
            Var::Perspective(ref $with) => $body,
            Var::Noise(ref $with) => $body,
            Var::Julian(ref $with) => $body,
            Var::Juliascope(ref $with) => $body,
            Var::Blur(ref $with) => $body,
            Var::GaussianBlur(ref $with) => $body,
            Var::RadialBlur(ref $with) => $body,
            Var::Pie(ref $with) => $body,
            Var::Ngon(ref $with) => $body,
            Var::Curl(ref $with) => $body,
            Var::Rectangles(ref $with) => $body,
            Var::Arch(ref $with) => $body,
            Var::Tangent(ref $with) => $body,
            Var::Square(ref $with) => $body,
            Var::Rays(ref $with) => $body,
            Var::Blade(ref $with) => $body,
            Var::Secant2(ref $with) => $body,
            Var::Twintrian(ref $with) => $body,
            Var::Cross(ref $with) => $body,
            Var::Disc2(ref $with) => $body,
            Var::SuperShape(ref $with) => $body,
            Var::Flower(ref $with) => $body,
            Var::Conic(ref $with) => $body,
            Var::Parabola(ref $with) => $body,
            Var::Bent2(ref $with) => $body,
            Var::Bipolar(ref $with) => $body,
            Var::Boarders(ref $with) => $body,
            Var::Butterfly(ref $with) => $body,
            Var::Cell(ref $with) => $body,
            Var::Cpow(ref $with) => $body,
            Var::Curve(ref $with) => $body,
            Var::Edisc(ref $with) => $body,
            Var::Elliptic(ref $with) => $body,
            Var::Escher(ref $with) => $body,
            Var::Foci(ref $with) => $body,
            Var::Lazysusan(ref $with) => $body,
            Var::Loonie(ref $with) => $body,
            Var::PreBlur(ref $with) => $body,
            Var::Modulus(ref $with) => $body,
            Var::Oscilloscope(ref $with) => $body,
            Var::Polar2(ref $with) => $body,
            Var::Popcorn2(ref $with) => $body,
            Var::Scry(ref $with) => $body,
            Var::Separation(ref $with) => $body,
            Var::Split(ref $with) => $body,
            Var::Splits(ref $with) => $body,
            Var::Stripes(ref $with) => $body,
            Var::Wedge(ref $with) => $body,
            Var::WedgeJulia(ref $with) => $body,
            Var::WedgeSph(ref $with) => $body,
            Var::Whorl(ref $with) => $body,
            Var::Waves2(ref $with) => $body,
            Var::Exp(ref $with) => $body,
            Var::Log(ref $with) => $body,
            Var::Sin(ref $with) => $body,
            Var::Cos(ref $with) => $body,
            Var::Tan(ref $with) => $body,
            Var::Sec(ref $with) => $body,
            Var::Csc(ref $with) => $body,
            Var::Cot(ref $with) => $body,
            Var::Sinh(ref $with) => $body,
            Var::Cosh(ref $with) => $body,
            Var::Tanh(ref $with) => $body,
            Var::Sech(ref $with) => $body,
            Var::Csch(ref $with) => $body,
            Var::Coth(ref $with) => $body,
            Var::Auger(ref $with) => $body,
            Var::Flux(ref $with) => $body,
            Var::Mobius(ref $with) => $body,
        }
    };
}

impl Variation for Var {
    fn weight(&self) -> f64 {
        with_var!(self, v, v.weight())
    }
}
