use std::f64::consts::{FRAC_1_PI, FRAC_2_PI, FRAC_PI_2, FRAC_PI_4, PI};

use lazy_static::lazy_static;

use crate::{
    acos, atan2, cos, cosh, exp, ln, log10, pow, sin, sincos, sinh, sqr, sqrt, tan,
    utils::PanicCast,
    variations::{self, Var, Variation},
    with_var, Affine, Coordinate, Transform,
};

use super::{adjust_percentage, rng::Flam3Rng};

const EPS: f64 = 1e-10;
const M_PI: f64 = PI;
const M_PI_2: f64 = FRAC_PI_2;
const M_1_PI: f64 = FRAC_1_PI;
const M_2_PI: f64 = FRAC_2_PI;
const M_PI_4: f64 = FRAC_PI_4;

lazy_static! {
    static ref BUTTERFLY_WEIGHT: f64 = 4.0 / sqrt!(3.0 * M_PI);
}

fn badvalue(x: f64) -> bool {
    !x.is_finite() || (x > 1e10) || (x < -1e10)
}

macro_rules! upsert {
    ($var:expr, $($body:tt)*) => {{
        if !$var.is_some() {
            $var = Some($($body)*);
        }
        $var.unwrap()
    }}
}

#[derive(Default)]
pub struct VariationPrecalculations {
    persp_vsin: Option<f64>,
    persp_vfcos: Option<f64>,
    julian_r_n: Option<f64>,
    julian_cn: Option<f64>,
    wedge_julia_cf: Option<f64>,
    wedge_julia_r_n: Option<f64>,
    wedge_julia_cn: Option<f64>,
    juliascope_r_n: Option<f64>,
    juliascope_cn: Option<f64>,
    radial_blur_spinvar: Option<f64>,
    radial_blur_zoomvar: Option<f64>,
    waves_dx2: Option<f64>,
    waves_dy2: Option<f64>,
    disc2_timespi: Option<f64>,
    disc2_cosadd: Option<f64>,
    disc2_sinadd: Option<f64>,
    super_shape_pm_4: Option<f64>,
    super_shape_pneg1_n1: Option<f64>,
}

impl VariationPrecalculations {
    pub fn new(transform: &Transform) -> Self {
        let mut precalcs = Self::default();

        for var in transform.variations.iter() {
            match var {
                Var::Perspective(ref v) => precalcs.perspective_precalc(v),
                Var::Julian(ref v) => precalcs.julian_precalc(v),
                Var::WedgeJulia(ref v) => precalcs.wedge_julia_precalc(v),
                Var::Juliascope(ref v) => precalcs.juliascope_precalc(v),
                Var::RadialBlur(ref v) => precalcs.radial_blur_precalc(v),
                Var::Waves(ref v) => precalcs.waves_precalc(transform, v),
                Var::Disc2(ref v) => precalcs.disc2_precalc(v),
                Var::SuperShape(ref v) => precalcs.super_shape_precalc(v),
                _ => {}
            }
        }

        precalcs
    }

    fn perspective_precalc(&mut self, var: &variations::Perspective) {
        let ang = var.angle * M_PI / 2.0;
        self.persp_vsin = Some(sin!(ang));
        self.persp_vfcos = Some(var.distance * cos!(ang));
    }

    fn julian_precalc(&mut self, var: &variations::Julian) {
        self.julian_r_n = Some(var.power.abs());
        self.julian_cn = Some(var.distance / var.power / 2.0);
    }

    fn wedge_julia_precalc(&mut self, var: &variations::WedgeJulia) {
        self.wedge_julia_cf = Some(1.0 - var.angle * var.count * M_1_PI * 0.5);
        self.wedge_julia_r_n = Some(var.power.abs());
        self.wedge_julia_cn = Some(var.dist / var.power / 2.0);
    }

    fn juliascope_precalc(&mut self, var: &variations::Juliascope) {
        self.juliascope_r_n = Some(var.power.abs());
        self.juliascope_cn = Some(var.distance / var.power / 2.0);
    }

    fn radial_blur_precalc(&mut self, var: &variations::RadialBlur) {
        let (spinvar, zoomvar) = sincos!(var.angle * M_PI / 2.0);
        self.radial_blur_spinvar = Some(spinvar);
        self.radial_blur_zoomvar = Some(zoomvar);
    }

    fn waves_precalc(&mut self, transform: &Transform, _var: &variations::Waves) {
        let dx = transform.coefficients[2][0];
        let dy = transform.coefficients[2][1];

        self.waves_dx2 = Some(1.0 / (dx * dx + EPS));
        self.waves_dy2 = Some(1.0 / (dy * dy + EPS));
    }

    fn disc2_precalc(&mut self, var: &variations::Disc2) {
        let add = var.twist;

        self.disc2_timespi = Some(var.rotate * M_PI);

        let (mut sinadd, mut cosadd) = sincos!(add);
        cosadd -= 1.0;

        if add > 2.0 * M_PI {
            let k = 1.0 + add - 2.0 * M_PI;
            cosadd *= k;
            sinadd *= k;
        }

        if add < -2.0 * M_PI {
            let k = 1.0 + add + 2.0 * M_PI;
            cosadd *= k;
            sinadd *= k;
        }

        self.disc2_cosadd = Some(cosadd);
        self.disc2_sinadd = Some(sinadd);
    }

    fn super_shape_precalc(&mut self, var: &variations::SuperShape) {
        self.super_shape_pm_4 = Some(var.m / 4.0);
        self.super_shape_pneg1_n1 = Some(-1.0 / var.n1);
    }
}

pub trait Flam3Variation: Variation {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    );

    fn is_pre(&self) -> bool {
        false
    }
}

pub struct Flam3IterHelper<'a> {
    p0: f64,
    p1: f64,
    tx: f64,
    ty: f64,

    rc: &'a mut Flam3Rng,

    sumsq: Option<f64>,
    sqrt: Option<f64>,
    atan: Option<f64>,
    sina: Option<f64>,
    cosa: Option<f64>,
    atanyx: Option<f64>,
}

impl<'a> Flam3IterHelper<'a> {
    fn new(coords: Coordinate, rc: &'a mut Flam3Rng) -> Self {
        Self {
            p0: 0.0,
            p1: 0.0,
            tx: coords.x,
            ty: coords.y,
            rc,
            sumsq: None,
            sqrt: None,
            atan: None,
            sina: None,
            cosa: None,
            atanyx: None,
        }
    }

    fn sumsq(&mut self) -> f64 {
        upsert!(self.sumsq, self.tx * self.tx + self.ty * self.ty)
    }

    fn sqrt(&mut self) -> f64 {
        upsert!(self.sqrt, sqrt!(self.sumsq()))
    }

    fn atan(&mut self) -> f64 {
        upsert!(self.atan, atan2!(self.tx, self.ty))
    }

    fn cosa(&mut self) -> f64 {
        upsert!(self.cosa, self.ty / self.sqrt())
    }

    fn sina(&mut self) -> f64 {
        upsert!(self.sina, self.tx / self.sqrt())
    }

    fn atanyx(&mut self) -> f64 {
        upsert!(self.atanyx, atan2!(self.ty, self.tx))
    }
}

impl Flam3Variation for Var {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        with_var!(self, v, v.apply(f, coeffs, precalc))
    }

    fn is_pre(&self) -> bool {
        with_var!(self, v, v.is_pre())
    }
}

impl Flam3Variation for variations::Linear {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        f.p0 += self.weight * f.tx;
        f.p1 += self.weight * f.ty;
    }
}

impl Flam3Variation for variations::Sinusoidal {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        f.p0 += self.weight * sin!(f.tx);
        f.p1 += self.weight * sin!(f.ty);
    }
}

impl Flam3Variation for variations::Spherical {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r2 = self.weight / (f.sumsq() + EPS);

        f.p0 += r2 * f.tx;
        f.p1 += r2 * f.ty;
    }
}

impl Flam3Variation for variations::Swirl {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r2 = f.sumsq();

        let (c1, c2) = sincos!(r2);
        let nx = c1 * f.tx - c2 * f.ty;
        let ny = c2 * f.tx + c1 * f.ty;

        f.p0 += self.weight * nx;
        f.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Horseshoe {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = self.weight / (f.sqrt() + EPS);

        f.p0 += (f.tx - f.ty) * (f.tx + f.ty) * r;
        f.p1 += 2.0 * f.tx * f.ty * r;
    }
}

impl Flam3Variation for variations::Polar {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let nx = f.atan() * M_1_PI;
        let ny = f.sqrt() - 1.0;

        f.p0 += self.weight * nx;
        f.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Handkerchief {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let a = f.atan();
        let r = f.sqrt();

        f.p0 += self.weight * r * sin!(a + r);
        f.p1 += self.weight * r * cos!(a - r);
    }
}

impl Flam3Variation for variations::Heart {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let a = f.sqrt() * f.atan();
        let r = self.weight * f.sqrt();

        let (sa, ca) = sincos!(a);

        f.p0 += r * sa;
        f.p1 += (-r) * ca;
    }
}

impl Flam3Variation for variations::Disc {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let a = f.atan() * M_1_PI;
        let r = M_PI * f.sqrt();
        let (sr, cr) = sincos!(r);

        f.p0 += self.weight * sr * a;
        f.p1 += self.weight * cr * a;
    }
}

impl Flam3Variation for variations::Spiral {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = f.sqrt() + EPS;
        let r1 = self.weight / r;
        let (sr, cr) = sincos!(r);

        f.p0 += r1 * (f.cosa() + sr);
        f.p1 += r1 * (f.sina() - cr);
    }
}

impl Flam3Variation for variations::Hyperbolic {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = f.sqrt() + EPS;

        f.p0 += self.weight * f.sina() / r;
        f.p1 += self.weight * f.cosa() * r;
    }
}

impl Flam3Variation for variations::Diamond {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = f.sqrt();
        let (sr, cr) = sincos!(r);

        f.p0 += self.weight * f.sina() * cr;
        f.p1 += self.weight * f.cosa() * sr;
    }
}

impl Flam3Variation for variations::Ex {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let a = f.atan();
        let r = f.sqrt();

        let n0 = sin!(a + r);
        let n1 = cos!(a - r);

        let m0 = n0 * n0 * n0 * r;
        let m1 = n1 * n1 * n1 * r;

        f.p0 += self.weight * (m0 + m1);
        f.p1 += self.weight * (m0 - m1);
    }
}

impl Flam3Variation for variations::Julia {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let mut a = 0.5 * f.atan();

        if f.rc.isaac_bit() {
            a += M_PI;
        }

        let r = self.weight * sqrt!(f.sqrt());

        let (sa, ca) = sincos!(a);

        f.p0 += r * ca;
        f.p1 += r * sa;
    }
}

impl Flam3Variation for variations::Bent {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let mut nx = f.tx;
        let mut ny = f.ty;

        if nx < 0.0 {
            nx *= 2.0;
        }
        if ny < 0.0 {
            ny /= 2.0;
        }

        f.p0 += self.weight * nx;
        f.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Waves {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        let c10 = coeffs[1][0];
        let c11 = coeffs[1][1];

        let nx = f.tx + c10 * sin!(f.ty * precalc.waves_dx2.unwrap());
        let ny = f.ty + c11 * sin!(f.tx * precalc.waves_dy2.unwrap());

        f.p0 += self.weight * nx;
        f.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Fisheye {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = 2.0 * self.weight / (f.sqrt() + 1.0);

        f.p0 += r * f.ty;
        f.p1 += r * f.tx;
    }
}

impl Flam3Variation for variations::Popcorn {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let dx = tan!(3.0 * f.ty);
        let dy = tan!(3.0 * f.tx);

        let nx = f.tx + coeffs[2][0] * sin!(dx);
        let ny = f.ty + coeffs[2][1] * sin!(dy);

        f.p0 += self.weight * nx;
        f.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Exponential {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let dx = self.weight * exp!(f.tx - 1.0);
        let dy = M_PI * f.ty;

        let (sdy, cdy) = sincos!(dy);

        f.p0 += dx * cdy;
        f.p1 += dx * sdy;
    }
}

impl Flam3Variation for variations::Power {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = self.weight * pow!(f.sqrt(), f.sina());

        f.p0 += r * f.cosa();
        f.p1 += r * f.sina();
    }
}

impl Flam3Variation for variations::Cosine {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let a = f.tx * M_PI;

        let (sa, ca) = sincos!(a);
        let nx = ca * cosh!(f.ty);
        let ny = -sa * sinh!(f.ty);

        f.p0 += self.weight * nx;
        f.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Rings {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let dx = coeffs[2][0] * coeffs[2][0] + EPS;
        let mut r = f.sqrt();
        r = self.weight * (((r + dx) % (2.0 * dx)) - dx + r * (1.0 - dx));

        f.p0 += r * f.cosa();
        f.p1 += r * f.sina();
    }
}

impl Flam3Variation for variations::Fan {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let dx = M_PI * (coeffs[2][0] * coeffs[2][0] + EPS);
        let dy = coeffs[2][1];
        let dx2 = 0.5 * dx;

        let mut a = f.atan();
        let r = self.weight * f.sqrt();

        a += if ((a + dy) % dx) > dx2 { -dx2 } else { dx2 };
        let (sa, ca) = sincos!(a);

        f.p0 += r * ca;
        f.p1 += r * sa;
    }
}

impl Flam3Variation for variations::Blob {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let mut r = f.sqrt();
        let a = f.atan();
        let bdiff = self.high - self.low;

        r *= self.low + bdiff * (0.5 + 0.5 * sin!(self.waves * a));

        f.p0 += self.weight * f.sina() * r;
        f.p1 += self.weight * f.cosa() * r;
    }
}

impl Flam3Variation for variations::Pdj {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let nx1 = cos!(self.b * f.tx);
        let nx2 = sin!(self.c * f.tx);
        let ny1 = sin!(self.a * f.ty);
        let ny2 = cos!(self.d * f.ty);

        f.p0 += self.weight * (ny1 - nx1);
        f.p1 += self.weight * (nx2 - ny2);
    }
}

impl Flam3Variation for variations::Fan2 {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let dy = self.y;
        let dx = M_PI * (self.x * self.x + EPS);
        let dx2 = 0.5 * dx;
        let mut a = f.atan();
        let r = self.weight * f.sqrt();

        let t = a + dy - dx * ((a + dy) / dx).trunc();

        if t > dx2 {
            a -= dx2;
        } else {
            a += dx2;
        }

        let (sa, ca) = sincos!(a);

        f.p0 += r * sa;
        f.p1 += r * ca;
    }
}

impl Flam3Variation for variations::Rings2 {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let mut r = f.sqrt();
        let dx = self.val * self.val + EPS;

        r += -2.0 * dx * ((r + dx) / (2.0 * dx)).trunc() + r * (1.0 - dx);

        f.p0 += self.weight * f.sina() * r;
        f.p1 += self.weight * f.cosa() * r;
    }
}

impl Flam3Variation for variations::Eyefish {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = (self.weight * 2.0) / (f.sqrt() + 1.0);

        f.p0 += r * f.tx;
        f.p1 += r * f.ty;
    }
}

impl Flam3Variation for variations::Bubble {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = self.weight / (0.25 * (f.sumsq()) + 1.0);

        f.p0 += r * f.tx;
        f.p1 += r * f.ty;
    }
}

impl Flam3Variation for variations::Cylinder {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        f.p0 += self.weight * sin!(f.tx);
        f.p1 += self.weight * f.ty;
    }
}

impl Flam3Variation for variations::Perspective {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        let t = 1.0 / (self.distance - f.ty * precalc.persp_vsin.unwrap());

        f.p0 += self.weight * self.distance * f.tx * t;
        f.p1 += self.weight * precalc.persp_vfcos.unwrap() * f.ty * t;
    }
}

impl Flam3Variation for variations::Noise {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let tmpr = f.rc.isaac_01() * 2.0 * M_PI;
        let (sinr, cosr) = sincos!(tmpr);

        let r = self.weight * f.rc.isaac_01();

        f.p0 += f.tx * r * cosr;
        f.p1 += f.ty * r * sinr;
    }
}

impl Flam3Variation for variations::Julian {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        let t_rnd = ((precalc.julian_r_n.unwrap()) * f.rc.isaac_01()).trunc();

        let tmpr = (f.atanyx() + 2.0 * M_PI * t_rnd) / self.power;

        let r = self.weight * pow!(f.sumsq(), precalc.julian_cn.unwrap());
        let (sina, cosa) = sincos!(tmpr);

        f.p0 += r * cosa;
        f.p1 += r * sina;
    }
}

impl Flam3Variation for variations::Juliascope {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        let t_rnd = ((precalc.juliascope_r_n.unwrap()) * f.rc.isaac_01()).trunc();

        let tmpr = if ((t_rnd.i32()) & 1) == 0 {
            (2.0 * M_PI * t_rnd + f.atanyx()) / self.power
        } else {
            (2.0 * M_PI * t_rnd - f.atanyx()) / self.power
        };

        let (sina, cosa) = sincos!(tmpr);

        let r = self.weight * pow!(f.sumsq(), precalc.juliascope_cn.unwrap());

        f.p0 += r * cosa;
        f.p1 += r * sina;
    }
}

impl Flam3Variation for variations::Blur {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let tmpr = f.rc.isaac_01() * 2.0 * M_PI;
        let (sinr, cosr) = sincos!(tmpr);

        let r = self.weight * f.rc.isaac_01();

        f.p0 += r * cosr;
        f.p1 += r * sinr;
    }
}

impl Flam3Variation for variations::GaussianBlur {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let ang = f.rc.isaac_01() * 2.0 * M_PI;
        let (sina, cosa) = sincos!(ang);

        let r = self.weight
            * (f.rc.isaac_01() + f.rc.isaac_01() + f.rc.isaac_01() + f.rc.isaac_01() - 2.0);

        f.p0 += r * cosa;
        f.p1 += r * sina;
    }
}

impl Flam3Variation for variations::RadialBlur {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        /* Get pseudo-gaussian */
        let rnd = self.weight
            * (f.rc.isaac_01() + f.rc.isaac_01() + f.rc.isaac_01() + f.rc.isaac_01() - 2.0);

        /* Calculate angle & zoom */
        let ra = f.sqrt();
        let tmpa = f.atanyx() + precalc.radial_blur_spinvar.unwrap() * rnd;
        let (sa, ca) = sincos!(tmpa);
        let rz = precalc.radial_blur_zoomvar.unwrap() * rnd - 1.0;

        f.p0 += ra * ca + rz * f.tx;
        f.p1 += ra * sa + rz * f.ty;
    }
}

impl Flam3Variation for variations::Pie {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let sl = (f.rc.isaac_01() * self.slices + 0.5).trunc();
        let a = self.rotation + 2.0 * M_PI * (sl + f.rc.isaac_01() * self.thickness) / self.slices;
        let r = self.weight * f.rc.isaac_01();
        let (sa, ca) = sincos!(a);

        f.p0 += r * ca;
        f.p1 += r * sa;
    }
}

impl Flam3Variation for variations::Ngon {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r_factor = pow!(f.sumsq(), self.power / 2.0);

        let theta = f.atanyx();
        let b = 2.0 * M_PI / self.sides;

        let mut phi = theta - (b * (theta / b).floor());
        if phi > b / 2.0 {
            phi -= b;
        }

        let mut amp = self.corners * (1.0 / (cos!(phi) + EPS) - 1.0) + self.circle;
        amp /= r_factor + EPS;

        f.p0 += self.weight * f.tx * amp;
        f.p1 += self.weight * f.ty * amp;
    }
}

impl Flam3Variation for variations::Curl {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let re = 1.0 + self.c1 * f.tx + self.c2 * (f.tx * f.tx - f.ty * f.ty);
        let im = self.c1 * f.ty + 2.0 * self.c2 * f.tx * f.ty;

        let r = self.weight / (re * re + im * im);

        f.p0 += (f.tx * re + f.ty * im) * r;
        f.p1 += (f.ty * re - f.tx * im) * r;
    }
}

impl Flam3Variation for variations::Rectangles {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        if self.x == 0.0 {
            f.p0 += self.weight * f.tx;
        } else {
            f.p0 += self.weight * ((2.0 * (f.tx / self.x).floor() + 1.0) * self.x - f.tx);
        }

        if self.y == 0.0 {
            f.p1 += self.weight * f.ty;
        } else {
            f.p1 += self.weight * ((2.0 * (f.ty / self.y).floor() + 1.0) * self.y - f.ty);
        }
    }
}

impl Flam3Variation for variations::Arch {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let ang = f.rc.isaac_01() * self.weight * M_PI;
        let (sinr, cosr) = sincos!(ang);

        f.p0 += self.weight * sinr;
        f.p1 += self.weight * (sinr * sinr) / cosr;
    }
}

impl Flam3Variation for variations::Tangent {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        f.p0 += self.weight * sin!(f.tx) / cos!(f.ty);
        f.p1 += self.weight * tan!(f.ty);
    }
}

impl Flam3Variation for variations::Square {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        f.p0 += self.weight * (f.rc.isaac_01() - 0.5);
        f.p1 += self.weight * (f.rc.isaac_01() - 0.5);
    }
}

impl Flam3Variation for variations::Rays {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let ang = self.weight * f.rc.isaac_01() * M_PI;
        let r = self.weight / (f.sumsq() + EPS);
        let tanr = self.weight * tan!(ang) * r;

        f.p0 += tanr * cos!(f.tx);
        f.p1 += tanr * sin!(f.ty);
    }
}

impl Flam3Variation for variations::Blade {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let r = f.rc.isaac_01() * self.weight * f.sqrt();

        let (sinr, cosr) = sincos!(r);

        f.p0 += self.weight * f.tx * (cosr + sinr);
        f.p1 += self.weight * f.tx * (cosr - sinr);
    }
}

impl Flam3Variation for variations::Secant2 {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /* Intended as a 'fixed' version of secant */

        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let r = self.weight * f.sqrt();
        let cr = cos!(r);
        let icr = 1.0 / cr;

        f.p0 += self.weight * f.tx;

        if cr < 0.0 {
            f.p1 += self.weight * (icr + 1.0);
        } else {
            f.p1 += self.weight * (icr - 1.0);
        }
    }
}

impl Flam3Variation for variations::Twintrian {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */
        let r = f.rc.isaac_01() * self.weight * f.sqrt();

        let (sinr, cosr) = sincos!(r);
        let mut diff = log10!(sinr * sinr) + cosr;

        if badvalue(diff) {
            diff = -30.0;
        }

        f.p0 += self.weight * f.tx * diff;
        f.p1 += self.weight * f.tx * (diff - sinr * M_PI);
    }
}

impl Flam3Variation for variations::Cross {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let s = f.tx * f.tx - f.ty * f.ty;
        let r = self.weight * sqrt!(1.0 / (sqr!(s) + EPS));

        f.p0 += f.tx * r;
        f.p1 += f.ty * r;
    }
}

impl Flam3Variation for variations::Disc2 {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        let t = precalc.disc2_timespi.unwrap() * (f.tx + f.ty);
        let (sinr, cosr) = sincos!(t);
        let r = self.weight * f.atan() / M_PI;

        f.p0 += (sinr + precalc.disc2_cosadd.unwrap()) * r;
        f.p1 += (cosr + precalc.disc2_sinadd.unwrap()) * r;
    }
}

impl Flam3Variation for variations::SuperShape {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        let theta = precalc.super_shape_pm_4.unwrap() * f.atanyx() + M_PI_4;

        let (st, ct) = sincos!(theta);

        let t1 = pow!(ct.abs(), self.n2);
        let t2 = pow!(st.abs(), self.n3);

        let myrnd = self.rnd;

        let r = self.weight
            * ((myrnd * f.rc.isaac_01() + (1.0 - myrnd) * f.sqrt()) - self.holes)
            * pow!(t1 + t2, precalc.super_shape_pneg1_n1.unwrap())
            / f.sqrt();

        f.p0 += r * f.tx;
        f.p1 += r * f.ty;
    }
}

impl Flam3Variation for variations::Flower {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let theta = f.atanyx();
        let r = self.weight * (f.rc.isaac_01() - self.holes) * cos!(self.petals * theta) / f.sqrt();

        f.p0 += r * f.tx;
        f.p1 += r * f.ty;
    }
}

impl Flam3Variation for variations::Conic {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let ct = f.tx / f.sqrt();
        let r = self.weight * (f.rc.isaac_01() - self.holes) * self.eccentricity
            / (1.0 + self.eccentricity * ct)
            / f.sqrt();

        f.p0 += r * f.tx;
        f.p1 += r * f.ty;
    }
}

impl Flam3Variation for variations::Parabola {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let r = f.sqrt();

        let (sr, cr) = sincos!(r);

        f.p0 += self.height * self.weight * sr * sr * f.rc.isaac_01();
        f.p1 += self.width * self.weight * cr * f.rc.isaac_01();
    }
}

impl Flam3Variation for variations::Bent2 {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let mut nx = f.tx;
        let mut ny = f.ty;

        if nx < 0.0 {
            nx *= self.x;
        }
        if ny < 0.0 {
            ny *= self.y;
        }

        f.p0 += self.weight * nx;
        f.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Bipolar {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let x2y2 = f.sumsq();
        let t = x2y2 + 1.0;
        let x2 = 2.0 * f.tx;
        let ps = -M_PI_2 * self.shift;
        let mut y = 0.5 * atan2!(2.0 * f.ty, x2y2 - 1.0) + ps;

        if y > M_PI_2 {
            y = -M_PI_2 + ((y + M_PI_2) % M_PI);
        } else if y < -M_PI_2 {
            y = M_PI_2 - ((M_PI_2 - y) % M_PI);
        }

        f.p0 += self.weight * 0.25 * M_2_PI * ln!((t + x2) / (t - x2));
        f.p1 += self.weight * M_2_PI * y;
    }
}

impl Flam3Variation for variations::Boarders {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let round_x = f.tx.round();
        let round_y = f.ty.round();
        let offset_x = f.tx - round_x;
        let offset_y = f.ty - round_y;

        if f.rc.isaac_01() >= 0.75 {
            f.p0 += self.weight * (offset_x * 0.5 + round_x);
            f.p1 += self.weight * (offset_y * 0.5 + round_y);
        } else if offset_x.abs() >= offset_y.abs() {
            if offset_x >= 0.0 {
                f.p0 += self.weight * (offset_x * 0.5 + round_x + 0.25);
                f.p1 += self.weight * (offset_y * 0.5 + round_y + 0.25 * offset_y / offset_x);
            } else {
                f.p0 += self.weight * (offset_x * 0.5 + round_x - 0.25);
                f.p1 += self.weight * (offset_y * 0.5 + round_y - 0.25 * offset_y / offset_x);
            }
        } else if offset_y >= 0.0 {
            f.p1 += self.weight * (offset_y * 0.5 + round_y + 0.25);
            f.p0 += self.weight * (offset_x * 0.5 + round_x + offset_x / offset_y * 0.25);
        } else {
            f.p1 += self.weight * (offset_y * 0.5 + round_y - 0.25);
            f.p0 += self.weight * (offset_x * 0.5 + round_x - offset_x / offset_y * 0.25);
        }
    }
}

impl Flam3Variation for variations::Butterfly {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /* wx is self.weight *4/sqrt(3*pi) */
        let wx = self.weight * *BUTTERFLY_WEIGHT;

        let y2 = f.ty * 2.0;
        let r = wx * sqrt!((f.ty * f.tx).abs() / (EPS + f.tx * f.tx + y2 * y2));

        f.p0 += r * f.tx;
        f.p1 += r * y2;
    }
}

impl Flam3Variation for variations::Cell {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let inv_cell_size = 1.0 / self.size;

        /* calculate input cell */
        let mut x = (f.tx * inv_cell_size).floor();
        let mut y = (f.ty * inv_cell_size).floor();

        /* Offset from cell origin */
        let dx = f.tx - x * self.size;
        let dy = f.ty - y * self.size;

        /* interleave cells */
        if y >= 0.0 {
            if x >= 0.0 {
                y *= 2.0;
                x *= 2.0;
            } else {
                y *= 2.0;
                x = -(2.0 * x + 1.0);
            }
        } else if x >= 0.0 {
            y = -(2.0 * y + 1.0);
            x *= 2.0;
        } else {
            y = -(2.0 * y + 1.0);
            x = -(2.0 * x + 1.0);
        }

        f.p0 += self.weight * (dx + x * self.size);
        f.p1 -= self.weight * (dy + y * self.size);
    }
}

impl Flam3Variation for variations::Cpow {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let a = f.atanyx();
        let lnr = 0.5 * ln!(f.sumsq());
        let va = 2.0 * M_PI / self.power;
        let vc = self.r / self.power;
        let vd = self.i / self.power;
        let ang = vc * a + vd * lnr + va * (self.power * f.rc.isaac_01()).trunc();

        let m = self.weight * exp!(vc * lnr - vd * a);

        let (sa, ca) = sincos!(ang);

        f.p0 += m * ca;
        f.p1 += m * sa;
    }
}

impl Flam3Variation for variations::Curve {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let mut pc_xlen = self.x_length * self.x_length;
        let mut pc_ylen = self.y_length * self.y_length;

        if pc_xlen < 1E-20 {
            pc_xlen = 1E-20;
        }

        if pc_ylen < 1E-20 {
            pc_ylen = 1E-20;
        }

        f.p0 += self.weight * (f.tx + self.x_amp * exp!(-f.ty * f.ty / pc_xlen));
        f.p1 += self.weight * (f.ty + self.y_amp * exp!(-f.tx * f.tx / pc_ylen));
    }
}

impl Flam3Variation for variations::Edisc {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let tmp = f.sumsq() + 1.0;
        let tmp2 = 2.0 * f.tx;
        let r1 = sqrt!(tmp + tmp2);
        let r2 = sqrt!(tmp - tmp2);
        let xmax = (r1 + r2) * 0.5;
        let a1 = ln!(xmax + sqrt!(xmax - 1.0));
        let a2 = -acos!(f.tx / xmax);
        let w = self.weight / 11.57034632;

        let (mut snv, csv) = sincos!(a1);

        let snhu = sinh!(a2);
        let cshu = cosh!(a2);

        if f.ty > 0.0 {
            snv = -snv;
        }

        f.p0 += w * cshu * csv;
        f.p1 += w * snhu * snv;
    }
}

impl Flam3Variation for variations::Elliptic {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let tmp = f.sumsq() + 1.0;
        let x2 = 2.0 * f.tx;
        let xmax = 0.5 * (sqrt!(tmp + x2) + sqrt!(tmp - x2));
        let a = f.tx / xmax;
        let mut b = 1.0 - sqr!(a);
        let mut ssx = xmax - 1.0;
        let w = self.weight / M_PI_2;

        if b < 0.0 {
            b = 0.0;
        } else {
            b = sqrt!(b);
        }

        if ssx < 0.0 {
            ssx = 0.0;
        } else {
            ssx = sqrt!(ssx);
        }

        f.p0 += w * atan2!(a, b);

        if f.ty > 0.0 {
            f.p1 += w * ln!(xmax + ssx);
        } else {
            f.p1 -= w * ln!(xmax + ssx);
        }
    }
}

impl Flam3Variation for variations::Escher {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let a = f.atanyx();
        let lnr = 0.5 * ln!(f.sumsq());

        let (seb, ceb) = sincos!(self.beta);

        let vc = 0.5 * (1.0 + ceb);
        let vd = 0.5 * seb;

        let m = self.weight * exp!(vc * lnr - vd * a);
        let n = vc * a + vd * lnr;

        let (sn, cn) = sincos!(n);

        f.p0 += m * cn;
        f.p1 += m * sn;
    }
}

impl Flam3Variation for variations::Foci {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let expx = exp!(f.tx) * 0.5;
        let expnx = 0.25 / expx;

        let (sn, cn) = sincos!(f.ty);
        let tmp = self.weight / (expx + expnx - cn);

        f.p0 += tmp * (expx - expnx);
        f.p1 += tmp * sn;
    }
}

impl Flam3Variation for variations::Lazysusan {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let x = f.tx - self.x;
        let y = f.ty + self.y;
        let mut r = sqrt!(sqr!(x) + sqr!(y));

        if r < self.weight {
            let a = atan2!(y, x) + self.spin + self.twist * (self.weight - r);
            let (sina, cosa) = sincos!(a);
            r *= self.weight;

            f.p0 += r * cosa + self.x;
            f.p1 += r * sina - self.y;
        } else {
            r = self.weight * (1.0 + self.space / r);

            f.p0 += r * x + self.x;
            f.p1 += r * y - self.y;
        }
    }
}

impl Flam3Variation for variations::Loonie {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let r2 = f.sumsq();
        let w2 = self.weight * self.weight;

        if r2 < w2 {
            let r = self.weight * sqrt!(w2 / r2 - 1.0);
            f.p0 += r * f.tx;
            f.p1 += r * f.ty;
        } else {
            f.p0 += self.weight * f.tx;
            f.p1 += self.weight * f.ty;
        }
    }
}

impl Flam3Variation for variations::PreBlur {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /* Get pseudo-gaussian */
        let rnd_g = self.weight
            * (f.rc.isaac_01() + f.rc.isaac_01() + f.rc.isaac_01() + f.rc.isaac_01() - 2.0);
        let rnd_a = f.rc.isaac_01() * 2.0 * M_PI;

        let (sin_a, cos_a) = sincos!(rnd_a);

        /* Note: original coordinate changed */
        f.tx += rnd_g * cos_a;
        f.ty += rnd_g * sin_a;
    }

    fn is_pre(&self) -> bool {
        true
    }
}

impl Flam3Variation for variations::Modulus {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let xr = 2.0 * self.x;
        let yr = 2.0 * self.y;

        if f.tx > self.x {
            f.p0 += self.weight * (-self.x + (f.tx + self.x) % xr);
        } else if f.tx < -self.x {
            f.p0 += self.weight * (self.x - (self.x - f.tx) % xr);
        } else {
            f.p0 += self.weight * f.tx;
        }

        if f.ty > self.y {
            f.p1 += self.weight * (-self.y + ((f.ty + self.y) % yr));
        } else if f.ty < -self.y {
            f.p1 += self.weight * (self.y - ((self.y - f.ty) % yr));
        } else {
            f.p1 += self.weight * f.ty;
        }
    }
}

impl Flam3Variation for variations::Oscilloscope {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let tpf = 2.0 * M_PI * self.frequency;

        let t = if self.damping == 0.0 {
            self.amplitude * cos!(tpf * f.tx) + self.separation
        } else {
            self.amplitude * exp!(-f.tx.abs() * self.damping) * cos!(tpf * f.tx) + self.separation
        };

        if f.ty.abs() <= t {
            f.p0 += self.weight * f.tx;
            f.p1 -= self.weight * f.ty;
        } else {
            f.p0 += self.weight * f.tx;
            f.p1 += self.weight * f.ty;
        }
    }
}

impl Flam3Variation for variations::Polar2 {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let p2v = self.weight / M_PI;

        f.p0 += p2v * f.atan();
        f.p1 += p2v / 2.0 * ln!(f.sumsq());
    }
}

impl Flam3Variation for variations::Popcorn2 {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        f.p0 += self.weight * (f.tx + self.x * sin!(tan!(f.ty * self.c)));
        f.p1 += self.weight * (f.ty + self.y * sin!(tan!(f.tx * self.c)));
    }
}

impl Flam3Variation for variations::Scry {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /*
         * Note that scry does not multiply by self.weight , but as the
         * values still approach 0 as the self.weight  approaches 0, it
         * should be ok
         */

        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let t = f.sumsq();
        let r = 1.0 / (f.sqrt() * (t + 1.0 / (self.weight + EPS)));

        f.p0 += f.tx * r;
        f.p1 += f.ty * r;
    }
}

impl Flam3Variation for variations::Separation {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let sx2 = self.x * self.x;
        let sy2 = self.y * self.y;

        if f.tx > 0.0 {
            f.p0 += self.weight * (sqrt!(f.tx * f.tx + sx2) - f.tx * self.x_inside);
        } else {
            f.p0 -= self.weight * (sqrt!(f.tx * f.tx + sx2) + f.tx * self.x_inside);
        }

        if f.ty > 0.0 {
            f.p1 += self.weight * (sqrt!(f.ty * f.ty + sy2) - f.ty * self.y_inside);
        } else {
            f.p1 -= self.weight * (sqrt!(f.ty * f.ty + sy2) + f.ty * self.y_inside);
        }
    }
}

impl Flam3Variation for variations::Split {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        if cos!(f.tx * self.x_size * M_PI) >= 0.0 {
            f.p1 += self.weight * f.ty;
        } else {
            f.p1 -= self.weight * f.ty;
        }

        if cos!(f.ty * self.y_size * M_PI) >= 0.0 {
            f.p0 += self.weight * f.tx;
        } else {
            f.p0 -= self.weight * f.tx;
        }
    }
}

impl Flam3Variation for variations::Splits {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        if f.tx >= 0.0 {
            f.p0 += self.weight * (f.tx + self.x);
        } else {
            f.p0 += self.weight * (f.tx - self.x);
        }

        if f.ty >= 0.0 {
            f.p1 += self.weight * (f.ty + self.y);
        } else {
            f.p1 += self.weight * (f.ty - self.y);
        }
    }
}

impl Flam3Variation for variations::Stripes {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let roundx = (f.tx + 0.5).floor();
        let offsetx = f.tx - roundx;

        f.p0 += self.weight * (offsetx * (1.0 - self.space) + roundx);
        f.p1 += self.weight * (f.ty + offsetx * offsetx * self.warp);
    }
}

impl Flam3Variation for variations::Wedge {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let mut r = f.sqrt();
        let mut a = f.atanyx() + self.swirl * r;
        let c = ((self.count * a + M_PI) * M_1_PI * 0.5).floor();

        let comp_fac = 1.0 - self.angle * self.count * M_1_PI * 0.5;

        a = a * comp_fac + c * self.angle;

        let (sa, ca) = sincos!(a);

        r = self.weight * (r + self.hole);

        f.p0 += r * ca;
        f.p1 += r * sa;
    }
}

impl Flam3Variation for variations::WedgeJulia {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        precalc: &mut VariationPrecalculations,
    ) {
        let r = self.weight * pow!(f.sumsq(), precalc.wedge_julia_cn.unwrap());
        let t_rnd = (precalc.wedge_julia_r_n.unwrap() * f.rc.isaac_01()).trunc();
        let mut a = (f.atanyx() + 2.0 * M_PI * t_rnd) / self.power;
        let c = ((self.count * a + M_PI) * M_1_PI * 0.5).floor();

        a = a * precalc.wedge_julia_cf.unwrap() + c * self.angle;

        let (sa, ca) = sincos!(a);

        f.p0 += r * ca;
        f.p1 += r * sa;
    }
}

impl Flam3Variation for variations::WedgeSph {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let mut r = 1.0 / (f.sqrt() + EPS);
        let mut a = f.atanyx() + self.swirl * r;
        let c = ((self.count * a + M_PI) * M_1_PI * 0.5).floor();

        let comp_fac = 1.0 - self.angle * self.count * M_1_PI * 0.5;

        a = a * comp_fac + c * self.angle;

        let (sa, ca) = sincos!(a);
        r = self.weight * (r + self.hole);

        f.p0 += r * ca;
        f.p1 += r * sa;
    }
}

impl Flam3Variation for variations::Whorl {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation weight in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let r = f.sqrt();

        let a = if r < self.weight {
            f.atanyx() + self.inside / (self.weight - r)
        } else {
            f.atanyx() + self.outside / (self.weight - r)
        };

        let (sa, ca) = sincos!(a);

        f.p0 += self.weight * r * ca;
        f.p1 += self.weight * r * sa;
    }
}

impl Flam3Variation for variations::Waves2 {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        f.p0 += self.weight * (f.tx + self.x_scale * sin!(f.ty * self.x_frequency));
        f.p1 += self.weight * (f.ty + self.y_scale * sin!(f.tx * self.y_frequency));
    }
}

impl Flam3Variation for variations::Exp {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let expe = exp!(f.tx);
        let (expsin, expcos) = sincos!(f.ty);
        f.p0 += self.weight * expe * expcos;
        f.p1 += self.weight * expe * expsin;
    }
}

impl Flam3Variation for variations::Log {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        f.p0 += self.weight * 0.5 * ln!(f.sumsq());
        f.p1 += self.weight * f.atanyx();
    }
}

impl Flam3Variation for variations::Sin {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (sinsin, sinacos) = sincos!(f.tx);
        let sinsinh = sinh!(f.ty);
        let sincosh = cosh!(f.ty);
        f.p0 += self.weight * sinsin * sincosh;
        f.p1 += self.weight * sinacos * sinsinh;
    }
}

impl Flam3Variation for variations::Cos {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (cossin, coscos) = sincos!(f.tx);
        let cossinh = sinh!(f.ty);
        let coscosh = cosh!(f.ty);
        f.p0 += self.weight * coscos * coscosh;
        f.p1 -= self.weight * cossin * cossinh;
    }
}

impl Flam3Variation for variations::Tan {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (tansin, tancos) = sincos!(2.0 * f.tx);
        let tansinh = sinh!(2.0 * f.ty);
        let tancosh = cosh!(2.0 * f.ty);
        let tanden = 1.0 / (tancos + tancosh);
        f.p0 += self.weight * tanden * tansin;
        f.p1 += self.weight * tanden * tansinh;
    }
}

impl Flam3Variation for variations::Sec {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (secsin, seccos) = sincos!(f.tx);
        let secsinh = sinh!(f.ty);
        let seccosh = cosh!(f.ty);
        let secden = 2.0 / (cos!(2.0 * f.tx) + cosh!(2.0 * f.ty));
        f.p0 += self.weight * secden * seccos * seccosh;
        f.p1 += self.weight * secden * secsin * secsinh;
    }
}

impl Flam3Variation for variations::Csc {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (cscsin, csccos) = sincos!(f.tx);
        let cscsinh = sinh!(f.ty);
        let csccosh = cosh!(f.ty);
        let cscden = 2.0 / (cosh!(2.0 * f.ty) - cos!(2.0 * f.tx));
        f.p0 += self.weight * cscden * cscsin * csccosh;
        f.p1 -= self.weight * cscden * csccos * cscsinh;
    }
}

impl Flam3Variation for variations::Cot {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (cotsin, cotcos) = sincos!(2.0 * f.tx);
        let cotsinh = sinh!(2.0 * f.ty);
        let cotcosh = cosh!(2.0 * f.ty);
        let cotden = 1.0 / (cotcosh - cotcos);
        f.p0 += self.weight * cotden * cotsin;
        f.p1 += self.weight * cotden * -1.0 * cotsinh;
    }
}

impl Flam3Variation for variations::Sinh {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (sinhsin, sinhcos) = sincos!(f.ty);
        let sinhsinh = sinh!(f.tx);
        let sinhcosh = cosh!(f.tx);
        f.p0 += self.weight * sinhsinh * sinhcos;
        f.p1 += self.weight * sinhcosh * sinhsin;
    }
}

impl Flam3Variation for variations::Cosh {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (coshsin, coshcos) = sincos!(f.ty);
        let coshsinh = sinh!(f.tx);
        let coshcosh = cosh!(f.tx);
        f.p0 += self.weight * coshcosh * coshcos;
        f.p1 += self.weight * coshsinh * coshsin;
    }
}

impl Flam3Variation for variations::Tanh {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (tanhsin, tanhcos) = sincos!(2.0 * f.ty);
        let tanhsinh = sinh!(2.0 * f.tx);
        let tanhcosh = cosh!(2.0 * f.tx);
        let tanhden = 1.0 / (tanhcos + tanhcosh);
        f.p0 += self.weight * tanhden * tanhsinh;
        f.p1 += self.weight * tanhden * tanhsin;
    }
}

impl Flam3Variation for variations::Sech {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (sechsin, sechcos) = sincos!(f.ty);
        let sechsinh = sinh!(f.tx);
        let sechcosh = cosh!(f.tx);
        let sechden = 2.0 / (cos!(2.0 * f.ty) + cosh!(2.0 * f.tx));
        f.p0 += self.weight * sechden * sechcos * sechcosh;
        f.p1 -= self.weight * sechden * sechsin * sechsinh;
    }
}

impl Flam3Variation for variations::Csch {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (cschsin, cschcos) = sincos!(f.ty);
        let cschsinh = sinh!(f.tx);
        let cschcosh = cosh!(f.tx);
        let cschden = 2.0 / (cosh!(2.0 * f.tx) - cos!(2.0 * f.ty));
        f.p0 += self.weight * cschden * cschsinh * cschcos;
        f.p1 -= self.weight * cschden * cschcosh * cschsin;
    }
}

impl Flam3Variation for variations::Coth {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let (cothsin, cothcos) = sincos!(2.0 * f.ty);
        let cothsinh = sinh!(2.0 * f.tx);
        let cothcosh = cosh!(2.0 * f.tx);
        let cothden = 1.0 / (cothcosh - cothcos);
        f.p0 += self.weight * cothden * cothsinh;
        f.p1 += self.weight * cothden * cothsin;
    }
}

impl Flam3Variation for variations::Auger {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let s = sin!(self.frequency * f.tx);
        let t = sin!(self.frequency * f.ty);
        let dy = f.ty + self.strength * (self.scale * s / 2.0 + f.ty.abs() * s);
        let dx = f.tx + self.strength * (self.scale * t / 2.0 + f.tx.abs() * t);

        f.p0 += self.weight * (f.tx + self.symmetry * (dx - f.tx));
        f.p1 += self.weight * dy;
    }
}

impl Flam3Variation for variations::Flux {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let xpw = f.tx + self.weight;
        let xmw = f.tx - self.weight;
        let avgr = self.weight
            * (2.0 + self.spread)
            * sqrt!(sqrt!(f.ty * f.ty + xpw * xpw) / sqrt!(f.ty * f.ty + xmw * xmw));
        let avga = (atan2!(f.ty, xmw) - atan2!(f.ty, xpw)) * 0.5;

        f.p0 += avgr * cos!(avga);
        f.p1 += avgr * sin!(avga);
    }
}

impl Flam3Variation for variations::Mobius {
    fn apply(
        &self,
        f: &mut Flam3IterHelper,
        _coeffs: &Affine,
        _precalc: &mut VariationPrecalculations,
    ) {
        let re_u = self.re_a * f.tx - self.im_a * f.ty + self.re_b;
        let im_u = self.re_a * f.ty + self.im_a * f.tx + self.im_b;
        let re_v = self.re_c * f.tx - self.im_c * f.ty + self.re_d;
        let im_v = self.re_c * f.ty + self.im_c * f.tx + self.im_d;

        let rad_v = self.weight / (re_v * re_v + im_v * im_v);

        f.p0 += rad_v * (re_u * re_v + im_u * im_v);
        f.p1 += rad_v * (im_u * re_v - re_u * im_v);
    }
}

pub fn apply_xform(
    xform: &Transform,
    p: &[f64; 4],
    q: &mut [f64; 4],
    precalc: &mut VariationPrecalculations,
    rc: &mut Flam3Rng,
) -> bool {
    let s1 = xform.color_speed;

    q[2] = s1 * xform.color + (1.0 - s1) * p[2];
    q[3] = adjust_percentage(xform.opacity);

    let transformed = xform.coefficients.transform(p.into());
    let mut f = Flam3IterHelper::new(transformed, rc);

    /* Pre-xforms go here, and modify the f.tx and f.ty values */
    for var in xform.variations.iter() {
        if var.is_pre() {
            var.apply(&mut f, &xform.coefficients, precalc);
        }
    }

    for var in xform.variations.iter() {
        if !var.is_pre() {
            var.apply(&mut f, &xform.coefficients, precalc);
        }
    }

    /* apply the post transform */
    (q[0], q[1]) = if !xform.post.is_identity() {
        xform.post.transform((&[f.p0, f.p1]).into()).into()
    } else {
        (f.p0, f.p1)
    };

    /* Check for badvalues and return randoms if bad */
    if badvalue(q[0]) || badvalue(q[1]) {
        q[0] = rc.isaac_11();
        q[1] = rc.isaac_11();
        false
    } else {
        true
    }
}
