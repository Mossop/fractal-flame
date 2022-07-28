use std::f64::consts::{FRAC_1_PI, FRAC_2_PI, FRAC_PI_2, FRAC_PI_4, PI};

use lazy_static::lazy_static;

use crate::math::{
    acos, atan2, cos, cosh, exp, ln, log10, pow, sin, sincos, sinh, sqr, sqrt, sum_sqr, tan,
};
use crate::{
    utils::PanicCast,
    variations::{self, Var, Variation},
    with_var, Affine, Coordinate, Transform,
};

use super::rng::Flam3Rng;
use super::{adjust_percentage, rng::IsaacRng};

const EPS: f64 = 1e-10;

lazy_static! {
    static ref BUTTERFLY_WEIGHT: f64 = 4.0 / sqrt!(3.0 * PI);
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

#[derive(Default, Clone)]
pub(crate) struct VariationPrecalculations {
    pub persp_vsin: f64,
    pub persp_vfcos: f64,
    pub julian_r_n: f64,
    pub julian_cn: f64,
    pub wedge_julia_cf: f64,
    pub wedge_julia_r_n: f64,
    pub wedge_julia_cn: f64,
    pub juliascope_r_n: f64,
    pub juliascope_cn: f64,
    pub radial_blur_spinvar: f64,
    pub radial_blur_zoomvar: f64,
    pub waves_dx2: f64,
    pub waves_dy2: f64,
    pub disc2_timespi: f64,
    pub disc2_cosadd: f64,
    pub disc2_sinadd: f64,
    pub super_shape_pm_4: f64,
    pub super_shape_pneg1_n1: f64,
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
                Var::Waves(_) => precalcs.waves_precalc(transform),
                Var::Disc2(ref v) => precalcs.disc2_precalc(v),
                Var::SuperShape(ref v) => precalcs.super_shape_precalc(v),
                _ => {}
            }
        }

        precalcs
    }

    fn perspective_precalc(&mut self, var: &variations::Perspective) {
        let ang = var.angle * PI / 2.0;
        self.persp_vsin = sin!(ang);
        self.persp_vfcos = var.distance * cos!(ang);
    }

    fn julian_precalc(&mut self, var: &variations::Julian) {
        self.julian_r_n = var.power.abs();
        self.julian_cn = var.distance / var.power / 2.0;
    }

    fn wedge_julia_precalc(&mut self, var: &variations::WedgeJulia) {
        self.wedge_julia_cf = 1.0 - var.angle * var.count * FRAC_1_PI * 0.5;
        self.wedge_julia_r_n = var.power.abs();
        self.wedge_julia_cn = var.dist / var.power / 2.0;
    }

    fn juliascope_precalc(&mut self, var: &variations::Juliascope) {
        self.juliascope_r_n = var.power.abs();
        self.juliascope_cn = var.distance / var.power / 2.0;
    }

    fn radial_blur_precalc(&mut self, var: &variations::RadialBlur) {
        let (spinvar, zoomvar) = sincos!(var.angle * PI / 2.0);
        self.radial_blur_spinvar = spinvar;
        self.radial_blur_zoomvar = zoomvar;
    }

    fn waves_precalc(&mut self, transform: &Transform) {
        let dx = transform.coefficients[2][0];
        let dy = transform.coefficients[2][1];

        self.waves_dx2 = 1.0 / (sqr!(dx) + EPS);
        self.waves_dy2 = 1.0 / (sqr!(dy) + EPS);
    }

    fn disc2_precalc(&mut self, var: &variations::Disc2) {
        let add = var.twist;

        self.disc2_timespi = var.rotate * PI;

        let (mut sinadd, mut cosadd) = sincos!(add);
        cosadd -= 1.0;

        if add > 2.0 * PI {
            let k = 1.0 + add - 2.0 * PI;
            cosadd *= k;
            sinadd *= k;
        }

        if add < -2.0 * PI {
            let k = 1.0 + add + 2.0 * PI;
            cosadd *= k;
            sinadd *= k;
        }

        self.disc2_cosadd = cosadd;
        self.disc2_sinadd = sinadd;
    }

    fn super_shape_precalc(&mut self, var: &variations::SuperShape) {
        self.super_shape_pm_4 = var.m / 4.0;
        self.super_shape_pneg1_n1 = -1.0 / var.n1;
    }
}

pub(crate) trait Flam3Variation: Variation {
    fn apply(
        &self,
        state: &mut IterationState,
        coeffs: &Affine,
        precalc: &VariationPrecalculations,
    );

    fn is_pre(&self) -> bool {
        false
    }
}

pub(crate) struct IterationState<'a> {
    p0: f64,
    p1: f64,
    tx: f64,
    ty: f64,

    rc: &'a mut IsaacRng,

    sumsq: Option<f64>,
    sqrt: Option<f64>,
    atan: Option<f64>,
    sina: Option<f64>,
    cosa: Option<f64>,
    atanyx: Option<f64>,
}

impl<'a> IterationState<'a> {
    fn new(coords: Coordinate<f64>, rc: &'a mut IsaacRng) -> Self {
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
        upsert!(self.sumsq, sum_sqr!(self.tx, self.ty))
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
        state: &mut IterationState,
        coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        with_var!(self, v, v.apply(state, coeffs, precalc));
    }

    fn is_pre(&self) -> bool {
        with_var!(self, v, v.is_pre())
    }
}

impl Flam3Variation for variations::Linear {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        state.p0 += self.weight * state.tx;
        state.p1 += self.weight * state.ty;
    }
}

impl Flam3Variation for variations::Sinusoidal {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        state.p0 += self.weight * sin!(state.tx);
        state.p1 += self.weight * sin!(state.ty);
    }
}

impl Flam3Variation for variations::Spherical {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r2 = self.weight / (state.sumsq() + EPS);

        state.p0 += r2 * state.tx;
        state.p1 += r2 * state.ty;
    }
}

impl Flam3Variation for variations::Swirl {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r2 = state.sumsq();

        let (c1, c2) = sincos!(r2);
        let nx = c1 * state.tx - c2 * state.ty;
        let ny = c2 * state.tx + c1 * state.ty;

        state.p0 += self.weight * nx;
        state.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Horseshoe {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = self.weight / (state.sqrt() + EPS);

        state.p0 += (state.tx - state.ty) * (state.tx + state.ty) * r;
        state.p1 += 2.0 * state.tx * state.ty * r;
    }
}

impl Flam3Variation for variations::Polar {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let nx = state.atan() * FRAC_1_PI;
        let ny = state.sqrt() - 1.0;

        state.p0 += self.weight * nx;
        state.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Handkerchief {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let a = state.atan();
        let r = state.sqrt();

        state.p0 += self.weight * r * sin!(a + r);
        state.p1 += self.weight * r * cos!(a - r);
    }
}

impl Flam3Variation for variations::Heart {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let a = state.sqrt() * state.atan();
        let r = self.weight * state.sqrt();

        let (sa, ca) = sincos!(a);

        state.p0 += r * sa;
        state.p1 += (-r) * ca;
    }
}

impl Flam3Variation for variations::Disc {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let a = state.atan() * FRAC_1_PI;
        let r = PI * state.sqrt();
        let (sr, cr) = sincos!(r);

        state.p0 += self.weight * sr * a;
        state.p1 += self.weight * cr * a;
    }
}

impl Flam3Variation for variations::Spiral {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = state.sqrt() + EPS;
        let r1 = self.weight / r;
        let (sr, cr) = sincos!(r);

        state.p0 += r1 * (state.cosa() + sr);
        state.p1 += r1 * (state.sina() - cr);
    }
}

impl Flam3Variation for variations::Hyperbolic {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = state.sqrt() + EPS;

        state.p0 += self.weight * state.sina() / r;
        state.p1 += self.weight * state.cosa() * r;
    }
}

impl Flam3Variation for variations::Diamond {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = state.sqrt();
        let (sr, cr) = sincos!(r);

        state.p0 += self.weight * state.sina() * cr;
        state.p1 += self.weight * state.cosa() * sr;
    }
}

impl Flam3Variation for variations::Ex {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let a = state.atan();
        let r = state.sqrt();

        let n0 = sin!(a + r);
        let n1 = cos!(a - r);

        let m0 = n0 * n0 * n0 * r;
        let m1 = n1 * n1 * n1 * r;

        state.p0 += self.weight * (m0 + m1);
        state.p1 += self.weight * (m0 - m1);
    }
}

impl Flam3Variation for variations::Julia {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let mut a = 0.5 * state.atan();

        if state.rc.next_bit() {
            a += PI;
        }

        let r = self.weight * sqrt!(state.sqrt());

        let (sa, ca) = sincos!(a);

        state.p0 += r * ca;
        state.p1 += r * sa;
    }
}

impl Flam3Variation for variations::Bent {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let mut nx = state.tx;
        let mut ny = state.ty;

        if nx < 0.0 {
            nx *= 2.0;
        }
        if ny < 0.0 {
            ny /= 2.0;
        }

        state.p0 += self.weight * nx;
        state.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Waves {
    fn apply(
        &self,
        state: &mut IterationState,
        coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        let c10 = coeffs[1][0];
        let c11 = coeffs[1][1];

        let nx = state.tx + c10 * sin!(state.ty * precalc.waves_dx2);
        let ny = state.ty + c11 * sin!(state.tx * precalc.waves_dy2);

        state.p0 += self.weight * nx;
        state.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Fisheye {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = 2.0 * self.weight / (state.sqrt() + 1.0);

        state.p0 += r * state.ty;
        state.p1 += r * state.tx;
    }
}

impl Flam3Variation for variations::Popcorn {
    fn apply(
        &self,
        state: &mut IterationState,
        coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let dx = tan!(3.0 * state.ty);
        let dy = tan!(3.0 * state.tx);

        let nx = state.tx + coeffs[2][0] * sin!(dx);
        let ny = state.ty + coeffs[2][1] * sin!(dy);

        state.p0 += self.weight * nx;
        state.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Exponential {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let dx = self.weight * exp!(state.tx - 1.0);
        let dy = PI * state.ty;

        let (sdy, cdy) = sincos!(dy);

        state.p0 += dx * cdy;
        state.p1 += dx * sdy;
    }
}

impl Flam3Variation for variations::Power {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = self.weight * pow!(state.sqrt(), state.sina());

        state.p0 += r * state.cosa();
        state.p1 += r * state.sina();
    }
}

impl Flam3Variation for variations::Cosine {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let a = state.tx * PI;

        let (sa, ca) = sincos!(a);
        let nx = ca * cosh!(state.ty);
        let ny = -sa * sinh!(state.ty);

        state.p0 += self.weight * nx;
        state.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Rings {
    fn apply(
        &self,
        state: &mut IterationState,
        coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let dx = coeffs[2][0] * coeffs[2][0] + EPS;
        let mut r = state.sqrt();
        r = self.weight * (((r + dx) % (2.0 * dx)) - dx + r * (1.0 - dx));

        state.p0 += r * state.cosa();
        state.p1 += r * state.sina();
    }
}

impl Flam3Variation for variations::Fan {
    fn apply(
        &self,
        state: &mut IterationState,
        coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let dx = PI * (coeffs[2][0] * coeffs[2][0] + EPS);
        let dy = coeffs[2][1];
        let dx2 = 0.5 * dx;

        let mut a = state.atan();
        let r = self.weight * state.sqrt();

        a += if ((a + dy) % dx) > dx2 { -dx2 } else { dx2 };
        let (sa, ca) = sincos!(a);

        state.p0 += r * ca;
        state.p1 += r * sa;
    }
}

impl Flam3Variation for variations::Blob {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let mut r = state.sqrt();
        let a = state.atan();
        let bdiff = self.high - self.low;

        r *= self.low + bdiff * (0.5 + 0.5 * sin!(self.waves * a));

        state.p0 += self.weight * state.sina() * r;
        state.p1 += self.weight * state.cosa() * r;
    }
}

impl Flam3Variation for variations::Pdj {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let nx1 = cos!(self.b * state.tx);
        let nx2 = sin!(self.c * state.tx);
        let ny1 = sin!(self.a * state.ty);
        let ny2 = cos!(self.d * state.ty);

        state.p0 += self.weight * (ny1 - nx1);
        state.p1 += self.weight * (nx2 - ny2);
    }
}

impl Flam3Variation for variations::Fan2 {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let dy = self.y;
        let dx = PI * (sqr!(self.x) + EPS);
        let dx2 = 0.5 * dx;
        let mut a = state.atan();
        let r = self.weight * state.sqrt();

        let t = a + dy - dx * ((a + dy) / dx).trunc();

        if t > dx2 {
            a -= dx2;
        } else {
            a += dx2;
        }

        let (sa, ca) = sincos!(a);

        state.p0 += r * sa;
        state.p1 += r * ca;
    }
}

impl Flam3Variation for variations::Rings2 {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let mut r = state.sqrt();
        let dx = self.val * self.val + EPS;

        r += -2.0 * dx * ((r + dx) / (2.0 * dx)).trunc() + r * (1.0 - dx);

        state.p0 += self.weight * state.sina() * r;
        state.p1 += self.weight * state.cosa() * r;
    }
}

impl Flam3Variation for variations::Eyefish {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = (self.weight * 2.0) / (state.sqrt() + 1.0);

        state.p0 += r * state.tx;
        state.p1 += r * state.ty;
    }
}

impl Flam3Variation for variations::Bubble {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = self.weight / (0.25 * (state.sumsq()) + 1.0);

        state.p0 += r * state.tx;
        state.p1 += r * state.ty;
    }
}

impl Flam3Variation for variations::Cylinder {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        state.p0 += self.weight * sin!(state.tx);
        state.p1 += self.weight * state.ty;
    }
}

impl Flam3Variation for variations::Perspective {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        let t = 1.0 / (self.distance - state.ty * precalc.persp_vsin);

        state.p0 += self.weight * self.distance * state.tx * t;
        state.p1 += self.weight * precalc.persp_vfcos * state.ty * t;
    }
}

impl Flam3Variation for variations::Noise {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let tmpr = state.rc.next_01() * 2.0 * PI;
        let (sinr, cosr) = sincos!(tmpr);

        let r = self.weight * state.rc.next_01();

        state.p0 += state.tx * r * cosr;
        state.p1 += state.ty * r * sinr;
    }
}

impl Flam3Variation for variations::Julian {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        let t_rnd = (precalc.julian_r_n * state.rc.next_01()).trunc();

        let tmpr = (state.atanyx() + 2.0 * PI * t_rnd) / self.power;

        let r = self.weight * pow!(state.sumsq(), precalc.julian_cn);
        let (sina, cosa) = sincos!(tmpr);

        state.p0 += r * cosa;
        state.p1 += r * sina;
    }
}

impl Flam3Variation for variations::Juliascope {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        let t_rnd = (precalc.juliascope_r_n * state.rc.next_01()).trunc();

        let tmpr = if ((t_rnd.i32()) & 1) == 0 {
            (2.0 * PI * t_rnd + state.atanyx()) / self.power
        } else {
            (2.0 * PI * t_rnd - state.atanyx()) / self.power
        };

        let (sina, cosa) = sincos!(tmpr);

        let r = self.weight * pow!(state.sumsq(), precalc.juliascope_cn);

        state.p0 += r * cosa;
        state.p1 += r * sina;
    }
}

impl Flam3Variation for variations::Blur {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let tmpr = state.rc.next_01() * 2.0 * PI;
        let (sinr, cosr) = sincos!(tmpr);

        let r = self.weight * state.rc.next_01();

        state.p0 += r * cosr;
        state.p1 += r * sinr;
    }
}

impl Flam3Variation for variations::GaussianBlur {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let ang = state.rc.next_01() * 2.0 * PI;
        let (sina, cosa) = sincos!(ang);

        let r = self.weight
            * (state.rc.next_01() + state.rc.next_01() + state.rc.next_01() + state.rc.next_01()
                - 2.0);

        state.p0 += r * cosa;
        state.p1 += r * sina;
    }
}

impl Flam3Variation for variations::RadialBlur {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        /* Get pseudo-gaussian */
        let rnd = self.weight
            * (state.rc.next_01() + state.rc.next_01() + state.rc.next_01() + state.rc.next_01()
                - 2.0);

        /* Calculate angle & zoom */
        let ra = state.sqrt();
        let tmpa = state.atanyx() + precalc.radial_blur_spinvar * rnd;
        let (sa, ca) = sincos!(tmpa);
        let rz = precalc.radial_blur_zoomvar * rnd - 1.0;

        state.p0 += ra * ca + rz * state.tx;
        state.p1 += ra * sa + rz * state.ty;
    }
}

impl Flam3Variation for variations::Pie {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let sl = (state.rc.next_01() * self.slices + 0.5).trunc();
        let a = self.rotation + 2.0 * PI * (sl + state.rc.next_01() * self.thickness) / self.slices;
        let r = self.weight * state.rc.next_01();
        let (sa, ca) = sincos!(a);

        state.p0 += r * ca;
        state.p1 += r * sa;
    }
}

impl Flam3Variation for variations::Ngon {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r_factor = pow!(state.sumsq(), self.power / 2.0);

        let theta = state.atanyx();
        let b = 2.0 * PI / self.sides;

        let mut phi = theta - (b * (theta / b).floor());
        if phi > b / 2.0 {
            phi -= b;
        }

        let mut amp = self.corners * (1.0 / (cos!(phi) + EPS) - 1.0) + self.circle;
        amp /= r_factor + EPS;

        state.p0 += self.weight * state.tx * amp;
        state.p1 += self.weight * state.ty * amp;
    }
}

impl Flam3Variation for variations::Curl {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let re = 1.0 + self.c1 * state.tx + self.c2 * (state.tx * state.tx - state.ty * state.ty);
        let im = self.c1 * state.ty + 2.0 * self.c2 * state.tx * state.ty;

        let r = self.weight / sum_sqr!(re, im);

        state.p0 += (state.tx * re + state.ty * im) * r;
        state.p1 += (state.ty * re - state.tx * im) * r;
    }
}

impl Flam3Variation for variations::Rectangles {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        if self.x == 0.0 {
            state.p0 += self.weight * state.tx;
        } else {
            state.p0 +=
                self.weight * ((2.0 * (state.tx / self.x).floor() + 1.0) * self.x - state.tx);
        }

        if self.y == 0.0 {
            state.p1 += self.weight * state.ty;
        } else {
            state.p1 +=
                self.weight * ((2.0 * (state.ty / self.y).floor() + 1.0) * self.y - state.ty);
        }
    }
}

impl Flam3Variation for variations::Arch {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let ang = state.rc.next_01() * self.weight * PI;
        let (sinr, cosr) = sincos!(ang);

        state.p0 += self.weight * sinr;
        state.p1 += self.weight * sqr!(sinr) / cosr;
    }
}

impl Flam3Variation for variations::Tangent {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        state.p0 += self.weight * sin!(state.tx) / cos!(state.ty);
        state.p1 += self.weight * tan!(state.ty);
    }
}

impl Flam3Variation for variations::Square {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        state.p0 += self.weight * (state.rc.next_01() - 0.5);
        state.p1 += self.weight * (state.rc.next_01() - 0.5);
    }
}

impl Flam3Variation for variations::Rays {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let ang = self.weight * state.rc.next_01() * PI;
        let r = self.weight / (state.sumsq() + EPS);
        let tanr = self.weight * tan!(ang) * r;

        state.p0 += tanr * cos!(state.tx);
        state.p1 += tanr * sin!(state.ty);
    }
}

impl Flam3Variation for variations::Blade {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let r = state.rc.next_01() * self.weight * state.sqrt();

        let (sinr, cosr) = sincos!(r);

        state.p0 += self.weight * state.tx * (cosr + sinr);
        state.p1 += self.weight * state.tx * (cosr - sinr);
    }
}

impl Flam3Variation for variations::Secant2 {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /* Intended as a 'fixed' version of secant */

        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let r = self.weight * state.sqrt();
        let cr = cos!(r);
        let icr = 1.0 / cr;

        state.p0 += self.weight * state.tx;

        if cr < 0.0 {
            state.p1 += self.weight * (icr + 1.0);
        } else {
            state.p1 += self.weight * (icr - 1.0);
        }
    }
}

impl Flam3Variation for variations::Twintrian {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */
        let r = state.rc.next_01() * self.weight * state.sqrt();

        let (sinr, cosr) = sincos!(r);
        let mut diff = log10!(sqr!(sinr)) + cosr;

        if badvalue(diff) {
            diff = -30.0;
        }

        state.p0 += self.weight * state.tx * diff;
        state.p1 += self.weight * state.tx * (diff - sinr * PI);
    }
}

impl Flam3Variation for variations::Cross {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let s = state.tx * state.tx - state.ty * state.ty;
        let r = self.weight * sqrt!(1.0 / (sqr!(s) + EPS));

        state.p0 += state.tx * r;
        state.p1 += state.ty * r;
    }
}

impl Flam3Variation for variations::Disc2 {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        let t = precalc.disc2_timespi * (state.tx + state.ty);
        let (sinr, cosr) = sincos!(t);
        let r = self.weight * state.atan() / PI;

        state.p0 += (sinr + precalc.disc2_cosadd) * r;
        state.p1 += (cosr + precalc.disc2_sinadd) * r;
    }
}

impl Flam3Variation for variations::SuperShape {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        let theta = precalc.super_shape_pm_4 * state.atanyx() + FRAC_PI_4;

        let (st, ct) = sincos!(theta);

        let t1 = pow!(ct.abs(), self.n2);
        let t2 = pow!(st.abs(), self.n3);

        let myrnd = self.rnd;

        let r = self.weight
            * ((myrnd * state.rc.next_01() + (1.0 - myrnd) * state.sqrt()) - self.holes)
            * pow!(t1 + t2, precalc.super_shape_pneg1_n1)
            / state.sqrt();

        state.p0 += r * state.tx;
        state.p1 += r * state.ty;
    }
}

impl Flam3Variation for variations::Flower {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let theta = state.atanyx();
        let r = self.weight * (state.rc.next_01() - self.holes) * cos!(self.petals * theta)
            / state.sqrt();

        state.p0 += r * state.tx;
        state.p1 += r * state.ty;
    }
}

impl Flam3Variation for variations::Conic {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let ct = state.tx / state.sqrt();
        let r = self.weight * (state.rc.next_01() - self.holes) * self.eccentricity
            / (1.0 + self.eccentricity * ct)
            / state.sqrt();

        state.p0 += r * state.tx;
        state.p1 += r * state.ty;
    }
}

impl Flam3Variation for variations::Parabola {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let r = state.sqrt();

        let (sr, cr) = sincos!(r);

        state.p0 += self.height * self.weight * sqr!(sr) * state.rc.next_01();
        state.p1 += self.width * self.weight * cr * state.rc.next_01();
    }
}

impl Flam3Variation for variations::Bent2 {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let mut nx = state.tx;
        let mut ny = state.ty;

        if nx < 0.0 {
            nx *= self.x;
        }
        if ny < 0.0 {
            ny *= self.y;
        }

        state.p0 += self.weight * nx;
        state.p1 += self.weight * ny;
    }
}

impl Flam3Variation for variations::Bipolar {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let x2y2 = state.sumsq();
        let t = x2y2 + 1.0;
        let x2 = 2.0 * state.tx;
        let ps = -FRAC_PI_2 * self.shift;
        let mut y = 0.5 * atan2!(2.0 * state.ty, x2y2 - 1.0) + ps;

        if y > FRAC_PI_2 {
            y = -FRAC_PI_2 + ((y + FRAC_PI_2) % PI);
        } else if y < -FRAC_PI_2 {
            y = FRAC_PI_2 - ((FRAC_PI_2 - y) % PI);
        }

        state.p0 += self.weight * 0.25 * FRAC_2_PI * ln!((t + x2) / (t - x2));
        state.p1 += self.weight * FRAC_2_PI * y;
    }
}

impl Flam3Variation for variations::Boarders {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let round_x = state.tx.round();
        let round_y = state.ty.round();
        let offset_x = state.tx - round_x;
        let offset_y = state.ty - round_y;

        if state.rc.next_01() >= 0.75 {
            state.p0 += self.weight * (offset_x * 0.5 + round_x);
            state.p1 += self.weight * (offset_y * 0.5 + round_y);
        } else if offset_x.abs() >= offset_y.abs() {
            if offset_x >= 0.0 {
                state.p0 += self.weight * (offset_x * 0.5 + round_x + 0.25);
                state.p1 += self.weight * (offset_y * 0.5 + round_y + 0.25 * offset_y / offset_x);
            } else {
                state.p0 += self.weight * (offset_x * 0.5 + round_x - 0.25);
                state.p1 += self.weight * (offset_y * 0.5 + round_y - 0.25 * offset_y / offset_x);
            }
        } else if offset_y >= 0.0 {
            state.p1 += self.weight * (offset_y * 0.5 + round_y + 0.25);
            state.p0 += self.weight * (offset_x * 0.5 + round_x + offset_x / offset_y * 0.25);
        } else {
            state.p1 += self.weight * (offset_y * 0.5 + round_y - 0.25);
            state.p0 += self.weight * (offset_x * 0.5 + round_x - offset_x / offset_y * 0.25);
        }
    }
}

impl Flam3Variation for variations::Butterfly {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /* wx is self.weight *4/sqrt(3*pi) */
        let wx = self.weight * *BUTTERFLY_WEIGHT;

        let y2 = state.ty * 2.0;
        let r = wx * sqrt!((state.ty * state.tx).abs() / (EPS + sum_sqr!(state.tx, y2)));

        state.p0 += r * state.tx;
        state.p1 += r * y2;
    }
}

impl Flam3Variation for variations::Cell {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let inv_cell_size = 1.0 / self.size;

        /* calculate input cell */
        let mut x = (state.tx * inv_cell_size).floor().i32();
        let mut y = (state.ty * inv_cell_size).floor().i32();

        /* Offset from cell origin */
        let dx = state.tx - x.f64() * self.size;
        let dy = state.ty - y.f64() * self.size;

        /* interleave cells */
        if y >= 0 {
            if x >= 0 {
                y *= 2;
                x *= 2;
            } else {
                y *= 2;
                x = -(2 * x + 1);
            }
        } else if x >= 0 {
            y = -(2 * y + 1);
            x *= 2;
        } else {
            y = -(2 * y + 1);
            x = -(2 * x + 1);
        }

        state.p0 += self.weight * (dx + x.f64() * self.size);
        state.p1 -= self.weight * (dy + y.f64() * self.size);
    }
}

impl Flam3Variation for variations::Cpow {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let a = state.atanyx();
        let lnr = 0.5 * ln!(state.sumsq());
        let va = 2.0 * PI / self.power;
        let vc = self.r / self.power;
        let vd = self.i / self.power;
        let ang = vc * a + vd * lnr + va * (self.power * state.rc.next_01()).trunc();

        let m = self.weight * exp!(vc * lnr - vd * a);

        let (sa, ca) = sincos!(ang);

        state.p0 += m * ca;
        state.p1 += m * sa;
    }
}

impl Flam3Variation for variations::Curve {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let mut pc_xlen = self.x_length * self.x_length;
        let mut pc_ylen = self.y_length * self.y_length;

        if pc_xlen < 1E-20 {
            pc_xlen = 1E-20;
        }

        if pc_ylen < 1E-20 {
            pc_ylen = 1E-20;
        }

        state.p0 += self.weight * (state.tx + self.x_amp * exp!(-state.ty * state.ty / pc_xlen));
        state.p1 += self.weight * (state.ty + self.y_amp * exp!(-state.tx * state.tx / pc_ylen));
    }
}

impl Flam3Variation for variations::Edisc {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let tmp = state.sumsq() + 1.0;
        let tmp2 = 2.0 * state.tx;
        let r1 = sqrt!(tmp + tmp2);
        let r2 = sqrt!(tmp - tmp2);
        let xmax = (r1 + r2) * 0.5;
        let a1 = ln!(xmax + sqrt!(xmax - 1.0));
        let a2 = -acos!(state.tx / xmax);
        let w = self.weight / 11.57034632;

        let (mut snv, csv) = sincos!(a1);

        let snhu = sinh!(a2);
        let cshu = cosh!(a2);

        if state.ty > 0.0 {
            snv = -snv;
        }

        state.p0 += w * cshu * csv;
        state.p1 += w * snhu * snv;
    }
}

impl Flam3Variation for variations::Elliptic {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let tmp = state.sumsq() + 1.0;
        let x2 = 2.0 * state.tx;
        let xmax = 0.5 * (sqrt!(tmp + x2) + sqrt!(tmp - x2));
        let a = state.tx / xmax;
        let mut b = 1.0 - sqr!(a);
        let mut ssx = xmax - 1.0;
        let w = self.weight / FRAC_PI_2;

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

        state.p0 += w * atan2!(a, b);

        if state.ty > 0.0 {
            state.p1 += w * ln!(xmax + ssx);
        } else {
            state.p1 -= w * ln!(xmax + ssx);
        }
    }
}

impl Flam3Variation for variations::Escher {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let a = state.atanyx();
        let lnr = 0.5 * ln!(state.sumsq());

        let (seb, ceb) = sincos!(self.beta);

        let vc = 0.5 * (1.0 + ceb);
        let vd = 0.5 * seb;

        let m = self.weight * exp!(vc * lnr - vd * a);
        let n = vc * a + vd * lnr;

        let (sn, cn) = sincos!(n);

        state.p0 += m * cn;
        state.p1 += m * sn;
    }
}

impl Flam3Variation for variations::Foci {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let expx = exp!(state.tx) * 0.5;
        let expnx = 0.25 / expx;

        let (sn, cn) = sincos!(state.ty);
        let tmp = self.weight / (expx + expnx - cn);

        state.p0 += tmp * (expx - expnx);
        state.p1 += tmp * sn;
    }
}

impl Flam3Variation for variations::Lazysusan {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let x = state.tx - self.x;
        let y = state.ty + self.y;
        let mut r = sqrt!(sum_sqr!(x, y));

        if r < self.weight {
            let a = atan2!(y, x) + self.spin + self.twist * (self.weight - r);
            let (sina, cosa) = sincos!(a);
            r *= self.weight;

            state.p0 += r * cosa + self.x;
            state.p1 += r * sina - self.y;
        } else {
            r = self.weight * (1.0 + self.space / r);

            state.p0 += r * x + self.x;
            state.p1 += r * y - self.y;
        }
    }
}

impl Flam3Variation for variations::Loonie {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation self.weight  in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let r2 = state.sumsq();
        let w2 = self.weight * self.weight;

        if r2 < w2 {
            let r = self.weight * sqrt!(w2 / r2 - 1.0);
            state.p0 += r * state.tx;
            state.p1 += r * state.ty;
        } else {
            state.p0 += self.weight * state.tx;
            state.p1 += self.weight * state.ty;
        }
    }
}

impl Flam3Variation for variations::PreBlur {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /* Get pseudo-gaussian */
        let rnd_g = self.weight
            * (state.rc.next_01() + state.rc.next_01() + state.rc.next_01() + state.rc.next_01()
                - 2.0);
        let rnd_a = state.rc.next_01() * 2.0 * PI;

        let (sin_a, cos_a) = sincos!(rnd_a);

        /* Note: original coordinate changed */
        state.tx += rnd_g * cos_a;
        state.ty += rnd_g * sin_a;
    }

    fn is_pre(&self) -> bool {
        true
    }
}

impl Flam3Variation for variations::Modulus {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let xr = 2.0 * self.x;
        let yr = 2.0 * self.y;

        if state.tx > self.x {
            state.p0 += self.weight * (-self.x + (state.tx + self.x) % xr);
        } else if state.tx < -self.x {
            state.p0 += self.weight * (self.x - (self.x - state.tx) % xr);
        } else {
            state.p0 += self.weight * state.tx;
        }

        if state.ty > self.y {
            state.p1 += self.weight * (-self.y + ((state.ty + self.y) % yr));
        } else if state.ty < -self.y {
            state.p1 += self.weight * (self.y - ((self.y - state.ty) % yr));
        } else {
            state.p1 += self.weight * state.ty;
        }
    }
}

impl Flam3Variation for variations::Oscilloscope {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let tpf = 2.0 * PI * self.frequency;

        let t = if self.damping == 0.0 {
            self.amplitude * cos!(tpf * state.tx) + self.separation
        } else {
            self.amplitude * exp!(-state.tx.abs() * self.damping) * cos!(tpf * state.tx)
                + self.separation
        };

        if state.ty.abs() <= t {
            state.p0 += self.weight * state.tx;
            state.p1 -= self.weight * state.ty;
        } else {
            state.p0 += self.weight * state.tx;
            state.p1 += self.weight * state.ty;
        }
    }
}

impl Flam3Variation for variations::Polar2 {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let p2v = self.weight / PI;

        state.p0 += p2v * state.atan();
        state.p1 += p2v / 2.0 * ln!(state.sumsq());
    }
}

impl Flam3Variation for variations::Popcorn2 {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        state.p0 += self.weight * (state.tx + self.x * sin!(tan!(state.ty * self.c)));
        state.p1 += self.weight * (state.ty + self.y * sin!(tan!(state.tx * self.c)));
    }
}

impl Flam3Variation for variations::Scry {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
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

        let t = state.sumsq();
        let r = 1.0 / (state.sqrt() * (t + 1.0 / (self.weight + EPS)));

        state.p0 += state.tx * r;
        state.p1 += state.ty * r;
    }
}

impl Flam3Variation for variations::Separation {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let sx2 = self.x * self.x;
        let sy2 = self.y * self.y;

        if state.tx > 0.0 {
            state.p0 += self.weight * (sqrt!(state.tx * state.tx + sx2) - state.tx * self.x_inside);
        } else {
            state.p0 -= self.weight * (sqrt!(state.tx * state.tx + sx2) + state.tx * self.x_inside);
        }

        if state.ty > 0.0 {
            state.p1 += self.weight * (sqrt!(state.ty * state.ty + sy2) - state.ty * self.y_inside);
        } else {
            state.p1 -= self.weight * (sqrt!(state.ty * state.ty + sy2) + state.ty * self.y_inside);
        }
    }
}

impl Flam3Variation for variations::Split {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        if cos!(state.tx * self.x_size * PI) >= 0.0 {
            state.p1 += self.weight * state.ty;
        } else {
            state.p1 -= self.weight * state.ty;
        }

        if cos!(state.ty * self.y_size * PI) >= 0.0 {
            state.p0 += self.weight * state.tx;
        } else {
            state.p0 -= self.weight * state.tx;
        }
    }
}

impl Flam3Variation for variations::Splits {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        if state.tx >= 0.0 {
            state.p0 += self.weight * (state.tx + self.x);
        } else {
            state.p0 += self.weight * (state.tx - self.x);
        }

        if state.ty >= 0.0 {
            state.p1 += self.weight * (state.ty + self.y);
        } else {
            state.p1 += self.weight * (state.ty - self.y);
        }
    }
}

impl Flam3Variation for variations::Stripes {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let roundx = (state.tx + 0.5).floor();
        let offsetx = state.tx - roundx;

        state.p0 += self.weight * (offsetx * (1.0 - self.space) + roundx);
        state.p1 += self.weight * (state.ty + sqr!(offsetx) * self.warp);
    }
}

impl Flam3Variation for variations::Wedge {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let mut r = state.sqrt();
        let mut a = state.atanyx() + self.swirl * r;
        let c = ((self.count * a + PI) * FRAC_1_PI * 0.5).floor();

        let comp_fac = 1.0 - self.angle * self.count * FRAC_1_PI * 0.5;

        a = a * comp_fac + c * self.angle;

        let (sa, ca) = sincos!(a);

        r = self.weight * (r + self.hole);

        state.p0 += r * ca;
        state.p1 += r * sa;
    }
}

impl Flam3Variation for variations::WedgeJulia {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        precalc: &VariationPrecalculations,
    ) {
        let r = self.weight * pow!(state.sumsq(), precalc.wedge_julia_cn);
        let t_rnd = (precalc.wedge_julia_r_n * state.rc.next_01()).trunc();
        let mut a = (state.atanyx() + 2.0 * PI * t_rnd) / self.power;
        let c = ((self.count * a + PI) * FRAC_1_PI * 0.5).floor();

        a = a * precalc.wedge_julia_cf + c * self.angle;

        let (sa, ca) = sincos!(a);

        state.p0 += r * ca;
        state.p1 += r * sa;
    }
}

impl Flam3Variation for variations::WedgeSph {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let mut r = 1.0 / (state.sqrt() + EPS);
        let mut a = state.atanyx() + self.swirl * r;
        let c = ((self.count * a + PI) * FRAC_1_PI * 0.5).floor();

        let comp_fac = 1.0 - self.angle * self.count * FRAC_1_PI * 0.5;

        a = a * comp_fac + c * self.angle;

        let (sa, ca) = sincos!(a);
        r = self.weight * (r + self.hole);

        state.p0 += r * ca;
        state.p1 += r * sa;
    }
}

impl Flam3Variation for variations::Whorl {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        /*
         * !!! Note !!!
         * This code uses the variation weight in a non-standard fashion, and
         * it may change or even be removed in future versions of flam3.
         */

        let r = state.sqrt();

        let a = if r < self.weight {
            state.atanyx() + self.inside / (self.weight - r)
        } else {
            state.atanyx() + self.outside / (self.weight - r)
        };

        let (sa, ca) = sincos!(a);

        state.p0 += self.weight * r * ca;
        state.p1 += self.weight * r * sa;
    }
}

impl Flam3Variation for variations::Waves2 {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        state.p0 += self.weight * (state.tx + self.x_scale * sin!(state.ty * self.x_frequency));
        state.p1 += self.weight * (state.ty + self.y_scale * sin!(state.tx * self.y_frequency));
    }
}

impl Flam3Variation for variations::Exp {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let expe = exp!(state.tx);
        let (expsin, expcos) = sincos!(state.ty);
        state.p0 += self.weight * expe * expcos;
        state.p1 += self.weight * expe * expsin;
    }
}

impl Flam3Variation for variations::Log {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        state.p0 += self.weight * 0.5 * ln!(state.sumsq());
        state.p1 += self.weight * state.atanyx();
    }
}

impl Flam3Variation for variations::Sin {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (sinsin, sinacos) = sincos!(state.tx);
        let sinsinh = sinh!(state.ty);
        let sincosh = cosh!(state.ty);
        state.p0 += self.weight * sinsin * sincosh;
        state.p1 += self.weight * sinacos * sinsinh;
    }
}

impl Flam3Variation for variations::Cos {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (cossin, coscos) = sincos!(state.tx);
        let cossinh = sinh!(state.ty);
        let coscosh = cosh!(state.ty);
        state.p0 += self.weight * coscos * coscosh;
        state.p1 -= self.weight * cossin * cossinh;
    }
}

impl Flam3Variation for variations::Tan {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (tansin, tancos) = sincos!(2.0 * state.tx);
        let tansinh = sinh!(2.0 * state.ty);
        let tancosh = cosh!(2.0 * state.ty);
        let tanden = 1.0 / (tancos + tancosh);
        state.p0 += self.weight * tanden * tansin;
        state.p1 += self.weight * tanden * tansinh;
    }
}

impl Flam3Variation for variations::Sec {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (secsin, seccos) = sincos!(state.tx);
        let secsinh = sinh!(state.ty);
        let seccosh = cosh!(state.ty);
        let secden = 2.0 / (cos!(2.0 * state.tx) + cosh!(2.0 * state.ty));
        state.p0 += self.weight * secden * seccos * seccosh;
        state.p1 += self.weight * secden * secsin * secsinh;
    }
}

impl Flam3Variation for variations::Csc {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (cscsin, csccos) = sincos!(state.tx);
        let cscsinh = sinh!(state.ty);
        let csccosh = cosh!(state.ty);
        let cscden = 2.0 / (cosh!(2.0 * state.ty) - cos!(2.0 * state.tx));
        state.p0 += self.weight * cscden * cscsin * csccosh;
        state.p1 -= self.weight * cscden * csccos * cscsinh;
    }
}

impl Flam3Variation for variations::Cot {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (cotsin, cotcos) = sincos!(2.0 * state.tx);
        let cotsinh = sinh!(2.0 * state.ty);
        let cotcosh = cosh!(2.0 * state.ty);
        let cotden = 1.0 / (cotcosh - cotcos);
        state.p0 += self.weight * cotden * cotsin;
        state.p1 += self.weight * cotden * -1.0 * cotsinh;
    }
}

impl Flam3Variation for variations::Sinh {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (sinhsin, sinhcos) = sincos!(state.ty);
        let sinhsinh = sinh!(state.tx);
        let sinhcosh = cosh!(state.tx);
        state.p0 += self.weight * sinhsinh * sinhcos;
        state.p1 += self.weight * sinhcosh * sinhsin;
    }
}

impl Flam3Variation for variations::Cosh {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (coshsin, coshcos) = sincos!(state.ty);
        let coshsinh = sinh!(state.tx);
        let coshcosh = cosh!(state.tx);
        state.p0 += self.weight * coshcosh * coshcos;
        state.p1 += self.weight * coshsinh * coshsin;
    }
}

impl Flam3Variation for variations::Tanh {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (tanhsin, tanhcos) = sincos!(2.0 * state.ty);
        let tanhsinh = sinh!(2.0 * state.tx);
        let tanhcosh = cosh!(2.0 * state.tx);
        let tanhden = 1.0 / (tanhcos + tanhcosh);
        state.p0 += self.weight * tanhden * tanhsinh;
        state.p1 += self.weight * tanhden * tanhsin;
    }
}

impl Flam3Variation for variations::Sech {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (sechsin, sechcos) = sincos!(state.ty);
        let sechsinh = sinh!(state.tx);
        let sechcosh = cosh!(state.tx);
        let sechden = 2.0 / (cos!(2.0 * state.ty) + cosh!(2.0 * state.tx));
        state.p0 += self.weight * sechden * sechcos * sechcosh;
        state.p1 -= self.weight * sechden * sechsin * sechsinh;
    }
}

impl Flam3Variation for variations::Csch {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (cschsin, cschcos) = sincos!(state.ty);
        let cschsinh = sinh!(state.tx);
        let cschcosh = cosh!(state.tx);
        let cschden = 2.0 / (cosh!(2.0 * state.tx) - cos!(2.0 * state.ty));
        state.p0 += self.weight * cschden * cschsinh * cschcos;
        state.p1 -= self.weight * cschden * cschcosh * cschsin;
    }
}

impl Flam3Variation for variations::Coth {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let (cothsin, cothcos) = sincos!(2.0 * state.ty);
        let cothsinh = sinh!(2.0 * state.tx);
        let cothcosh = cosh!(2.0 * state.tx);
        let cothden = 1.0 / (cothcosh - cothcos);
        state.p0 += self.weight * cothden * cothsinh;
        state.p1 += self.weight * cothden * cothsin;
    }
}

impl Flam3Variation for variations::Auger {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let s = sin!(self.frequency * state.tx);
        let t = sin!(self.frequency * state.ty);
        let dy = state.ty + self.strength * (self.scale * s / 2.0 + state.ty.abs() * s);
        let dx = state.tx + self.strength * (self.scale * t / 2.0 + state.tx.abs() * t);

        state.p0 += self.weight * (state.tx + self.symmetry * (dx - state.tx));
        state.p1 += self.weight * dy;
    }
}

impl Flam3Variation for variations::Flux {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let xpw = state.tx + self.weight;
        let xmw = state.tx - self.weight;
        let avgr = self.weight
            * (2.0 + self.spread)
            * sqrt!(
                sqrt!(state.ty * state.ty + sqr!(xpw)) / sqrt!(state.ty * state.ty + sqr!(xmw))
            );
        let avga = (atan2!(state.ty, xmw) - atan2!(state.ty, xpw)) * 0.5;

        state.p0 += avgr * cos!(avga);
        state.p1 += avgr * sin!(avga);
    }
}

impl Flam3Variation for variations::Mobius {
    fn apply(
        &self,
        state: &mut IterationState,
        _coeffs: &Affine,
        _precalc: &VariationPrecalculations,
    ) {
        let re_u = self.re_a * state.tx - self.im_a * state.ty + self.re_b;
        let im_u = self.re_a * state.ty + self.im_a * state.tx + self.im_b;
        let re_v = self.re_c * state.tx - self.im_c * state.ty + self.re_d;
        let im_v = self.re_c * state.ty + self.im_c * state.tx + self.im_d;

        let rad_v = self.weight / sum_sqr!(re_v, im_v);

        state.p0 += rad_v * (re_u * re_v + im_u * im_v);
        state.p1 += rad_v * (im_u * re_v - re_u * im_v);
    }
}

pub(crate) fn apply_xform(
    xform: &Transform,
    p: &[f64; 4],
    q: &mut [f64; 4],
    precalc: &VariationPrecalculations,
    rc: &mut IsaacRng,
) -> bool {
    let s1 = xform.color_speed;

    q[2] = s1 * xform.color + (1.0 - s1) * p[2];
    q[3] = adjust_percentage(xform.opacity);

    let transformed = xform.coefficients.transform(p.into());
    let mut state = IterationState::new(transformed, rc);

    /* Pre-xforms go here, and modify the state.tx and state.ty values */
    for var in xform.variations.iter() {
        if var.is_pre() {
            var.apply(&mut state, &xform.coefficients, precalc);
        }
    }

    for var in xform.variations.iter() {
        if !var.is_pre() {
            var.apply(&mut state, &xform.coefficients, precalc);
        }
    }

    /* apply the post transform */
    (q[0], q[1]) = if !xform.post.is_identity() {
        xform.post.transform((&[state.p0, state.p1]).into()).into()
    } else {
        (state.p0, state.p1)
    };

    /* Check for badvalues and return randoms if bad */
    if badvalue(q[0]) || badvalue(q[1]) {
        q[0] = rc.next_11();
        q[1] = rc.next_11();
        false
    } else {
        true
    }
}
