use std::ops::{Index, IndexMut, Mul, MulAssign};

use palette::Srgba;

use crate::utils::PanicCast;

#[derive(Clone)]
pub(crate) struct Bucket<T> {
    pub red: T,
    pub green: T,
    pub blue: T,
    pub alpha: T,
    pub density: T,
}

impl<T> Default for Bucket<T>
where
    T: Default,
{
    fn default() -> Bucket<T> {
        Bucket {
            red: Default::default(),
            green: Default::default(),
            blue: Default::default(),
            alpha: Default::default(),
            density: Default::default(),
        }
    }
}

impl<T> Bucket<T>
where
    T: Default + Copy + PanicCast,
{
    pub fn accumulator(&self) -> Accumulator<f64> {
        Accumulator {
            red: self.red.f64(),
            green: self.green.f64(),
            blue: self.blue.f64(),
            alpha: self.alpha.f64(),
        }
    }
}

#[derive(Clone)]
pub struct Accumulator<T> {
    pub red: T,
    pub green: T,
    pub blue: T,
    pub alpha: T,
}

impl<T> Default for Accumulator<T>
where
    T: Default,
{
    fn default() -> Accumulator<T> {
        Accumulator {
            red: Default::default(),
            green: Default::default(),
            blue: Default::default(),
            alpha: Default::default(),
        }
    }
}

impl MulAssign<f64> for Accumulator<f64> {
    fn mul_assign(&mut self, rhs: f64) {
        self.red *= rhs;
        self.green *= rhs;
        self.blue *= rhs;
        self.alpha *= rhs;
    }
}

impl Mul<f64> for &Accumulator<f64> {
    type Output = Accumulator<f64>;

    fn mul(self, rhs: f64) -> Accumulator<f64> {
        Accumulator {
            red: self.red * rhs,
            green: self.green * rhs,
            blue: self.blue * rhs,
            alpha: self.alpha * rhs,
        }
    }
}

impl<T> Index<usize> for Accumulator<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.red,
            1 => &self.green,
            2 => &self.blue,
            3 => &self.alpha,
            _ => panic!("Out of bounds"),
        }
    }
}

impl<T> IndexMut<usize> for Accumulator<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.red,
            1 => &mut self.green,
            2 => &mut self.blue,
            3 => &mut self.alpha,
            _ => panic!("Out of bounds"),
        }
    }
}

pub(crate) trait RenderStorage {
    type BucketField: PanicCast + Default + Copy;
    type AccumulatorField: PanicCast + Default + Copy;

    fn increase_bucket_field(bucket: &mut Self::BucketField, delta: f64);
    fn increase_accumulator_field(accumulator: &mut Self::AccumulatorField, delta: f64);
    fn into_accumulator(val: f64) -> Self::AccumulatorField;

    fn increase_bucket(dest: &mut Bucket<Self::BucketField>, color: Srgba<f64>, logvis: f64) {
        Self::increase_bucket_field(&mut dest.red, color.red * logvis);
        Self::increase_bucket_field(&mut dest.green, color.green * logvis);
        Self::increase_bucket_field(&mut dest.blue, color.blue * logvis);
        Self::increase_bucket_field(&mut dest.alpha, color.alpha * logvis);
        Self::increase_bucket_field(&mut dest.density, logvis);
    }

    fn increase_accumulator(
        dest: &mut Accumulator<Self::AccumulatorField>,
        pixel_delta: &Accumulator<f64>,
    ) {
        Self::increase_accumulator_field(&mut dest.red, pixel_delta.red);
        Self::increase_accumulator_field(&mut dest.green, pixel_delta.green);
        Self::increase_accumulator_field(&mut dest.blue, pixel_delta.blue);
        Self::increase_accumulator_field(&mut dest.alpha, pixel_delta.alpha);
    }
}

pub struct RenderStorageAtomicFloat {}

impl RenderStorage for RenderStorageAtomicFloat {
    type BucketField = u32;
    type AccumulatorField = f32;

    fn increase_bucket_field(bucket: &mut Self::BucketField, delta: f64) {
        *bucket = bucket.saturating_add(delta as u32)
    }

    fn increase_accumulator_field(acc: &mut Self::AccumulatorField, delta: f64) {
        *acc += delta.f32();
    }

    fn into_accumulator(val: f64) -> Self::AccumulatorField {
        val.f32()
    }
}
