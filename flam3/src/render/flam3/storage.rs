use std::{
    cmp::max,
    ops::{Index, IndexMut, Mul, MulAssign},
};

use palette::Srgba;

use crate::{rect::Rect, utils::PanicCast, Genome};

use super::{
    filters::DensityEstimatorFilters,
    rng::IsaacRng,
    thread::{de_thread, empty_de_thread, iter_thread, DensityEstimationStorage, IterationStorage},
    variations::VariationPrecalculations,
    DensityEstimationContext, IterationContext, TransformSelector,
};

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

pub(super) trait RenderStorage {
    type AccumulatorField: PanicCast + Default + Copy;

    fn new(width: usize, height: usize) -> Self;
    fn into_accumulator(value: f64) -> Self::AccumulatorField;

    fn width(&self) -> usize;
    fn height(&self) -> usize;
    fn run_iteration_threads(
        &mut self,
        cp: &Genome,
        context: &IterationContext,
        rng: &mut IsaacRng,
        thread_count: usize,
    ) -> Result<(), String>;
    fn run_de_threads(
        &mut self,
        de_filters: DensityEstimatorFilters,
        k1: f64,
        k2: f64,
        supersample: u32,
        thread_count: usize,
    );

    fn accumulators(self) -> Rect<Accumulator<Self::AccumulatorField>>;
}

pub(super) struct RenderStorageAtomicFloat {
    width: usize,
    height: usize,
    buckets: Rect<Bucket<u32>>,
    accumulators: Rect<Accumulator<f32>>,
}

impl RenderStorageAtomicFloat {
    fn increase_bucket_field(field: &mut u32, delta: f64) {
        *field = field.saturating_add(delta as u32)
    }

    fn increase_accumulator_field(field: &mut f32, delta: f64) {
        *field += delta.f32();
    }
}

impl RenderStorage for RenderStorageAtomicFloat {
    type AccumulatorField = f32;

    fn new(width: usize, height: usize) -> Self {
        RenderStorageAtomicFloat {
            width,
            height,
            buckets: Rect::rectangle(width, height),
            accumulators: Rect::rectangle(width, height),
        }
    }

    fn into_accumulator(value: f64) -> Self::AccumulatorField {
        value.f32()
    }

    fn width(&self) -> usize {
        self.width
    }

    fn height(&self) -> usize {
        self.height
    }

    fn run_iteration_threads(
        &mut self,
        cp: &Genome,
        context: &IterationContext,
        rng: &mut IsaacRng,
        thread_count: usize,
    ) -> Result<(), String> {
        let xform_distrib = TransformSelector::new(cp)?;
        let final_xform = cp
            .final_transform
            .as_ref()
            .map(|xf| (xf.clone(), VariationPrecalculations::new(xf)));

        for _i in 0..thread_count {
            iter_thread(
                context,
                &xform_distrib,
                &final_xform,
                &mut IsaacRng::from_rng(rng),
                self,
            )?;
        }

        Ok(())
    }

    fn run_de_threads(
        &mut self,
        de_filters: DensityEstimatorFilters,
        k1: f64,
        k2: f64,
        supersample: u32,
        thread_count: usize,
    ) {
        if de_filters.filters.is_empty() {
            empty_de_thread(k1, k2, self);
        } else {
            let myspan = self.height.u32() - 2 * (supersample - 1) + 1;
            let thread_count = max(thread_count, myspan.usize());
            let swath = myspan / thread_count.u32();

            for i in 0..thread_count {
                let start_row = i.u32() * swath;
                let end_row = if i == thread_count - 1 {
                    self.height.u32() //myspan.i32();
                } else {
                    ((i.u32() + 1) * swath).u32()
                };

                //  Set up the contents of the helper structure
                let context = DensityEstimationContext {
                    supersample,
                    de_filters: de_filters.clone(),
                    k1,
                    k2,
                    curve: de_filters.curve,
                    start_row,
                    end_row,
                };

                de_thread(context, self);
            }
        }
    }

    fn accumulators(self) -> Rect<Accumulator<f32>> {
        self.accumulators
    }
}

impl IterationStorage for RenderStorageAtomicFloat {
    fn increase_bucket(&mut self, x: usize, y: usize, color: Srgba<f64>, logvis: f64) {
        Self::increase_bucket_field(&mut self.buckets[(x, y)].red, color.red * logvis);
        Self::increase_bucket_field(&mut self.buckets[(x, y)].green, color.green * logvis);
        Self::increase_bucket_field(&mut self.buckets[(x, y)].blue, color.blue * logvis);
        Self::increase_bucket_field(&mut self.buckets[(x, y)].alpha, color.alpha * logvis);
        Self::increase_bucket_field(&mut self.buckets[(x, y)].density, logvis);
    }
}

impl DensityEstimationStorage for RenderStorageAtomicFloat {
    fn width(&self) -> usize {
        self.width
    }

    fn height(&self) -> usize {
        self.height
    }

    fn density(&self, x: usize, y: usize) -> f64 {
        self.buckets[(x, y)].density.f64()
    }

    fn accumulate(&self, x: usize, y: usize) -> Accumulator<f64> {
        self.buckets[(x, y)].accumulator()
    }

    fn increase_accumulator(&mut self, x: usize, y: usize, pixel_delta: &Accumulator<f64>) {
        Self::increase_accumulator_field(&mut self.accumulators[(x, y)].red, pixel_delta.red);
        Self::increase_accumulator_field(&mut self.accumulators[(x, y)].green, pixel_delta.green);
        Self::increase_accumulator_field(&mut self.accumulators[(x, y)].blue, pixel_delta.blue);
        Self::increase_accumulator_field(&mut self.accumulators[(x, y)].alpha, pixel_delta.alpha);
    }
}
