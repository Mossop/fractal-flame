use std::{
    cmp::min,
    ops::{Index, IndexMut, Mul, MulAssign},
};

use crate::{rect::Rect, render::ThreadingMode, utils::PanicCast, Transform};

use super::{
    filters::DensityEstimatorFilters, rng::IsaacRng, variations::VariationPrecalculations,
    DensityEstimationContext, IterationContext, TransformSelector,
};

pub mod atomic;
pub mod sync;

trait Field: Default {
    fn as_f64(&self) -> f64;
}

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

#[derive(Clone)]
pub struct Accumulator<T> {
    pub red: T,
    pub green: T,
    pub blue: T,
    pub alpha: T,
}

impl<T> Accumulator<T> {
    pub fn map<R, F>(&self, mut mapper: F) -> Accumulator<R>
    where
        F: FnMut(&T) -> R,
    {
        Accumulator {
            red: mapper(&self.red),
            green: mapper(&self.green),
            blue: mapper(&self.blue),
            alpha: mapper(&self.alpha),
        }
    }
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
    fn threading_model() -> ThreadingMode;
    fn new(width: usize, height: usize) -> Self;
    fn height(&self) -> u32;

    fn reset_buckets(&mut self);
    fn run_iteration_threads(
        &mut self,
        transforms: TransformSelector,
        final_transform: Option<(Transform, VariationPrecalculations)>,
        context: &IterationContext,
        rng: &mut IsaacRng,
        thread_count: usize,
    ) -> Result<(), String>;

    fn accumulators(self) -> Rect<Accumulator<f64>>;

    fn run_empty_de_thread(&mut self, k1: f64, k2: f64);
    fn run_all_de_threads(&mut self, contexts: Vec<DensityEstimationContext>);

    fn run_de_threads(
        &mut self,
        de_filters: DensityEstimatorFilters,
        k1: f64,
        k2: f64,
        supersample: u32,
        thread_count: usize,
    ) {
        if de_filters.filters.is_empty() {
            self.run_empty_de_thread(k1, k2);
        } else {
            let myspan = self.height() - 2 * (supersample - 1) + 1;
            let thread_count = min(thread_count, myspan.usize());
            let swath = myspan / thread_count.u32();

            let mut contexts = Vec::new();
            contexts.reserve(thread_count);

            for i in 0..thread_count {
                let start_row = i.u32() * swath;
                let end_row = if i == thread_count - 1 {
                    self.height()
                } else {
                    ((i.u32() + 1) * swath).u32()
                };

                //  Set up the contents of the helper structure
                contexts.push(DensityEstimationContext {
                    supersample,
                    de_filters: de_filters.clone(),
                    k1,
                    k2,
                    curve: de_filters.curve,
                    start_row,
                    end_row,
                });
            }

            self.run_all_de_threads(contexts);
        }
    }
}
