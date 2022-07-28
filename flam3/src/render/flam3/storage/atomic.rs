use std::{
    panic::resume_unwind,
    sync::{
        atomic::{AtomicU32, Ordering},
        Arc,
    },
    thread,
};

use atomic_float::{AtomicF32, AtomicF64};

use crate::{
    rect::Rect,
    render::{
        flam3::{
            rng::IsaacRng,
            thread::{
                de_thread, empty_de_thread, iter_thread, DensityEstimationThreadStorage,
                IterationThreadStorage,
            },
            variations::VariationPrecalculations,
            DensityEstimationContext, IterationContext, TransformSelector,
        },
        ThreadingMode,
    },
    utils::PanicCast,
    Transform,
};

use super::{Accumulator, Bucket, Field, RenderStorage};

trait AtomicField: Field {
    fn increase(&self, delta: f64);
}

impl AtomicField for AtomicU32 {
    fn increase(&self, delta: f64) {
        #![allow(unused_must_use)]
        self.fetch_update(Ordering::SeqCst, Ordering::SeqCst, |v| {
            Some(v.saturating_add(delta as u32))
        });
    }
}

impl Field for AtomicU32 {
    fn as_f64(&self) -> f64 {
        self.load(Ordering::Relaxed) as f64
    }
}

impl AtomicField for AtomicF32 {
    fn increase(&self, delta: f64) {
        #![allow(unused_must_use)]
        self.fetch_update(Ordering::SeqCst, Ordering::SeqCst, |v| {
            Some(v + delta as f32)
        });
    }
}

impl Field for AtomicF32 {
    fn as_f64(&self) -> f64 {
        self.load(Ordering::Relaxed) as f64
    }
}

impl AtomicField for AtomicF64 {
    fn increase(&self, delta: f64) {
        #![allow(unused_must_use)]
        self.fetch_update(Ordering::SeqCst, Ordering::SeqCst, |v| Some(v + delta));
    }
}

impl Field for AtomicF64 {
    fn as_f64(&self) -> f64 {
        self.load(Ordering::Relaxed)
    }
}

pub(crate) struct AtomicRenderStorage<BField, AField> {
    buckets: Arc<Rect<Bucket<BField>>>,
    accumulators: Arc<Rect<Accumulator<AField>>>,
}

impl<BField, AField> AtomicRenderStorage<BField, AField> {
    fn de_storage(&self) -> AtomicDensityEstimationThreadStorage<BField, AField> {
        AtomicDensityEstimationThreadStorage {
            buckets: self.buckets.clone(),
            accumulators: self.accumulators.clone(),
        }
    }
}

impl<BField, AField> RenderStorage for AtomicRenderStorage<BField, AField>
where
    AField: AtomicField + Sync + Send + 'static,
    BField: AtomicField + Sync + Send + 'static,
{
    fn threading_model() -> ThreadingMode {
        ThreadingMode::Atomic
    }

    fn new(width: usize, height: usize) -> Self {
        AtomicRenderStorage {
            buckets: Arc::new(Rect::default(width, height)),
            accumulators: Arc::new(Rect::default(width, height)),
        }
    }

    fn height(&self) -> u32 {
        self.accumulators.height().u32()
    }

    fn reset_buckets(&mut self) {
        self.buckets = Arc::new(Rect::default(
            self.accumulators.width(),
            self.accumulators.height(),
        ));
    }

    fn run_iteration_threads(
        &mut self,
        transforms: TransformSelector,
        final_transform: Option<(Transform, VariationPrecalculations)>,
        context: &IterationContext,
        rng: &mut IsaacRng,
        thread_count: usize,
    ) -> Result<(), String> {
        let mut handles = Vec::new();
        handles.reserve(thread_count);

        for _i in 0..thread_count {
            let mut storage = AtomicIterationThreadStorage {
                buckets: self.buckets.clone(),
            };

            let context = context.clone();
            let transforms = transforms.clone();
            let final_transform = final_transform.clone();
            let mut rng = IsaacRng::from_rng(rng);

            handles.push(thread::spawn(move || {
                iter_thread(
                    &context,
                    &transforms,
                    &final_transform,
                    &mut rng,
                    &mut storage,
                )
            }));
        }

        for handle in handles {
            match handle.join() {
                Ok(r) => r?,
                Err(e) => resume_unwind(e),
            }
        }

        Ok(())
    }

    fn run_empty_de_thread(&mut self, k1: f64, k2: f64) {
        empty_de_thread(k1, k2, &mut self.de_storage());
    }

    fn run_all_de_threads(&mut self, contexts: Vec<DensityEstimationContext>) {
        let mut handles = Vec::new();
        handles.reserve(contexts.len());

        for context in contexts {
            let mut storage = self.de_storage();
            handles.push(thread::spawn(move || de_thread(context, &mut storage)));
        }

        for handle in handles {
            if let Err(e) = handle.join() {
                resume_unwind(e);
            }
        }
    }

    fn accumulators(self) -> Rect<Accumulator<f64>> {
        assert_eq!(Arc::strong_count(&self.accumulators), 1);
        self.accumulators.map(|a| a.map(|f| f.as_f64()))
    }
}

struct AtomicIterationThreadStorage<BField> {
    buckets: Arc<Rect<Bucket<BField>>>,
}

impl<BField> IterationThreadStorage for AtomicIterationThreadStorage<BField>
where
    BField: AtomicField,
{
    fn increase_bucket(&mut self, x: usize, y: usize, color: palette::Srgba<f64>, logvis: f64) {
        let bucket = &self.buckets[(x, y)];
        bucket.red.increase(color.red * logvis);
        bucket.green.increase(color.green * logvis);
        bucket.blue.increase(color.blue * logvis);
        bucket.alpha.increase(color.alpha * logvis);
        bucket.density.increase(logvis);
    }
}

struct AtomicDensityEstimationThreadStorage<BField, AField> {
    buckets: Arc<Rect<Bucket<BField>>>,
    accumulators: Arc<Rect<Accumulator<AField>>>,
}

impl<BField, AField> DensityEstimationThreadStorage
    for AtomicDensityEstimationThreadStorage<BField, AField>
where
    AField: AtomicField,
    BField: AtomicField,
{
    fn width(&self) -> usize {
        self.accumulators.width()
    }

    fn height(&self) -> usize {
        self.accumulators.height()
    }

    fn density(&self, x: usize, y: usize) -> f64 {
        self.buckets[(x, y)].density.as_f64()
    }

    fn accumulate(&self, x: usize, y: usize) -> Accumulator<f64> {
        let bucket = &self.buckets[(x, y)];
        Accumulator {
            red: bucket.red.as_f64(),
            green: bucket.green.as_f64(),
            blue: bucket.blue.as_f64(),
            alpha: bucket.alpha.as_f64(),
        }
    }

    fn increase_accumulator(&mut self, x: usize, y: usize, pixel: &Accumulator<f64>) {
        let accumulator = &self.accumulators[(x, y)];
        accumulator.red.increase(pixel.red);
        accumulator.green.increase(pixel.green);
        accumulator.blue.increase(pixel.blue);
        accumulator.alpha.increase(pixel.alpha);
    }
}

pub(crate) type RenderStorageAtomicInt = AtomicRenderStorage<AtomicU32, AtomicU32>;
pub(crate) type RenderStorageAtomicFloat = AtomicRenderStorage<AtomicU32, AtomicF32>;
pub(crate) type RenderStorageAtomicDouble = AtomicRenderStorage<AtomicF64, AtomicF64>;
