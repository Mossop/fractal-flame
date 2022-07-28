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

trait SyncField: Field {
    fn increase(&mut self, delta: f64);
}

impl SyncField for u32 {
    fn increase(&mut self, delta: f64) {
        *self = self.saturating_add(delta as u32);
    }
}

impl Field for u32 {
    fn as_f64(&self) -> f64 {
        *self as f64
    }
}

impl SyncField for f32 {
    fn increase(&mut self, delta: f64) {
        *self += delta as f32;
    }
}

impl Field for f32 {
    fn as_f64(&self) -> f64 {
        *self as f64
    }
}

impl SyncField for f64 {
    fn increase(&mut self, delta: f64) {
        *self += delta;
    }
}

impl Field for f64 {
    fn as_f64(&self) -> f64 {
        *self as f64
    }
}

pub(crate) struct SyncRenderStorage<BField, AField> {
    buckets: Rect<Bucket<BField>>,
    accumulators: Rect<Accumulator<AField>>,
}

impl<BField, AField> RenderStorage for SyncRenderStorage<BField, AField>
where
    AField: SyncField,
    BField: SyncField,
{
    fn threading_model() -> ThreadingMode {
        ThreadingMode::Sync
    }

    fn new(width: usize, height: usize) -> Self {
        SyncRenderStorage {
            buckets: Rect::default(width, height),
            accumulators: Rect::default(width, height),
        }
    }

    fn height(&self) -> u32 {
        self.accumulators.height().u32()
    }

    fn reset_buckets(&mut self) {
        self.buckets = Rect::default(self.accumulators.width(), self.accumulators.height());
    }

    fn run_iteration_threads(
        &mut self,
        transforms: TransformSelector,
        final_transform: Option<(Transform, VariationPrecalculations)>,
        context: &IterationContext,
        rng: &mut IsaacRng,
        thread_count: usize,
    ) -> Result<(), String> {
        for _i in 0..thread_count {
            iter_thread(
                context,
                &transforms,
                &final_transform,
                &mut IsaacRng::from_rng(rng),
                self,
            )?;
        }

        Ok(())
    }

    fn run_empty_de_thread(&mut self, k1: f64, k2: f64) {
        empty_de_thread(k1, k2, self);
    }

    fn run_all_de_threads(&mut self, contexts: Vec<DensityEstimationContext>) {
        for context in contexts {
            de_thread(context, self);
        }
    }

    fn accumulators(self) -> Rect<Accumulator<f64>> {
        self.accumulators.map(|a| a.map(|f| f.as_f64()))
    }
}

impl<BField, AField> IterationThreadStorage for SyncRenderStorage<BField, AField>
where
    BField: SyncField,
{
    fn increase_bucket(&mut self, x: usize, y: usize, color: palette::Srgba<f64>, logvis: f64) {
        let bucket = &mut self.buckets[(x, y)];
        bucket.red.increase(color.red * logvis);
        bucket.green.increase(color.green * logvis);
        bucket.blue.increase(color.blue * logvis);
        bucket.alpha.increase(color.alpha * logvis);
        bucket.density.increase(logvis);
    }
}

impl<BField, AField> DensityEstimationThreadStorage for SyncRenderStorage<BField, AField>
where
    AField: SyncField,
    BField: SyncField,
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
        let accumulator = &mut self.accumulators[(x, y)];
        accumulator.red.increase(pixel.red);
        accumulator.green.increase(pixel.green);
        accumulator.blue.increase(pixel.blue);
        accumulator.alpha.increase(pixel.alpha);
    }
}

pub(crate) type RenderStorageSyncInt = SyncRenderStorage<u32, u32>;
pub(crate) type RenderStorageSyncFloat = SyncRenderStorage<u32, f32>;
pub(crate) type RenderStorageSyncDouble = SyncRenderStorage<f64, f64>;
