pub(crate) trait PanicCast {
    fn i32(self) -> i32;
    fn u8(self) -> u8;
    fn u16(self) -> u16;
    fn u32(self) -> u32;
    fn u64(self) -> u64;
    fn usize(self) -> usize;
    fn f32(self) -> f32;
    fn f64(self) -> f64;
}

impl PanicCast for usize {
    fn i32(self) -> i32 {
        self.try_into().unwrap()
    }

    fn u8(self) -> u8 {
        self.try_into().unwrap()
    }

    fn u16(self) -> u16 {
        self.try_into().unwrap()
    }

    fn u32(self) -> u32 {
        self.try_into().unwrap()
    }

    fn u64(self) -> u64 {
        self.try_into().unwrap()
    }

    fn usize(self) -> usize {
        self
    }

    fn f32(self) -> f32 {
        self as f32
    }

    fn f64(self) -> f64 {
        self as f64
    }
}

impl PanicCast for i32 {
    fn i32(self) -> i32 {
        self
    }

    fn u8(self) -> u8 {
        self.try_into().unwrap()
    }

    fn u16(self) -> u16 {
        self.try_into().unwrap()
    }

    fn u32(self) -> u32 {
        self.try_into().unwrap()
    }

    fn u64(self) -> u64 {
        self.try_into().unwrap()
    }

    fn usize(self) -> usize {
        self.try_into().unwrap()
    }

    fn f32(self) -> f32 {
        self as f32
    }

    fn f64(self) -> f64 {
        self.into()
    }
}

impl PanicCast for u32 {
    fn i32(self) -> i32 {
        self.try_into().unwrap()
    }

    fn u8(self) -> u8 {
        self as u8
    }

    fn u16(self) -> u16 {
        self as u16
    }

    fn u32(self) -> u32 {
        self
    }

    fn u64(self) -> u64 {
        self as u64
    }

    fn usize(self) -> usize {
        self.try_into().unwrap()
    }

    fn f32(self) -> f32 {
        self as f32
    }

    fn f64(self) -> f64 {
        self.into()
    }
}

impl PanicCast for f32 {
    fn i32(self) -> i32 {
        if self > i32::MAX as f32 {
            panic!("Out of bounds conversion.");
        }
        if self < i32::MIN as f32 {
            panic!("Out of bounds conversion.");
        }

        self as i32
    }

    fn u8(self) -> u8 {
        self.u64().try_into().unwrap()
    }

    fn u16(self) -> u16 {
        self.u64().try_into().unwrap()
    }

    fn u32(self) -> u32 {
        self.u64().try_into().unwrap()
    }

    fn u64(self) -> u64 {
        if self < u64::MIN as f32 {
            panic!("Out of bounds conversion.");
        }
        if self > u64::MAX as f32 {
            panic!("Out of bounds conversion.");
        }

        self as u64
    }

    fn usize(self) -> usize {
        self.u64().try_into().unwrap()
    }

    fn f32(self) -> f32 {
        self
    }

    fn f64(self) -> f64 {
        self as f64
    }
}

impl PanicCast for f64 {
    fn i32(self) -> i32 {
        if self > i32::MAX as f64 {
            panic!("Out of bounds conversion.");
        }
        if self < i32::MIN as f64 {
            panic!("Out of bounds conversion.");
        }

        self as i32
    }

    fn u8(self) -> u8 {
        self.u64().try_into().unwrap()
    }

    fn u16(self) -> u16 {
        self.u64().try_into().unwrap()
    }

    fn u32(self) -> u32 {
        self.u64().try_into().unwrap()
    }

    fn u64(self) -> u64 {
        if self < u64::MIN as f64 {
            panic!("Out of bounds conversion.");
        }
        if self > u64::MAX as f64 {
            panic!("Out of bounds conversion.");
        }

        self as u64
    }

    fn usize(self) -> usize {
        self.u64().try_into().unwrap()
    }

    fn f32(self) -> f32 {
        self as f32
    }

    fn f64(self) -> f64 {
        self
    }
}

/// Simulates the -freciprocal-math optimisation when enabled.
#[macro_export]
macro_rules! fastdiv {
    ($num:expr, $den:expr) => {
        if cfg!(fastmath) {
            $num * (1.0 / $den)
        } else {
            $num / $den
        }
    };
}
