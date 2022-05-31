#[cfg(feature = "libm")]
macro_rules! sin {
    ($e:expr) => {
        libm::sin($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! sin {
    ($e:expr) => {
        $e.sin()
    };
}

#[cfg(feature = "libm")]
macro_rules! sinh {
    ($e:expr) => {
        libm::sinh($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! sinh {
    ($e:expr) => {
        $e.sinh()
    };
}

#[cfg(feature = "libm")]
macro_rules! cos {
    ($e:expr) => {
        libm::cos($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! cos {
    ($e:expr) => {
        $e.cos()
    };
}

#[cfg(feature = "libm")]
macro_rules! acos {
    ($e:expr) => {
        libm::acos($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! acos {
    ($e:expr) => {
        $e.acos()
    };
}

#[cfg(feature = "libm")]
macro_rules! cosh {
    ($e:expr) => {
        libm::cosh($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! cosh {
    ($e:expr) => {
        $e.cosh()
    };
}

#[cfg(feature = "libm")]
macro_rules! tan {
    ($e:expr) => {
        libm::tan($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! tan {
    ($e:expr) => {
        $e.tan()
    };
}

#[cfg(feature = "libm")]
macro_rules! atan2 {
    ($a:expr, $b:expr) => {
        libm::atan2($a, $b)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! atan2 {
    ($a:expr, $b:expr) => {
        $a.atan2($b)
    };
}

#[cfg(feature = "libm")]
macro_rules! sincos {
    ($e:expr) => {
        libm::sincos($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! sincos {
    ($e:expr) => {
        $e.sin_cos()
    };
}

#[cfg(feature = "libm")]
macro_rules! sqrt {
    ($e:expr) => {
        libm::sqrt($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! sqrt {
    ($e:expr) => {
        $e.sqrt()
    };
}

#[cfg(feature = "libm")]
macro_rules! exp {
    ($e:expr) => {
        libm::exp($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! exp {
    ($e:expr) => {
        $e.exp()
    };
}

#[cfg(feature = "libm")]
macro_rules! pow {
    ($a:expr, $b:expr) => {
        libm::pow($a, $b)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! pow {
    ($a:expr, $b:expr) => {
        $a.powf($b)
    };
}

macro_rules! sqr {
    ($a:expr) => {
        $a * $a
    };
}

macro_rules! sum_sqr {
    ($a:expr, $b:expr) => {
        $a * $a + $b * $b
    };
}

#[cfg(feature = "libm")]
macro_rules! log10 {
    ($e:expr) => {
        libm::log10($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! log10 {
    ($e:expr) => {
        $e.log10()
    };
}

#[cfg(feature = "libm")]
macro_rules! ln {
    ($e:expr) => {
        libm::log($e)
    };
}
#[cfg(not(feature = "libm"))]
macro_rules! ln {
    ($e:expr) => {
        $e.ln()
    };
}

pub(crate) use {
    acos, atan2, cos, cosh, exp, ln, log10, pow, sin, sincos, sinh, sqr, sqrt, sum_sqr, tan,
};
