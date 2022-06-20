#[cfg(feature = "libm")]
mod internal {
    macro_rules! sin {
        ($e:expr) => {
            libm::sin($e)
        };
    }

    macro_rules! sinh {
        ($e:expr) => {
            libm::sinh($e)
        };
    }

    macro_rules! cos {
        ($e:expr) => {
            libm::cos($e)
        };
    }

    macro_rules! acos {
        ($e:expr) => {
            libm::acos($e)
        };
    }

    macro_rules! cosh {
        ($e:expr) => {
            libm::cosh($e)
        };
    }

    macro_rules! tan {
        ($e:expr) => {
            libm::tan($e)
        };
    }

    macro_rules! atan2 {
        ($a:expr, $b:expr) => {
            libm::atan2($a, $b)
        };
    }

    macro_rules! sincos {
        ($e:expr) => {
            libm::sincos($e)
        };
    }

    macro_rules! sqrt {
        ($e:expr) => {
            libm::sqrt($e)
        };
    }

    macro_rules! exp {
        ($e:expr) => {
            libm::exp($e)
        };
    }

    macro_rules! pow {
        ($a:expr, $b:expr) => {
            libm::pow($a, $b)
        };
    }

    macro_rules! log10 {
        ($e:expr) => {
            libm::log10($e)
        };
    }

    macro_rules! ln {
        ($e:expr) => {
            libm::log($e)
        };
    }

    pub(crate) use {acos, atan2, cos, cosh, exp, ln, log10, pow, sin, sincos, sinh, sqrt, tan};
}

#[cfg(not(feature = "libm"))]
mod internal {
    macro_rules! sin {
        ($e:expr) => {
            $e.sin()
        };
    }

    macro_rules! sinh {
        ($e:expr) => {
            $e.sinh()
        };
    }

    macro_rules! cos {
        ($e:expr) => {
            $e.cos()
        };
    }

    macro_rules! acos {
        ($e:expr) => {
            $e.acos()
        };
    }

    macro_rules! cosh {
        ($e:expr) => {
            $e.cosh()
        };
    }

    macro_rules! tan {
        ($e:expr) => {
            $e.tan()
        };
    }

    macro_rules! atan2 {
        ($a:expr, $b:expr) => {
            $a.atan2($b)
        };
    }

    macro_rules! sincos {
        ($e:expr) => {
            $e.sin_cos()
        };
    }

    macro_rules! sqrt {
        ($e:expr) => {
            $e.sqrt()
        };
    }

    macro_rules! exp {
        ($e:expr) => {
            $e.exp()
        };
    }

    macro_rules! pow {
        ($a:expr, $b:expr) => {
            $a.powf($b)
        };
    }

    macro_rules! ln {
        ($e:expr) => {
            $e.ln()
        };
    }

    macro_rules! log10 {
        ($e:expr) => {
            $e.log10()
        };
    }

    pub(crate) use {acos, atan2, cos, cosh, exp, ln, log10, pow, sin, sincos, sinh, sqrt, tan};
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

pub(crate) use internal::{
    acos, atan2, cos, cosh, exp, ln, log10, pow, sin, sincos, sinh, sqrt, tan,
};
pub(crate) use {sqr, sum_sqr};
