#[macro_export]
macro_rules! sin {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::sin($e)
        } else {
            $e.sin()
        }
    };
}

#[macro_export]
macro_rules! sinh {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::sinh($e)
        } else {
            $e.sinh()
        }
    };
}

#[macro_export]
macro_rules! cos {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::cos($e)
        } else {
            $e.cos()
        }
    };
}

#[macro_export]
macro_rules! acos {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::acos($e)
        } else {
            $e.acos()
        }
    };
}

#[macro_export]
macro_rules! cosh {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::cosh($e)
        } else {
            $e.cosh()
        }
    };
}

#[macro_export]
macro_rules! tan {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::tan($e)
        } else {
            $e.tan()
        }
    };
}

#[macro_export]
macro_rules! atan2 {
    ($a:expr, $b:expr) => {
        if cfg!(feature = "libm") {
            libm::atan2($a, $b)
        } else {
            $a.atan2($b)
        }
    };
}

#[macro_export]
macro_rules! sincos {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::sincos($e)
        } else {
            $e.sin_cos()
        }
    };
}

#[macro_export]
macro_rules! sqrt {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::sqrt($e)
        } else {
            $e.sqrt()
        }
    };
}

#[macro_export]
macro_rules! exp {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::exp($e)
        } else {
            $e.exp()
        }
    };
}

#[macro_export]
macro_rules! pow {
    ($a:expr, $b:expr) => {
        if cfg!(feature = "libm") {
            libm::pow($a, $b)
        } else {
            $a.powf($b)
        }
    };
}

#[macro_export]
macro_rules! sqr {
    ($a:expr) => {
        $a * $a
    };
}

#[macro_export]
macro_rules! log10 {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::log10($e)
        } else {
            $e.log10()
        }
    };
}

#[macro_export]
macro_rules! ln {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::log($e)
        } else {
            $e.ln()
        }
    };
}
