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
macro_rules! sincos {
    ($e:expr) => {
        if cfg!(feature = "libm") {
            libm::sincos($e)
        } else {
            $e.sin_cos()
        }
    };
}
