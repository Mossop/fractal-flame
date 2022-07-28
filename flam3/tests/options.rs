mod utils;

// These tests only pass with the libm feature enabled.
#[cfg(feature = "libm")]
mod options {
    use flam3::RenderOptions;

    use super::utils::run_test;

    fn do_test(name: &str, options: RenderOptions) {
        run_test("options", name, options);
    }

    macro_rules! render_test {
        ($name:ident, {
            $(
                $field_name:ident : $field_value:expr
            ),*$(,)+
        }) => {
            #[test]
            fn $name() {
                do_test(
                    stringify!($name),
                    RenderOptions {
                         $(
                             $field_name: $field_value,
                        )*
                        ..Default::default()
                    },
                );
            }
        };
    }

    render_test!(earlyclip, {
        threads: Some(1),
        earlyclip: true,
    });
}
