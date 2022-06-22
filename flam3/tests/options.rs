// These tests only pass with the libm feature enabled.
#[cfg(feature = "libm")]
mod options {
    use std::{fs::File, io::BufReader, path::PathBuf};

    use flam3::{flam3_from_reader, render, RenderOptions};
    use image::{buffer::ConvertBuffer, open, RgbImage};
    use image_compare::rgb_hybrid_compare;

    const TEST_SEED: &str = "foobar";
    const MAX_ERROR: f64 = 0.00001;

    fn do_test(name: &str, options: RenderOptions) {
        let mut test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        test_dir.push("tests");
        test_dir.push("options");

        println!("Testing {}", name);

        let mut flame = test_dir.to_owned();
        flame.push(format!("{}.flam3", name));

        let genome = flam3_from_reader(BufReader::new(File::open(flame).unwrap()))
            .unwrap()
            .remove(0);

        let generated_image: RgbImage = render(genome, options).unwrap().convert();

        let mut image = test_dir.to_owned();
        image.push(format!("{}.png", name));

        let expected_image = open(image).unwrap().into_rgb8();
        assert_eq!(generated_image.width(), expected_image.width());
        assert_eq!(generated_image.height(), expected_image.height());

        let diff = rgb_hybrid_compare(&generated_image, &expected_image).unwrap();
        let score = 1.0 - diff.score;

        if score > MAX_ERROR {
            panic!("Result differed by {:.6}", score);
        }
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
                        isaac_seed: Some(TEST_SEED.to_string()),
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
        earlyclip: true,
    });
}
