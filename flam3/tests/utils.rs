use std::{fs::File, io::BufReader, path::PathBuf};

use flam3::{flam3_from_reader, render, RenderOptions};
use image::{buffer::ConvertBuffer, open, RgbImage};
use image_compare::rgb_hybrid_compare;

const TEST_SEED: &str = "foobar";
const MAX_ERROR: f64 = 0.00001;

pub fn run_test(path: &str, name: &str, mut options: RenderOptions) {
    let mut test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    test_dir.push("tests");
    test_dir.push(path);

    println!("Testing {}", name);

    if options.isaac_seed.is_none() {
        options.isaac_seed = Some(TEST_SEED.to_string())
    }

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
