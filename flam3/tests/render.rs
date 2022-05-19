use std::{
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use flam3::{flam3_from_reader, render, RenderOptions};
use png::{BitDepth, ColorType};

const TEST_SEED: &str = "foobar";

fn do_test(test_dir: &Path, name: &str) {
    println!("Testing {}", name);

    let mut flame = test_dir.to_owned();
    flame.push(format!("{}.flam3", name));

    let genome = flam3_from_reader(BufReader::new(File::open(flame).unwrap()))
        .unwrap()
        .remove(0);

    let data = render(
        genome,
        RenderOptions {
            isaac_seed: Some(TEST_SEED.to_string()),
            ..Default::default()
        },
    )
    .unwrap();

    let mut image = test_dir.to_owned();
    image.push(format!("{}.png", name));

    let decoder = png::Decoder::new(File::open(image).unwrap());
    let mut reader = decoder.read_info().unwrap();
    let mut buf = vec![0; reader.output_buffer_size()];
    let info = reader.next_frame(&mut buf).unwrap();

    assert_eq!(info.color_type, ColorType::Rgba);
    assert_eq!(info.bit_depth, BitDepth::Eight);

    let bytes = &buf[..info.buffer_size()];
    assert_eq!(bytes.len(), data.len());

    let mut index = 0;
    let mut diff_sum = 0;
    for y in 0..info.height {
        for x in 0..info.width {
            let diff = if data[index] > bytes[index] {
                data[index] - bytes[index]
            } else {
                bytes[index] - data[index]
            };

            if diff > 0 {
                println!("Rendered data differed by {} at x={}, y={}.", diff, x, y);
                diff_sum += diff;
            }

            index += 1;
        }
    }

    if diff_sum > 10 {
        panic!("Total image difference surpassed the limit.");
    }
}

#[test]
fn render_tests() {
    let mut test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    test_dir.push("tests");

    do_test(&test_dir, "test1");
    do_test(&test_dir, "test2");
}
