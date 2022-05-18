use std::{
    cmp::{max, min},
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use flam3::{flam3_from_reader, render, RenderOptions};

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
    let bytes = &buf[..info.buffer_size()];

    assert_eq!(bytes.len(), data.len());

    for (index, data) in data.iter().enumerate() {
        if data != &bytes[index] {
            let diff = max(data, &bytes[index]) - min(data, &bytes[index]);
            panic!("Rendered data differed by {} at offset {}.", diff, index);
        }
    }
}

#[test]
fn render_tests() {
    let mut test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    test_dir.push("tests");

    do_test(&test_dir, "test1");
    do_test(&test_dir, "test2");
}
