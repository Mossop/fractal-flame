use std::{fs::File, io::BufReader, path::PathBuf};

use flam3::{flam3_from_reader, render, RenderOptions};
use png::{BitDepth, ColorType};

const TEST_SEED: &str = "foobar";

fn do_test(name: &str) {
    let mut test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    test_dir.push("tests");

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

            diff_sum += diff;
            index += 1;
        }
    }

    if diff_sum > 10 {
        panic!(
            "Total image difference ({}) surpassed the limit ({}).",
            diff_sum, 10
        );
    }
}

macro_rules! render_test {
    ($name:ident) => {
        mod $name {
            #[test]
            fn render_test() {
                super::do_test(stringify!($name));
            }
        }
    };
}

render_test!(test000);
render_test!(test001);
render_test!(test002);
render_test!(test003);
render_test!(test004);
render_test!(test005);
render_test!(test006);
render_test!(test007);
render_test!(test008);
render_test!(test009);
render_test!(test010);
render_test!(test011);
render_test!(test012);
render_test!(test013);
render_test!(test014);
render_test!(test015);
render_test!(test016);
render_test!(test017);
render_test!(test018);
render_test!(test019);
render_test!(test020);
render_test!(test021);
render_test!(test022);
render_test!(test023);
render_test!(test024);
render_test!(test025);
render_test!(test026);
render_test!(test027);
render_test!(test028);
render_test!(test029);
render_test!(test030);
render_test!(test031);
render_test!(test032);
render_test!(test033);
render_test!(test034);
render_test!(test035);
render_test!(test036);
render_test!(test037);
render_test!(test038);
render_test!(test039);
render_test!(test040);
render_test!(test041);
render_test!(test042);
render_test!(test043);
render_test!(test044);
render_test!(test045);
render_test!(test046);
render_test!(test047);
render_test!(test048);
render_test!(test049);
render_test!(test050);
render_test!(test051);
render_test!(test052);
render_test!(test053);
render_test!(test054);
render_test!(test055);
render_test!(test056);
render_test!(test057);
render_test!(test058);
render_test!(test059);
render_test!(test060);
render_test!(test061);
render_test!(test062);
render_test!(test063);
render_test!(test064);
render_test!(test065);
render_test!(test066);
render_test!(test067);
render_test!(test068);
render_test!(test069);
render_test!(test070);
render_test!(test071);
render_test!(test072);
render_test!(test073);
render_test!(test074);
render_test!(test075);
render_test!(test076);
render_test!(test077);
render_test!(test078);
render_test!(test079);
render_test!(test080);
render_test!(test081);
render_test!(test082);
render_test!(test083);
render_test!(test084);
render_test!(test085);
render_test!(test086);
render_test!(test087);
render_test!(test088);
render_test!(test089);
render_test!(test090);
render_test!(test091);
render_test!(test092);
render_test!(test093);
render_test!(test094);
render_test!(test095);
render_test!(test096);
render_test!(test097);
render_test!(test098);
render_test!(test099);
render_test!(test100);
