mod utils;

// These tests only pass with the libm feature enabled.
#[cfg(feature = "libm")]
mod render {
    use super::utils::run_test;

    fn do_test(name: &str) {
        run_test("render", name, Default::default());
    }

    macro_rules! render_test {
        ($name:ident) => {
            #[test]
            fn $name() {
                do_test(stringify!($name));
            }
        };
    }

    render_test!(render000);
    render_test!(render001);
    render_test!(render002);
    render_test!(render003);
    render_test!(render004);
    render_test!(render005);
    render_test!(render006);
    render_test!(render007);
    render_test!(render008);
    render_test!(render009);
    render_test!(render010);
    render_test!(render011);
    render_test!(render012);
    render_test!(render013);
    render_test!(render014);
    render_test!(render015);
    render_test!(render016);
    render_test!(render017);
    render_test!(render018);
    render_test!(render019);
    render_test!(render020);
    render_test!(render021);
    render_test!(render022);
    render_test!(render023);
    render_test!(render024);
    render_test!(render025);
    render_test!(render026);
    render_test!(render027);
    render_test!(render028);
    render_test!(render029);
    render_test!(render030);
    render_test!(render031);
    render_test!(render032);
    render_test!(render033);
    render_test!(render034);
    render_test!(render035);
    render_test!(render036);
    render_test!(render037);
    render_test!(render038);
    render_test!(render039);
    render_test!(render040);
    render_test!(render041);
    render_test!(render042);
    render_test!(render043);
    render_test!(render044);
    render_test!(render045);
    render_test!(render046);
    render_test!(render047);
    render_test!(render048);
    render_test!(render049);
    render_test!(render050);
    render_test!(render051);
    render_test!(render052);
    render_test!(render053);
    render_test!(render054);
    render_test!(render055);
    render_test!(render056);
    render_test!(render057);
    render_test!(render058);
    render_test!(render059);
    render_test!(render060);
    render_test!(render061);
    render_test!(render062);
    render_test!(render063);
    render_test!(render064);
    render_test!(render065);
    render_test!(render066);
    render_test!(render067);
    render_test!(render068);
    render_test!(render069);
    render_test!(render070);
    render_test!(render071);
    render_test!(render072);
    render_test!(render073);
    render_test!(render074);
    render_test!(render075);
    render_test!(render076);
    render_test!(render077);
    render_test!(render078);
    render_test!(render079);
    render_test!(render080);
    render_test!(render081);
    render_test!(render082);
    render_test!(render083);
    render_test!(render084);
    render_test!(render085);
    render_test!(render086);
    render_test!(render087);
    render_test!(render088);
    render_test!(render089);
    render_test!(render090);
    render_test!(render091);
    render_test!(render092);
    render_test!(render093);
    render_test!(render094);
    render_test!(render095);
    render_test!(render096);
    render_test!(render097);
    render_test!(render098);
    render_test!(render099);
    render_test!(render100);
    render_test!(render101);
    render_test!(render102);
    render_test!(render103);
    render_test!(render104);
}
