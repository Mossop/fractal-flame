use libm::{cos, sin, tan};

#[test]
fn math() {
    println!("{:.40}", sin(0.1_f64));
    println!("{:.40}", 0.1_f64.sin());
    println!("{:.40}", cos(0.1_f64));
    println!("{:.40}", 0.1_f64.cos());
    println!("{:.40}", tan(0.1_f64));
    println!("{:.40}", 0.1_f64.tan());
}
