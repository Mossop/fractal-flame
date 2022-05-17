use rand::{thread_rng, RngCore};

use crate::utils::PanicCast;

const RANDSIZL: usize = 4;
const RANDSIZ: usize = 1 << RANDSIZL;

fn ind(mm: &[u64; RANDSIZ], x: usize) -> u64 {
    mm[(x >> 2) & (RANDSIZ - 1)]
}

macro_rules! increment {
    ($a:expr, $b:expr) => {
        $a = $a.wrapping_add($b)
    };
}

macro_rules! mix {
    ($a:expr, $b: expr, $c: expr, $d: expr, $e: expr, $f: expr, $g: expr, $h: expr) => {
        $a ^= $b << 11;
        increment!($d, $a);
        increment!($b, $c);
        $b ^= $c >> 2;
        increment!($e, $b);
        increment!($c, $d);
        $c ^= $d << 8;
        increment!($f, $c);
        increment!($d, $e);
        $d ^= $e >> 16;
        increment!($g, $d);
        increment!($e, $f);
        $e ^= $f << 10;
        increment!($h, $e);
        increment!($f, $g);
        $f ^= $g >> 4;
        increment!($a, $f);
        increment!($g, $h);
        $g ^= $h << 8;
        increment!($b, $g);
        increment!($h, $a);
        $h ^= $a >> 9;
        increment!($c, $h);
        increment!($a, $b);
    };
}

macro_rules! rngstep {
    ($randmem:expr, $randrsl:expr, $mix:expr, $a:expr, $b:expr, $m:expr, $m2:expr, $r:expr) => {
        let x = $randmem[$m];
        $a = ($a ^ $mix).wrapping_add($randmem[$m2]) & 0xffffffff;
        $m2 += 1;
        let y = ind(&$randmem, x as usize).wrapping_add($a).wrapping_add($b) & 0xffffffff;
        $randmem[$m] = y;
        $m += 1;
        $b = ind(&$randmem, (y >> RANDSIZL) as usize).wrapping_add(x) & 0xffffffff;
        $randrsl[$r] = $b;
        $r += 1;
    };
}

#[derive(Clone)]
pub struct Flam3Rng {
    randcnt: usize,
    randrsl: [u64; RANDSIZ],
    randmem: [u64; RANDSIZ],
    randa: u64,
    randb: u64,
    randc: u64,
}

impl Default for Flam3Rng {
    fn default() -> Self {
        let mut seed = [0; RANDSIZ * 8];
        thread_rng().fill_bytes(&mut seed);

        Self::irandinit(Some(seed))
    }
}

impl Flam3Rng {
    pub fn from_seed(str: &str) -> Flam3Rng {
        let mut seed = [0; RANDSIZ * 8];
        let mut bytes = str.as_bytes();
        if bytes.len() > RANDSIZ * 8 {
            bytes = &bytes[0..RANDSIZ * 8]
        }
        seed[0..bytes.len()].copy_from_slice(str.as_bytes());

        Self::irandinit(Some(seed))
    }

    pub fn from_rng(rng: &mut Flam3Rng) -> Flam3Rng {
        let mut seed = [0; RANDSIZ * 8];

        for i in (0..RANDSIZ * 8).step_by(8) {
            (&mut seed[i..i + 4]).copy_from_slice(&rng.irand().to_le_bytes());
        }

        Self::irandinit(Some(seed))
    }

    fn isaac(&mut self) {
        let mut a = self.randa;
        self.randc += 1;
        let mut b = self.randb + self.randc;

        let mut r = 0;
        let mut m = 0;
        let mut m2 = RANDSIZ / 2;
        let mend = RANDSIZ / 2;

        while m < mend {
            rngstep!(self.randmem, self.randrsl, a << 13, a, b, m, m2, r);
            rngstep!(self.randmem, self.randrsl, a >> 6, a, b, m, m2, r);
            rngstep!(self.randmem, self.randrsl, a << 2, a, b, m, m2, r);
            rngstep!(self.randmem, self.randrsl, a >> 16, a, b, m, m2, r);
        }

        m2 = 0;
        while m2 < mend {
            rngstep!(self.randmem, self.randrsl, a << 13, a, b, m, m2, r);
            rngstep!(self.randmem, self.randrsl, a >> 6, a, b, m, m2, r);
            rngstep!(self.randmem, self.randrsl, a << 2, a, b, m, m2, r);
            rngstep!(self.randmem, self.randrsl, a >> 16, a, b, m, m2, r);
        }

        self.randa = a;
        self.randb = b;
    }

    fn irandinit(seed: Option<[u8; RANDSIZ * 8]>) -> Self {
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = [0x9e3779b9_u64; 8];
        let mut m = [0_u64; RANDSIZ];
        let r = [0_u64; RANDSIZ];

        for _i in 0..4 {
            mix!(a, b, c, d, e, f, g, h);
        }

        if let Some(bytes) = seed {
            let mut seed = [0_u64; RANDSIZ];
            for (i, seed_val) in seed.iter_mut().enumerate().take(RANDSIZ) {
                let pos = i * 8;
                let mut le = [0_u8; 8];
                le.copy_from_slice(&bytes[pos..(pos + 8)]);
                *seed_val = u64::from_le_bytes(le);
            }

            /* initialize using the contents of r[] as the seed */
            for i in (0..RANDSIZ).step_by(8) {
                increment!(a, seed[i]);
                increment!(b, seed[i + 1]);
                increment!(c, seed[i + 2]);
                increment!(d, seed[i + 3]);
                increment!(e, seed[i + 4]);
                increment!(f, seed[i + 5]);
                increment!(g, seed[i + 6]);
                increment!(h, seed[i + 7]);
                mix!(a, b, c, d, e, f, g, h);
                m[i] = a;
                m[i + 1] = b;
                m[i + 2] = c;
                m[i + 3] = d;
                m[i + 4] = e;
                m[i + 5] = f;
                m[i + 6] = g;
                m[i + 7] = h;
            }
            /* do a second pass to make all of the seed affect all of m */
            for i in (0..RANDSIZ).step_by(8) {
                increment!(a, m[i]);
                increment!(b, m[i + 1]);
                increment!(c, m[i + 2]);
                increment!(d, m[i + 3]);
                increment!(e, m[i + 4]);
                increment!(f, m[i + 5]);
                increment!(g, m[i + 6]);
                increment!(h, m[i + 7]);
                mix!(a, b, c, d, e, f, g, h);
                m[i] = a;
                m[i + 1] = b;
                m[i + 2] = c;
                m[i + 3] = d;
                m[i + 4] = e;
                m[i + 5] = f;
                m[i + 6] = g;
                m[i + 7] = h;
            }
        } else {
            for i in (0..RANDSIZ).step_by(8) {
                /* fill in mm[] with messy stuff */
                mix!(a, b, c, d, e, f, g, h);
                m[i] = a;
                m[i + 1] = b;
                m[i + 2] = c;
                m[i + 3] = d;
                m[i + 4] = e;
                m[i + 5] = f;
                m[i + 6] = g;
                m[i + 7] = h;
            }
        }

        let mut rng = Self {
            randa: 0,
            randb: 0,
            randc: 0,
            randmem: m,
            randrsl: r,
            randcnt: RANDSIZ,
        };

        rng.isaac();
        rng
    }

    /// Returns an unsigned 32-bit random number.
    pub fn irand(&mut self) -> u32 {
        if self.randcnt == 0 {
            self.isaac();
            self.randcnt = RANDSIZ - 1;
        } else {
            self.randcnt -= 1;
        }

        self.randrsl[self.randcnt] as u32
    }

    /// Returns a random bit.
    pub fn isaac_bit(&mut self) -> bool {
        (self.irand() & 1) != 0
    }

    /// Returns a random number in the range 0..1.
    pub fn isaac_01(&mut self) -> f64 {
        (self.irand() as i32 & 0xfffffff).f64() / 0xfffffff.f64()
    }

    /// Returns a random number in the range -1..1.
    pub fn isaac_11(&mut self) -> f64 {
        ((self.irand() as i32 & 0xfffffff) - 0x7ffffff).f64() / 0x7ffffff.f64()
    }
}

#[cfg(test)]
mod tests {
    use crate::render::flam3::rng::Flam3Rng;

    #[test]
    fn empty_seed() {
        let mut rng = Flam3Rng::from_seed("");

        assert_eq!(rng.irand(), 2452851915);
        assert_eq!(rng.irand(), 1340468986);
        assert_eq!(rng.irand(), 3346877890);
        assert_eq!(rng.irand(), 3366115673);
        assert_eq!(rng.irand(), 61449112);
        assert_eq!(rng.irand(), 1440946401);
        assert_eq!(rng.irand(), 2756852084);
    }

    #[test]
    fn irand() {
        let mut rng = Flam3Rng::from_seed("foobar");

        assert_eq!(rng.irand(), 358156358);
        assert_eq!(rng.irand(), 3974535195);
        assert_eq!(rng.irand(), 916074978);
        assert_eq!(rng.irand(), 1958166261);
        assert_eq!(rng.irand(), 96173260);
        assert_eq!(rng.irand(), 224373016);
        assert_eq!(rng.irand(), 3800469598);
    }

    #[test]
    fn skipped() {
        let mut rng = Flam3Rng::from_seed("foobar");

        for _ in 0..500 {
            rng.irand();
        }

        rng = rng.clone();

        for _ in 0..500 {
            rng.irand();
        }

        assert_eq!(rng.irand(), 451929136);
        assert_eq!(rng.irand(), 3706281000);
        assert_eq!(rng.irand(), 891491316);
        assert_eq!(rng.irand(), 3771320315);
        assert_eq!(rng.irand(), 1241496118);
        assert_eq!(rng.irand(), 1321834101);
        assert_eq!(rng.irand(), 4179049759);
        assert_eq!(rng.irand(), 800856678);
        assert_eq!(rng.irand(), 567422021);
        assert_eq!(rng.irand(), 1711913236);
    }

    fn rounded(val: f64) -> String {
        format!("{:.15}", val)
    }

    #[test]
    fn isaac_01() {
        let mut rng = Flam3Rng::from_seed("foobar");

        assert_eq!(rounded(rng.isaac_01()), "0.334236407034980");
        assert_eq!(rounded(rng.isaac_01()), "0.806297405832624");
        assert_eq!(rounded(rng.isaac_01()), "0.412645229744335");
        assert_eq!(rounded(rng.isaac_01()), "0.294737775976724");
        assert_eq!(rounded(rng.isaac_01()), "0.358273313784127");
        assert_eq!(rounded(rng.isaac_01()), "0.835854622855241");
        assert_eq!(rounded(rng.isaac_01()), "0.157852523616897");
        assert_eq!(rounded(rng.isaac_01()), "0.329166286919885");
        assert_eq!(rounded(rng.isaac_01()), "0.222156201385544");
        assert_eq!(rounded(rng.isaac_01()), "0.804080962404910");
    }

    #[test]
    fn isaac_11() {
        let mut rng = Flam3Rng::from_seed("foobar");

        assert_eq!(rounded(rng.isaac_11()), "-0.331527183439785");
        assert_eq!(rounded(rng.isaac_11()), "0.612594817672631");
        assert_eq!(rounded(rng.isaac_11()), "-0.174709537436884");
        assert_eq!(rounded(rng.isaac_11()), "-0.410524445850584");
        assert_eq!(rounded(rng.isaac_11()), "-0.283453369762401");
        assert_eq!(rounded(rng.isaac_11()), "0.671709251938084");
        assert_eq!(rounded(rng.isaac_11()), "-0.684294951590113");
        assert_eq!(rounded(rng.isaac_11()), "-0.341667423707749");
        assert_eq!(rounded(rng.isaac_11()), "-0.555687595573720");
        assert_eq!(rounded(rng.isaac_11()), "0.608161930800691");
    }

    #[test]
    fn isaac_bit() {
        let mut rng = Flam3Rng::from_seed("foobar");

        assert!(!rng.isaac_bit());
        assert!(rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(rng.isaac_bit());
        assert!(rng.isaac_bit());
        assert!(rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(rng.isaac_bit());
        assert!(!rng.isaac_bit());
        assert!(!rng.isaac_bit());
    }

    #[test]
    fn from_rng() {
        let mut rng = Flam3Rng::from_seed("foobar");

        assert_eq!(rng.irand(), 358156358);
        assert_eq!(rng.irand(), 3974535195);
        assert_eq!(rng.irand(), 916074978);
        assert_eq!(rng.irand(), 1958166261);

        rng = Flam3Rng::from_rng(&mut rng);

        assert_eq!(rng.irand(), 4279217379);
        assert_eq!(rng.irand(), 1406626236);
        assert_eq!(rng.irand(), 2627953629);
        assert_eq!(rng.irand(), 3924363182);
        assert_eq!(rng.irand(), 1953805307);
        assert_eq!(rng.irand(), 2271928233);
        assert_eq!(rng.irand(), 4227309479);
        assert_eq!(rng.irand(), 1961718375);
        assert_eq!(rng.irand(), 2344761049);
        assert_eq!(rng.irand(), 1051335410);
        assert_eq!(rng.irand(), 2622661706);
    }
}
