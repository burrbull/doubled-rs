use super::*;

impl Upper for f32 {
    #[inline]
    fn upper(self) -> Self {
        f32::from_bits(self.to_bits() & 0x_ffff_f000)
    }
}

impl FromMask for Doubled<f32> {
    type Mask = u32;
    fn from_mask(u0: Self::Mask, u1: Self::Mask) -> Self {
        Self::new(f32::from_bits(u0), f32::from_bits(u1))
    }
}

impl Check for f32 {
    fn check(self) -> bool {
        self.is_infinite() || self.is_nan()
    }
}

#[inline]
fn fabsfk(x: f32) -> f32 {
    f32::from_bits(0x_7fff_ffff & x.to_bits())
}

impl core::convert::From<f32> for Doubled<f32> {
    #[inline]
    fn from(f: f32) -> Self {
        Self::new(f, 0.)
    }
}

impl core::convert::From<Doubled<f32>> for f32 {
    #[inline]
    fn from(f: Doubled<f32>) -> Self {
        f.0 + f.1
    }
}

impl core::convert::From<f64> for Doubled<f32> {
    #[inline]
    fn from(f: f64) -> Self {
        Self::from_f64(f)
    }
}

impl Doubled<f32> {
    #[inline]
    pub fn abs(self) -> Self {
        if self.0 < 0. {
            Self::new(-self.0, -self.1)
        } else {
            self
        }
    }

    #[inline]
    pub fn square(self) -> Self {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        let r0 = self.0 * self.0;
        Self::new(
            r0,
            xh * xh - r0 + (xh + xh) * xl + xl * xl + self.0 * (self.1 + self.1),
        )
    }

    #[inline]
    pub fn square_as_f(self) -> f32 {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        xh * self.1 + xh * self.1 + xl * xl + (xh * xl + xh * xl) + xh * xh
    }

    #[inline]
    pub fn mul_as_f(self, other: Self) -> f32 {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        let yh = other.0.upper();
        let yl = other.0 - yh;
        self.1 * yh + xh * other.1 + xl * yl + xh * yl + xl * yh + xh * yh
    }

    #[inline]
    pub fn from_f64(f: f64) -> Self {
        let x = f as f32;
        Self::new(x, (f - (x as f64)) as f32)
    }

    pub fn recip(self) -> Self {
        let t = 1. / self.0;
        let dh = self.0.upper();
        let dl = self.0 - dh;
        let th = t.upper();
        let tl = t - th;
        Doubled::new(
            t,
            t * (1. - dh * th - dh * tl - dl * th - dl * tl - self.1 * t),
        )
    }
}

impl CheckOrder for Doubled<f32> {
    #[inline]
    fn check_order(self, other: Self) {
        self.0.check_order(other.0)
    }
}

impl CheckOrder<f32> for Doubled<f32> {
    #[inline]
    fn check_order(self, other: f32) {
        self.0.check_order(other)
    }
}

impl CheckOrder<Doubled<f32>> for f32 {
    #[inline]
    fn check_order(self, other: Doubled<f32>) {
        self.check_order(other.0)
    }
}

impl CheckOrder for f32 {
    #[inline]
    fn check_order(self, other: Self) {
        debug_assert!(self.check() || other.check() || fabsfk(self) >= fabsfk(other));
    }
}

impl core::ops::Add<Doubled<f32>> for f32 {
    type Output = Doubled<f32>;
    #[inline]
    fn add(self, other: Doubled<f32>) -> Self::Output {
        let r0 = self + other.0;
        let v = r0 - self; // == other.0
        Doubled::new(r0, self - (r0 - v) + (other.0 - v) + other.1) // [other.0+self, other.1]
    }
}

impl core::ops::Mul for Doubled<f32> {
    type Output = Self;
    #[inline]
    fn mul(self, other: Self) -> Self {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        let yh = other.0.upper();
        let yl = other.0 - yh;
        let r0 = self.0 * other.0;
        Self::new(
            r0,
            xh * yh - r0 + xl * yh + xh * yl + xl * yl + self.0 * other.1 + self.1 * other.0,
        )
    }
}

impl core::ops::Mul<f32> for Doubled<f32> {
    type Output = Self;
    #[inline]
    fn mul(self, other: f32) -> Self {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        let yh = other.upper();
        let yl = other - yh;
        let r0 = self.0 * other;
        Self::new(
            r0,
            xh * yh - r0 + xl * yh + xh * yl + xl * yl + self.1 * other,
        )
    }
}

impl core::ops::Div for Doubled<f32> {
    type Output = Self;
    #[inline]
    fn div(self, other: Self) -> Self {
        let t = 1. / other.0;
        let dh = other.0.upper();
        let dl = other.0 - dh;
        let th = t.upper();
        let tl = t - th;
        let nhh = self.0.upper();
        let nhl = self.0 - nhh;

        let q0 = self.0 * t;

        let u = -q0
            + nhh * th
            + nhh * tl
            + nhl * th
            + nhl * tl
            + q0 * (1. - dh * th - dh * tl - dl * th - dl * tl);

        Self::new(q0, t * (self.1 - q0 * other.1) + u)
    }
}

impl AsDoubled for f32 {
    #[inline]
    fn as_doubled(self) -> Doubled<Self> {
        Doubled::new(self, 0.)
    }
}

impl MulAsDoubled for f32 {
    #[inline]
    fn mul_as_doubled(self, other: Self) -> Doubled<Self> {
        let xh = self.upper();
        let xl = self - xh;
        let yh = other.upper();
        let yl = other - yh;
        let r0 = self * other;
        Doubled::new(r0, xh * yh - r0 + xl * yh + xh * yl + xl * yl)
    }
}

impl RecipAsDoubled for f32 {
    fn recip_as_doubled(self) -> Doubled<f32> {
        let t = 1. / self;
        let dh = self.upper();
        let dl = self - dh;
        let th = t.upper();
        let tl = t - th;
        Doubled::new(t, t * (1. - dh * th - dh * tl - dl * th - dl * tl))
    }
}
