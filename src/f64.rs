use super::*;

impl Upper for f64 {
    #[inline]
    fn upper(self) -> Self {
        f64::from_bits(self.to_bits() & 0x_ffff_ffff_f800_0000)
    }
}

impl FromMask for Doubled<f64> {
    type Mask = u64;
    fn from_mask(u0: Self::Mask, u1: Self::Mask) -> Self {
        Self::new(f64::from_bits(u0), f64::from_bits(u1))
    }
}

impl Check for f64 {
    fn check(self) -> bool {
        self.is_infinite() || self.is_nan()
    }
}

#[inline]
fn fabsk(x: f64) -> f64 {
    f64::from_bits(0x7fff_ffff_ffff_ffff & x.to_bits())
}

impl core::convert::From<f64> for Doubled<f64> {
    #[inline]
    fn from(f: f64) -> Self {
        Self::new(f, 0.)
    }
}

impl core::convert::From<Doubled<f64>> for f64 {
    #[inline]
    fn from(f: Doubled<f64>) -> Self {
        f.0 + f.1
    }
}

impl Doubled<f64> {
    #[inline]
    pub fn abs(self) -> Self {
        Self::new(
            if self.0 < 0. { -self.0 } else { self.0 },
            if self.0 < 0. { -self.1 } else { self.1 },
        )
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
    pub fn square_as_f(self) -> f64 {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        xh * self.1 + xh * self.1 + xl * xl + (xh * xl + xh * xl) + xh * xh
    }

    #[inline]
    pub fn mul_as_f(self, other: Self) -> f64 {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        let yh = other.0.upper();
        let yl = other.0 - yh;
        self.1 * yh + xh * other.1 + xl * yl + xh * yl + xl * yh + xh * yh
    }
}

impl CheckOrder for Doubled<f64> {
    fn check_order(self, other: Self) {
        debug_assert!(
            self.0 == 0.
                || self.0.check()
                || other.0.check()
                || fabsk(self.0) >= fabsk(other.0)
                || ((fabsk(self.0 + other.0) <= fabsk(self.0))
                    && (fabsk(self.0 + other.0) <= fabsk(other.0))),
            "[ddadd_d2_d2_d2 : {:e} {:e}]\n",
            self.0,
            other.0
        );
    }
}

impl CheckOrder<f64> for Doubled<f64> {
    fn check_order(self, other: f64) {
        debug_assert!(
            self.0.check()
                || other.check()
                || fabsk(self.0) >= fabsk(other)
                || ((fabsk(self.0 + other) <= fabsk(self.0))
                    && (fabsk(self.0 + other) <= fabsk(other))),
            "[ddadd_d2_d2_d : {:e} {:e}]\n",
            self.0,
            other
        );
    }
}

impl CheckOrder<Doubled<f64>> for f64 {
    fn check_order(self, other: Doubled<f64>) {
        debug_assert!(
            self.check()
                || other.0.check()
                || fabsk(self) >= fabsk(other.0)
                || ((fabsk(self + other.0) <= fabsk(self))
                    && (fabsk(self + other.0) <= fabsk(other.0))),
            "[ddadd_d2_d_d2 : {:e} {:e}]\n",
            self,
            other.0
        );
    }
}

impl CheckOrder for f64 {
    fn check_order(self, other: Self) {
        debug_assert!(
            self.check()
                || other.check()
                || fabsk(self) >= fabsk(other)
                || ((fabsk(self + other) <= fabsk(self)) && (fabsk(self + other) <= fabsk(other))),
            "[ddadd_d2_d_d : {:e}, {:e}]\n",
            self,
            other
        );
    }
}

impl core::ops::Add<Doubled<f64>> for f64 {
    type Output = Doubled<f64>;
    #[inline]
    fn add(self, other: Doubled<f64>) -> Self::Output {
        let r0 = self + other.0;
        let v = r0 - self; // == other.0
        Doubled::new(r0, self - (r0 - v) + (other.0 - v) + other.1) // [other.0+self, other.1]
    }
}

impl core::ops::Mul for Doubled<f64> {
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

impl core::ops::Mul<f64> for Doubled<f64> {
    type Output = Self;
    #[inline]
    fn mul(self, other: f64) -> Self {
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

impl core::ops::Div for Doubled<f64> {
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

impl AsDoubled for f64 {
    #[inline]
    fn as_doubled(self) -> Doubled<Self> {
        Doubled::new(self, 0.)
    }
}

impl MulAsDoubled for f64 {
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

impl RecPre for Doubled<f64> {
    fn recpre(self) -> Doubled<f64> {
        let t = 1. / self.0;
        let dh = self.0.upper();
        let dl = self.0 - dh;
        let th = t.upper();
        let tl = t - th;
        let q0 = t;
        Self::new(
            q0,
            t * (1. - dh * th - dh * tl - dl * th - dl * tl - self.1 * t),
        )
    }
}

impl RecPre<Doubled<f64>> for f64 {
    fn recpre(self) -> Doubled<f64> {
        let t = 1. / self;
        let dh = self.upper();
        let dl = self - dh;
        let th = t.upper();
        let tl = t - th;
        let q0 = t;
        Doubled::new(q0, t * (1. - dh * th - dh * tl - dl * th - dl * tl))
    }
}
