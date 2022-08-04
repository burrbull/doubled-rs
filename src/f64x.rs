macro_rules! impl_doubled_f64 {
    ($size:literal) => {
        use crate::*;

        type F64x = packed_simd::Simd<[f64; $size]>;
        type U64x = packed_simd::Simd<[u64; $size]>;

        impl Upper for F64x {
            #[inline]
            fn upper(self) -> Self {
                F64x::from_bits(
                    U64x::from_bits(self) & U64x::splat((0x_ffff_ffff << 32) + 0x_f800_0000),
                )
            }
        }

        impl FromMask for Doubled<F64x> {
            type Mask = U64x;
            fn from_mask(u0: Self::Mask, u1: Self::Mask) -> Self {
                Self::new(F64x::from_bits(u0), F64x::from_bits(u1))
            }
        }

        impl core::convert::From<F64x> for Doubled<F64x> {
            #[inline]
            fn from(f: F64x) -> Self {
                Self::new(f, F64x::splat(0.))
            }
        }

        impl core::convert::From<Doubled<F64x>> for F64x {
            #[inline]
            fn from(f: Doubled<F64x>) -> Self {
                f.0 + f.1
            }
        }

        impl Doubled<F64x> {
            #[inline]
            pub const fn splat(value: Doubled<f64>) -> Self {
                Self::new(F64x::splat(value.0), F64x::splat(value.1))
            }

            #[inline]
            pub fn abs(self) -> Self {
                Self::new(
                    self.0.abs(),
                    F64x::from_bits(
                        U64x::from_bits(self.1)
                            ^ (U64x::from_bits(self.0) & U64x::from_bits(F64x::splat(-0.))),
                    ),
                )
            }

            #[cfg(target_feature = "fma")]
            #[inline]
            pub fn square(self) -> Self {
                let r0 = self.0 * self.0;
                Self::new(
                    r0,
                    (self.0 + self.0).mul_adde(self.1, self.0.mul_sube(self.0, r0)),
                )
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            pub fn square(self) -> Self {
                let xh = self.0.upper();
                let xl = self.0 - xh;
                let r0 = self.0 * self.0;
                Self::new(
                    r0,
                    xh * xh + (-r0) + (xh + xh) * xl + xl * xl + self.0 * (self.1 + self.1),
                )
            }

            #[cfg(target_feature = "fma")]
            #[inline]
            pub fn square_as_f(self) -> F64x {
                self.0.mul_adde(self.0, self.0 * self.1 + self.0 * self.1)
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            pub fn square_as_f(self) -> F64x {
                let xh = self.0.upper();
                let xl = self.0 - xh;

                xh * self.1 + xh * self.1 + xl * xl + (xh * xl + xh * xl) + xh * xh
            }

            #[inline]
            pub fn sqrt(self) -> Self {
                let t = (self.0 + self.1).sqrt();
                ((self + t.mul_as_doubled(t)) * t.recpre_as_doubled()).scale(F64x::splat(0.5))
            }

            #[cfg(target_feature = "fma")]
            #[inline]
            pub fn mul_as_f(self, other: Self) -> F64x {
                self.0
                    .mul_adde(other.0, self.1.mul_adde(other.0, self.0 * other.1))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            pub fn mul_as_f(self, other: Self) -> F64x {
                let xh = self.0.upper();
                let xl = self.0 - xh;
                let yh = other.0.upper();
                let yl = other.0 - yh;

                self.1 * yh + xh * other.1 + xl * yl + xh * yl + xl * yh + xh * yh
            }
        }

        impl core::ops::Add<Doubled<F64x>> for F64x {
            type Output = Doubled<Self>;
            #[inline]
            fn add(self, other: Doubled<Self>) -> Self::Output {
                let r0 = self + other.0;
                let v = r0 - self;
                Doubled::new(r0, self - (r0 - v) + (other.0 - v) + other.1)
            }
        }

        impl core::ops::Mul for Doubled<F64x> {
            type Output = Self;
            #[cfg(target_feature = "fma")]
            #[inline]
            fn mul(self, other: Self) -> Self {
                let r0 = self.0 * other.0;
                Self::new(
                    r0,
                    self.0.mul_adde(
                        other.1,
                        self.1.mul_adde(other.0, self.0.mul_sube(other.0, r0)),
                    ),
                )
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn mul(self, other: Self) -> Self {
                let xh = self.0.upper();
                let xl = self.0 - xh;
                let yh = other.0.upper();
                let yl = other.0 - yh;
                let r0 = self.0 * other.0;
                Self::new(
                    r0,
                    xh * yh
                        + (-r0)
                        + xl * yh
                        + xh * yl
                        + xl * yl
                        + self.0 * other.1
                        + self.1 * other.0,
                )
            }
        }

        impl core::ops::Mul<F64x> for Doubled<F64x> {
            type Output = Self;
            #[cfg(target_feature = "fma")]
            #[inline]
            fn mul(self, other: F64x) -> Self {
                let r0 = self.0 * other;
                Self::new(r0, self.1.mul_adde(other, self.0.mul_sube(other, r0)))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn mul(self, other: F64x) -> Self {
                let xh = self.0.upper();
                let xl = self.0 - xh;
                let yh = other.upper();
                let yl = other - yh;
                let r0 = self.0 * other;
                Self::new(
                    r0,
                    xh * yh + (-r0) + xl * yh + xh * yl + xl * yl + self.1 * other,
                )
            }
        }

        impl core::ops::Div for Doubled<F64x> {
            type Output = Self;
            #[cfg(target_feature = "fma")]
            #[inline]
            fn div(self, other: Self) -> Self {
                let t = other.0.recpre();

                let q0 = self.0 * t;
                let u = t.mul_sube(self.0, q0);
                let mut q1 = other.1.fmanp(t, other.0.fmanp(t, F64x::splat(1.)));
                q1 = q0.mul_adde(q1, self.1.mul_adde(t, u));

                Self::new(q0, q1)
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn div(self, other: Self) -> Self {
                let t = other.0.recpre();
                let dh = other.0.upper();
                let dl = other.0 - dh;
                let th = t.upper();
                let tl = t - th;
                let nhh = self.0.upper();
                let nhl = self.0 - nhh;

                let q0 = self.0 * t;

                let u = nhh * th - q0
                    + nhh * tl
                    + nhl * th
                    + nhl * tl
                    + q0 * (F64x::splat(1.) - dh * th - dh * tl - dl * th - dl * tl);

                Self::new(q0, t.mul_add(self.1 - q0 * other.1, u))
            }
        }

        impl CheckOrder for Doubled<F64x> {
            fn check_order(self, _other: Self) {}
        }

        impl CheckOrder<F64x> for Doubled<F64x> {
            fn check_order(self, _other: F64x) {}
        }

        impl CheckOrder<Doubled<F64x>> for F64x {
            fn check_order(self, _other: Doubled<F64x>) {}
        }

        impl CheckOrder for F64x {
            fn check_order(self, _other: Self) {}
        }

        impl AsDoubled for F64x {
            #[inline]
            fn as_doubled(self) -> Doubled<Self> {
                Doubled::new(self, Self::splat(0.))
            }
        }

        impl MulAsDoubled for F64x {
            #[cfg(target_feature = "fma")]
            #[inline]
            fn mul_as_doubled(self, other: Self) -> Doubled<Self> {
                let r0 = self * other;
                Doubled::new(r0, self.mul_sube(other, r0))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn mul_as_doubled(self, other: Self) -> Doubled<Self> {
                let xh = self.upper();
                let xl = self - xh;
                let yh = other.upper();
                let yl = other - yh;
                let r0 = self * other;
                Doubled::new(r0, xh * yh + (-r0) + xl * yh + xh * yl + xl * yl)
            }
        }

        impl RecPre for Doubled<F64x> {
            #[cfg(target_feature = "fma")]
            #[inline]
            fn recpre(self) -> Self {
                let q0 = self.0.recpre();
                Self::new(q0, q0 * self.1.fmanp(q0, self.0.fmanp(q0, F64x::splat(1))))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn recpre(self) -> Doubled<F64x> {
                let t = self.0.recpre();
                let dh = self.0.upper();
                let dl = self.0 - dh;
                let th = t.upper();
                let tl = t - th;
                let q0 = t;
                Self::new(
                    q0,
                    t * (F64x::splat(1.) - dh * th - dh * tl - dl * th - dl * tl - self.1 * t),
                )
            }
        }

        impl RecPreAsDoubled for F64x {
            #[cfg(target_feature = "fma")]
            #[inline]
            fn recpre_as_doubled(self) -> Doubled<Self> {
                let q0 = self.recpre();
                Doubled::new(q0, q0 * (self, q0, Self::splat(1.)).fmanp())
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn recpre_as_doubled(self) -> Doubled<Self> {
                let t = self.recpre();
                let dh = self.upper();
                let dl = self - dh;
                let th = t.upper();
                let tl = t - th;
                let q0 = t;
                Doubled::new(
                    q0,
                    t * (Self::splat(1.) - dh * th - dh * tl - dl * th - dl * tl),
                )
            }
        }
    };
}

pub mod f64x2 {
    use packed_simd::*;
    impl_doubled_f64!(2);
}

pub mod f64x4 {
    use packed_simd::*;
    impl_doubled_f64!(4);
}

pub mod f64x8 {
    use packed_simd::*;
    impl_doubled_f64!(8);
}
