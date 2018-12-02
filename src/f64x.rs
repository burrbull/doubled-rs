macro_rules! impl_f2_f64 {
    ($f64x:ident, $u64x:ident, $m64x:ident) => {
        use crate::*;

        impl FromU32 for $u64x {
            fn from_u32(i: (u32, u32)) -> Self {
                let mut a = [0_u32; $u64x::lanes() * 2];
                for j in 0..$u64x::lanes() {
                    a[2 * j] = i.0;
                    a[2 * j + 1] = i.1;
                }
                $u64x::from_bits(Simd::<[u32; $u64x::lanes() * 2]>::from_slice_aligned(&a))
            }
        }

        #[inline]
        fn vupper_vd_vd(d: $f64x) -> $f64x {
            $f64x::from_bits($u64x::from_bits(d) & $u64x::from_u32((0x_ffff_ffff, 0x_f800_0000)))
        }
        impl FloatConsts for $f64x {
            const INFINITY: Self = Self::splat(core::f64::INFINITY);
            const NEG_INFINITY: Self = Self::splat(core::f64::NEG_INFINITY);
            const NAN: Self = Self::splat(core::f64::NAN);
        }

        impl IsInf for $f64x {
            type Mask = $m64x;
            #[inline]
            fn is_infinite(self) -> Self::Mask {
                self.eq(Self::INFINITY) | self.eq(Self::NEG_INFINITY)
            }
        }

        impl IsNan for $f64x {
            type Mask = $m64x;
            #[inline]
            fn is_nan(self) -> Self::Mask {
                self.ne(self)
            }
        }

        impl FromMask for Doubled<$f64x> {
            type Mask = $u64x;
            fn from_mask(u0: Self::Mask, u1: Self::Mask) -> Self {
                Self::new($f64x::from_bits(u0), $f64x::from_bits(u1))
            }
        }

        impl core::convert::From<$f64x> for Doubled<$f64x> {
            #[inline]
            fn from(f: $f64x) -> Self {
                Self::new(f, $f64x::splat(0.))
            }
        }

        impl core::convert::From<Doubled<$f64x>> for $f64x {
            #[inline]
            fn from(f: Doubled<$f64x>) -> Self {
                f.0 + f.1
            }
        }

        impl Doubled<$f64x> {
            #[inline]
            pub fn abs(self) -> Self {
                Self::new(
                    self.0.abs(),
                    $f64x::from_bits(
                        $u64x::from_bits(self.1)
                            ^ ($u64x::from_bits(self.0) & $u64x::from_bits($f64x::splat(-0.))),
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
                let xh = vupper_vd_vd(self.0);
                let xl = self.0 - xh;
                let r0 = self.0 * self.0;
                Self::new(
                    r0,
                    xh * xh + (-r0) + (xh + xh) * xl + xl * xl + self.0 * (self.1 + self.1),
                )
            }

            #[cfg(target_feature = "fma")]
            #[inline]
            pub fn square_as_f(self) -> $f64x {
                self.0.mul_adde(self.0, self.0 * self.1 + self.0 * self.1)
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            pub fn square_as_f(self) -> $f64x {
                let xh = vupper_vd_vd(self.0);
                let xl = self.0 - xh;

                xh * self.1 + xh * self.1 + xl * xl + (xh * xl + xh * xl) + xh * xh
            }

            #[inline]
            pub fn sqrt(self) -> Self {
                let t = (self.0 + self.1).sqrt();
                ((self + t.mul_as_doubled(t)) * t.recpre_as_doubled()).scale($f64x::splat(0.5))
            }

            #[cfg(target_feature = "fma")]
            #[inline]
            pub fn mul_as_f(self, other: Self) -> $f64x {
                self.0
                    .mul_adde(other.0, self.1.mul_adde(other.0, self.0 * other.1))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            pub fn mul_as_f(self, other: Self) -> $f64x {
                let xh = vupper_vd_vd(self.0);
                let xl = self.0 - xh;
                let yh = vupper_vd_vd(other.0);
                let yl = other.0 - yh;

                self.1 * yh + xh * other.1 + xl * yl + xh * yl + xl * yh + xh * yh
            }
        }

        impl core::convert::From<(f64, f64)> for Doubled<$f64x> {
            fn from(f: (f64, f64)) -> Self {
                Self::new($f64x::splat(f.0), $f64x::splat(f.1))
            }
        }

        impl core::ops::Add<Doubled<$f64x>> for $f64x {
            type Output = Doubled<Self>;
            #[inline]
            fn add(self, other: Doubled<Self>) -> Self::Output {
                let r0 = self + other.0;
                let v = r0 - self;
                Doubled::new(r0, self - (r0 - v) + (other.0 - v) + other.1)
            }
        }

        impl core::ops::Mul for Doubled<$f64x> {
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
                let xh = vupper_vd_vd(self.0);
                let xl = self.0 - xh;
                let yh = vupper_vd_vd(other.0);
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

        impl core::ops::Mul<$f64x> for Doubled<$f64x> {
            type Output = Self;
            #[cfg(target_feature = "fma")]
            #[inline]
            fn mul(self, other: $f64x) -> Self {
                let r0 = self.0 * other;
                Self::new(r0, self.1.mul_adde(other, self.0.mul_sube(other, r0)))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn mul(self, other: $f64x) -> Self {
                let xh = vupper_vd_vd(self.0);
                let xl = self.0 - xh;
                let yh = vupper_vd_vd(other);
                let yl = other - yh;
                let r0 = self.0 * other;
                Self::new(
                    r0,
                    xh * yh + (-r0) + xl * yh + xh * yl + xl * yl + self.1 * other,
                )
            }
        }

        impl core::ops::Div for Doubled<$f64x> {
            type Output = Self;
            #[cfg(target_feature = "fma")]
            #[inline]
            fn div(self, other: Self) -> Self {
                let t = other.0.recpre();

                let q0 = self.0 * t;
                let u = t.mul_sube(self.0, q0);
                let mut q1 = other.1.fmanp(t, other.0.fmanp(t, $f64x::splat(1.)));
                q1 = q0.mul_adde(q1, self.1.mul_adde(t, u));

                Self::new(q0, q1)
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn div(self, other: Self) -> Self {
                let t = other.0.recpre();
                let dh = vupper_vd_vd(other.0);
                let dl = other.0 - dh;
                let th = vupper_vd_vd(t);
                let tl = t - th;
                let nhh = vupper_vd_vd(self.0);
                let nhl = self.0 - nhh;

                let q0 = self.0 * t;

                let u = nhh * th - q0
                    + nhh * tl
                    + nhl * th
                    + nhl * tl
                    + q0 * ($f64x::splat(1.) - dh * th - dh * tl - dl * th - dl * tl);

                Self::new(q0, t.mul_add(self.1 - q0 * other.1, u))
            }
        }

        impl CheckOrder for Doubled<$f64x> {
            fn check_order(self, _other: Self) {}
        }

        impl CheckOrder<$f64x> for Doubled<$f64x> {
            fn check_order(self, _other: $f64x) {}
        }

        impl CheckOrder<Doubled<$f64x>> for $f64x {
            fn check_order(self, _other: Doubled<$f64x>) {}
        }

        impl CheckOrder for $f64x {
            fn check_order(self, _other: Self) {}
        }

        impl AsDoubled for $f64x {
            #[inline]
            fn as_doubled(self) -> Doubled<Self> {
                Doubled::new(self, Self::splat(0.))
            }
        }

        impl MulAsDoubled for $f64x {
            #[cfg(target_feature = "fma")]
            #[inline]
            fn mul_as_doubled(self, other: Self) -> Doubled<Self> {
                let r0 = self * other;
                Doubled::new(r0, self.mul_sube(other, r0))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn mul_as_doubled(self, other: Self) -> Doubled<Self> {
                let xh = vupper_vd_vd(self);
                let xl = self - xh;
                let yh = vupper_vd_vd(other);
                let yl = other - yh;
                let r0 = self * other;
                Doubled::new(r0, xh * yh + (-r0) + xl * yh + xh * yl + xl * yl)
            }
        }

        impl RecPre for Doubled<$f64x> {
            #[cfg(target_feature = "fma")]
            #[inline]
            fn recpre(self) -> Self {
                let q0 = self.0.recpre();
                Self::new(q0, q0 * self.1.fmanp(q0, self.0.fmanp(q0, $f64x::splat(1))))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn recpre(self) -> Doubled<$f64x> {
                let t = self.0.recpre();
                let dh = vupper_vd_vd(self.0);
                let dl = self.0 - dh;
                let th = vupper_vd_vd(t);
                let tl = t - th;
                let q0 = t;
                Self::new(
                    q0,
                    t * ($f64x::splat(1.) - dh * th - dh * tl - dl * th - dl * tl - self.1 * t),
                )
            }
        }

        impl RecPreAsDoubled for $f64x {
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
                let dh = vupper_vd_vd(self);
                let dl = self - dh;
                let th = vupper_vd_vd(t);
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
    impl_f2_f64!(f64x2, u64x2, m64x2);
}

pub mod f64x4 {
    use packed_simd::*;
    impl_f2_f64!(f64x4, u64x4, m64x4);
}

pub mod f64x8 {
    use packed_simd::*;
    impl_f2_f64!(f64x8, u64x8, m64x8);
}
