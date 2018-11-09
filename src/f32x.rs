macro_rules! impl_f2_f32 {
    ($f32x:ident, $u32x:ident, $m32x:ident) => {
        use crate::*;

        #[inline]
        fn vupper_vf_vf(d: $f32x) -> $f32x {
            ($u32x::from_bits(d) & $u32x::splat(0xfffff000)).into_bits()
        }
        
        impl IsInf for $f32x {
          type Mask = $m32x;
          #[inline]
          fn isinf(self) -> Self::Mask {
             self.abs().eq(Self::splat(SLEEF_INFINITY_F))
          }
          #[inline]
          fn ispinf(self) -> Self::Mask {
            self.eq(Self::splat(SLEEF_INFINITY_F))
          }
        }
        impl IsNan for $f32x {
          type Mask = $m32x;
          #[inline]
          fn isnan(self) -> Self::Mask {
            self.ne(self)
          }
        }

        impl Doubled<$f32x> {
            #[inline]
            pub fn abs(self) -> Self {
                Self::new(
                    self.0.abs(),
                    ($u32x::from_bits(self.1)
                            ^ ($u32x::from_bits(self.0) & $u32x::from_bits($f32x::splat(-0.)))
                    ).into_bits()
                )
            }

            #[cfg(target_feature = "fma")]
            #[inline]
            pub fn square(self) -> Self {
                let r0 = self.0 * self.0;
                Self::new(r0, (self.0 + self.0).mul_adde(self.1, self.0.mul_sube(self.0, r0)))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            pub fn square(self) -> Self {
                let xh = vupper_vf_vf(self.0);
                let xl = self.0 - xh;
                let r0 = self.0 * self.0;

                let mut t = xh.mul_add(xh, -r0);
                t = (xh + xh).mul_add(xl, t);
                t = xl.mul_add(xl, t);
                t = self.0.mul_add(self.1 + self.1, t);
                Self::new(r0, t)
            }

            #[cfg(target_feature = "fma")]
            #[inline]
            pub fn square_as_f(self) -> $f32x {
                self.0.mul_adde(self.0, self.0 * self.1 + self.0 * self.1)
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            pub fn square_as_f(self) -> $f32x {
                let xh = vupper_vf_vf(self.0);
                let xl = self.0 - xh;
                xh * self.1 + xh * self.1 + xl * xl + (xh * xl + xh * xl) + xh * xh
            }

            #[cfg(feature = "enable_recsqrt_sp")]
            #[inline]
            pub fn sqrt(self) -> Self {
                let x = (self.0 + self.1).rsqrte();
                let r = self * x;
                (r * (r * x + $f32x::splat(-3.0))).scale($f32x::splat(-0.5))
            }
            #[cfg(not(feature = "enable_recsqrt_sp"))]
            #[inline]
            pub fn sqrt(self) -> Self {
                let t = (self.0 + self.1).sqrt();
                ((self + t.mul_as_doubled(t)) * t.recpre()).scale($f32x::splat(0.5))
            }

            #[cfg(target_feature = "fma")]
            #[inline]
            pub fn mul_as_f(self, other: Self) -> $f32x {
                self.0.mul_adde(other.0, self.1.mul_adde(other.0, self.0 * other.1))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            pub fn mul_as_f(self, other: Self) -> $f32x {
                let xh = vupper_vf_vf(self.0);
                let xl = self.0 - xh;
                let yh = vupper_vf_vf(other.0);
                let yl = other.0 - yh;
                self.1 * yh + xh * other.1 + xl * yl + xh * yl + xl * yh + xh * yh
            }
        }

        impl core::convert::From<f64> for Doubled<$f32x> {
            fn from(d: f64) -> Self {
                Self::new(
                    $f32x::splat(d as f32),
                    $f32x::splat((d as f32) - (d as f32)),
                )
            }
        }

        impl core::convert::From<(f32, f32)> for Doubled<$f32x> {
            fn from(f: (f32, f32)) -> Self {
                Self::new($f32x::splat(f.0), $f32x::splat(f.1))
            }
        }

        impl core::ops::Add<Doubled<$f32x>> for $f32x {
            type Output = Doubled<$f32x> ;
            #[inline]
            fn add(self, other: Doubled<$f32x>) -> Self::Output {
                let r0 = self + other.0;
                let v = r0 - self;
                Doubled::new(r0, self - (r0 - v) + (other.0 - v) + other.1)
            }
        }

        impl core::ops::Mul for Doubled<$f32x> {
            type Output = Self;
            #[cfg(target_feature = "fma")]
            #[inline]
            fn mul(self, other: Self) -> Self {
                let r0 = self.0 * other.0;
                Self::new(
                    r0,
                    self.0
                        .mul_adde(other.1, self.1.mul_adde(other.0, self.0.mul_sube(other.0, r0))),
                )
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn mul(self, other: Self) -> Self {
                let xh = vupper_vf_vf(self.0);
                let xl = self.0 - xh;
                let yh = vupper_vf_vf(other.0);
                let yl = other.0 - yh;
                let r0 = self.0 * other.0;

                let mut t = xh.mul_add(yh, -r0);
                t = xl.mul_add(yh, t);
                t = xh.mul_add(yl, t);
                t = xl.mul_add(yl, t);
                t = self.0.mul_add(other.1, t);
                t = self.1.mul_add(other.0, t);
                Self::new(r0, t)
            }
        }

        impl core::ops::Mul<$f32x> for Doubled<$f32x> {
            type Output = Self;
            #[cfg(target_feature = "fma")]
            #[inline]
            fn mul(self, other: $f32x) -> Self {
                let r0 = self.0 * other;
                Self::new(r0, self.1.mul_adde(other, self.0.fmanp(other, r0)))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn mul(self, other: $f32x) -> Self {
                let xh = vupper_vf_vf(self.0);
                let xl = self.0 - xh;
                let yh = vupper_vf_vf(other);
                let yl = other - yh;
                let r0 = self.0 * other;

                let mut t = xh.mul_add(yh, -r0);
                t = xl.mul_add(yh, t);
                t = xh.mul_add(yl, t);
                t = xl.mul_add(yl, t);
                t = self.1.mul_add(other, t);
                Self::new(r0, t)
            }
        }

        impl core::ops::Div for Doubled<$f32x> {
            type Output = Self;
            #[cfg(target_feature = "fma")]
            #[inline]
            fn div(self, other: Self) -> Self {
                let t = other.0.recpre();

                let q0 = self.0 * t;
                let u = t.mul_sube(self.0, q0);
                let mut q1 = other.1.fmanp(t, other.0.fmanp(t, $f32x::splat(1)));
                q1 = q0.mul_adde(q1, self.1.mul_adde(t, u));

                Self::new(q0, q1)
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn div(self, other: Self) -> Self {
                let t = other.0.recpre();
                let dh = vupper_vf_vf(other.0);
                let dl = other.0 - dh;
                let th = vupper_vf_vf(t);
                let tl = t - th;
                let nhh = vupper_vf_vf(self.0);
                let nhl = self.0 - nhh;

                let q0 = self.0 * t;

                let mut w = $f32x::splat(-1.);
                w = dh.mul_add(th, w);
                w = dh.mul_add(tl, w);
                w = dl.mul_add(th, w);
                w = dl.mul_add(tl, w);
                w = -w;

                let mut u = nhh.mul_add(th, -q0);
                u = nhh.mul_add(tl, u);
                u = nhl.mul_add(th, u);
                u = nhl.mul_add(tl, u);
                u = q0.mul_add(w, u);

                Self::new(q0, t.mul_add(self.1 - q0 * other.1, u))
            }
        }

        impl CheckOrder for Doubled<$f32x> {
            fn check_order(self, _other: Self) {
            }
        }

        impl CheckOrder<$f32x> for Doubled<$f32x> {
            fn check_order(self, _other: $f32x) {
            }
        }

        impl CheckOrder<Doubled<$f32x>> for $f32x {
            fn check_order(self, _other: Doubled<$f32x>) {
            }
        }

        impl CheckOrder for $f32x {
            fn check_order(self, _other: Self) {
            }
        }
        impl MulAsDoubled for $f32x {
            #[cfg(target_feature = "fma")]
            #[inline]
            fn mul_as_doubled(self, other: Self) -> Doubled<Self> {
                let r0 = self * other;
                Self::new(r0, self.mul_sube(other, r0))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn mul_as_doubled(self, other: Self) -> Doubled<Self> {
                let xh = vupper_vf_vf(self);
                let xl = self - xh;
                let yh = vupper_vf_vf(other);
                let yl = other - yh;
                let r0 = self * other;

                let mut t = xh.mul_add(yh, -r0);
                t = xl.mul_add(yh, t);
                t = xh.mul_add(yl, t);
                t = xl.mul_add(yl, t);
                Doubled::new(r0, t)
            }
        }

        impl RecPre for Doubled<$f32x> {
            #[cfg(target_feature = "fma")]
            #[inline]
            fn recpre(self) -> Self {
                let q0 = self.0.recpre();
                Self::new(q0, q0 * self.1.fmanp(q0, self.0.fmanp(q0, $f32x::splat(1))))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn recpre(self) -> Self {
                let t = self.0.recpre();
                let dh = vupper_vf_vf(self.0);
                let dl = self.0 - dh;
                let th = vupper_vf_vf(t);
                let tl = t - th;
                let q0 = t;

                let mut u = $f32x::splat(-1.);
                u = dh.mul_add(th, u);
                u = dh.mul_add(tl, u);
                u = dl.mul_add(th, u);
                u = dl.mul_add(tl, u);
                u = self.1.mul_add(t, u);
                Self::new(q0, (-t) * u)
            }
        }

        impl RecPreAsDoubled for $f32x {
            #[cfg(target_feature = "fma")]
            #[inline]
            fn recpre_as_doubled(self) -> Doubled<Self>  {
                let q0 = self.recpre();
                Doubled::new(q0, q0 * self.fmanp(q0, Self::splat(1)))
            }
            #[cfg(not(target_feature = "fma"))]
            #[inline]
            fn recpre_as_doubled(self) -> Doubled<Self> {
                let t = self.recpre();
                let dh = vupper_vf_vf(self);
                let dl = self - dh;
                let th = vupper_vf_vf(t);
                let tl = t - th;
                let q0 = t;

                let mut u = Self::splat(-1.);
                u = dh.mul_add(th, u);
                u = dh.mul_add(tl, u);
                u = dl.mul_add(th, u);
                u = dl.mul_add(tl, u);
                Doubled::new(q0, (-t) * u)
            }
        }
    };
}

pub mod f32x2 {
    use packed_simd::*;
    impl_f2_f32!(f32x2, u32x2, m32x2);
}

pub mod f32x4 {
    use packed_simd::*;
    impl_f2_f32!(f32x4, u32x4, m32x4);
}

pub mod f32x8 {
    use packed_simd::*;
    impl_f2_f32!(f32x8, u32x8, m32x8);
}

pub mod f32x16 {
    use packed_simd::*;
    impl_f2_f32!(f32x16, u32x16, m32x16);
}
