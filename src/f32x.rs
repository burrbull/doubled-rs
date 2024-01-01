use crate::*;
use core::simd::{num::SimdFloat, LaneCount, Simd, SupportedLaneCount};
use std::simd::StdFloat;

type F32x<const N: usize> = Simd<f32, N>;
type U32x<const N: usize> = Simd<u32, N>;

impl<const N: usize> Upper for F32x<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn upper(self) -> Self {
        Self::from_bits(self.to_bits() & U32x::splat(0x_fff_ff000))
    }
}

impl<const N: usize> FromMask for Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    type Mask = U32x<N>;
    fn from_mask(u0: Self::Mask, u1: Self::Mask) -> Self {
        Self::new(F32x::from_bits(u0), F32x::from_bits(u1))
    }
}

impl<const N: usize> core::convert::From<F32x<N>> for Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn from(f: F32x<N>) -> Self {
        Self::new(f, F32x::splat(0.))
    }
}

impl<const N: usize> core::convert::From<Doubled<F32x<N>>> for F32x<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn from(f: Doubled<F32x<N>>) -> Self {
        f.0 + f.1
    }
}

impl<const N: usize> Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    pub const fn splat(value: Doubled<f32>) -> Self {
        Self::new(
            F32x::from_array([value.0; N]),
            F32x::from_array([value.1; N]),
        )
    }

    #[inline]
    pub fn abs(self) -> Self {
        Self::new(
            self.0.abs(),
            F32x::from_bits(self.1.to_bits() ^ (self.0.to_bits() & F32x::splat(-0.).to_bits())),
        )
    }

    #[cfg(target_feature = "fma")]
    #[inline]
    pub fn square(self) -> Self {
        let r0 = self.0 * self.0;
        Self::new(
            r0,
            (self.0 + self.0).mul_add(self.1, self.0.mul_add(self.0, -r0)),
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
    pub fn square_as_f(self) -> F32x<N> {
        self.0.mul_add(self.0, self.0 * self.1 + self.0 * self.1)
    }
    #[cfg(not(target_feature = "fma"))]
    #[inline]
    pub fn square_as_f(self) -> F32x<N> {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        xh * self.1 + xh * self.1 + xl * xl + (xh * xl + xh * xl) + xh * xh
    }

    #[cfg(feature = "enable_recsqrt_sp")]
    #[inline]
    pub fn sqrt(self) -> Self {
        let x = (self.0 + self.1).rsqrte();
        let r = self * x;
        (r * (r * x + F32x::splat(-3.0))).scale(F32x::splat(-0.5))
    }
    #[cfg(not(feature = "enable_recsqrt_sp"))]
    #[inline]
    pub fn sqrt(self) -> Self {
        let t = (self.0 + self.1).sqrt();
        ((self + t.mul_as_doubled(t)) * t.recip_as_doubled()).scale(F32x::splat(0.5))
    }

    #[cfg(target_feature = "fma")]
    #[inline]
    pub fn mul_as_f(self, other: Self) -> F32x<N> {
        self.0
            .mul_add(other.0, self.1.mul_add(other.0, self.0 * other.1))
    }
    #[cfg(not(target_feature = "fma"))]
    #[inline]
    pub fn mul_as_f(self, other: Self) -> F32x<N> {
        let xh = self.0.upper();
        let xl = self.0 - xh;
        let yh = other.0.upper();
        let yl = other.0 - yh;
        self.1 * yh + xh * other.1 + xl * yl + xh * yl + xl * yh + xh * yh
    }
    #[cfg(target_feature = "fma")]
    #[inline]
    pub fn recip(self) -> Self {
        let q0 = self.0.recip();
        Self::new(
            q0,
            q0 * (-self.1).mul_add(q0, (-self.0).mul_add(q0, F32x::splat(1.))),
        )
    }
    #[cfg(not(target_feature = "fma"))]
    #[inline]
    pub fn recip(self) -> Self {
        let t = self.0.recip();
        let dh = self.0.upper();
        let dl = self.0 - dh;
        let th = t.upper();
        let tl = t - th;
        let q0 = t;
        Self::new(
            q0,
            t * (F32x::splat(1.) - dh * th - dh * tl - dl * th - dl * tl - self.1 * t),
        )
    }
}

impl<const N: usize> core::convert::From<f64> for Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn from(d: f64) -> Self {
        Self::splat(Doubled::<f32>::from_f64(d))
    }
}

impl<const N: usize> core::ops::Add<Doubled<F32x<N>>> for F32x<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    type Output = Doubled<F32x<N>>;
    #[inline]
    fn add(self, other: Doubled<F32x<N>>) -> Self::Output {
        let r0 = self + other.0;
        let v = r0 - self;
        Doubled::new(r0, self - (r0 - v) + (other.0 - v) + other.1)
    }
}

impl<const N: usize> core::ops::Mul for Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    type Output = Self;
    #[cfg(target_feature = "fma")]
    #[inline]
    fn mul(self, other: Self) -> Self {
        let r0 = self.0 * other.0;
        Self::new(
            r0,
            self.0.mul_add(
                other.1,
                self.1.mul_add(other.0, self.0.mul_add(other.0, -r0)),
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
            xh * yh + (-r0) + xl * yh + xh * yl + xl * yl + self.0 * other.1 + self.1 * other.0,
        )
    }
}

impl<const N: usize> core::ops::Mul<F32x<N>> for Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    type Output = Self;
    #[cfg(target_feature = "fma")]
    #[inline]
    fn mul(self, other: F32x<N>) -> Self {
        let r0 = self.0 * other;
        Self::new(r0, self.1.mul_add(other, self.0.mul_add(other, -r0)))
    }
    #[cfg(not(target_feature = "fma"))]
    #[inline]
    fn mul(self, other: F32x<N>) -> Self {
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

impl<const N: usize> core::ops::Div for Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    type Output = Self;
    #[cfg(target_feature = "fma")]
    #[inline]
    fn div(self, other: Self) -> Self {
        let t = other.0.recip();

        let q0 = self.0 * t;
        let u = t.mul_add(self.0, -q0);
        let mut q1 = (-other.1).mul_add(t, (-other.0).mul_add(t, F32x::splat(1.)));
        q1 = q0.mul_add(q1, self.1.mul_add(t, u));

        Self::new(q0, q1)
    }
    #[cfg(not(target_feature = "fma"))]
    #[inline]
    fn div(self, other: Self) -> Self {
        let t = other.0.recip();
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
            + q0 * (F32x::splat(1.) - dh * th - dh * tl - dl * th - dl * tl);

        Self::new(q0, t * (self.1 - q0 * other.1) + u)
    }
}

impl<const N: usize> CheckOrder for Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn check_order(self, _other: Self) {}
}

impl<const N: usize> CheckOrder<F32x<N>> for Doubled<F32x<N>>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn check_order(self, _other: F32x<N>) {}
}

impl<const N: usize> CheckOrder<Doubled<F32x<N>>> for F32x<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn check_order(self, _other: Doubled<F32x<N>>) {}
}

impl<const N: usize> CheckOrder for F32x<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn check_order(self, _other: Self) {}
}

impl<const N: usize> AsDoubled for F32x<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn as_doubled(self) -> Doubled<Self> {
        Doubled::new(self, Self::splat(0.))
    }
}

impl<const N: usize> MulAsDoubled for F32x<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[cfg(target_feature = "fma")]
    #[inline]
    fn mul_as_doubled(self, other: Self) -> Doubled<Self> {
        let r0 = self * other;
        Doubled::new(r0, self.mul_add(other, -r0))
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

impl<const N: usize> RecipAsDoubled for F32x<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[cfg(target_feature = "fma")]
    #[inline]
    fn recip_as_doubled(self) -> Doubled<Self> {
        let q0 = self.recip();
        Doubled::new(q0, q0 * (-self).mul_add(q0, Self::splat(1.)))
    }
    #[cfg(not(target_feature = "fma"))]
    #[inline]
    fn recip_as_doubled(self) -> Doubled<Self> {
        let t = self.recip();
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
