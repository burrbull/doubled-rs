#![deny(warnings)]
#![no_std]

pub const SLEEF_INFINITY: f64 = (1e+300 * 1e+300);
pub const SLEEF_NAN: f64 = SLEEF_INFINITY - SLEEF_INFINITY;
pub const SLEEF_INFINITY_F: f32 = SLEEF_INFINITY as f32;
pub const SLEEF_NAN_F: f32 = SLEEF_NAN as f32;

mod f32;
mod f64;

pub trait IsInf {
    type Mask;
    fn isinf(self) -> Self::Mask;
    fn ispinf(self) -> Self::Mask;
}
pub trait IsNan {
    type Mask;
    fn isnan(self) -> Self::Mask;
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Doubled<T>(pub T, pub T);

impl<T> Doubled<T> where T: Sized {
    #[inline]
    pub fn new(x0: T, x1: T) -> Self {
        Doubled(x0, x1)
    }

}
pub trait Normalize {
    fn normalize(self) -> Self;
}

impl<T> Normalize for Doubled<T>
where
    T: Copy + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    fn normalize(self) -> Self {
        let s0 = self.0 + self.1;
        Self::new(s0, self.0 - s0 + self.1)
    }
}
pub trait Scale<T> {
    fn scale(self, other: T) -> Self;
}

impl<T> Scale<T> for Doubled<T>
where
    T: Copy + core::ops::Mul<Output=T>
{
    #[inline]
    fn scale(self, other: T) -> Self {
        Self::new(self.0 * other, self.1 * other)
    }
}

pub trait Check {
    fn check(self) -> bool;
}

pub trait CheckOrder<T=Self> {
    fn check_order(self, other: T);
}

pub trait AddAsDoubled where Self: Sized {
    fn add_as_doubled(self, other: Self) -> Doubled<Self>;
}
pub trait AddCheckedAsDoubled : CheckOrder<Self> where Self: Sized {
    fn add_checked_as_doubled(self, other: Self) -> Doubled<Self>;
}
pub trait MulAsDoubled where Self: Sized {
    fn mul_as_doubled(self, other: Self) -> Doubled<Self>;
}

pub trait RecPre<T = Self> {
    fn recpre(self) -> T;
}

pub trait RecPreAsDoubled
where
    Self: Sized,
{
    fn recpre_as_doubled(self) -> Doubled<Self>;
}

pub trait AddChecked<T = Self> : CheckOrder<T> {
    type Output;
    fn add_checked(self, other: T) -> Self::Output;
}

pub trait AddCheckedAssign<T> {
    fn add_checked_assign(&mut self, other: T);
}

pub trait SubChecked<T = Self> : CheckOrder<T> {
    type Output;
    fn sub_checked(self, other: T) -> Self::Output;
}

impl<T> AddChecked for Doubled<T>
where
    Doubled<T>: CheckOrder,
    T: Copy + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    type Output = Self;
    #[inline]
    fn add_checked(self, other: Self) -> Self::Output {
        // |self| >= |other|
        self.check_order(other);
        let r0 = self.0 + other.0;
        Self::new(r0, self.0 - r0 + other.0 + self.1 + other.1)
    }
}


impl<T> AddChecked<T> for Doubled<T>
where
    Doubled<T>: CheckOrder<T>,
    T: Copy + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    type Output = Self;
    #[inline]
    fn add_checked(self, other: T) -> Self::Output {
        // |self| >= |other|
        self.check_order(other);
        let r0 = self.0 + other;
        Self::new(r0, self.0 - r0 + other + self.1)
    }
}

impl<T> AddChecked<Doubled<T>> for T
where
    T: Copy + CheckOrder<Doubled<T>> + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    type Output = Doubled<T> ;
    #[inline]
    fn add_checked(self, other: Doubled<T>) -> Self::Output {
        self.check_order(other);
        let r0 = self + other.0;
        Doubled::new(r0, self - r0 + other.0 + other.1)
    }
}

impl<T> AddCheckedAssign<Self> for Doubled<T>
where
    Self: Copy + AddChecked<Output=Self>
{
    #[inline]
    fn add_checked_assign(&mut self, other: Self) {
        *self = (*self).add_checked(other);
    }
}

impl<T> AddCheckedAssign<T> for Doubled<T>
where
    Self: AddChecked<T, Output=Self>,
    T: Copy
{
    #[inline]
    fn add_checked_assign(&mut self, other: T) {
        *self = (*self).add_checked(other);
    }
}

impl<T> SubChecked for Doubled<T>
where
    Doubled<T>: CheckOrder,
    T: Copy + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    type Output = Self;
    #[inline]
    fn sub_checked(self, other: Self) -> Self::Output {
        // |self| >= |other|
        self.check_order(other);
        let r0 = self.0 - other.0;
        Self::new(r0, self.0 - r0 - other.0 + self.1 - other.1)
    }
}

impl<T> SubChecked<T> for Doubled<T>
where
    Doubled<T>: CheckOrder<T>,
    T: Copy + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    type Output = Self;
    #[inline]
    fn sub_checked(self, other: T) -> Self::Output {
        // |self| >= |other|
        self.check_order(other);
        let r0 = self.0 - other;
        Self::new(r0, self.0 - r0 - other + self.1)
    }
}


impl<T> core::ops::Neg for Doubled<T>
where
    T: Copy + core::ops::Neg<Output=T>
{
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.0, -self.1)
    }
}

impl<T> core::ops::Add for Doubled<T>
where
    T: Copy + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    type Output = Self;
    #[inline]
    fn add(self, other: Self) -> Self {
        let r0 = self.0 + other.0;
        let v = r0 - self.0;
        Self::new(r0, self.0 - (r0 - v) + (other.0 - v) + (self.1 + other.1))
    }
}

impl<T> core::ops::AddAssign for Doubled<T>
where
    Self: Copy + core::ops::Add<Output=Self>
{
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl<T> core::ops::Add<T> for Doubled<T>
where
    T: Copy + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    type Output = Self;
    #[inline]
    fn add(self, other: T) -> Self {
        let r0 = self.0 + other;
        let v = r0 - self.0;
        Self::new(r0, self.0 - (r0 - v) + (other - v) + self.1)
    }
}


impl<T> core::ops::AddAssign<T> for Doubled<T>
where
    Self: Copy + core::ops::Add<T, Output=Self>
{
    #[inline]
    fn add_assign(&mut self, other: T) {
        *self = *self + other;
    }
}


impl<T> core::ops::MulAssign for Doubled<T>
where
    Self: Copy + core::ops::Mul<Output=Self>
{
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl<T> core::ops::MulAssign<T> for Doubled<T>
where
    Self: Copy + core::ops::Mul<T, Output=Self>
{
    #[inline]
    fn mul_assign(&mut self, other: T) {
        *self = *self * other;
    }
}


impl<T> AddAsDoubled for T
where
    Self: Copy + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    #[inline]
    fn add_as_doubled(self, other: Self) -> Doubled<Self> {
        let r0 = self + other;
        let v = r0 - self;
        Doubled::new(r0, self - (r0 - v) + (other - v))
    }
}

impl<T> AddCheckedAsDoubled for T
where
    Self: Copy + CheckOrder + core::ops::Add<Output=T> + core::ops::Sub<Output=T>
{
    #[inline]
    fn add_checked_as_doubled(self, other: Self) -> Doubled<Self> {
        self.check_order(other);
        let r0 = self + other;
        Doubled::new(r0, self - r0 + other)
    }
}

pub trait FromU32 {
    fn from_u32(i: (u32, u32)) -> Self;
}

#[cfg(feature = "packed_simd")]
pub mod f32x;
#[cfg(feature = "packed_simd")]
pub mod f64x;