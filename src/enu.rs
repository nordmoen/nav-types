use crate::ned::NED;
use crate::utils::RealFieldCopy;
use crate::Access;
use na::Vector3;
use core::convert::Into;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// East North Up vector
///
/// This struct represents a vector in the ENU coordinate system.
/// See: [ENU](https://en.wikipedia.org/wiki/Axes_conventions) for a general
/// description.
///
/// # Note
/// ENU implements `Into` on all operations. This means that any vector that
/// can be converted to ENU with a `From` implementation will automatically
/// be able to work with ENU at the cost of the `Into` conversion.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct ENU<N: RealFieldCopy>(Vector3<N>);

impl<N: RealFieldCopy> ENU<N> {
    /// Create a new ENU vector
    pub fn new(e: N, n: N, u: N) -> ENU<N> {
        ENU(Vector3::new(e, n, u))
    }

    /// Computes the L2 (Euclidean) norm of this vector
    pub fn norm(&self) -> N {
        self.0.norm()
    }
}

impl<N: RealFieldCopy> ENU<N> {
    /// Get the East component of this vector
    pub fn east(&self) -> N {
        self.0.x
    }

    /// Get the North component of this vector
    pub fn north(&self) -> N {
        self.0.y
    }

    /// Get the Up component of this vector
    pub fn up(&self) -> N {
        self.0.z
    }
}

impl<N: RealFieldCopy + Neg<Output = N>> From<ENU<N>> for NED<N> {
    /// Convert `ENU` vectors into `NED`
    fn from(e: ENU<N>) -> Self {
        NED::new(e.north(), e.east(), -e.up())
    }
}


impl<N: RealFieldCopy + Add<N, Output = N>, T: Into<ENU<N>>> Add<T> for ENU<N> {
    type Output = ENU<N>;
    fn add(self, right: T) -> Self::Output {
        ENU(self.0 + right.into().0)
    }
}

impl<N: RealFieldCopy + AddAssign<N>, T: Into<ENU<N>>> AddAssign<T> for ENU<N> {
    fn add_assign(&mut self, right: T) {
        self.0 += right.into().0
    }
}

impl<N: RealFieldCopy + Sub<N, Output = N>, T: Into<ENU<N>>> Sub<T> for ENU<N> {
    type Output = ENU<N>;
    fn sub(self, right: T) -> Self::Output {
        ENU(self.0 - right.into().0)
    }
}

impl<N: RealFieldCopy + SubAssign<N>, T: Into<ENU<N>>> SubAssign<T> for ENU<N> {
    fn sub_assign(&mut self, right: T) {
        self.0 -= right.into().0
    }
}

impl<N: RealFieldCopy + Mul<N, Output = N>> Mul<N> for ENU<N> {
    type Output = ENU<N>;
    fn mul(self, right: N) -> Self::Output {
        ENU(self.0 * right)
    }
}

impl<N: RealFieldCopy + MulAssign<N>> MulAssign<N> for ENU<N> {
    fn mul_assign(&mut self, right: N) {
        self.0 *= right
    }
}

impl<N: RealFieldCopy + Div<N, Output = N>> Div<N> for ENU<N> {
    type Output = ENU<N>;
    fn div(self, right: N) -> Self::Output {
        ENU(self.0 / right)
    }
}

impl<N: RealFieldCopy + DivAssign<N>> DivAssign<N> for ENU<N> {
    fn div_assign(&mut self, right: N) {
        self.0 /= right
    }
}

impl<N: RealFieldCopy> Access<Vector3<N>> for ENU<N> {
    fn access(self) -> Vector3<N> {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ned::NED;

    quickcheck! {
        fn create_enu(e: f32, n: f32, u: f32) -> () {
            ENU::new(e, n, u);
        }

        fn get_components(e: f32, n: f32, u: f32) -> () {
            let vec = ENU::new(e, n, u);
            assert_eq!(vec.east(), e);
            assert_eq!(vec.north(), n);
            assert_eq!(vec.up(), u);
        }

        fn from_ned(n: f32, e: f32, d: f32) -> () {
            let ned = NED::new(n, e, d);
            let enu = ENU::from(ned);
            assert_eq!(enu.east(), e);
            assert_eq!(enu.north(), n);
            assert_eq!(enu.up(), -d);
        }
    }
}
