use crate::Access;
use na::{BaseFloat, Norm, Vector3};
use std::convert::Into;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

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
pub struct ENU<N>(Vector3<N>);

impl<N> ENU<N> {
    /// Crate a new ENU vector
    pub fn new(e: N, n: N, u: N) -> ENU<N> {
        ENU(Vector3::new(e, n, u))
    }
}

impl<N: Copy> ENU<N> {
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

impl<N: Copy + Add<N, Output = N>, T: Into<ENU<N>>> Add<T> for ENU<N> {
    type Output = ENU<N>;
    fn add(self, right: T) -> Self::Output {
        ENU(self.0 + right.into().0)
    }
}

impl<N: Copy + AddAssign<N>, T: Into<ENU<N>>> AddAssign<T> for ENU<N> {
    fn add_assign(&mut self, right: T) {
        self.0 += right.into().0
    }
}

impl<N: Copy + Sub<N, Output = N>, T: Into<ENU<N>>> Sub<T> for ENU<N> {
    type Output = ENU<N>;
    fn sub(self, right: T) -> Self::Output {
        ENU(self.0 - right.into().0)
    }
}

impl<N: Copy + SubAssign<N>, T: Into<ENU<N>>> SubAssign<T> for ENU<N> {
    fn sub_assign(&mut self, right: T) {
        self.0 -= right.into().0
    }
}

impl<N: Copy + Mul<N, Output = N>> Mul<N> for ENU<N> {
    type Output = ENU<N>;
    fn mul(self, right: N) -> Self::Output {
        ENU(self.0 * right)
    }
}

impl<N: Copy + MulAssign<N>> MulAssign<N> for ENU<N> {
    fn mul_assign(&mut self, right: N) {
        self.0 *= right
    }
}

impl<N: Copy + Div<N, Output = N>> Div<N> for ENU<N> {
    type Output = ENU<N>;
    fn div(self, right: N) -> Self::Output {
        ENU(self.0 / right)
    }
}

impl<N: Copy + DivAssign<N>> DivAssign<N> for ENU<N> {
    fn div_assign(&mut self, right: N) {
        self.0 /= right
    }
}

impl<N: BaseFloat> Norm for ENU<N> {
    type NormType = N;

    fn norm_squared(&self) -> N {
        self.0.norm_squared()
    }
    fn normalize(&self) -> Self {
        ENU(self.0.normalize())
    }
    fn normalize_mut(&mut self) -> N {
        self.0.normalize_mut()
    }
    fn try_normalize(&self, min_norm: Self::NormType) -> Option<Self> {
        self.0.try_normalize(min_norm).map(ENU)
    }
    fn try_normalize_mut(&mut self, min_norm: Self::NormType) -> Option<Self::NormType> {
        self.0.try_normalize_mut(min_norm)
    }
}

impl<N> Access<Vector3<N>> for ENU<N> {
    fn access(self) -> Vector3<N> {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ned::NED;

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
