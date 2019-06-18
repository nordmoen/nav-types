use crate::enu::ENU;
use crate::Access;
use na::{BaseFloat, Norm, Vector3};
use std::convert::From;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// North East Down vector
///
/// This struct represents a vector in the NED coordinate system.
/// See: [NED](https://en.wikipedia.org/wiki/North_east_down) for a general
/// description.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct NED<N>(Vector3<N>);

impl<N> NED<N> {
    /// Create a new NED vector
    pub fn new(n: N, e: N, d: N) -> NED<N> {
        NED(Vector3::new(n, e, d))
    }
}

impl<N: Copy> NED<N> {
    /// Get the North component of this vector
    pub fn north(&self) -> N {
        self.0.x
    }

    /// Get the East component of this vector
    pub fn east(&self) -> N {
        self.0.y
    }

    /// Get the Down component of this vector
    pub fn down(&self) -> N {
        self.0.z
    }
}

impl<N: Copy + Neg<Output = N>> From<NED<N>> for ENU<N> {
    /// Convert `NED` vectors into `ENU`
    fn from(e: NED<N>) -> Self {
        ENU::new(e.east(), e.north(), -e.down())
    }
}

impl<N: Copy + Add<N, Output = N>> Add<NED<N>> for NED<N> {
    type Output = NED<N>;
    fn add(self, right: NED<N>) -> NED<N> {
        NED(self.0 + right.0)
    }
}

impl<N: Copy + AddAssign<N>> AddAssign<NED<N>> for NED<N> {
    fn add_assign(&mut self, right: NED<N>) {
        self.0 += right.0
    }
}

impl<N: Copy + Sub<N, Output = N>> Sub<NED<N>> for NED<N> {
    type Output = NED<N>;
    fn sub(self, right: NED<N>) -> NED<N> {
        NED(self.0 - right.0)
    }
}

impl<N: Copy + SubAssign<N>> SubAssign<NED<N>> for NED<N> {
    fn sub_assign(&mut self, right: NED<N>) {
        self.0 -= right.0
    }
}

impl<N: Copy + Mul<N, Output = N>> Mul<N> for NED<N> {
    type Output = NED<N>;
    fn mul(self, right: N) -> NED<N> {
        NED(self.0 * right)
    }
}

impl<N: Copy + MulAssign<N>> MulAssign<N> for NED<N> {
    fn mul_assign(&mut self, right: N) {
        self.0 *= right
    }
}

impl<N: Copy + Div<N, Output = N>> Div<N> for NED<N> {
    type Output = NED<N>;
    fn div(self, right: N) -> NED<N> {
        NED(self.0 / right)
    }
}

impl<N: Copy + DivAssign<N>> DivAssign<N> for NED<N> {
    fn div_assign(&mut self, right: N) {
        self.0 /= right
    }
}

impl<N: BaseFloat> Norm for NED<N> {
    type NormType = N;

    fn norm_squared(&self) -> N {
        self.0.norm_squared()
    }
    fn normalize(&self) -> Self {
        NED(self.0.normalize())
    }
    fn normalize_mut(&mut self) -> N {
        self.0.normalize_mut()
    }
    fn try_normalize(&self, min_norm: Self::NormType) -> Option<Self> {
        self.0.try_normalize(min_norm).map(NED)
    }
    fn try_normalize_mut(&mut self, min_norm: Self::NormType) -> Option<Self::NormType> {
        self.0.try_normalize_mut(min_norm)
    }
}

impl<N> Access<Vector3<N>> for NED<N> {
    fn access(self) -> Vector3<N> {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enu::ENU;

    quickcheck! {
        fn create_ned(n: f32, e: f32, d: f32) -> () {
            NED::new(n, e, d);
        }

        fn get_components(n: f32, e: f32, d: f32) -> () {
            let vec = NED::new(n, e, d);
            assert_eq!(vec.north(), n);
            assert_eq!(vec.east(), e);
            assert_eq!(vec.down(), d);
        }

        fn into_enu(n: f32, e: f32, d: f32) -> () {
            let ned = NED::new(n, e, d);
            let enu: ENU<_> = ned.into();
            assert_eq!(n, enu.north());
            assert_eq!(e, enu.east());
            assert_eq!(d, -enu.up());
        }

        fn add_enu(n: f32, e: f32, d: f32) -> () {
            let ned = NED::new(n, e, d);
            let enu = ENU::new(e, n, -d);
            let sum = enu + ned;
            let twi = ned * 2.0;
            assert_eq!(sum.north(), twi.north());
            assert_eq!(sum.east(), twi.east());
            assert_eq!(sum.up(), -twi.down());
        }
    }
}
