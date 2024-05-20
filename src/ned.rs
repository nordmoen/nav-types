use crate::enu::ENU;
use crate::utils::RealFieldCopy;
use crate::Access;
use na::Vector3;
use core::convert::From;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// North East Down vector
///
/// This struct represents a vector in the NED coordinate system.
/// See: [NED](https://en.wikipedia.org/wiki/North_east_down) for a general
/// description.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct NED<N: RealFieldCopy>(Vector3<N>);

impl<N: RealFieldCopy> NED<N> {
    /// Create a new NED vector
    pub fn new(n: N, e: N, d: N) -> NED<N> {
        NED(Vector3::new(n, e, d))
    }

    /// Computes the L2 (Euclidean) norm of this vector
    pub fn norm(&self) -> N {
        self.0.norm()
    }
}

impl<N: RealFieldCopy> NED<N> {
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

impl<N: RealFieldCopy + Neg<Output = N>> From<NED<N>> for ENU<N> {
    /// Convert `NED` vectors into `ENU`
    fn from(e: NED<N>) -> Self {
        ENU::new(e.east(), e.north(), -e.down())
    }
}

impl<N: RealFieldCopy + Add<N, Output = N>> Add<NED<N>> for NED<N> {
    type Output = NED<N>;
    fn add(self, right: NED<N>) -> NED<N> {
        NED(self.0 + right.0)
    }
}

impl<N: RealFieldCopy + AddAssign<N>> AddAssign<NED<N>> for NED<N> {
    fn add_assign(&mut self, right: NED<N>) {
        self.0 += right.0
    }
}

impl<N: RealFieldCopy + Sub<N, Output = N>> Sub<NED<N>> for NED<N> {
    type Output = NED<N>;
    fn sub(self, right: NED<N>) -> NED<N> {
        NED(self.0 - right.0)
    }
}

impl<N: RealFieldCopy + SubAssign<N>> SubAssign<NED<N>> for NED<N> {
    fn sub_assign(&mut self, right: NED<N>) {
        self.0 -= right.0
    }
}

impl<N: RealFieldCopy + Mul<N, Output = N>> Mul<N> for NED<N> {
    type Output = NED<N>;
    fn mul(self, right: N) -> NED<N> {
        NED(self.0 * right)
    }
}

impl<N: RealFieldCopy + MulAssign<N>> MulAssign<N> for NED<N> {
    fn mul_assign(&mut self, right: N) {
        self.0 *= right
    }
}

impl<N: RealFieldCopy + Div<N, Output = N>> Div<N> for NED<N> {
    type Output = NED<N>;
    fn div(self, right: N) -> NED<N> {
        NED(self.0 / right)
    }
}

impl<N: RealFieldCopy + DivAssign<N>> DivAssign<N> for NED<N> {
    fn div_assign(&mut self, right: N) {
        self.0 /= right
    }
}

impl<N: RealFieldCopy> Access<Vector3<N>> for NED<N> {
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

        fn from_enu(n: f32, e: f32, u: f32) -> () {
            let enu = ENU::new(e, n, u);
            let ned = NED::from(enu);
            assert_eq!(ned.east(), e);
            assert_eq!(ned.north(), n);
            assert_eq!(ned.down(), -u);
        }
    }
}
