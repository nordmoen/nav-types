use ::Access;
use ::enu::ENU;
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

impl<N: Copy + Neg<Output = N>> From<ENU<N>> for NED<N> {
    fn from(e: ENU<N>) -> Self {
        NED::new(e.north(), e.east(), -e.up())
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

impl<N: BaseFloat> Norm<N> for NED<N> {
    fn norm_squared(&self) -> N {
        self.0.norm_squared()
    }
    fn normalize(&self) -> Self {
        NED(self.0.normalize())
    }
    fn normalize_mut(&mut self) -> N {
        self.0.normalize_mut()
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
    use ::enu::ENU;

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

        fn from_enu(e: f32, n: f32, u: f32) -> () {
            let enu = ENU::new(e, n, u);
            let ned = NED::from(enu);
            assert_eq!(ned.north(), n);
            assert_eq!(ned.east(), e);
            assert_eq!(ned.down(), -u);
        }
    }
}
