use ::Access;
use ::ned::NED;
use na::{Vector3, Norm, BaseFloat};
use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign,
    Neg};

/// East North Up vector
///
/// This struct represents a vector in the ENU coordinate system.
/// See: [ENU](https://en.wikipedia.org/wiki/Axes_conventions) for a general
/// description.
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

impl<N: Copy + Neg<Output=N>> From<NED<N>> for ENU<N> {
    fn from(n: NED<N>) -> ENU<N> {
        ENU::new(n.east(), n.north(), -n.down())
    }
}

impl<N: Copy + Add<N, Output=N>> Add<ENU<N>> for ENU<N> {
    type Output = ENU<N>;
    fn add(self, right: ENU<N>) -> ENU<N> {
        ENU(self.0 + right.0)
    }
}

impl<N: Copy + AddAssign<N>> AddAssign<ENU<N>> for ENU<N> {
    fn add_assign(&mut self, right: ENU<N>) {
        self.0 += right.0
    }
}

impl<N: Copy + Sub<N, Output=N>> Sub<ENU<N>> for ENU<N> {
    type Output = ENU<N>;
    fn sub(self, right: ENU<N>) -> ENU<N> {
        ENU(self.0 - right.0)
    }
}

impl<N: Copy + SubAssign<N>> SubAssign<ENU<N>> for ENU<N> {
    fn sub_assign(&mut self, right: ENU<N>) {
        self.0 -= right.0
    }
}

impl<N: Copy + Mul<N, Output=N>> Mul<N> for ENU<N> {
    type Output = ENU<N>;
    fn mul(self, right: N) -> ENU<N> {
        ENU(self.0 * right)
    }
}

impl<N: Copy + MulAssign<N>> MulAssign<N> for ENU<N> {
    fn mul_assign(&mut self, right: N) {
        self.0 *= right
    }
}

impl<N: Copy + Div<N, Output=N>> Div<N> for ENU<N> {
    type Output = ENU<N>;
    fn div(self, right: N) -> Self {
        ENU(self.0 / right)
    }
}

impl<N: Copy + DivAssign<N>> DivAssign<N> for ENU<N> {
    fn div_assign(&mut self, right: N) {
        self.0 /= right
    }
}

impl<N: BaseFloat> Norm<N> for ENU<N> {
    fn norm_squared(&self) -> N {
        self.0.norm_squared()
    }
    fn normalize(&self) -> Self {
        ENU(self.0.normalize())
    }
    fn normalize_mut(&mut self) -> N {
        self.0.normalize_mut()
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
    use ::ned::NED;

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
