use ::Access;
use ::enu::ENU;
use ::nvector::NVector;
use ::wgs84::{ECCENTRICITY_SQ, SEMI_MAJOR_AXIS, SEMI_MINOR_AXIS, WGS84};
use na::{Matrix3, Point3, Transpose};
use num_traits::Float;
use std::convert::From;
use std::ops::{Add, AddAssign, Sub, SubAssign};

/// Earth Centered Earth Fixed position
///
/// This struct represents a position in the ECEF coordinate system.
/// See: [ECEF](https://en.wikipedia.org/wiki/ECEF) for general description.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct ECEF<N>(Point3<N>);

impl<N> ECEF<N> {
    /// Create a new ECEF position
    pub fn new(x: N, y: N, z: N) -> ECEF<N> {
        ECEF(Point3::new(x, y, z))
    }
}

impl<N: Copy> ECEF<N> {
    /// Get the X component of this position
    pub fn x(&self) -> N {
        self.0.x
    }

    /// Get the Y component of this position
    pub fn y(&self) -> N {
        self.0.y
    }

    /// Get the Z component of this position
    pub fn z(&self) -> N {
        self.0.z
    }
}

impl<N: Float> ECEF<N> {
    /// Create a rotation matrix from ECEF frame to ENU frame
    fn r_en(self) -> Matrix3<N> {
        let wgs = WGS84::from(self);
        Matrix3::new(-wgs.longitude().sin(),
                     wgs.longitude().cos(),
                     N::zero(),
                     -wgs.latitude().sin() * wgs.longitude().cos(),
                     -wgs.latitude().sin() * wgs.longitude().sin(),
                     wgs.latitude().cos(),
                     wgs.latitude().cos() * wgs.longitude().cos(),
                     wgs.latitude().cos() * wgs.longitude().sin(),
                     wgs.latitude().sin())

    }

    /// Create a rotation matrix from ENU frame to ECEF frame
    fn r_ne(self) -> Matrix3<N> {
        self.r_en().transpose()
    }
}

impl<N: Float> Add<ENU<N>> for ECEF<N> {
    type Output = ECEF<N>;
    fn add(self, right: ENU<N>) -> ECEF<N> {
        let pos = self.0 + self.r_en() * right.access();
        ECEF(pos)
    }
}

impl<N: Float + AddAssign> AddAssign<ENU<N>> for ECEF<N> {
    fn add_assign(&mut self, right: ENU<N>) {
        self.0 += self.r_en() * right.access();
    }
}

impl<N: Float> Sub<ECEF<N>> for ECEF<N> {
    type Output = ENU<N>;
    fn sub(self, right: ECEF<N>) -> ENU<N> {
        let vec = self.0 - right.0;
        let mat = right.r_ne();
        let enu = mat * vec;
        ENU::new(enu.x, enu.y, enu.z)
    }
}

impl<N: Float> Sub<ENU<N>> for ECEF<N> {
    type Output = ECEF<N>;
    fn sub(self, right: ENU<N>) -> ECEF<N> {
        let vec_e = self.r_en() * right.access();
        ECEF(self.0 - vec_e)
    }
}

impl<N: Float + SubAssign> SubAssign<ENU<N>> for ECEF<N> {
    fn sub_assign(&mut self, right: ENU<N>) {
        let vec_e = self.r_en() * right.access();
        self.0 -= vec_e;
    }
}

// This macro implements most standard operations for position types that
// can be converted to ECEF
macro_rules! ecef_impl {
    ($T:ident) => (
        impl<N: Float> Add<ENU<N>> for $T<N> {
            type Output = $T<N>;
            fn add(self, right: ENU<N>) -> Self {
                $T::from(ECEF::from(self) + right)
            }
        }

        impl<N: Float> AddAssign<ENU<N>> for $T<N> {
            fn add_assign(&mut self, right: ENU<N>) {
                *self = $T::from(ECEF::from(*self) + right);
            }
        }

        impl<N: Float> Sub<ENU<N>> for $T<N> {
            type Output = $T<N>;
            fn sub(self, right: ENU<N>) -> $T<N> {
                $T::from(ECEF::from(self) - right)
            }
        }

        impl<N: Float> Sub<$T<N>> for $T<N> {
            type Output = ENU<N>;
            fn sub(self, right: $T<N>) -> ENU<N> {
                ECEF::from(self) - ECEF::from(right)
            }
        }

        impl<N: Float> SubAssign<ENU<N>> for $T<N> {
            fn sub_assign(&mut self, right: ENU<N>) {
                *self = $T::from(ECEF::from(*self) - right);
            }
        }
    )
}

// In accordance with Gade(2010) all N-Vector operations works through
// converting to ECEF and then back
ecef_impl!(NVector);
// The only way to work with Latitude/Longitude is to convert to ECEF
ecef_impl!(WGS84);

impl<N: Float> From<WGS84<N>> for ECEF<N> {
    fn from(f: WGS84<N>) -> ECEF<N> {
        // Conversion from:
        // http://download.springer.com/static/pdf/723/art%253A10.1007%252Fs00190-004-0375-4
        // .pdf?originUrl=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1007%2Fs00190-004
        // -0375-4&token2=exp=1470729285~acl=%2Fstatic%2Fpdf%2F723%2Fart%25253A10.1007
        // %25252Fs00190-004-0375-4.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com
        // %252Farticle%252F10.1007%252Fs00190-004-0375-4*~hmac=4aaedb6b71f13fc9a9ce5175b
        // 4538c3ec38ddf11b77531d3bd2af75ee1fc2061
        let a = N::from(SEMI_MAJOR_AXIS).unwrap();
        let ecc_part = N::from(ECCENTRICITY_SQ).unwrap();
        let sin_part = N::from(0.5).unwrap() *
                       (N::from(1.0).unwrap() - (N::from(2.0).unwrap() * f.latitude()).cos());

        let n = a / (N::from(1.0).unwrap() - ecc_part * sin_part).sqrt();
        let h = f.altitude();

        let x = (h + n) * f.latitude().cos() * f.longitude().cos();
        let y = (h + n) * f.latitude().cos() * f.longitude().sin();
        let z = (h + n - N::from(ECCENTRICITY_SQ).unwrap() * n) * f.latitude().sin();

        ECEF::new(x, y, z)
    }
}

impl<N: Float> From<NVector<N>> for ECEF<N> {
    fn from(f: NVector<N>) -> ECEF<N> {
        // Constants used for calculation
        let x = f.vector().z;
        let y = f.vector().y;
        let z = -f.vector().x;
        let a_over_b = N::from(SEMI_MAJOR_AXIS).unwrap().powi(2) /
                       N::from(SEMI_MINOR_AXIS).unwrap().powi(2);
        // Multiplication part
        let mul = N::from(SEMI_MINOR_AXIS).unwrap() /
                  (x.powi(2) + a_over_b * y.powi(2) + a_over_b * z.powi(2)).sqrt();
        // NOTE: The following has been rearranged to follow ECEF convention
        // that Z points towards the north pole
        ECEF::new(-(mul * a_over_b * z + z * f.altitude()),
                  mul * a_over_b * y + y * f.altitude(),
                  mul * x + x * f.altitude())
    }
}

#[cfg(test)]
mod tests {
    use ::Access;
    use ::nvector::NVector;
    use ::wgs84::WGS84;
    use assert::close;
    use na::Norm;
    use super::*;

    fn ecef_close(a: ECEF<f64>, b: ECEF<f64>) {
        close(a.x(), b.x(), 0.000001);
        close(a.y(), b.y(), 0.000001);
        close(a.z(), b.z(), 0.000001);
    }

    quickcheck! {
        fn from_wgs84(wgs: WGS84<f64>) -> () {
            let test = WGS84::from(ECEF::from(wgs));

            close(wgs.latitude(), test.latitude(), 0.000001);
            close(wgs.longitude(), test.longitude(), 0.000001);
            close(wgs.altitude(), test.altitude(), 0.000001);
        }

        fn from_nvector(wgs: WGS84<f64>) -> () {
            let ans = ECEF::from(wgs);
            let test = ECEF::from(NVector::from(wgs));

            ecef_close(ans, test);
        }

        fn rotation(a: WGS84<f64>, b: WGS84<f64>) -> () {
            // Convert position into ECEF, we use WGS84 for simplicity
            let ecef_a = ECEF::from(a);
            let ecef_b = ECEF::from(b);

            // Calculate the vector between the two positions in ENU frame
            // and in Earth frame
            let vec_enu = ecef_b - ecef_a;
            let vec_e = ecef_b.0 - ecef_a.0;

            // Retrieve the rotation matrix for A converting ENU -> Earth
            let r_en = ecef_a.r_en();
            let vec_e2 = r_en * vec_enu.access();

            // These should be equivalent:
            close(vec_e.as_ref(), vec_e2.as_ref(), 0.0000001);
        }

        fn add_vector(a: WGS84<f64>, b: WGS84<f64>) -> () {
            let ecef_a = ECEF::from(a);
            let ecef_b = ECEF::from(b);

            let vec_ab = ecef_b - ecef_a;
            let ecef_b2 = ecef_a + vec_ab;

            ecef_close(ecef_b, ecef_b2);
        }

        fn distance(a: WGS84<f64>, b: WGS84<f64>) -> () {
            let dist_ab = (ECEF::from(b) - ECEF::from(a)).norm();
            let dist_ba = (ECEF::from(a) - ECEF::from(b)).norm();

            close(dist_ab, dist_ba, 0.000001);
        }
    }
}
