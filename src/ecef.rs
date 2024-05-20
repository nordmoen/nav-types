use crate::enu::ENU;
use crate::nvector::NVector;
use crate::utils::RealFieldCopy;
use crate::wgs84::{ECCENTRICITY_SQ, SEMI_MAJOR_AXIS, SEMI_MINOR_AXIS, WGS84};
use crate::Access;
use na::{Matrix3, Point3};
use core::convert::{From, Into};
use core::ops::{Add, AddAssign, Sub, SubAssign};

/// Earth Centered Earth Fixed position
///
/// This struct represents a position in the ECEF coordinate system.
/// See: [ECEF](https://en.wikipedia.org/wiki/ECEF) for general description.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct ECEF<N: RealFieldCopy>(Point3<N>);

impl<N: RealFieldCopy> ECEF<N> {
    /// Create a new ECEF position
    pub fn new(x: N, y: N, z: N) -> ECEF<N> {
        ECEF(Point3::new(x, y, z))
    }
}

impl<N: RealFieldCopy> ECEF<N> {
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

impl<N: RealFieldCopy> ECEF<N> {
    /// Create a rotation matrix from ECEF frame to ENU frame
    fn r_enu_from_ecef(self) -> Matrix3<N> {
        let wgs = WGS84::from(self);
        Matrix3::new(
            -wgs.longitude_radians().sin(),
            wgs.longitude_radians().cos(),
            N::zero(),
            -wgs.latitude_radians().sin() * wgs.longitude_radians().cos(),
            -wgs.latitude_radians().sin() * wgs.longitude_radians().sin(),
            wgs.latitude_radians().cos(),
            wgs.latitude_radians().cos() * wgs.longitude_radians().cos(),
            wgs.latitude_radians().cos() * wgs.longitude_radians().sin(),
            wgs.latitude_radians().sin(),
        )
    }

    /// Create a rotation matrix from ENU frame to ECEF frame
    fn r_ecef_from_enu(self) -> Matrix3<N> {
        self.r_enu_from_ecef().transpose()
    }

    /// Euclidean distance between two ECEF positions
    pub fn distance(&self, other: &ECEF<N>) -> N
    where
        N: RealFieldCopy,
    {
        (other.0 - self.0).norm()
    }
}

impl<N, T> Add<T> for ECEF<N>
where
    N: RealFieldCopy,
    T: Into<ENU<N>>,
{
    type Output = ECEF<N>;
    fn add(self, right: T) -> ECEF<N> {
        let enu = T::into(right);
        let pos = self.0 + self.r_ecef_from_enu() * enu.access();
        ECEF(pos)
    }
}

impl<N, T> AddAssign<T> for ECEF<N>
where
    N: RealFieldCopy + AddAssign,
    T: Into<ENU<N>>,
{
    fn add_assign(&mut self, right: T) {
        let enu = T::into(right);
        self.0 = self.0 + (self.r_ecef_from_enu() * enu.access());
    }
}

impl<N: RealFieldCopy> Sub<ECEF<N>> for ECEF<N> {
    type Output = ENU<N>;
    fn sub(self, right: ECEF<N>) -> ENU<N> {
        let enu = right.r_enu_from_ecef() * (self.0 - right.0);
        ENU::new(enu.x, enu.y, enu.z)
    }
}

impl<N, T> Sub<T> for ECEF<N>
where
    N: RealFieldCopy,
    T: Into<ENU<N>>,
{
    type Output = ECEF<N>;
    fn sub(self, right: T) -> ECEF<N> {
        let enu = T::into(right);
        ECEF(self.0 - (self.r_ecef_from_enu() * enu.access()))
    }
}

impl<N, T> SubAssign<T> for ECEF<N>
where
    N: RealFieldCopy + SubAssign,
    T: Into<ENU<N>>,
{
    fn sub_assign(&mut self, right: T) {
        let enu = T::into(right);
        self.0 = self.0 - (self.r_ecef_from_enu() * enu.access());
    }
}

// This macro implements most standard operations for position types that
// can be converted to ECEF
macro_rules! ecef_impl {
    ($T:ident) => {
        impl<N, T> Add<T> for $T<N>
        where
            N: RealFieldCopy,
            T: Into<ENU<N>>,
        {
            type Output = $T<N>;
            fn add(self, right: T) -> Self {
                $T::from(ECEF::from(self) + right)
            }
        }

        impl<N, T> AddAssign<T> for $T<N>
        where
            N: RealFieldCopy + AddAssign,
            T: Into<ENU<N>>,
        {
            fn add_assign(&mut self, right: T) {
                *self = $T::from(ECEF::from(*self) + right);
            }
        }

        impl<N, T> Sub<T> for $T<N>
        where
            N: RealFieldCopy,
            T: Into<ENU<N>>,
        {
            type Output = $T<N>;
            fn sub(self, right: T) -> $T<N> {
                $T::from(ECEF::from(self) - right)
            }
        }

        impl<N: RealFieldCopy> Sub<$T<N>> for $T<N> {
            type Output = ENU<N>;
            fn sub(self, right: $T<N>) -> ENU<N> {
                ECEF::from(self) - ECEF::from(right)
            }
        }

        impl<N, T> SubAssign<T> for $T<N>
        where
            N: RealFieldCopy + SubAssign,
            T: Into<ENU<N>>,
        {
            fn sub_assign(&mut self, right: T) {
                *self = $T::from(ECEF::from(*self) - right);
            }
        }
    };
}

// In accordance with Gade(2010) all N-Vector operations works through
// converting to ECEF and then back
ecef_impl!(NVector);
// The only way to work with Latitude/Longitude is to convert to ECEF
ecef_impl!(WGS84);

impl<N: RealFieldCopy> From<WGS84<N>> for ECEF<N> {
    fn from(wgs: WGS84<N>) -> ECEF<N> {
        // Conversion from:
        // https://doi.org/10.1007/s00190-004-0375-4
        let semi_major_axis = N::from_f64(SEMI_MAJOR_AXIS).unwrap();
        let ecc_part = N::from_f64(ECCENTRICITY_SQ).unwrap();
        let sin_part = N::from_f64(0.5).unwrap()
            * (N::one() - (N::from_f64(2.0).unwrap() * wgs.latitude_radians()).cos());

        let n = semi_major_axis / (N::one() - ecc_part * sin_part).sqrt();
        let altitude = wgs.altitude();

        let x = (altitude + n) * wgs.latitude_radians().cos() * wgs.longitude_radians().cos();
        let y = (altitude + n) * wgs.latitude_radians().cos() * wgs.longitude_radians().sin();
        let z = (altitude + n - N::from_f64(ECCENTRICITY_SQ).unwrap() * n)
            * wgs.latitude_radians().sin();

        ECEF::new(x, y, z)
    }
}

impl<N: RealFieldCopy> From<NVector<N>> for ECEF<N> {
    fn from(f: NVector<N>) -> ECEF<N> {
        // Constants used for calculation
        let x = f.vector().z;
        let y = f.vector().y;
        let z = -f.vector().x;
        let a_over_b = N::from_f64(SEMI_MAJOR_AXIS).unwrap().powi(2)
            / N::from_f64(SEMI_MINOR_AXIS).unwrap().powi(2);
        // Multiplication part
        let mul = N::from_f64(SEMI_MINOR_AXIS).unwrap()
            / (x.powi(2) + a_over_b * y.powi(2) + a_over_b * z.powi(2)).sqrt();
        // NOTE: The following has been rearranged to follow ECEF convention
        // that Z points towards the north pole
        ECEF::new(
            -(mul * a_over_b * z + z * f.altitude()),
            mul * a_over_b * y + y * f.altitude(),
            mul * x + x * f.altitude(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enu::ENU;
    use crate::ned::NED;
    use crate::nvector::NVector;
    use crate::wgs84::WGS84;
    use crate::Access;
    use assert::close;

    // Helper method to check that two ECEF positions are equal
    fn ecef_close(a: ECEF<f64>, b: ECEF<f64>) {
        close(a.x(), b.x(), 0.000001);
        close(a.y(), b.y(), 0.000001);
        close(a.z(), b.z(), 0.000001);
    }

    quickcheck! {
        fn from_wgs84(wgs: WGS84<f64>) -> () {
            let test = WGS84::from(ECEF::from(wgs));

            close(wgs.latitude_radians(), test.latitude_radians(), 0.000001);
            close(wgs.longitude_radians(), test.longitude_radians(), 0.000001);
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
            let r_ecef_from_enu = ecef_a.r_ecef_from_enu();
            let vec_e2 = r_ecef_from_enu * vec_enu.access();

            // These should be equivalent:
            close(vec_e.as_ref(), vec_e2.as_ref(), 0.0000001);
        }

        fn add_vector(a: WGS84<f64>, b: WGS84<f64>) -> () {
            // Test that the difference between two positions added to
            // the first position leads to the second position
            let ecef_a = ECEF::from(a);
            let ecef_b = ECEF::from(b);

            let vec_ab = ecef_b - ecef_a;
            let ecef_b2 = ecef_a + vec_ab;

            ecef_close(ecef_b, ecef_b2);
        }

        // fn distance(a: WGS84<f64>, b: WGS84<f64>) -> () {
        //     // Test that the vector A->B and B->A is of equal length
        //     let dist_ab = (ECEF::from(b) - ECEF::from(a)).norm();
        //     let dist_ba = (ECEF::from(a) - ECEF::from(b)).norm();

        //     close(dist_ab, dist_ba, 0.000001);
        // }

        fn add_ned(a: WGS84<f64>, n: f64, e: f64, d: f64) -> () {
            // Test that adding a random NED vector to an ECEF position
            // is the same as adding the ENU version of the vector
            let ecef = ECEF::from(a);
            let enu_vec = ENU::new(e, n, -d);
            let ned_vec = NED::new(n, e, d);

            let ecef_enu = ecef + enu_vec;
            let ecef_ned = ecef + ned_vec;
            ecef_close(ecef_enu, ecef_ned);
        }

        fn sub_ned(a: WGS84<f64>, n: f64, e: f64, d: f64) -> () {
            // See `add_ned`
            let ecef = ECEF::from(a);
            let enu_vec = ENU::new(e, n, -d);
            let ned_vec = NED::new(n, e, d);

            let ecef_enu = ecef - enu_vec;
            let ecef_ned = ecef - ned_vec;
            ecef_close(ecef_enu, ecef_ned);
        }

        fn add_enu(a: WGS84<f64>, e: f64, n: f64, u: f64) -> () {
            // Test that adding a random ENU vector to an ECEF position
            // than taking the difference between those two positions
            // result in the original ENU vector
            let ecef = ECEF::from(a);
            let enu = ENU::new(e, n, u);

            let ecef_2 = ecef + enu;
            let enu_2 = ecef_2 - ecef;
            close(enu.access().as_ref(), enu_2.access().as_ref(), 0.0000001);
        }
    }
}
