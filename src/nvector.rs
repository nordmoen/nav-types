use ::ecef::ECEF;
use ::enu::ENU;
use ::wgs84::{ECCENTRICITY_SQ, SEMI_MAJOR_AXIS, WGS84};
use na::Vector3;
use num_traits::Float;
use std::convert::From;
use std::ops::{Add, AddAssign, Sub, SubAssign};

/// N-Vector position
///
/// The N-Vector represents unique points on the earth's surface.
/// The advantage of N-Vectors is that they have no inconsistencies around
/// the poles compared to WGS84 Latitude, Longitude format.
/// See: [nvector](http://www.navlab.net/nvector/) for detailed information.
#[derive(PartialEq, Clone, Copy, Debug)]
pub struct NVector<N> {
    vec: Vector3<N>,
    alt: N,
}

impl<N> NVector<N> {
    /// Create a new NVector
    pub fn new(vec: Vector3<N>, altitude: N) -> NVector<N> {
        NVector {
            vec: vec,
            alt: altitude,
        }
    }
}

impl<N: Copy> NVector<N> {
    /// Get the vector component of this position
    pub fn vector(&self) -> Vector3<N> {
        self.vec
    }

    /// Get the altitude of this position
    pub fn altitude(&self) -> N {
        self.alt
    }
}

impl<N: Float> Add<ENU<N>> for NVector<N> {
    type Output = NVector<N>;
    fn add(self, right: ENU<N>) -> NVector<N> {
        // In accordance with Gade(2010), first convert to ECEF position
        // add vector in ECEF space and convert back
        NVector::from(ECEF::from(self) + right)
    }
}

impl<N: Float> AddAssign<ENU<N>> for NVector<N> {
    fn add_assign(&mut self, right: ENU<N>) {
        // In accordance with Gade(2010), first convert to ECEF position
        // add vector in ECEF space and convert back
        let new = NVector::from(ECEF::from(*self) + right);
        self.vec = new.vec;
        self.alt = new.alt;
    }
}

impl<N: Float> Sub<ENU<N>> for NVector<N> {
    type Output = NVector<N>;
    fn sub(self, right: ENU<N>) -> NVector<N> {
        // In accordance with Gade(2010), first convert to ECEF position
        // add vector in ECEF space and convert back
        NVector::from(ECEF::from(self) - right)
    }
}

impl<N: Float> Sub<NVector<N>> for NVector<N> {
    type Output = ENU<N>;
    fn sub(self, right: NVector<N>) -> ENU<N> {
        // In accordance with Gade(2010), first convert to ECEF position
        // add vector in ECEF space and convert back
        ECEF::from(self) - ECEF::from(right)
    }
}

impl<N: Float> SubAssign<ENU<N>> for NVector<N> {
    fn sub_assign(&mut self, right: ENU<N>) {
        // In accordance with Gade(2010), first convert to ECEF position
        // add vector in ECEF space and convert back
        let new = NVector::from(ECEF::from(*self) - right);
        self.vec = new.vec;
        self.alt = new.alt;
    }
}

impl<N: Float> From<WGS84<N>> for NVector<N> {
    fn from(f: WGS84<N>) -> NVector<N> {
        // This implementation defines the ECEF coordinate system to have the Z
        // axes point directly north, this affects the way which N-vectors are
        // defined. See: Table 2 in Gade(2010).
        // NOTE: This is consistent with the ECEF implementation in this crate
        let vec = Vector3::new(f.longitude().cos() * f.latitude().cos(),
                               f.longitude().sin() * f.latitude().cos(),
                               f.latitude().sin());
        NVector::new(vec, f.altitude())
    }
}

impl<N: Float> From<ECEF<N>> for NVector<N> {
    fn from(f: ECEF<N>) -> NVector<N> {
        // These are often used constants below:
        // a²
        let a_sq = N::from(SEMI_MAJOR_AXIS).unwrap().powi(2);
        // e²
        let e_2 = N::from(ECCENTRICITY_SQ).unwrap();
        // e⁴
        let e_4 = e_2.powi(2);

        // Select correct axis form ECEF vector
        let x = f.z();
        let y = f.y();
        let z = -f.x();

        let p = (y.powi(2) + z.powi(2)) / a_sq;
        let q = ((N::one() - e_2) / a_sq) * x.powi(2);
        let r = (p + q - e_4) / N::from(6.0).unwrap();
        let s = (e_4 * p * q) / (N::from(4.0).unwrap() * r.powi(3));
        let t = (N::one() + s + (s * (N::from(2.0).unwrap() + s)).sqrt()).cbrt();
        let u = r * (N::one() + t + t.recip());
        let v = (u.powi(2) + e_4 * q).sqrt();
        let w = e_2 * ((u + v - q) / (N::from(2.0).unwrap() * v));
        let k = (u + v + w.powi(2)).sqrt() - w;
        let d = (k * (y.powi(2) + z.powi(2)).sqrt()) / (k + e_2);

        let altitude = ((k + e_2 - N::one()) / k) * (d.powi(2) + x.powi(2)).sqrt();

        let denom = (d.powi(2) + x.powi(2)).sqrt().recip();
        let mul = k / (k + e_2);
        let vec = Vector3::new(-mul * z * denom, mul * y * denom, x * denom);

        NVector::new(vec, altitude)
    }
}

#[cfg(test)]
mod tests {
    use ::ecef::ECEF;
    use ::wgs84::WGS84;
    use assert::close;
    use super::*;

    quickcheck! {
        fn from_wgs84(wgs: WGS84<f64>) -> () {
            let test = WGS84::from(NVector::from(wgs));

            close(wgs.latitude(), test.latitude(), 0.000001);
            close(wgs.longitude(), test.longitude(), 0.000001);
            close(wgs.altitude(), test.altitude(), 0.000001);
        }

        fn from_ecef(wgs: WGS84<f64>) -> () {
            let ans = NVector::from(wgs);
            let test = NVector::from(ECEF::from(wgs));

            close(ans.altitude(), test.altitude(), 0.000001);
            close(ans.vector().as_ref(), test.vector().as_ref(), 0.000001)
        }
    }
}
