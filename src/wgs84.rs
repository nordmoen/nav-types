use crate::ecef::ECEF;
use crate::nvector::NVector;
use num_traits::Float;
use std::convert::From;
use std::f32::consts::FRAC_PI_2;

#[cfg(test)]
use quickcheck::{Arbitrary, Gen};

pub const INVERSE_FLATTENING: f64 = 298.257_223_563;
pub const FLATTENING: f64 = 1.0 / INVERSE_FLATTENING;
/// Semi-major axis of WGS84 ellipsoid in meters. This represents the
/// radius of the WGS84 ellipsoid.
pub const SEMI_MAJOR_AXIS: f64 = 6_378_137.0;
pub const SEMI_MINOR_AXIS: f64 = SEMI_MAJOR_AXIS * (1.0 - FLATTENING);
pub const ECCENTRICITY_SQ: f64 = 2.0 * FLATTENING - FLATTENING * FLATTENING;

/// Geodetic position
///
/// This struct represents a position in the geodetic system on the WGS84
/// ellipsoid.
/// See: [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System) for
/// more information.
#[derive(PartialEq, Clone, Copy, Debug)]
pub struct WGS84<N> {
    // Represented as radians
    lat: N,
    // Represented as radians
    lon: N,
    alt: N,
}

impl<N: Float> WGS84<N> {
    /// Create a new WGS84 position
    ///
    /// # Arguments
    /// - `latitude` in degrees
    /// - `longitude` in degrees
    /// - `altitude` in meters
    ///
    /// # Panics
    /// This will panic if `latitude` or `longitude` are not defined on the
    /// WGS84 ellipsoid.
    pub fn new(latitude: N, longitude: N, altitude: N) -> WGS84<N> {
        assert!(
            latitude.abs() <= N::from(90.0).unwrap(),
            "Latitude must be in the range [-90, 90]"
        );
        assert!(
            longitude.abs() <= N::from(180.0).unwrap(),
            "Longitude must be in the range [-180, 180]"
        );
        WGS84 {
            lat: latitude.to_radians(),
            lon: longitude.to_radians(),
            alt: altitude,
        }
    }

    /// Try to create a new WGS84 position
    ///
    /// # Arguments
    /// - `latitude` in degrees
    /// - `longitude` in degrees
    /// - `altitude` in meters
    pub fn try_new(latitude: N, longitude: N, altitude: N) -> Option<WGS84<N>> {
        if latitude.abs() <= N::from(90.0).unwrap() && longitude.abs() <= N::from(180.0).unwrap() {
            Some(WGS84::new(latitude, longitude, altitude))
        } else {
            None
        }
    }

    /// Get latitude of position, in degrees
    pub fn latitude_degrees(&self) -> N {
        self.lat.to_degrees()
    }

    /// Get longitude of position, in degrees
    pub fn longitude_degrees(&self) -> N {
        self.lon.to_degrees()
    }

    /// Distance between two WGS84 positions
    ///
    /// This function uses the haversin formula to calculate the distance
    /// between two positions. For more control convert to ECEF and calculate
    /// the difference.
    ///
    /// # Examples
    /// ```rust
    /// use nav_types::WGS84;
    ///
    /// let oslo = WGS84::new(59.95, 10.75, 0.0);
    /// let stockholm = WGS84::new(59.329444, 18.068611, 0.0);
    ///
    /// println!("Great circle distance between Oslo and Stockholm: {:?}",
    ///     oslo.distance(&stockholm));
    /// ```
    pub fn distance(&self, other: &WGS84<N>) -> N {
        let delta_lat = other.latitude() - self.latitude();
        let delta_lon = other.longitude() - self.longitude();

        let a = (delta_lat / N::from(2.0).unwrap()).sin().powi(2)
            + self.latitude().cos()
                * other.latitude().cos()
                * (delta_lon / N::from(2.0).unwrap()).sin().powi(2);
        let c = N::from(2.0).unwrap() * a.sqrt().asin();

        N::from(SEMI_MAJOR_AXIS).unwrap() * c + (self.altitude() - other.altitude()).abs()
    }
}

impl<N: Copy> WGS84<N> {
    /// Get altitude of position
    pub fn altitude(&self) -> N {
        self.alt
    }
    /// Get latitude in radians
    pub fn latitude(&self) -> N {
        self.lat
    }
    /// Get longitude in radians
    pub fn longitude(&self) -> N {
        self.lon
    }
}

impl<N: Float> From<NVector<N>> for WGS84<N> {
    fn from(f: NVector<N>) -> WGS84<N> {
        // This implementation defines the ECEF coordinate system to have the Z
        // axes point directly north, this affects the way which N-vectors are
        // defined. See: Table 2 in Gade(2010).
        // NOTE: This is consistent with the ECEF implementation in this crate
        let x = f.vector().z;
        let y = f.vector().y;
        let z = -f.vector().x;

        // See equation (5) in Gade(2010)
        let longitude = y.atan2(-z);
        // Equatorial component, see equation (6) Gade(2010)
        let equat = (y.powi(2) + z.powi(2)).sqrt();
        // See equation (6) Gade(2010)
        let latitude = x.atan2(equat);

        WGS84 {
            lat: latitude,
            lon: longitude,
            alt: f.altitude(),
        }
    }
}

impl<N: Float> From<ECEF<N>> for WGS84<N> {
    fn from(f: ECEF<N>) -> WGS84<N> {
        // Conversion from:
        // http://download.springer.com/static/pdf/723/art%253A10.1007%252Fs00190-004-0375-4
        // .pdf?originUrl=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1007%2Fs00190-004-0375-4&
        // token
        // 2=exp=1470729285~acl=%2Fstatic%2Fpdf%2F723%2Fart%25253A10.1007%25252Fs00190-004-0375-4.
        // pdf%
        // 3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Farticle%252F10.1007%252Fs00190-
        // 004-0375-4*~hmac=4aaedb6b71f13fc9a9ce5175b4538c3ec38ddf11b77531d3bd2af75ee1fc2061
        //
        // These are often used constants below:
        // a²
        let a_sq = N::from(SEMI_MAJOR_AXIS).unwrap().powi(2);
        // e²
        let e_2 = N::from(ECCENTRICITY_SQ).unwrap();
        // e⁴
        let e_4 = N::from(ECCENTRICITY_SQ).unwrap().powi(2);

        let p = (f.x().powi(2) + f.y().powi(2)) / a_sq;
        let q = ((N::one() - e_2) / a_sq) * f.z().powi(2);
        let r = (p + q - e_4) / N::from(6.0).unwrap();
        let s = e_4 * ((p * q) / (N::from(4.0).unwrap() * r.powi(3)));
        let t = (N::one() + s + (s * (N::from(2.0).unwrap() + s)).sqrt()).cbrt();
        let u = r * (N::one() + t + t.recip());
        let v = (u.powi(2) + e_4 * q).sqrt();
        let w = e_2 * ((u + v - q) / (N::from(2.0).unwrap() * v));
        let k = (u + v + w.powi(2)).sqrt() - w;
        let d = (k * (f.x().powi(2) + f.y().powi(2)).sqrt()) / (k + e_2);
        let pi_half = N::from(FRAC_PI_2).unwrap();

        let altitude = ((k + e_2 - N::one()) / k) * (d.powi(2) + f.z().powi(2)).sqrt();
        let latitude = N::from(2.0).unwrap() * f.z().atan2(d + (d.powi(2) + f.z().powi(2)).sqrt());
        let longitude = if f.y() >= N::zero() {
            pi_half
                - N::from(2.0).unwrap()
                    * f.x().atan2((f.x().powi(2) + f.y().powi(2)).sqrt() + f.y())
        } else {
            -pi_half
                + N::from(2.0).unwrap()
                    * f.x().atan2((f.x().powi(2) + f.y().powi(2)).sqrt() - f.y())
        };

        WGS84 {
            lat: latitude,
            lon: longitude,
            alt: altitude,
        }
    }
}

#[cfg(test)]
impl Arbitrary for WGS84<f64> {
    fn arbitrary<G: Gen>(g: &mut G) -> WGS84<f64> {
        let lat = g.gen_range(-90.0, 90.0);
        let lon = g.gen_range(-180.0, 180.0);
        // NOTE: Minimum altitude is radius of the earth, however, due to
        // float precision we have shrunk that down a bit, the values
        // below might still cause problems
        let alt = g.gen_range(-6300000.0, 10000000.0);

        WGS84::new(lat, lon, alt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert::close;
    use enu::ENU;
    use quickcheck::{quickcheck, TestResult};

    fn create_wgs84(latitude: f32, longitude: f32, altitude: f32) -> TestResult {
        // This function is used to check that illegal latitude and longitude
        // values panics
        if latitude.abs() <= 90.0 && longitude.abs() <= 180.0 {
            // If both latitude and longitude are within acceptable ranges
            // we tell quickcheck to discard the test so that it will
            // re-generate a test with different parameters hopefully
            // testing other parameters which will fail
            TestResult::discard()
        } else {
            // If either latitude or longitude is outside acceptable range
            // the test must fail
            TestResult::must_fail(move || {
                WGS84::new(latitude, longitude, altitude);
            })
        }
    }

    #[test]
    fn test_create_wgs84() {
        quickcheck(create_wgs84 as fn(f32, f32, f32) -> TestResult);
    }

    #[test]
    fn dateline() {
        // This test ensures that when moving west from the dateline
        // the longitude resets to -180
        let a = WGS84::new(20.0, 180.0, 0.0);
        let vec = ENU::new(-10.0, 0.0, 0.0);

        let ans = a + vec;
        // NOTE: Precision here is rather arbitrary and depends on the
        // length of the vector used and so on, we are mostly interested
        // in seeing that the longitude is reset to negative
        close(ans.longitude_degrees(), -180.0, 0.001);
        close(ans.latitude_degrees(), 20.0, 0.001);
    }

    quickcheck! {
        fn distance_haversin(a: WGS84<f64>, b: WGS84<f64>) -> () {
            // Trivially true:
            close(a.distance(&b), b.distance(&a), 0.0000001);
        }
    }
}
