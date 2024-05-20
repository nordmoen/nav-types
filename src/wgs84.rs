use crate::ecef::ECEF;
use crate::nvector::NVector;
use crate::utils::RealFieldCopy;
use core::convert::From;

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
/// The `serde` feature allows this to be Serialized / Deserialized.
/// If serialized into json, it will look like this. Enabled thought
/// the `serde` feature
/// ```json
/// {
///    "latitude": 0.0,
///    "longitude": 0.0,
///    "altitude": 0.0
/// }
/// ```
/// Note: latitude and longitude values will be in radians
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Clone, Copy, Debug)]
pub struct WGS84<N> {
    // Represented as radians
    #[cfg_attr(feature = "serde", serde(rename = "latitude"))]
    lat: N,
    // Represented as radians
    #[cfg_attr(feature = "serde", serde(rename = "longitude"))]
    lon: N,
    #[cfg_attr(feature = "serde", serde(rename = "altitude"))]
    alt: N,
}

impl<N: RealFieldCopy> WGS84<N>
where
    f64: From<N>,
{
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
    pub fn from_degrees_and_meters(latitude: N, longitude: N, altitude: N) -> WGS84<N> {
        assert!(
            latitude.abs() <= N::from_f64(90.0).unwrap(),
            "Latitude must be in the range [-90, 90]"
        );
        assert!(
            longitude.abs() <= N::from_f64(180.0).unwrap(),
            "Longitude must be in the range [-180, 180]"
        );
        WGS84 {
            lat: N::from_f64(f64::from(latitude).to_radians()).unwrap(),
            lon: N::from_f64(f64::from(longitude).to_radians()).unwrap(),
            alt: altitude,
        }
    }

    /// Try to create a new WGS84 position
    ///
    /// # Arguments
    /// - `latitude` in degrees
    /// - `longitude` in degrees
    /// - `altitude` in meters
    pub fn try_from_degrees_and_meters(latitude: N, longitude: N, altitude: N) -> Option<WGS84<N>> {
        if latitude.abs() <= N::from_f64(90.0).unwrap()
            && longitude.abs() <= N::from_f64(180.0).unwrap()
        {
            Some(WGS84::from_degrees_and_meters(
                latitude, longitude, altitude,
            ))
        } else {
            None
        }
    }

    /// Create a new WGS84 position
    ///
    /// # Arguments
    /// - `latitude` in radians
    /// - `longitude` in radians
    /// - `altitude` in meters
    ///
    /// # Panics
    /// This will panic if `latitude` or `longitude` are not defined on the
    /// WGS84 ellipsoid.
    pub fn from_radians_and_meters(latitude: N, longitude: N, altitude: N) -> WGS84<N> {
        assert!(
            latitude.abs() <= N::from_f64(core::f64::consts::FRAC_PI_2).unwrap(),
            "Latitude must be in the range [-π/2, π/2]"
        );
        assert!(
            longitude.abs() <= N::from_f64(core::f64::consts::PI).unwrap(),
            "Longitude must be in the range [-π, π]"
        );
        WGS84 {
            lat: latitude,
            lon: longitude,
            alt: altitude,
        }
    }

    /// Try to create a new WGS84 position
    ///
    /// # Arguments
    /// - `latitude` in radians
    /// - `longitude` in radians
    /// - `altitude` in meters
    pub fn try_from_radians_and_meters(latitude: N, longitude: N, altitude: N) -> Option<WGS84<N>> {
        if latitude.abs() <= N::from_f64(core::f64::consts::FRAC_PI_2).unwrap()
            && longitude.abs() <= N::from_f64(core::f64::consts::PI).unwrap()
        {
            Some(WGS84 {
                lat: latitude,
                lon: longitude,
                alt: altitude,
            })
        } else {
            None
        }
    }

    /// Get latitude of position, in degrees
    pub fn latitude_degrees(&self) -> N {
        N::from_f64(f64::from(self.lat).to_degrees()).unwrap()
    }

    /// Get longitude of position, in degrees
    pub fn longitude_degrees(&self) -> N {
        N::from_f64(f64::from(self.lon).to_degrees()).unwrap()
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
    /// let oslo = WGS84::from_degrees_and_meters(59.95, 10.75, 0.0);
    /// let stockholm = WGS84::from_degrees_and_meters(59.329444, 18.068611, 0.0);
    ///
    /// println!("Great circle distance between Oslo and Stockholm: {:?}",
    ///     oslo.distance(&stockholm));
    /// ```
    pub fn distance(&self, other: &WGS84<N>) -> N {
        let delta_lat = other.latitude_radians() - self.latitude_radians();
        let delta_lon = other.longitude_radians() - self.longitude_radians();

        let a = (delta_lat / N::from_f64(2.0).unwrap()).sin().powi(2)
            + self.latitude_radians().cos()
                * other.latitude_radians().cos()
                * (delta_lon / N::from_f64(2.0).unwrap()).sin().powi(2);
        let c = N::from_f64(2.0).unwrap() * a.sqrt().asin();

        N::from_f64(SEMI_MAJOR_AXIS).unwrap() * c + (self.altitude() - other.altitude()).abs()
    }
}

impl<N: Copy> WGS84<N> {
    /// Get altitude of position, in meters
    pub fn altitude(&self) -> N {
        self.alt
    }
    /// Get latitude of position, in radians
    pub fn latitude_radians(&self) -> N {
        self.lat
    }
    /// Get longitude of position, in radians
    pub fn longitude_radians(&self) -> N {
        self.lon
    }
}

impl<N: RealFieldCopy> From<NVector<N>> for WGS84<N> {
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

impl<N: RealFieldCopy> From<ECEF<N>> for WGS84<N> {
    #![allow(clippy::many_single_char_names)]
    fn from(ecef: ECEF<N>) -> WGS84<N> {
        // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#The_application_of_Ferrari's_solution
        let a = N::from_f64(SEMI_MAJOR_AXIS).unwrap();
        let b = N::from_f64(SEMI_MINOR_AXIS).unwrap();
        let r_squared = (ecef.x() * ecef.x()) + (ecef.y() * ecef.y());
        let r = r_squared.sqrt();
        let z_squared = ecef.z() * ecef.z();
        let z = z_squared.sqrt();
        let a_squared = a * a;
        let b_squared = b * b;
        let e_squared = N::one() - (b_squared / a_squared);
        let e_dot_squared = (a_squared - b_squared) / b_squared;
        let f = N::from_f64(54.0).unwrap() * b_squared * z_squared;
        let g = r_squared + ((N::one() - e_squared) * z_squared)
            - (e_squared * (a_squared - b_squared));
        let g_squared = g * g;
        let c = (e_squared * e_squared * f * r_squared) / (g * g_squared);
        let s = (N::one() + c + ((c * c) + c + c).sqrt()).cbrt();
        let s_plus_one_over_s_plus_one = s + (N::one() / s) + N::one();
        let p = f
            / (N::from_f64(3.0).unwrap()
                * s_plus_one_over_s_plus_one
                * s_plus_one_over_s_plus_one
                * g_squared);
        let q = (N::one() + (N::from_f64(2.0).unwrap() * e_squared * e_squared * p)).sqrt();
        let r_0 = (-(p * e_squared * r) / (N::one() + q))
            + (a_squared / N::from_f64(2.0).unwrap() * (N::one() + (N::one() / q))
                - ((p * (N::one() - e_squared) * z_squared) / (q * (N::one() + q)))
                - (p * r_squared / N::from_f64(2.0).unwrap()))
            .sqrt();
        let r_minus_e_squared_r0 = r - (e_squared * r_0);
        let r_minus_e_squared_r0_squared = r_minus_e_squared_r0 * r_minus_e_squared_r0;
        let u = (r_minus_e_squared_r0_squared + z_squared).sqrt();
        let v = (r_minus_e_squared_r0_squared + ((N::one() - e_squared) * z_squared)).sqrt();
        let z_0 = (b_squared * z) / (a * v);
        let h = u * (N::one() - (b_squared / (a * v)));
        let phi = ecef.z().signum() * ((z + (e_dot_squared * z_0)) / r).atan().abs();
        let lambda = ecef.y().atan2(ecef.x());
        WGS84 {
            lat: phi,
            lon: lambda,
            alt: h,
        }
    }
}

#[cfg(test)]
impl Arbitrary for WGS84<f64> {
    fn arbitrary<G: Gen>(g: &mut G) -> WGS84<f64> {
        use rand::Rng;
        let lat = g.gen_range(-90.0, 90.0);
        let lon = g.gen_range(-180.0, 180.0);
        // NOTE: Minimum altitude is radius of the earth, however, due to
        // float precision we have shrunk that down a bit, the values
        // below might still cause problems
        let alt = g.gen_range(-6300000.0, 10000000.0);

        WGS84::from_degrees_and_meters(lat, lon, alt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enu::ENU;
    use assert::close;
    use quickcheck::{quickcheck, TestResult};

    #[test]
    #[cfg_attr(not(feature = "serde"), ignore)]
    fn test_ser_de() {
        #[cfg(feature = "serde")]
        {
            use serde_test::{assert_tokens, Token};
            let oslo: WGS84<f64> = WGS84 {
                lat: 1.0463,
                lon: 0.1876,
                alt: 0.0,
            };
            assert_tokens(
                &oslo,
                &[
                    Token::Struct {
                        name: "WGS84",
                        len: 3,
                    },
                    Token::Str("latitude"),
                    Token::F64(1.0463),
                    Token::Str("longitude"),
                    Token::F64(0.1876),
                    Token::Str("altitude"),
                    Token::F64(0.0),
                    Token::StructEnd,
                ],
            );
        }
        #[cfg(not(feature = "serde"))]
        {
            panic!("This test requires the serde feature to be enabled");
        }
    }

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
                WGS84::from_degrees_and_meters(latitude, longitude, altitude);
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
        // the longitude resets to 180
        let a = WGS84::from_degrees_and_meters(20.0, 180.0, 0.0);
        let vec = ENU::new(-10.0, 0.0, 0.0);

        let ans = a + vec;
        // NOTE: Precision here is rather arbitrary and depends on the
        // length of the vector used and so on, we are mostly interested
        // in seeing that the longitude is reset to positive
        close(ans.longitude_degrees(), 180.0, 0.001);
        close(ans.latitude_degrees(), 20.0, 0.001);
    }

    #[test]
    fn conversion_inversion_ecef() {
        let oslo: WGS84<f64> = WGS84::from_degrees_and_meters(59.95, 10.75, 0.0);
        let stockholm: WGS84<f64> = WGS84::from_degrees_and_meters(59.329444, 18.068611, 0.0);

        for &place in [oslo, stockholm].iter() {
            let distance = place.distance(&WGS84::from(ECEF::from(place)));
            close(distance, 0.0, 0.00000001);
        }
    }

    #[test]
    fn conversion_ecef() {
        let oslo_wgs84: WGS84<f64> = WGS84::from_degrees_and_meters(59.95, 10.75, 0.0);
        let oslo_ecef: ECEF<f64> = ECEF::new(3145735.0, 597236.0, 5497690.0);

        for &(place_wgs84, place_ecef) in [(oslo_wgs84, oslo_ecef)].iter() {
            let distance = place_wgs84.distance(&WGS84::from(place_ecef));
            close(distance, 0.0, 1.0);
        }
    }

    #[test]
    fn add_enu() {
        let oslo: WGS84<f64> = WGS84::from_degrees_and_meters(59.95, 10.75, 0.0);
        let oslo_high: WGS84<f64> = WGS84::from_degrees_and_meters(59.95, 10.75, 100.0);

        let stockholm: WGS84<f64> = WGS84::from_degrees_and_meters(59.329444, 18.068611, 0.0);
        let stockholm_high: WGS84<f64> =
            WGS84::from_degrees_and_meters(59.329444, 18.068611, 100.0);

        for &(place, place_high) in [(oslo, oslo_high), (stockholm, stockholm_high)].iter() {
            let distance =
                ECEF::from(place_high).distance(&(ECEF::from(place) + ENU::new(0.0, 0.0, 100.0)));
            close(distance, 0.0, 0.00001);
        }
    }

    quickcheck! {
        fn distance_haversin(a: WGS84<f64>, b: WGS84<f64>) -> () {
            // Trivially true:
            close(a.distance(&b), b.distance(&a), 0.0000001);
        }
    }
}
