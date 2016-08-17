extern crate nav_types;

use nav_types::{WGS84, ECEF};

fn main() {
    // We start with our position in ECEF format
    let position = ECEF::new(299621.0, -5149462.0, 3738963.0);

    // We can now easily convert this position into different formats
    // for this concrete example we convert to latitude and longitude
    // using the WGS-84 ellipsoid
    let lat_lon = WGS84::from(position);

    // NOTE: To do the exact same as the example from Gade, we could do
    // `WGS84::from(NVector::from(position))`

    println!("{:?} in latitude and longitude: {:?}", position, lat_lon);
}
