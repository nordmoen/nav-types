extern crate nav_types;

use nav_types::WGS84;

fn main() {
    // First we define two positions in Latitude/Longitude,
    // this library only implement operations with WGS-84 ellipsoid
    // so we use that type
    let pos_a = WGS84::new(36.12, -86.67, 0.0);
    let pos_b = WGS84::new(33.94, -118.40, 0.0);

    // We could convert each of the positions to n-vectors with
    // `NVector::from(pos_a)`, but for the simple operations
    // we need it is not necessary

    // Calculate the delta vector. This is already rotated into
    // ENU format so no further work is needed
    let vec = pos_b - pos_a;

    // Next we calculate the azimuth given by atan2(east, north)
    let azimuth = f64::atan2(vec.east(), vec.north());

    println!("Vector: {:?}, azimuth: {:?}", vec, azimuth);
}
