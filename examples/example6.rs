extern crate nav_types;

use nav_types::{AER, WGS84};

fn main() {
    // We start with the position of our "vehicle"
    let pos_a = WGS84::from_degrees_and_meters(50.0, 10.0, 100.0);
    // Next we have the vector from the "vehicle" to the desired object,
    // given by azimuth, elevation and range
    let vec = AER::from_degrees_and_meters(45.0, -10.0, 100.0);

    // We then calculate the position of our mystery object
    let object = pos_a + vec;

    println!("Position of object is: {:?}", object);
}
