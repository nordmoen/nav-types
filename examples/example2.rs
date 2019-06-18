extern crate nav_types;

use nav_types::{ENU, WGS84};

fn main() {
    // We start with the position of our "vehicle"
    let pos_a = WGS84::new(36.12, -86.67, 0.0);
    // Next we have the vector from the "vehicle" to the desired object,
    // it is assumed that the vector is rotated into proper ENU or NED
    // format
    let vec = ENU::new(10.0, 11.0, 0.0);

    // We then calculate the position of our mystery object
    let object = pos_a + vec;

    println!("Position of object is: {:?}", object);
}
