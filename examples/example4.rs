extern crate nav_types;

use nav_types::{ECEF, WGS84};

fn main() {
    // We start with a position given in latitude, longitude and altitude
    let position = WGS84::new_deg(36.12, -86.67, 0.0);

    // This can then easily be converted into ECEF by changing the type
    let ecef = ECEF::from(position);

    // NOTE: To do the same as Gade we could do:
    // `ECEF::from(NVector::from(position)`, however all types in
    // this library can be converted from one another and so this is
    // not strictly necessary

    println!("Position: {:?} in ECEF: {:?}", position, ecef);
}
