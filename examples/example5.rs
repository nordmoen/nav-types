extern crate nalgebra as na;
extern crate nav_types;

use na::{Cross, Dot, Norm};
use nav_types::{NVector, WGS84};

static EARTH_RADIUS: f64 = 6371008.8;

fn main() {
    // First we define two positions in Latitude/Longitude,
    // this library only implement operations with WGS-84 ellipsoid
    // so we use that type
    let pos_a = WGS84::new(36.12, -86.67, 0.0);
    let pos_b = WGS84::new(33.94, -118.40, 0.0);

    // We then convert these positions into n-vectors:
    let n_a = NVector::from(pos_a);
    let n_b = NVector::from(pos_b);

    // To calculate the surface distance we use cross and dot
    // product
    let surface_distance = f64::atan2(
        n_a.vector().cross(&n_b.vector()).norm(),
        n_a.vector().dot(&n_b.vector()),
    ) * EARTH_RADIUS;

    // Euclidean distance can be calculated by subtracting the positions
    // and taking the norm
    let euclid_distance = (n_b - n_a).norm();

    println!(
        "Surface distance {:?}m, euclidean distance {:?}m",
        surface_distance, euclid_distance
    );
}
