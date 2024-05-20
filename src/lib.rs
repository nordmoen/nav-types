#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(feature = "dev", allow(unstable_features))]

//! Easily work with global positions and vectors.
//!
//! This crate is intended to allow easy and safe interoperability between
//! different global position and vectors in different coordinate systems.
//! The create is built on top of [nalgebra](http://nalgebra.org/doc/nalgebra/)
//! and strives to be compatible with types from that crate.
//!
//! # Example
//! ```rust
//! use nav_types::{ECEF, WGS84, ENU, NED};
//!
//! // Define positions using latitude and longitude on the WGS84 ellipsoid
//! let oslo = WGS84::from_degrees_and_meters(59.95, 10.75, 0.0);
//! let stockholm = WGS84::from_degrees_and_meters(59.329444, 18.068611, 0.0);
//!
//! println!("Great circle distance between Oslo and Stockholm: {:?}",
//!     oslo.distance(&stockholm));
//!
//! // Calculate vectors between positions
//! // This is equivalent of doint `ECEF::from(stockholm) - ECEF::from(oslo)`
//! let vec = stockholm - oslo;
//! println!("Vector between Oslo and Stockholm: {:?}", vec);
//!
//! // Easily convert between ENU and NED vectors
//! let ned_vec = NED::new(1f32, 0.0, 1.0);
//! let new_vec = vec + ned_vec;
//!
//! // Add vectors to positions
//! let stockholm_2 = ECEF::from(oslo) + new_vec;
//! ```
//!
//! # Where to begin
//! If you are unsure where to begin I would recommend you to research which
//! type of date you have available. Most of the time this is where you
//! should start. If no source is available or you are free to choose I recommend
//! that you start with `ECEF` and `ENU`. `ECEF` is efficient to calculate with
//! and `ENU` is quite straight forward and used in many contexts.

extern crate nalgebra as na;

#[cfg(test)]
#[macro_use]
extern crate quickcheck;
#[cfg(test)]
extern crate assert;

mod ecef;
mod enu;
mod ned;
mod nvector;
mod utils;
mod wgs84;

pub use self::utils::RealFieldCopy;
pub use self::ecef::ECEF;
pub use self::enu::ENU;
pub use self::ned::NED;
pub use self::nvector::NVector;
pub use self::wgs84::WGS84;

// This is a private trait to access the underlying structure, this is used
// so that this crate can access the implementation details without users
// of this crate knowing the underlying details.
trait Access<N> {
    // Transfer access of underlying vector
    fn access(self) -> N;
}
