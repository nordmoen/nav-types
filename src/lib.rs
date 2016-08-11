#![cfg_attr(feature = "dev", allow(unstable_features))]
#![cfg_attr(feature = "dev", feature(plugin))]
#![cfg_attr(feature = "dev", plugin(clippy))]

//! Easily work with global positions and vectors.
//!
//! This crate is intended to allow easy and safe interoperability between
//! different global position and vectors in different coordinate systems.
//! The create is built on top of [nalgebra](http://nalgebra.org/doc/nalgebra/)
//! and strives to be compatible with types from that crate.
//!
//! # Example
//! ```rust
//! use nav_types::ECEF;
//! use nav_types::ENU;
//! use nav_types::NED;
//! use nav_types::WGS84;
//!
//! // Define positions using latitude and longitude on the WGS84 ellipsoid
//! let oslo = WGS84::new(59.95, 10.75, 0.0);
//! let stockholm = WGS84::new(59.329444, 18.068611, 0.0);
//!
//! println!("Great circle distance between Oslo and Stockholm: {:?}",
//!     oslo.distance(&stockholm));
//!
//! // Calculate vectors between positions
//! let vec = ECEF::from(stockholm) - ECEF::from(oslo);
//! println!("Vector between Oslo and Stockholm: {:?}", vec);
//!
//! // Easily convert between ENU and NED vectors
//! let ned_vec = NED::from(vec);
//!
//! // Add vectors to positions
//! let stockholm_2 = ECEF::from(oslo) + vec + ENU::new(0.0, 10.0, 0.0);
//! ```

extern crate nalgebra as na;
extern crate num_traits;

#[cfg(test)]
#[macro_use]
extern crate quickcheck;
#[cfg(test)]
extern crate assert;

mod ecef;
mod enu;
mod ned;
mod nvector;
mod wgs84;

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
