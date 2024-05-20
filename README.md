# nav-types

[![Build
Status](https://github.com/nordmoen/nav-types/actions/workflows/rust.yml/badge.svg)](https://github.com/nordmoen/nav-types/actions/workflows/rust.yml) [![Crates.io](https://img.shields.io/crates/v/nav-types.svg?maxAge=2592000)](https://crates.io/crates/nav-types)

This crate is designed to make it easy to work with positions and vectors that
have meaning as geographic entities.

[Documentation](https://nordmoen.github.io/nav-types)

## Example

Place this in your `Cargo.toml`

```toml
nav-types = "0.5.2"
```

and use it to calculate vectors and position:

```rust
extern crate nav_types;

use nav_types::WGS84;

let pos_a = WGS84::from_degrees_and_meters(36.12, -86.67, 0.0);
let pos_b = WGS84::from_degrees_and_meters(33.94, -118.40, 0.0);

println!("Distance between a and b: {:.2}m", a.distance(&b));
```

All position formats can work with vectors as long as the vectors are defined in
some coordinate system.

```rust
use nav_types::{WGS84, ENU};

let pos_a = WGS84::from_degrees_and_meters(36.12, -86.67, 0.0);

let vec = ENU::new(0.0, 0.0, 10.0);
let pos_a_10m_up = pos_a + vec;

// Or with `NED` vector

let ned_vec = NED::new(0.0, 0.0, -10.0);
let pos_a_10m_up_2 = pos_a + ned_vec;
```

## Note about usage

The source is based on [`nalgebra`](http://nalgebra.org) and some methods are
only available if importing traits from `nalgebra`.

Library works in `no_std` when `libm` feature is be enabled. This feature
replaces `std` math functions in `nalgebra` with `libm`.

### Performance

Currently the only way to calculate vectors between latitude and longitude
positions is to convert to `ECEF` format and calculate the difference there.
This conversion happens behind the scenes and for this reason could be a source
of some surprise. It is therefor advised to try and use `ECEF` for as long as
possible only converting to and from at the beginning and end.

On my laptop (absolute numbers will differ on different machines, but relative
differences should be similar)

```bash
running 10 tests
test ecef::add_vector    ... bench:       1,321 ns/iter (+/- 21)
test ecef::difference    ... bench:       1,306 ns/iter (+/- 16)
test ecef::from_nvector  ... bench:          21 ns/iter (+/- 1)
test ecef::from_wgs84    ... bench:         492 ns/iter (+/- 23)
test nvector::add_vector ... bench:       1,486 ns/iter (+/- 61)
test nvector::difference ... bench:       1,352 ns/iter (+/- 230)
test nvector::from_ecef  ... bench:         144 ns/iter (+/- 3)
test nvector::from_wgs84 ... bench:         382 ns/iter (+/- 16)
test wgs84::add_vector   ... bench:       2,154 ns/iter (+/- 302)
test wgs84::difference   ... bench:       2,305 ns/iter (+/- 304)
```
