# nav-types
This crate is designed to make it easy to work with positions and vectors that
have meaning as geographic entities.

## Note about usage
The source is based on [`nalgebra`](http://nalgebra.org) and some methods are
only available if importing traits from `nalgebra`.

### Performance
Currently the only way to calculate vectors between latitude and longitude
positions is to convert to `ECEF` format and calculate the difference there.
This conversion happens behind the scenes and for this reason could be a source
of some surprise. It is therefor advised to try and use `ECEF` for as long as
possible only converting to and from at the beginning and end.

To see these result in action run the benchmark tests which shows the difference
between the formats.

# Example
```rust
extern crate nav_types;

use nav_types::WGS84;

let pos_a = WGS84::new(36.12, -86.67, 0.0);
let pos_b = WGS84::new(33.94, -118.40, 0.0);

println!("Distance between a and b: {:.2}m", a.distance(&b));
```

All position formats can work with vectors as long as the vectors are defined in
some coordinate system.
```rust
use nav_types::{WGS84, ENU};

let pos_a = WGS84::new(36.12, -86.67, 0.0);

let vec = ENU::new(0.0, 0.0, 10.0);
let pos_a_10m_up = pos_a + vec;
```
