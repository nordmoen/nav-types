#![feature(test)]
extern crate nav_types;
extern crate test;

#[cfg(test)]
mod ecef {
    use nav_types::*;
    use test::Bencher;
    #[bench]
    fn from_wgs84(b: &mut Bencher) {
        let a = WGS84::from_degrees_and_meters(36.12, -86.67, 12.3);
        b.iter(|| ECEF::from(a));
    }
    #[bench]
    fn from_nvector(b: &mut Bencher) {
        let a = NVector::from(WGS84::from_degrees_and_meters(36.12, -86.67, 12.3));
        b.iter(|| ECEF::from(a));
    }
    #[bench]
    fn difference(b: &mut Bencher) {
        let pos_a = ECEF::from(WGS84::from_degrees_and_meters(36.12, -86.67, 0.0));
        let pos_b = ECEF::from(WGS84::from_degrees_and_meters(33.94, -118.40, 0.0));
        b.iter(|| pos_b - pos_a);
    }
    #[bench]
    fn add_vector(b: &mut Bencher) {
        let pos_a = ECEF::from(WGS84::from_degrees_and_meters(36.12, -86.67, 0.0));
        let pos_b = ECEF::from(WGS84::from_degrees_and_meters(33.94, -118.40, 0.0));

        let vec = pos_b - pos_a;
        b.iter(|| pos_a + vec);
    }
}
#[cfg(test)]
mod nvector {
    use nav_types::*;
    use test::Bencher;
    #[bench]
    fn from_wgs84(b: &mut Bencher) {
        let a = WGS84::from_degrees_and_meters(36.12, -86.67, 12.3);
        b.iter(|| NVector::from(a));
    }
    #[bench]
    fn from_ecef(b: &mut Bencher) {
        let a = ECEF::from(WGS84::from_degrees_and_meters(36.12, -86.67, 12.3));
        b.iter(|| NVector::from(a));
    }
    #[bench]
    fn difference(b: &mut Bencher) {
        let pos_a = NVector::from(WGS84::from_degrees_and_meters(36.12, -86.67, 0.0));
        let pos_b = NVector::from(WGS84::from_degrees_and_meters(33.94, -118.40, 0.0));
        b.iter(|| pos_b - pos_a);
    }
    #[bench]
    fn add_vector(b: &mut Bencher) {
        let pos_a = NVector::from(WGS84::from_degrees_and_meters(36.12, -86.67, 0.0));
        let pos_b = NVector::from(WGS84::from_degrees_and_meters(33.94, -118.40, 0.0));

        let vec = pos_b - pos_a;
        b.iter(|| pos_a + vec);
    }
}

#[cfg(test)]
mod wgs84 {
    use nav_types::*;
    use test::Bencher;

    #[bench]
    fn difference(b: &mut Bencher) {
        let pos_a = WGS84::from_degrees_and_meters(36.12, -86.67, 0.0);
        let pos_b = WGS84::from_degrees_and_meters(33.94, -118.40, 0.0);

        b.iter(|| pos_b - pos_a);
    }

    #[bench]
    fn add_vector(b: &mut Bencher) {
        let pos_a = WGS84::from_degrees_and_meters(36.12, -86.67, 0.0);
        let pos_b = WGS84::from_degrees_and_meters(33.94, -118.40, 0.0);

        let vec = pos_b - pos_a;
        b.iter(|| pos_a + vec);
    }

    #[bench]
    fn from_ecef(b: &mut Bencher) {
        let ecef = ECEF::from(WGS84::from_degrees_and_meters(36.12, -86.67, 0.0));

        b.iter(|| WGS84::from(ecef));
    }

    #[bench]
    fn from_nvector(b: &mut Bencher) {
        let nvec = NVector::from(WGS84::from_degrees_and_meters(36.12, -86.67, 0.0));

        b.iter(|| WGS84::from(nvec));
    }
}
