use na::RealField;

/// Trait that combined RealField and Copy. Necessary to preserve simba v0.5.1 behavior
pub trait RealFieldCopy: RealField + Copy {}
impl<N: RealField + Copy> RealFieldCopy for N {}
