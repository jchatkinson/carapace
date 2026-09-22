/// Coordinate-transformation strategy for frame elements: relates an
/// element's local (basic) stiffness/forces to the global system. Closed
/// enum (§2.1) — `Linear` and `PDelta` are both "cheap" per the
/// implementation plan §3.4; `Corotational` (large-displacement) is its own
/// milestone (M9) if/when actually needed, not a third variant here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomTransf {
    /// Small-displacement: local/global relationship fixed at the element's
    /// original (undeformed) orientation, no correction for the effect of
    /// axial force on bending stiffness.
    Linear,
    /// Same small-displacement geometry as `Linear`, plus a linearized
    /// geometric-stiffness correction from the current axial force (a
    /// first-order P-Delta effect) — not full corotational tracking.
    PDelta,
}
