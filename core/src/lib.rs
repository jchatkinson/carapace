pub mod analysis;
pub mod model;
pub mod testkit;

/// Closed-form check value only — `δ = PL/(AE)`. Kept around as the M0
/// oracle constant; M1 onward computes this through the real
/// `Domain`/`Element`/`Analysis` architecture instead (see
/// `tests/m1_truss.rs`).
pub fn axial_displacement(load: f64, length: f64, area: f64, modulus: f64) -> f64 {
    load * length / (area * modulus)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_xara_poc_truss_case() {
        // Same inputs as xara/scratchpad/poc/src/main.cpp's Stage-0 check:
        // L=100, A=2, E=30000, P=50 -> dx = 0.0833333333 (validated against
        // Xara's own native build, not just hand-derived).
        let dx = axial_displacement(50.0, 100.0, 2.0, 30000.0);
        assert!((dx - 0.0833333333).abs() < 1e-9);
    }
}
