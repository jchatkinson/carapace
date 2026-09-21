/// Uniaxial material catalog. Closed enum, not a trait object (§3.1 of the
/// implementation plan) — the set is fixed at compile time. Leaf variants
/// only through M2; recursive composites (Parallel/Series/MinMax) land at
/// M7.
#[derive(Debug, Clone, Copy)]
pub enum Material {
    Elastic {
        e: f64,
    },
    /// Elastic-perfectly-plastic *envelope*: `stress = clamp(e*strain, -e*eyp, e*eyp)`.
    /// Reversible (no permanent-set tracking) — a simplified, stateless
    /// stand-in for OpenSees' `ElasticPPMaterial`, which stores plastic
    /// strain and unloads elastically from wherever it last yielded rather
    /// than back through the origin. Exact for monotonic loading that never
    /// reverses direction; diverges from the real material on unload/
    /// reload. Revisit alongside the other stateful leaves (Steel01,
    /// Concrete01, ...) at M7, once `Analysis` has real Newton iteration
    /// (M4) and materials gain trial/commit state.
    ElasticPP {
        e: f64,
        eyp: f64,
    },
    /// One-sided gap: zero stiffness until the strain closes the gap
    /// (`strain + gap < 0`, i.e. compressive strain past the opening), then
    /// linear elastic. Stateless — a function of current strain only, same
    /// caveat as `ElasticPP` above.
    Gap {
        e: f64,
        gap: f64,
    },
    /// "No tension": zero stiffness in tension (`strain >= 0`), linear
    /// elastic in compression. Equivalent to `Gap { e, gap: 0.0 }`, kept as
    /// its own variant since OpenSees exposes it separately and it's the
    /// simplest possible nonlinear material (useful as the base case for M2
    /// wiring).
    Ent {
        e: f64,
    },
}

impl Material {
    /// Given the current strain, return `(stress, tangent)`. Leaf materials
    /// are stateless for now (a pure function of current strain, no
    /// history); stateful leaves (Steel01, Concrete01, ...) will need
    /// `&mut self` and internal history variables when they land at M7.
    pub fn stress_tangent(&self, strain: f64) -> (f64, f64) {
        match self {
            Material::Elastic { e } => (e * strain, *e),
            Material::ElasticPP { e, eyp } => {
                let yield_stress = e * eyp;
                let trial_stress = e * strain;
                if trial_stress.abs() <= yield_stress {
                    (trial_stress, *e)
                } else {
                    (yield_stress.copysign(trial_stress), 0.0)
                }
            }
            Material::Gap { e, gap } => {
                // `<=`, not `<`: at exactly zero closure the tangent must
                // stay nonzero (matching the just-closed regime) or a
                // system starting from zero displacement forms a singular
                // initial tangent — there is no elastic-vs-open distinction
                // to make at a single point anyway.
                let closure = strain + gap;
                if closure <= 0.0 {
                    (e * closure, *e)
                } else {
                    (0.0, 0.0)
                }
            }
            Material::Ent { e } => {
                // See `Gap` above for why this is `<=`, not `<`.
                if strain <= 0.0 {
                    (e * strain, *e)
                } else {
                    (0.0, 0.0)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elastic_pp_clamps_symmetrically() {
        let m = Material::ElasticPP { e: 100.0, eyp: 0.01 };
        assert_eq!(m.stress_tangent(0.005), (0.5, 100.0));
        assert_eq!(m.stress_tangent(0.02), (1.0, 0.0));
        assert_eq!(m.stress_tangent(-0.02), (-1.0, 0.0));
    }

    #[test]
    fn gap_only_engages_past_closure() {
        let m = Material::Gap { e: 100.0, gap: 0.01 };
        assert_eq!(m.stress_tangent(0.0), (0.0, 0.0));
        assert_eq!(m.stress_tangent(-0.005), (0.0, 0.0));
        assert_eq!(m.stress_tangent(-0.02), (100.0 * (-0.02 + 0.01), 100.0));
    }

    #[test]
    fn ent_carries_no_tension() {
        let m = Material::Ent { e: 100.0 };
        assert_eq!(m.stress_tangent(0.01), (0.0, 0.0));
        assert_eq!(m.stress_tangent(-0.01), (-1.0, 100.0));
    }
}
