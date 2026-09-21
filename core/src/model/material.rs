/// Uniaxial material catalog. Closed enum, not a trait object (§3.1 of the
/// implementation plan) — the set is fixed at compile time. Leaf variants
/// only through M6; recursive composites (Parallel/Series/MinMax) land at
/// M7.
///
/// # Trial vs. commit
///
/// A `Material` value holds only *committed* state (e.g. `ElasticPP`'s
/// plastic strain `ep`) — there is no OpenSees-style parallel `T`/`C`
/// field pair. Instead there are two pure (`&self`, non-mutating) views
/// over that committed state:
///
/// - [`trial_stress_tangent`](Material::trial_stress_tangent) — called
///   every Newton iteration, as many times as needed, always relative to
///   the same fixed committed state. Since it can't mutate anything,
///   there's no way for one iteration's (possibly discarded) trial to leak
///   into the next — a real correctness requirement for path-dependent
///   materials, not just a convenience: `setTrialStrain` in OpenSees is
///   `&mut self`, but its `T*` fields are always recomputed from the `C*`
///   fields, i.e. it's *this* pure-function property, not the mutation
///   itself, that OpenSees's parallel field pair exists to guarantee.
/// - [`commit`](Material::commit) — called once, only after a step's
///   Newton iteration actually converges (or, for `Algorithm::Linear`,
///   after its one unconditional solve), producing the new committed
///   `Material`. `Domain::commit` is what calls this and writes the
///   result back into the model; nothing else ever mutates a `Material`.
///
/// This gets rollback "for free": since nothing is mutated until
/// `commit`, a step that fails partway through a Newton loop has left no
/// material in an inconsistent state to revert — only nodal displacement
/// (mutated eagerly every iteration, which is correct Newton behavior)
/// needs `Analysis`/`TransientAnalysis`'s snapshot-and-restore-on-failure.
#[derive(Debug, Clone, Copy)]
pub enum Material {
    Elastic {
        e: f64,
    },
    /// Elastic-perfectly-plastic, with real permanent-set tracking via the
    /// committed plastic strain `ep` — standard 1D return-mapping
    /// plasticity. Construct via [`Material::elastic_pp`], not the struct
    /// literal: `ep` is internal history, not a material parameter, and
    /// must start at zero.
    ElasticPP {
        e: f64,
        eyp: f64,
        ep: f64,
    },
    /// One-sided gap: zero stiffness until the strain closes the gap
    /// (`strain + gap < 0`, i.e. compressive strain past the opening), then
    /// linear elastic. Stateless — a genuine function of current strain
    /// only (no permanent set to track), unlike `ElasticPP`.
    Gap {
        e: f64,
        gap: f64,
    },
    /// "No tension": zero stiffness in tension (`strain >= 0`), linear
    /// elastic in compression. Equivalent to `Gap { e, gap: 0.0 }`, kept as
    /// its own variant since OpenSees exposes it separately and it's the
    /// simplest possible nonlinear material (useful as the base case for M2
    /// wiring). Stateless, same reasoning as `Gap`.
    Ent {
        e: f64,
    },
}

impl Material {
    pub fn elastic_pp(e: f64, eyp: f64) -> Self {
        Material::ElasticPP { e, eyp, ep: 0.0 }
    }

    /// `(stress, tangent)` for a trial evaluation at `strain` — read-only,
    /// safe to call any number of times from any Newton iteration. See the
    /// type-level doc comment for why this must never mutate `self`.
    pub fn trial_stress_tangent(&self, strain: f64) -> (f64, f64) {
        let (stress, tangent, _next) = self.evaluate(strain);
        (stress, tangent)
    }

    /// The new committed `Material` if `strain` were accepted as the
    /// converged state. Called once per converged step, by `Domain::commit`
    /// — see the type-level doc comment.
    pub fn commit(&self, strain: f64) -> Material {
        let (_stress, _tangent, next) = self.evaluate(strain);
        next
    }

    /// Shared core: `(stress, tangent, next_committed_state)`. Leaf
    /// materials without history (`Elastic`, `Gap`, `Ent`) return `*self`
    /// unchanged as `next_committed_state` — committing them is a no-op,
    /// correctly, since they have nothing to remember.
    fn evaluate(&self, strain: f64) -> (f64, f64, Material) {
        match *self {
            Material::Elastic { e } => (e * strain, e, *self),
            Material::ElasticPP { e, eyp, ep } => {
                let yield_stress = e * eyp;
                // Elastic predictor relative to the *committed* plastic
                // strain — standard 1D return mapping.
                let trial_stress = e * (strain - ep);
                if trial_stress.abs() <= yield_stress {
                    (trial_stress, e, Material::ElasticPP { e, eyp, ep })
                } else {
                    let stress = yield_stress.copysign(trial_stress);
                    // Return to the yield surface: the plastic strain
                    // absorbs whatever additional strain the elastic
                    // predictor overshot by.
                    let next_ep = strain - stress / e;
                    (
                        stress,
                        0.0,
                        Material::ElasticPP {
                            e,
                            eyp,
                            ep: next_ep,
                        },
                    )
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
                    (e * closure, e, *self)
                } else {
                    (0.0, 0.0, *self)
                }
            }
            Material::Ent { e } => {
                // See `Gap` above for why this is `<=`, not `<`.
                if strain <= 0.0 {
                    (e * strain, e, *self)
                } else {
                    (0.0, 0.0, *self)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elastic_pp_clamps_symmetrically_on_first_loading() {
        // First loading from ep=0: identical to the old reversible-envelope
        // behavior, since nothing has been committed yet to make history
        // matter.
        let m = Material::elastic_pp(100.0, 0.01);
        assert_eq!(m.trial_stress_tangent(0.005), (0.5, 100.0));
        assert_eq!(m.trial_stress_tangent(0.02), (1.0, 0.0));
        assert_eq!(m.trial_stress_tangent(-0.02), (-1.0, 0.0));
    }

    #[test]
    fn elastic_pp_remembers_permanent_set_after_commit() {
        let m = Material::elastic_pp(100.0, 0.01);

        // Load past yield (strain=0.02 well beyond eyp=0.01) and commit.
        let m = m.commit(0.02);
        let Material::ElasticPP { ep, .. } = m else { panic!() };
        assert!((ep - 0.01).abs() < 1e-12, "expected ep=0.01, got {ep}");

        // Partially unload to strain=0.015 (still positive, but now less
        // than the peak). A stateless envelope would still clamp this to
        // the yield force (|100*0.015|=1.5 > fy=1.0), but real plasticity
        // unloads *elastically* from the committed ep, not back through
        // the origin: trial_stress = e*(strain-ep) = 100*(0.015-0.01)=0.5,
        // well inside the yield surface.
        let (stress, tangent) = m.trial_stress_tangent(0.015);
        assert!((stress - 0.5).abs() < 1e-12, "expected elastic unload stress=0.5, got {stress}");
        assert_eq!(tangent, 100.0, "should be elastic (unloading), not the plastic tangent");
    }

    #[test]
    fn gap_only_engages_past_closure() {
        let m = Material::Gap { e: 100.0, gap: 0.01 };
        assert_eq!(m.trial_stress_tangent(0.0), (0.0, 0.0));
        assert_eq!(m.trial_stress_tangent(-0.005), (0.0, 0.0));
        assert_eq!(m.trial_stress_tangent(-0.02), (100.0 * (-0.02 + 0.01), 100.0));
    }

    #[test]
    fn ent_carries_no_tension() {
        let m = Material::Ent { e: 100.0 };
        assert_eq!(m.trial_stress_tangent(0.01), (0.0, 0.0));
        assert_eq!(m.trial_stress_tangent(-0.01), (-1.0, 100.0));
    }
}
