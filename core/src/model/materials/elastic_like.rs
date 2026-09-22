use super::Material;

impl Material {
    pub fn elastic_pp(e: f64, eyp: f64) -> Self {
        Material::ElasticPP { e, eyp, ep: 0.0 }
    }
}

pub(super) fn evaluate_elastic(e: f64, strain: f64) -> (f64, f64, Material) {
    (e * strain, e, Material::Elastic { e })
}

pub(super) fn evaluate_elastic_pp(e: f64, eyp: f64, ep: f64, strain: f64) -> (f64, f64, Material) {
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

pub(super) fn evaluate_gap(e: f64, gap: f64, strain: f64) -> (f64, f64, Material) {
    // `<=`, not `<`: at exactly zero closure the tangent must
    // stay nonzero (matching the just-closed regime) or a
    // system starting from zero displacement forms a singular
    // initial tangent — there is no elastic-vs-open distinction
    // to make at a single point anyway.
    let closure = strain + gap;
    if closure <= 0.0 {
        (e * closure, e, Material::Gap { e, gap })
    } else {
        (0.0, 0.0, Material::Gap { e, gap })
    }
}

pub(super) fn evaluate_ent(e: f64, strain: f64) -> (f64, f64, Material) {
    // See `evaluate_gap` above for why this is `<=`, not `<`.
    if strain <= 0.0 {
        (e * strain, e, Material::Ent { e })
    } else {
        (0.0, 0.0, Material::Ent { e })
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
