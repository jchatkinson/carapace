use super::Material;

impl Material {
    /// Equal-factor (`1.0` each) parallel composite. For custom factors,
    /// build the `Material::Parallel(Vec<(Material, f64)>)` variant
    /// directly.
    pub fn parallel(materials: Vec<Material>) -> Self {
        Material::Parallel(materials.into_iter().map(|m| (m, 1.0)).collect())
    }

    /// Series composite with OpenSees' default tolerance (`1e-10`) and a
    /// generous internal iteration cap (`30`) — since, unlike OpenSees'
    /// `SeriesMaterial`, `evaluate` can't lean on many external Newton
    /// calls to slowly converge a warm-started trial state (see the type's
    /// doc), the internal loop needs to reach convergence by itself.
    pub fn series(children: Vec<Material>) -> Self {
        let n = children.len();
        Material::Series {
            children,
            child_strains: vec![0.0; n],
            stress: 0.0,
            max_iter: 30,
            tol: 1e-10,
        }
    }

    pub fn min_max(inner: Material, min_strain: f64, max_strain: f64) -> Self {
        Material::MinMax {
            inner: Box::new(inner),
            min_strain,
            max_strain,
            failed: false,
        }
    }
}

pub(super) fn evaluate_parallel(children: &[(Material, f64)], strain: f64) -> (f64, f64, Material) {
    let mut stress = 0.0;
    let mut tangent = 0.0;
    let mut next_children = Vec::with_capacity(children.len());
    for (child, factor) in children {
        let (s, t, next) = child.evaluate(strain);
        stress += factor * s;
        tangent += factor * t;
        next_children.push((next, *factor));
    }
    (stress, tangent, Material::Parallel(next_children))
}

pub(super) fn evaluate_min_max(
    inner: &Material,
    min_strain: f64,
    max_strain: f64,
    failed: bool,
    strain: f64,
) -> (f64, f64, Material) {
    if failed || strain <= min_strain || strain >= max_strain {
        let tangent = 1.0e-8 * inner.initial_tangent();
        (
            0.0,
            tangent,
            Material::MinMax {
                inner: Box::new(inner.clone()),
                min_strain,
                max_strain,
                failed: true,
            },
        )
    } else {
        let (s, t, next_inner) = inner.evaluate(strain);
        (
            s,
            t,
            Material::MinMax {
                inner: Box::new(next_inner),
                min_strain,
                max_strain,
                failed: false,
            },
        )
    }
}

/// Ported from `SeriesMaterial::setTrialStrain` (`SeriesMaterial.cpp`)
/// — flexibility-based iteration, but self-contained: applies the
/// stress correction `dq` every inner iteration (a proper Newton
/// update) rather than only once after the loop, since — unlike
/// OpenSees, which carries `Tstress`/`Ttangent` across many external
/// `setTrialStrain` calls to slowly converge — this must converge
/// within one `evaluate` call. See the type's doc comment and the M7
/// stage-2 handoff notes.
pub(super) fn evaluate_series(m: &Material, total_strain: f64) -> (f64, f64, Material) {
    let Material::Series {
        children,
        child_strains,
        stress: cstress,
        max_iter,
        tol,
    } = m
    else {
        unreachable!()
    };

    let n = children.len();
    let mut e = child_strains.clone();
    let mut child_stress = vec![0.0; n];
    let mut flex = vec![0.0; n];
    let mut tstress = *cstress;
    let mut series_tangent = 0.0;

    let safe_inv = |t: f64| {
        if t.abs() > 1.0e-12 {
            1.0 / t
        } else if t < 0.0 {
            -1.0e12
        } else {
            1.0e12
        }
    };

    for _ in 0..*max_iter {
        let mut f = 0.0;
        let mut vr = 0.0;
        for i in 0..n {
            let ds = tstress - child_stress[i];
            e[i] += flex[i] * ds;
            let (s, t) = children[i].trial_stress_tangent(e[i]);
            child_stress[i] = s;
            flex[i] = safe_inv(t);
            let ds = tstress - child_stress[i];
            let de = flex[i] * ds;
            f += flex[i];
            vr += e[i] + de;
        }
        series_tangent = safe_inv(f);
        let dv = total_strain - vr;
        let dq = series_tangent * dv;
        tstress += dq;
        if (dq * dv).abs() < *tol {
            break;
        }
    }

    let next_children = children
        .iter()
        .zip(e.iter())
        .map(|(child, &strain)| child.evaluate(strain).2)
        .collect();

    (
        tstress,
        series_tangent,
        Material::Series {
            children: next_children,
            child_strains: e,
            stress: tstress,
            max_iter: *max_iter,
            tol: *tol,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parallel_of_two_elastics_equals_one_elastic_with_summed_stiffness() {
        let combined = Material::parallel(vec![Material::Elastic { e: 100.0 }, Material::Elastic { e: 50.0 }]);
        let single = Material::Elastic { e: 150.0 };
        assert_eq!(combined.trial_stress_tangent(0.01), single.trial_stress_tangent(0.01));
    }

    #[test]
    fn parallel_history_is_independent_per_child() {
        // Two ElasticPP with different yield points in parallel: combined
        // response is the sum of each independently-clamped branch.
        let combined = Material::parallel(vec![Material::elastic_pp(100.0, 0.01), Material::elastic_pp(50.0, 0.02)]);
        let (stress, _tangent) = combined.trial_stress_tangent(0.05);
        // Both branches are well past yield: 100*0.01 + 50*0.02 = 2.0.
        assert!((stress - 2.0).abs() < 1e-9, "got {stress}");
    }

    #[test]
    fn series_of_two_equal_elastics_halves_the_stiffness() {
        let series = Material::series(vec![Material::Elastic { e: 100.0 }, Material::Elastic { e: 100.0 }]);
        let (stress, tangent) = series.trial_stress_tangent(0.02);
        // Two equal springs in series: combined stiffness = E/2, so for a
        // given total strain the stress equals a single spring of
        // stiffness 50 at that same total strain.
        assert!((tangent - 50.0).abs() < 1e-6, "expected tangent 50, got {tangent}");
        assert!((stress - 50.0 * 0.02).abs() < 1e-6, "got {stress}");
    }

    #[test]
    fn series_splits_strain_so_each_child_sees_equal_stress() {
        let series = Material::series(vec![Material::Elastic { e: 100.0 }, Material::Elastic { e: 300.0 }]);
        let Material::Series { child_strains, .. } = series.commit(0.04) else {
            panic!()
        };
        // Equilibrium requires equal stress: 100*e0 = 300*e1, and
        // compatibility requires e0+e1=0.04 => e0=0.03, e1=0.01.
        assert!((child_strains[0] - 0.03).abs() < 1e-6, "got {child_strains:?}");
        assert!((child_strains[1] - 0.01).abs() < 1e-6, "got {child_strains:?}");
    }

    #[test]
    fn min_max_fails_permanently_once_strain_exits_bounds() {
        let m = Material::min_max(Material::Elastic { e: 100.0 }, -0.01, 0.01);
        let m = m.commit(0.02); // exceeds max_strain
        let Material::MinMax { failed, .. } = m else { panic!() };
        assert!(failed);

        let (stress, tangent) = m.trial_stress_tangent(0.0); // back within bounds
        assert_eq!(stress, 0.0, "failed material must read zero stress even if strain re-enters bounds");
        assert!((tangent - 1.0e-8 * 100.0).abs() < 1e-12, "got {tangent}");
    }

    #[test]
    fn min_max_passes_through_inner_material_before_failure() {
        let m = Material::min_max(Material::Elastic { e: 100.0 }, -0.01, 0.01);
        assert_eq!(m.trial_stress_tangent(0.005), (0.5, 100.0));
    }
}
