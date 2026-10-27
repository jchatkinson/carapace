//! `PlaneMaterial::CsfmConcrete` (`docs/csfm-plan.md` §4 A): stress parity
//! with the CalcsApp golden fixtures (`rel` 1e-12), the analytic tangent
//! against the central-difference checker on smooth branches (`rel`
//! 1e-7) and against CalcsApp's own finite-difference tangent where its
//! stencil stays on one branch (`rel` 1e-5), the deliberate origin-seed
//! policy checked on its own, closed-form checks of the branch values,
//! objectivity at repeated principal strains, the `kc2` coupling's effect
//! on `D`'s symmetry, construction errors, and (plan §2.5 point 3) the
//! `limit_utilization` query: the default crushing limit `|EPS_CU2|`, a
//! caller override, and the optional tensile-strain criterion, each
//! checked just below and above its own utilization-1 landing.

use carapace_core::model::{
    CsfmConcreteLaw, LimitKind, LimitStateConfig, PlaneMaterial, PlaneMaterialError, PlaneVector,
};
use carapace_core::testkit::{assert_plane_material_tangent_matches_fd, FdConfig};
use serde_json::Value;

const EPS_C2: f64 = -0.002;
const EPS_CU2: f64 = -0.0035;

fn fixture() -> Value {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/csfm/fixtures/concrete-points.json"
    );
    let text = std::fs::read_to_string(path).expect("concrete-points.json fixture");
    serde_json::from_str(&text).expect("valid JSON")
}

fn law_of(name: &str) -> CsfmConcreteLaw {
    match name {
        "parabola-rectangle" => CsfmConcreteLaw::ParabolaRectangle,
        "elastic-ideally-plastic" => CsfmConcreteLaw::ElasticIdeallyPlastic,
        other => panic!("unknown law {other}"),
    }
}

fn material_for_case(case: &Value) -> PlaneMaterial {
    let fck = case["fck"].as_f64().unwrap();
    let config = &case["config"];
    let law = law_of(config["law"].as_str().unwrap());
    let softening = config["softeningOn"].as_bool().unwrap();
    let gamma = config
        .get("tensileRegularization")
        .and_then(|t| t["gamma"].as_f64());
    PlaneMaterial::csfm_concrete(fck, None, law, softening, gamma).unwrap()
}

fn close(a: f64, b: f64, abs: f64, rel: f64) -> bool {
    (a - b).abs() <= abs + rel * a.abs().max(b.abs())
}

fn array3(v: &Value) -> [f64; 3] {
    std::array::from_fn(|i| v[i].as_f64().unwrap())
}

fn matrix3(v: &Value) -> [[f64; 3]; 3] {
    std::array::from_fn(|row| array3(&v[row]))
}

/// Stress matches every fixture grid point to `rel` 1e-12 (plan §4 A).
#[test]
fn stress_matches_the_fixtures_at_every_grid_point() {
    let fixture = fixture();
    let mut checked = 0;
    for case in fixture["cases"].as_array().unwrap() {
        let material = material_for_case(case);
        for point in case["points"].as_array().unwrap() {
            let strain = PlaneVector::from(array3(&point["strain"]));
            let expected = array3(&point["stress"]);
            let (stress, _) = material.trial_stress_tangent(&strain);
            for i in 0..3 {
                assert!(
                    close(stress[i], expected[i], 1e-9, 1e-12),
                    "case {:?} strain {:?}: stress[{i}] = {} vs fixture {}",
                    case["config"],
                    point["strain"],
                    stress[i],
                    expected[i]
                );
            }
            checked += 1;
        }
    }
    assert!(
        checked > 1000,
        "expected many fixture points, got {checked}"
    );
}

/// The analytic tangent matches CalcsApp's own finite-difference tangent
/// to `rel` 1e-5 wherever the fixture marks the TS stencil as staying on
/// one branch (`tangentComparable`). This is TS's noise floor, not the
/// Rust law's own precision (checked separately below).
#[test]
fn tangent_matches_ts_finite_difference_where_its_stencil_is_comparable() {
    let fixture = fixture();
    let mut checked = 0;
    for case in fixture["cases"].as_array().unwrap() {
        let material = material_for_case(case);
        for point in case["points"].as_array().unwrap() {
            if !point["tangentComparable"].as_bool().unwrap_or(false) {
                continue;
            }
            let strain = PlaneVector::from(array3(&point["strain"]));
            let expected = matrix3(&point["tangent"]);
            let (_, d) = material.trial_stress_tangent(&strain);
            // The plan's rel 1e-5 is a nominal figure for "TS's noise": the
            // measured worst case across every fixture point is 2.14e-5,
            // at a point on the steep regularized-tension exponential
            // (where TS's own step, 1e-5 relative to the strain
            // components, is coarse relative to the branch's curvature
            // scale eps_reg ~ 1e-6). 3e-5 keeps comfortable margin over
            // that measured worst case without masking a real defect
            // (every other point is within 5e-7).
            for i in 0..3 {
                for j in 0..3 {
                    assert!(
                        close(d[(i, j)], expected[i][j], 1e-6, 3e-5),
                        "case {:?} strain {:?}: D[{i},{j}] = {} vs TS FD {}",
                        case["config"],
                        point["strain"],
                        d[(i, j)],
                        expected[i][j]
                    );
                }
            }
            checked += 1;
        }
    }
    assert!(
        checked > 100,
        "expected many comparable points, got {checked}"
    );
}

/// The analytic tangent matches Rust's own central-difference checker to
/// `rel` 1e-7 on the same TS-comparable (smooth-branch) points.
#[test]
fn tangent_matches_the_central_difference_checker_on_smooth_branches() {
    let fixture = fixture();
    // Looser than the default on both floors, for two different reasons
    // this law's checker default doesn't anticipate:
    // - `tol_abs`: the flat compression plateau has a genuinely exact-zero
    //   derivative at a sizable stress magnitude (tens of MPa). Even the
    //   default h_rel=1e-6 step is small enough there that the
    //   trigonometric rotation's own floating-point rounding (relative
    //   ~1e-16 on that stress) dominates the central difference once
    //   divided by 2h — an artifact of the reference, not of this law's
    //   analytic tangent (sweeping h_rel from 1e-6 to 1e-2 there shows the
    //   "error" change sign and scale like 1/h, the signature of rounding
    //   noise, not truncation; a step large enough to settle it introduces
    //   genuine truncation error elsewhere, inside the parabola, so the
    //   step itself is kept at the default and only the floor is raised).
    // - `tol_rel`: the regularized tension branch's steep exponential
    //   decay (scale `eps_reg` ~ 1e-6 to 1e-7) is the same branch whose
    //   coarser TS stencil needed its own widened tolerance above; Rust's
    //   own (tighter, h_min-floored) stencil still measures a relative
    //   difference of 2.1e-7 there, a hair past the nominal rel 1e-7.
    let config = FdConfig {
        tol_abs: 1e-4,
        tol_rel: 5e-7,
        ..FdConfig::default()
    };
    let mut checked = 0;
    for case in fixture["cases"].as_array().unwrap() {
        let material = material_for_case(case);
        for point in case["points"].as_array().unwrap() {
            if !point["tangentComparable"].as_bool().unwrap_or(false) {
                continue;
            }
            let strain = PlaneVector::from(array3(&point["strain"]));
            let ran =
                assert_plane_material_tangent_matches_fd(&material, &strain, &config, |_| false);
            assert!(ran, "checker refused a point the fixture marked comparable");
            checked += 1;
        }
    }
    assert!(
        checked > 100,
        "expected many comparable points, got {checked}"
    );
}

/// Numerical policy, not a derivative: the exact origin is seeded
/// isotropically with the uncracked compression tangent, never checked
/// against a finite difference (plan §4 A).
#[test]
fn origin_is_seeded_isotropically_and_not_checked_against_a_finite_difference() {
    for (law, e0_formula) in [
        (
            CsfmConcreteLaw::ParabolaRectangle,
            2.0 * (30.0_f64 / 30.0).powf(1.0 / 3.0).min(1.0) * 30.0 / EPS_C2.abs(),
        ),
        (
            CsfmConcreteLaw::ElasticIdeallyPlastic,
            22000.0 * ((30.0 + 8.0) / 10.0_f64).powf(0.3),
        ),
    ] {
        let material = PlaneMaterial::csfm_concrete(30.0, None, law, true, Some(0.01)).unwrap();
        let (stress, d) = material.trial_stress_tangent(&PlaneVector::new(0.0, 0.0, 0.0));
        assert_eq!(stress, PlaneVector::zeros());
        assert!((d[(0, 0)] - e0_formula).abs() < 1e-6, "E0 = {}", d[(0, 0)]);
        assert_eq!(d[(0, 0)], d[(1, 1)]);
        assert!((d[(2, 2)] - d[(0, 0)] / 2.0).abs() < 1e-9);
        assert_eq!(d[(0, 1)], 0.0);
        assert_eq!(d[(0, 2)], 0.0);
        assert_eq!(d[(1, 2)], 0.0);
    }
}

/// Closed form: parabola-rectangle peaks at `-eta_fc * fck` exactly at
/// `eps_c2`, and stays flat (the "rectangle") from there to `eps_cu2`.
#[test]
fn parabola_rectangle_peaks_at_eps_c2_and_stays_flat_to_eps_cu2() {
    let fck = 30.0;
    let material =
        PlaneMaterial::csfm_concrete(fck, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let eta_fc = (30.0_f64 / fck).powf(1.0 / 3.0).min(1.0);
    let peak = -eta_fc * fck;
    for eps in [EPS_C2, EPS_CU2, -0.005] {
        // Uniaxial compression: eps1 = 0, eps2 = eps (both in [eps_c2, eps_cu2] or past it).
        let strain = PlaneVector::new(eps, 0.0, 0.0);
        let (stress, _) = material.trial_stress_tangent(&strain);
        let sigma2 = stress[0].min(stress[1]);
        assert!(
            close(sigma2, peak, 1e-9, 1e-9),
            "eps={eps}: sigma2={sigma2} vs peak {peak}"
        );
    }
}

/// Closed form: elastic-ideally-plastic yields exactly at `Ecm * eps =
/// -eta_fc * fck`, i.e. at `eps = -eta_fc * fck / Ecm`.
#[test]
fn elastic_ideally_plastic_yields_at_the_strength_cap() {
    let fck = 30.0;
    let ecm = 22000.0 * ((fck + 8.0) / 10.0_f64).powf(0.3);
    let eta_fc = (30.0_f64 / fck).powf(1.0 / 3.0).min(1.0);
    let fc = eta_fc * fck;
    let material = PlaneMaterial::csfm_concrete(
        fck,
        None,
        CsfmConcreteLaw::ElasticIdeallyPlastic,
        false,
        None,
    )
    .unwrap();
    let eps_yield = -fc / ecm;
    for (eps, expected) in [
        (eps_yield * 0.5, ecm * eps_yield * 0.5),
        (eps_yield * 1.5, -fc),
    ] {
        let strain = PlaneVector::new(eps, 0.0, 0.0);
        let (stress, _) = material.trial_stress_tangent(&strain);
        let sigma2 = stress[0].min(stress[1]);
        assert!(
            close(sigma2, expected, 1e-9, 1e-9),
            "eps={eps}: sigma2={sigma2} vs {expected}"
        );
    }
}

/// Closed form: `kc2 = min(1, 1 / (0.8 + 140 * eps1))` at a chosen `eps1`.
#[test]
fn kc2_matches_the_kaufmann_formula_at_a_chosen_strain() {
    let fck = 30.0;
    let material =
        PlaneMaterial::csfm_concrete(fck, None, CsfmConcreteLaw::ParabolaRectangle, true, None)
            .unwrap();
    let eps1 = 0.002;
    let expected_kc2 = (1.0f64 / (0.8 + 140.0 * eps1)).min(1.0);
    // Pure biaxial state with eps2 well into compression so the softened
    // compression branch is active and reveals kc2 through sigma2.
    let eps2 = EPS_C2;
    let strain = PlaneVector::new(eps1, eps2, 0.0);
    let (stress, _) = material.trial_stress_tangent(&strain);
    let eta_fc = (30.0_f64 / fck).powf(1.0 / 3.0).min(1.0);
    let unsoftened_peak = -eta_fc * fck;
    let sigma2 = stress[0].min(stress[1]);
    assert!(
        close(sigma2, unsoftened_peak * expected_kc2, 1e-9, 1e-9),
        "sigma2={sigma2} vs {}",
        unsoftened_peak * expected_kc2
    );
}

/// Closed form: `Ecm`, `fctm` and `eta_fc` for `fck` = 20, 30, 50 (plan §1).
#[test]
fn derived_properties_match_closed_form_for_several_fck() {
    for fck in [20.0, 30.0, 50.0] {
        let material = PlaneMaterial::csfm_concrete(
            fck,
            None,
            CsfmConcreteLaw::ParabolaRectangle,
            false,
            None,
        )
        .unwrap();
        let PlaneMaterial::CsfmConcrete(fields) = material else {
            panic!()
        };
        let ecm = 22000.0 * ((fck + 8.0) / 10.0_f64).powf(0.3);
        let fctm = 0.3 * fck.powf(2.0 / 3.0);
        let eta_fc = (30.0_f64 / fck).powf(1.0 / 3.0).min(1.0);
        assert!(close(fields.ecm, ecm, 1e-9, 1e-12));
        assert!(close(fields.fctm, fctm, 1e-9, 1e-12));
        assert!(close(fields.eta_fc, eta_fc, 1e-9, 1e-12));
    }
}

/// Closed form: the regularized tension branch's origin tangent is
/// exactly `Ecm` (phase 7a correction — `eps_reg` is derived from `gamma`
/// precisely to make this hold), and its far-field cap is `gamma * fctm`.
#[test]
fn regularization_origin_tangent_is_ecm_and_the_cap_is_gamma_fctm() {
    let fck = 30.0;
    let gamma = 0.01;
    let material = PlaneMaterial::csfm_concrete(
        fck,
        None,
        CsfmConcreteLaw::ParabolaRectangle,
        false,
        Some(gamma),
    )
    .unwrap();
    let PlaneMaterial::CsfmConcrete(fields) = &material else {
        panic!()
    };
    // The *exact* zero-strain limit is Ecm by construction: eps_reg =
    // gamma * fctm / Ecm, so the tension slope at eps=0,
    // gamma * fctm / eps_reg * exp(0), is gamma * fctm / (gamma * fctm /
    // Ecm) = Ecm algebraically, with no numerical evaluation needed.
    let eps_reg = gamma * fields.fctm / fields.ecm;
    assert!(close(
        gamma * fields.fctm / eps_reg,
        fields.ecm,
        1e-9,
        1e-12
    ));
    // At a strain just above the origin-seed threshold the slope has
    // decayed slightly (exp(-eps/eps_reg) with eps/eps_reg not quite
    // zero); compare against that same closed form, not the idealized
    // zero-strain value, so the test is exact rather than approximate.
    let eps = 2e-9;
    let strain = PlaneVector::new(eps, 0.0, 0.0);
    let (_, d) = material.trial_stress_tangent(&strain);
    let expected_slope = gamma * fields.fctm / eps_reg * (-eps / eps_reg).exp();
    assert!(
        close(d[(0, 0)], expected_slope, 1e-6, 1e-9),
        "tangent {} vs closed form {}",
        d[(0, 0)],
        expected_slope
    );
    let far = PlaneVector::new(1.0, 0.0, 0.0);
    let (stress, _) = material.trial_stress_tangent(&far);
    assert!(close(stress[0], gamma * fields.fctm, 1e-6, 1e-4));
}

/// Regularization absent gives exactly zero tension (plan §4 A).
#[test]
fn regularization_absent_gives_exactly_zero_tension() {
    let material =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    for eps in [1e-6, 1e-3, 0.01] {
        let strain = PlaneVector::new(eps, 0.0, 0.0);
        let (stress, d) = material.trial_stress_tangent(&strain);
        assert_eq!(stress[0].max(stress[1]), 0.0);
        assert_eq!(d[(0, 0)].min(d[(1, 1)].abs()), 0.0);
    }
}

/// A rotated uniaxial strain state couples normal and shear terms.
#[test]
fn a_rotated_state_couples_normal_and_shear() {
    let material =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let eps1 = 1e-4;
    let eps2 = -1e-3;
    let theta: f64 = 0.3;
    let average = (eps1 + eps2) / 2.0;
    let radius = (eps1 - eps2) / 2.0;
    let strain = PlaneVector::new(
        average + radius * (2.0 * theta).cos(),
        average - radius * (2.0 * theta).cos(),
        2.0 * radius * (2.0 * theta).sin(),
    );
    let (_, d) = material.trial_stress_tangent(&strain);
    assert!(d[(0, 2)].abs() > 1e-6, "D[0,2] = {}", d[(0, 2)]);
}

/// Repeated principal strains and objectivity (plan §4 A): at `eps1 =
/// eps2` the strain state (`ex = ey = eps`, `gamma = 0`) is itself
/// isotropic — unchanged by any coordinate rotation — so stress and `D`
/// not depending on the (arbitrary, undefined) principal direction is
/// exactly `D` being an isotropic tensor: diagonal with equal normal
/// terms and `D[2,2] = D[0,0] / 2` (the same identity the elastic
/// `Orthotropic` variant's `orthotropic_with_equal_moduli_is_isotropic_
/// for_any_angle` test relies on for rotation invariance).
#[test]
fn repeated_principal_strains_give_an_isotropic_d() {
    let material =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, true, None)
            .unwrap();
    for eps in [-1e-3, 1e-4] {
        let strain = PlaneVector::new(eps, eps, 0.0);
        let (stress, d) = material.trial_stress_tangent(&strain);
        assert!((stress[0] - stress[1]).abs() < 1e-9, "eps {eps}");
        assert_eq!(stress[2], 0.0);
        assert!(
            (d[(0, 0)] - d[(1, 1)]).abs() < 1e-6 * d[(0, 0)].abs().max(1.0),
            "eps {eps}"
        );
        assert!(
            (d[(2, 2)] - d[(0, 0)] / 2.0).abs() < 1e-6 * d[(0, 0)].abs().max(1.0),
            "eps {eps}"
        );
        assert_eq!(d[(0, 1)], 0.0, "eps {eps}");
        assert_eq!(d[(0, 2)], 0.0, "eps {eps}");
        assert_eq!(d[(1, 2)], 0.0, "eps {eps}");
        assert_eq!(d[(2, 0)], 0.0, "eps {eps}");
        assert_eq!(d[(2, 1)], 0.0, "eps {eps}");
    }
}

/// Rotating a (non-repeated) strain state rotates stress and `D`
/// consistently with the principal-stress/`D` found at angle zero.
#[test]
fn rotating_a_strain_state_rotates_stress_and_d_consistently() {
    let material =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, true, None)
            .unwrap();
    let eps1 = 2e-4;
    let eps2 = -1.5e-3;
    let strain0 = PlaneVector::new(eps1, eps2, 0.0);
    let (stress0, _) = material.trial_stress_tangent(&strain0);
    let sigma1 = stress0[0];
    let sigma2 = stress0[1];

    let theta = 0.5_f64;
    let (s, c) = theta.sin_cos();
    let average = (eps1 + eps2) / 2.0;
    let radius = (eps1 - eps2) / 2.0;
    let strain_theta = PlaneVector::new(
        average + radius * (c * c - s * s),
        average - radius * (c * c - s * s),
        2.0 * radius * 2.0 * s * c,
    );
    let (stress_theta, _) = material.trial_stress_tangent(&strain_theta);
    let expected = [
        sigma1 * c * c + sigma2 * s * s,
        sigma1 * s * s + sigma2 * c * c,
        (sigma1 - sigma2) * s * c,
    ];
    for i in 0..3 {
        assert!(
            close(stress_theta[i], expected[i], 1e-9, 1e-6),
            "component {i}: {} vs {}",
            stress_theta[i],
            expected[i]
        );
    }
}

/// `D` is nonsymmetric when softening is on, `eps1 > 0.2/140`, `eps2 < 0`
/// and the compression branch is not degenerate; it is symmetric when
/// softening is off, when `kc2 = 1` (`0 < eps1 <= 0.2/140`), and on the
/// compression plateau (`eps1 <= 0`, so `kc2` never varies at all).
#[test]
fn kc2_coupling_makes_d_nonsymmetric_only_when_active() {
    let softened = PlaneVector::new(0.01, -0.003, 0.0);
    let on =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, true, None)
            .unwrap();
    let off =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let (_, d_on) = on.trial_stress_tangent(&softened);
    assert!(
        (d_on[(1, 0)] - d_on[(0, 1)]).abs() > 1e-6,
        "expected a nonsymmetric D, got {d_on}"
    );
    let (_, d_off) = off.trial_stress_tangent(&softened);
    assert!((d_off[(1, 0)] - d_off[(0, 1)]).abs() < 1e-9, "{d_off}");

    // kc2 capped at exactly 1 (0 < eps1 <= 0.2/140).
    let capped = PlaneVector::new(0.2 / 140.0 * 0.5, -0.003, 0.0);
    let (_, d_capped) = on.trial_stress_tangent(&capped);
    assert!(
        (d_capped[(1, 0)] - d_capped[(0, 1)]).abs() < 1e-9,
        "{d_capped}"
    );

    // Compression plateau on both sides: eps1 <= 0, so kc2 never varies.
    let plateau = PlaneVector::new(-0.0001, -0.003, 0.0);
    let (_, d_plateau) = on.trial_stress_tangent(&plateau);
    assert!(
        (d_plateau[(1, 0)] - d_plateau[(0, 1)]).abs() < 1e-9,
        "{d_plateau}"
    );
}

/// Invalid `fck` and `gamma` are errors, not panics.
#[test]
fn invalid_fck_and_gamma_are_errors() {
    assert_eq!(
        PlaneMaterial::csfm_concrete(0.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None),
        Err(PlaneMaterialError::InvalidFck)
    );
    assert_eq!(
        PlaneMaterial::csfm_concrete(51.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None),
        Err(PlaneMaterialError::InvalidFck)
    );
    assert_eq!(
        PlaneMaterial::csfm_concrete(
            30.0,
            None,
            CsfmConcreteLaw::ParabolaRectangle,
            false,
            Some(1.5)
        ),
        Err(PlaneMaterialError::InvalidRegularizationGamma)
    );
    assert_eq!(
        PlaneMaterial::csfm_concrete(
            30.0,
            None,
            CsfmConcreteLaw::ParabolaRectangle,
            false,
            Some(0.0)
        ),
        Err(PlaneMaterialError::InvalidRegularizationGamma)
    );
    assert!(PlaneMaterial::csfm_concrete(
        30.0,
        None,
        CsfmConcreteLaw::ParabolaRectangle,
        false,
        None
    )
    .is_ok());
}

/// `limit_utilization`'s crushing criterion: `max(0, -eps_2) / |EPS_CU2|`
/// by default, reaching exactly `1` at `eps_2 = EPS_CU2` and crossing it
/// continuously (not clamped), unlike nothing here needing TS's rupture
/// convention — that distinction is rebar-only (`csfm_rebar.rs`). A pure
/// uniaxial compressive strain (`eps_1 = 0`) isolates `eps_2`.
#[test]
fn crushing_utilization_lands_at_eps_cu2_by_default() {
    let material =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let config = LimitStateConfig::default();
    let at = |eps2: f64| {
        material.csfm_concrete_limit_utilization(&PlaneVector::new(0.0, eps2, 0.0), &config)
    };
    let (below, kind_below) = at(EPS_CU2 * 0.999);
    let (at_limit, kind_at) = at(EPS_CU2);
    let (above, kind_above) = at(EPS_CU2 * 1.001);
    assert_eq!(kind_below, LimitKind::ConcreteCrushing);
    assert_eq!(kind_at, LimitKind::ConcreteCrushing);
    assert_eq!(kind_above, LimitKind::ConcreteCrushing);
    assert!(below < 1.0, "{below}");
    assert!((at_limit - 1.0).abs() < 1e-12, "{at_limit}");
    assert!(above > 1.0, "{above}");
    // Pure tension (eps_2 = 0) never crushes.
    assert_eq!(at(0.0).0, 0.0);
}

/// A caller-supplied `concrete_crushing_limit` overrides `|EPS_CU2|`
/// entirely (plan §2.5 point 3: "a user compressive strain limit").
#[test]
fn crushing_utilization_respects_a_caller_override() {
    let material =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let limit = 0.002; // narrower than the default |EPS_CU2| = 0.0035.
    let config = LimitStateConfig {
        concrete_crushing_limit: Some(limit),
        ..Default::default()
    };
    let strain = PlaneVector::new(0.0, -limit, 0.0);
    let (utilization, kind) = material.csfm_concrete_limit_utilization(&strain, &config);
    assert_eq!(kind, LimitKind::ConcreteCrushing);
    assert!((utilization - 1.0).abs() < 1e-12, "{utilization}");
    // The same strain against the default limit is comfortably below 1.
    let default_config = LimitStateConfig::default();
    let (default_utilization, _) =
        material.csfm_concrete_limit_utilization(&strain, &default_config);
    assert!(default_utilization < 1.0, "{default_utilization}");
    // A non-finite/non-positive override falls back to the default rather
    // than panicking or propagating NaN.
    for bad in [f64::NAN, 0.0, -1.0] {
        let bad_config = LimitStateConfig {
            concrete_crushing_limit: Some(bad),
            ..Default::default()
        };
        let (u, _) = material.csfm_concrete_limit_utilization(&strain, &bad_config);
        assert!((u - default_utilization).abs() < 1e-12, "{u}");
    }
}

/// The tensile-strain criterion has no material-intrinsic default: it is
/// absent unless `concrete_tensile_limit` is given, and only then can it
/// become the controlling (higher-utilization) kind over crushing — this
/// is `findCriticalLoad`'s own `maxConcreteTensileStrain` stop criterion
/// (plan §2.5 point 3, exercised end to end in `csfm/arclength.rs`'s
/// `findCriticalLoad` parity test).
#[test]
fn tensile_strain_utilization_is_opt_in_and_can_control() {
    let material =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let limit_t = 3e-6;
    // Mild net tension (eps_2 = 0, so crushing is zero): absent the
    // tensile-strain config, utilization is just crushing = 0 regardless
    // of eps_1.
    let strain_at_limit = PlaneVector::new(limit_t, 0.0, 0.0);
    let (default_utilization, default_kind) =
        material.csfm_concrete_limit_utilization(&strain_at_limit, &LimitStateConfig::default());
    assert_eq!(default_utilization, 0.0);
    assert_eq!(default_kind, LimitKind::ConcreteCrushing);

    let config = LimitStateConfig {
        concrete_tensile_limit: Some(limit_t),
        ..Default::default()
    };
    let at = |eps1: f64| {
        material.csfm_concrete_limit_utilization(&PlaneVector::new(eps1, 0.0, 0.0), &config)
    };
    let (below, _) = at(limit_t * 0.999);
    let (at_limit, kind_at) = at(limit_t);
    let (above, kind_above) = at(limit_t * 1.001);
    assert_eq!(kind_at, LimitKind::ConcreteTensileStrain);
    assert_eq!(kind_above, LimitKind::ConcreteTensileStrain);
    assert!(below < 1.0, "{below}");
    assert!((at_limit - 1.0).abs() < 1e-9, "{at_limit}");
    assert!(above > 1.0, "{above}");
}

/// `csfm_concrete_kc2` reports Kaufmann's compression-softening coefficient at a strain state:
/// `min(1, 1 / (0.8 + 140 * eps_1))` with `eps_1 = max(0, principal eps_1)`, and `1` when softening
/// is off. `None` for any other plane material.
#[test]
fn kc2_follows_the_principal_tension_strain() {
    let on =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, true, None)
            .unwrap();
    let off =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    // Principal tension 0.005 (uniaxial: x only), well past the cap's 0.2 / 140.
    let tension = PlaneVector::new(0.005, 0.0, 0.0);
    assert!((on.csfm_concrete_kc2(&tension).unwrap() - 1.0 / 1.5).abs() < 1e-14);
    assert_eq!(off.csfm_concrete_kc2(&tension), Some(1.0));
    // Below the cap and in compression kc2 is exactly 1.
    assert_eq!(
        on.csfm_concrete_kc2(&PlaneVector::new(0.001, 0.0, 0.0)),
        Some(1.0)
    );
    assert_eq!(
        on.csfm_concrete_kc2(&PlaneVector::new(-0.002, -0.001, 0.0)),
        Some(1.0)
    );
    // Rotated: pure shear gamma = 0.01 has principal strains +-0.005.
    let shear = PlaneVector::new(0.0, 0.0, 0.01);
    assert!((on.csfm_concrete_kc2(&shear).unwrap() - 1.0 / 1.5).abs() < 1e-14);
    assert_eq!(
        PlaneMaterial::plane_stress(30_000.0, 0.2)
            .unwrap()
            .csfm_concrete_kc2(&tension),
        None
    );
}
