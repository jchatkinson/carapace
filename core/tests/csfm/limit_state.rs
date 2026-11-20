//! The `material_limit` stop criterion (`docs/csfm-plan.md` §2.5 point 3)
//! on `Integrator::ArcLength`, and the equivalent per-step diagnostic
//! (`StepResult::material_limit`) on plain `Integrator::LoadControl`. A
//! minimal one-`Tri3`/one-`EmbeddedBar` fixture, pulled in tension until
//! the rebar's rupture utilization reaches `1`. The thorough exit-criterion
//! coverage (one model per `LimitKind`, exact landing checked just below
//! and above) is `arclength.rs`'s own material-limit tests; this file
//! covers the mechanism itself.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, AnalysisError, ArcDirection, ArcLength, ArcScales, ArcStopCriteria,
    ConstraintHandler, ConvergenceTest, Integrator, MaterialLimitTarget, StopReason,
    TangentStrategy,
};
use carapace_core::model::{
    CsfmConcreteLaw, Domain, Element, EmbeddedBar, LimitKind, LimitStateConfig, Material, Node,
    NodeId, PlaneMaterial, Tri3,
};

const LEG_LENGTH: f64 = 100.0;

/// A `Tri3` (`CsfmConcrete`, `fck = 30`, exact zero tension: `gamma =
/// None`) with legs of length `LEG_LENGTH` along `x` and `y`, one
/// `EmbeddedBar` (`CsfmRebar`, tension stiffening on, area `1.0`) along
/// the `x` leg. Every DOF but node 1's `ux` is fixed (node 0 fully, node
/// 1's `uy`, node 2 fully), so this reduces to one free DOF: a nonlinear
/// spring combining the concrete's (CST) and the bar's axial stiffness in
/// parallel along `x`, with a unit reference load at node 1's `ux`.
///
/// Past the tiny origin-seed strain, the concrete's principal direction
/// aligns exactly with `x` (`eps_y = gamma_xy = 0` throughout — this
/// one-DOF strain state has no shear or transverse component), so its
/// zero-tension law gives it **exactly** zero stress and tangent there:
/// the bar alone carries the tension past cracking, and (since the bar's
/// own tangent is clamped to `[Es/1e4, 10 Es]`, always strictly positive)
/// the one-DOF system is never singular regardless of the concrete's
/// state — no tensile regularization is needed to keep this fixture
/// well-posed.
///
/// Closed form (used only to calibrate the test schedules below, not
/// asserted as a tolerance): the bar's own coefficient already folds in
/// `1 / LEG_LENGTH` (`EmbeddedBar`'s `coeffs`), and its resistance is
/// `stress * area * LEG_LENGTH * coeffs`, so that `LEG_LENGTH` cancels
/// exactly and (with `area = 1`, `direction = [1, 0]`) internal resistance
/// at `ux1` is just `stress_bar` (the Tri3 contributes exactly `0` there
/// past cracking). Equilibrium with the unit reference load therefore
/// gives `stress_bar = lambda` directly, so rupture (`stress_bar = ft =
/// 540`) is at `lambda = 540`.
fn build() -> (Domain, f64, f64) {
    let coords = [[0.0, 0.0], [LEG_LENGTH, 0.0], [0.0, LEG_LENGTH]];
    let mut domain: Domain = Domain::new();
    let ids: [NodeId; 3] = std::array::from_fn(|i| {
        let mut node = Node::new(coords[i]);
        node = match i {
            0 => node.fix(0).fix(1),
            1 => node.fix(1),
            _ => node.fix(0).fix(1),
        };
        domain.add_node(node)
    });
    domain.add_element(Element::Tri3(Tri3::new(
        ids[0],
        ids[1],
        ids[2],
        1.0,
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap(),
    )));
    let rebar = Material::csfm_rebar(
        200_000.0, 500.0, 540.0, 0.05, 2.9, 32000.0, 16.0, 0.02, 0.67, true,
    )
    .unwrap();
    let Material::CsfmRebar(fields) = &rebar else {
        panic!()
    };
    let eps_m_at_ft = fields.eps_m_at_ft;
    domain.add_element(Element::EmbeddedBar(EmbeddedBar::in_tri3(
        ids,
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.0],
        LEG_LENGTH,
        1.0,
        rebar,
    )));
    domain.load_node(ids[1], 0, 1.0);
    (domain, eps_m_at_ft, LEG_LENGTH)
}

fn combined() -> ConvergenceTest {
    ConvergenceTest::Combined {
        force_tol: 1e-9,
        moment_tol: 1e-9,
        relative_tol: 1e-9,
        displacement_tol: None,
        max_iter: 30,
    }
}

/// `Integrator::ArcLength` with `material_limit: Some({ exact: true })`
/// stops exactly when the rebar's rupture utilization reaches `1`
/// (closed-form at `lambda = 540`, per `build`'s doc comment).
#[test]
fn arc_length_lands_exactly_on_the_rebar_rupture_utilization() {
    let (domain, eps_m_at_ft, length) = build();
    let displacement_scale = eps_m_at_ft * length; // u1x at rupture.
    let mut arc = ArcLength::fixed(
        0.1,
        ArcScales::Explicit {
            displacement: displacement_scale,
            rotation: None,
            load: 540.0,
        },
    );
    arc.direction = ArcDirection::Increasing;
    arc.stop = ArcStopCriteria {
        material_limit: Some(Box::new(MaterialLimitTarget {
            limits: LimitStateConfig::default(),
            exact: true,
        })),
        ..Default::default()
    };
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::ArcLength(arc))
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(combined())
        .build(domain);

    let mut stopped = false;
    for _ in 0..200 {
        let result = analysis.step().expect("arc-length step should converge");
        let arc_info = result.arc.expect("ArcLength reports continuation info");
        let (utilization, kind) = result
            .material_limit
            .expect("the model has a CSFM limit state throughout");
        assert_eq!(kind, LimitKind::RebarRupture);
        if let Some(stop) = arc_info.stop {
            assert_eq!(stop.reason, StopReason::MaterialLimit);
            assert!(stop.landed_exactly, "exact landing was requested");
            assert!(
                (utilization - 1.0).abs() < 1e-6,
                "landed utilization {utilization}"
            );
            assert!(
                (result.load_factor - 540.0).abs() < 1.0,
                "landed lambda {} should be near the closed-form 540",
                result.load_factor
            );
            stopped = true;
            break;
        }
    }
    assert!(stopped, "the phase must stop once utilization reaches 1");
}

/// Plain `Integrator::LoadControl` reports the same controlling
/// utilization and kind after every accepted step (plan §2.5 point 3's
/// "equivalent check after each accepted load-control step"), rising
/// monotonically with the load factor and approaching (but, with no
/// bordered corrector, never exactly landing on) utilization `1` — and
/// once the schedule asks for more than the bar's rupture capacity (its
/// rupture plateau pins stress at `ft`, plan §1, so this one-DOF system's
/// resistance is bounded above by `ft` exactly), load control simply
/// fails to converge, exactly as `arclength_softening.rs`'s baseline
/// load-control test fails past an ordinary softening peak. This is the
/// mechanism §2.5 point 4's arc-length-past-the-peak criterion exists for.
#[test]
fn load_control_approaches_but_cannot_land_on_or_pass_rupture() {
    let (domain, _eps_m_at_ft, _length) = build();
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 5.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(combined())
        .build(domain);

    let mut last_utilization = 0.0;
    let mut failed = false;
    for _ in 0..150 {
        match analysis.step() {
            Ok(result) => {
                let (utilization, kind) = result
                    .material_limit
                    .expect("the model has a CSFM limit state throughout");
                assert_eq!(kind, LimitKind::RebarRupture);
                assert!(
                    utilization >= last_utilization - 1e-9,
                    "utilization should rise monotonically under monotonically increasing load: {utilization} < {last_utilization}"
                );
                assert!(
                    utilization <= 1.0 + 1e-6,
                    "this one-DOF system's resistance cannot exceed the bar's rupture capacity: {utilization}"
                );
                last_utilization = utilization;
            }
            Err(AnalysisError::FailedToConverge { .. }) => {
                failed = true;
                break;
            }
            Err(other) => panic!("unexpected error: {other:?}"),
        }
    }
    assert!(
        failed,
        "load control must eventually fail once the schedule exceeds the bar's capacity"
    );
    assert!(
        last_utilization > 0.99,
        "the last converged step should be right at the rupture limit: {last_utilization}"
    );
}
