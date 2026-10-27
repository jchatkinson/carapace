//! §4 G "Material limits with exact landing" (plan `docs/csfm-plan.md`):
//! one small, single-purpose Rust-only model per `LimitKind` sub-case (no
//! TS parity needed, per the plan) — a lightly reinforced tie (tension
//! rupture, TCM/POM routing), an over-reinforced strut (concrete
//! crushing), and a bare-steel case (`|eps| = eps_u`). Each runs under
//! `Integrator::ArcLength` with `material_limit: Some({ exact: true })`
//! and must land on utilization `1` within the landing tolerance;
//! utilization is then checked just below and just above that landing to
//! distinguish exact landing (`== 1`, not clamped) from TS's strict `>`
//! rupture semantics (`Material::csfm_rebar_ruptured`).
//!
//! All three share one geometry (see `one_dof_tie_or_strut`): a `Tri3`
//! with legs along `x`/`y`, one `EmbeddedBar` along the `x` leg, every DOF
//! but node 1's `ux` fixed. This reduces every model to one free DOF, so
//! each one's own closed form (worked out in each test's doc comment) is
//! exact, not approximate — the same reduction `limit_state.rs`'s
//! mechanism tests use.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ArcDirection, ArcLength, ArcScales, ArcStopCriteria,
    ConstraintHandler, ConvergenceTest, Integrator, MaterialLimitTarget, StopReason,
    TangentStrategy,
};
use carapace_core::model::{
    CsfmConcreteLaw, Domain, Element, EmbeddedBar, LimitKind, LimitStateConfig, Material, Node,
    NodeId, PlaneMaterial, Tri3,
};

const LEG_LENGTH: f64 = 100.0;

/// One `Tri3` (`CsfmConcrete`) with legs of `LEG_LENGTH` along `x`/`y` and
/// one `EmbeddedBar` (`material`, area `bar_area`) along the `x` leg.
/// Every DOF but node 1's `ux` is fixed, so this reduces to one free DOF
/// (see `limit_state.rs::build`'s doc comment for why that makes the
/// concrete contribute **exactly** zero resistance past cracking in pure
/// tension, and why the reduced system is never singular regardless). A
/// unit reference load acts at node 1's `ux`; a positive load pulls the
/// bar in tension, negative pushes the leg into compression.
fn one_dof_tie_or_strut(
    concrete: PlaneMaterial,
    material: Material,
    bar_area: f64,
) -> (
    Domain,
    carapace_core::model::ElementId,
    carapace_core::model::ElementId,
) {
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
    let concrete_id = domain.add_element(Element::Tri3(Tri3::new(
        ids[0], ids[1], ids[2], 1.0, concrete,
    )));
    let bar_id = domain.add_element(Element::EmbeddedBar(EmbeddedBar::in_tri3(
        ids,
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.0],
        LEG_LENGTH,
        bar_area,
        material,
    )));
    domain.load_node(ids[1], 0, 1.0);
    (domain, concrete_id, bar_id)
}

fn combined() -> ConvergenceTest {
    ConvergenceTest::Combined {
        force_tol: 1e-8,
        moment_tol: 1e-8,
        relative_tol: 1e-9,
        displacement_tol: None,
        max_iter: 30,
    }
}

/// Runs `domain` under `Integrator::ArcLength` with `material_limit: Some
/// ({ exact: true })` (increasing direction) until it lands, returning the
/// landed `(load_factor, bar_element, controlling_kind)`. Panics if it
/// doesn't land within `step_budget` steps.
fn land_on_material_limit(
    domain: Domain,
    scales: ArcScales,
    radius: (f64, f64, f64),
    step_budget: usize,
) -> (
    carapace_core::analysis::Analysis,
    f64,
    carapace_core::model::ElementId,
    LimitKind,
) {
    let mut arc = ArcLength::adaptive(radius.0, radius.1, radius.2, scales);
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
        .try_build(domain)
        .unwrap_or_else(|e| panic!("invalid model: {e:?}"));

    for _ in 0..step_budget {
        let result = analysis
            .step()
            .unwrap_or_else(|e| panic!("arc-length step failed before landing: {e:?}"));
        let arc_info = result.arc.expect("ArcLength reports continuation info");
        if let Some(stop) = arc_info.stop {
            assert_eq!(
                stop.reason,
                StopReason::MaterialLimit,
                "unexpected stop reason"
            );
            assert!(stop.landed_exactly, "exact landing was requested");
            let (element, utilization, kind) = analysis
                .domain()
                .controlling_limit_utilization(&LimitStateConfig::default())
                .expect("a limit state must still be present at the landed step");
            assert!(
                (utilization - 1.0).abs() < 1e-6,
                "landed utilization should be 1: {utilization}"
            );
            return (analysis, result.load_factor, element, kind);
        }
    }
    panic!("did not land on the material limit within {step_budget} steps");
}

fn rebar_fields(material: &Material) -> &carapace_core::model::CsfmRebarFields {
    match material {
        Material::CsfmRebar(fields) => fields,
        _ => panic!("expected CsfmRebar"),
    }
}

/// **The lightly reinforced tie** (tension rupture, TCM routing since
/// `rho_eff = 0.02 > rho_crit`). Closed form (concrete contributes
/// exactly `0` resistance past cracking in this pure-tension state, same
/// as `limit_state.rs::build`): internal resistance at `ux1` is `stress_bar
/// * area_bar` undivided by `LEG_LENGTH` (the bar's own `1/LEG_LENGTH`
/// coefficient cancels the `* LEG_LENGTH` in its resistance formula), so
/// with `area_bar = 1`, equilibrium gives `stress_bar = lambda`, landing
/// at `lambda = ft = 540`.
#[test]
fn lightly_reinforced_tie_lands_on_rebar_rupture() {
    let rebar = Material::csfm_rebar(
        200_000.0, 500.0, 540.0, 0.05, 2.9, 32000.0, 16.0, 0.02, 0.67, true,
    )
    .unwrap();
    let eps_m_at_ft = rebar_fields(&rebar).eps_m_at_ft;
    let concrete =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let (domain, _concrete_id, bar_id) = one_dof_tie_or_strut(concrete, rebar, 1.0);

    let (analysis, lambda, landed_id, kind) = land_on_material_limit(
        domain,
        ArcScales::Explicit {
            displacement: eps_m_at_ft * LEG_LENGTH,
            rotation: None,
            load: 540.0,
        },
        (0.1, 0.01, 0.5),
        200,
    );
    assert_eq!(kind, LimitKind::RebarRupture);
    assert_eq!(landed_id, bar_id, "the bar must be the controlling element");
    assert!((lambda - 540.0).abs() < 1.0, "lambda {lambda}");

    let strain = analysis.domain().element_gauss_responses(bar_id).unwrap()[0].strain[0];
    assert!(
        (strain - eps_m_at_ft).abs() < 1e-6 * eps_m_at_ft,
        "strain {strain}"
    );
    let carapace_core::model::Element::EmbeddedBar(bar) = analysis.domain().element(bar_id) else {
        panic!()
    };
    assert!(
        !bar.material().csfm_rebar_ruptured(strain),
        "landed utilization 1 must not read as ruptured under TS's strict > threshold"
    );
    assert!(
        bar.material().csfm_rebar_ruptured(strain * (1.0 + 1e-6)),
        "just past the landing, ruptured() must flip true"
    );
    let (just_below, _) = bar
        .material()
        .csfm_rebar_limit_utilization(strain * (1.0 - 1e-6), None);
    let (just_above, _) = bar
        .material()
        .csfm_rebar_limit_utilization(strain * (1.0 + 1e-6), None);
    assert!(just_below < 1.0, "{just_below}");
    assert!(just_above > 1.0, "{just_above}");
}

/// **The bare-steel case** (`tension_stiffening: false`): the same tie
/// geometry, pulled in tension until `|eps| = eps_u` — the bilinear bare
/// branch rather than TCM/POM. Same closed form as the tie (concrete
/// contributes `0` past cracking): landing at `lambda = ft = 540`.
#[test]
fn bare_steel_case_lands_at_eps_u() {
    let rebar = Material::csfm_rebar(
        200_000.0, 500.0, 540.0, 0.05, 2.9, 32000.0, 16.0, 0.02, 0.67, false,
    )
    .unwrap();
    let eps_u = rebar_fields(&rebar).eps_u;
    let concrete =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let (domain, _concrete_id, bar_id) = one_dof_tie_or_strut(concrete, rebar, 1.0);

    let (analysis, lambda, landed_id, kind) = land_on_material_limit(
        domain,
        ArcScales::Explicit {
            displacement: eps_u * LEG_LENGTH,
            rotation: None,
            load: 540.0,
        },
        (0.1, 0.01, 0.5),
        200,
    );
    assert_eq!(kind, LimitKind::RebarRupture);
    assert_eq!(landed_id, bar_id, "the bar must be the controlling element");
    assert!((lambda - 540.0).abs() < 1.0, "lambda {lambda}");

    let strain = analysis.domain().element_gauss_responses(bar_id).unwrap()[0].strain[0];
    assert!((strain - eps_u).abs() < 1e-6 * eps_u, "strain {strain}");
    let carapace_core::model::Element::EmbeddedBar(bar) = analysis.domain().element(bar_id) else {
        panic!()
    };
    assert!(
        !bar.material().csfm_rebar_ruptured(strain),
        "landed utilization 1 must not read as ruptured under TS's strict > threshold"
    );
    assert!(bar.material().csfm_rebar_ruptured(strain * (1.0 + 1e-6)));
    let (just_below, _) = bar
        .material()
        .csfm_rebar_limit_utilization(strain * (1.0 - 1e-6), None);
    let (just_above, _) = bar
        .material()
        .csfm_rebar_limit_utilization(strain * (1.0 + 1e-6), None);
    assert!(just_below < 1.0, "{just_below}");
    assert!(just_above > 1.0, "{just_above}");
}

/// **The over-reinforced strut** (concrete crushing): the same
/// one-free-DOF geometry, pushed into compression. A generously sized bar
/// (`50 mm^2`, well inside the bare-steel elastic/hardening range at the
/// strain that crushes the concrete, so it never controls) makes this
/// "over-reinforced" in the ordinary RC design sense: the section fails by
/// the concrete crushing, not by the steel.
///
/// Unlike the tension cases, the concrete's own stiffness contributes
/// materially here (it isn't driven to exactly zero resistance the way
/// zero-tension concrete is past cracking), so there's no exact one-line
/// closed form to pick `ArcScales::Explicit` from the way the two tests
/// above do. `displacement`/`load` below are just order-of-magnitude
/// starting scales (a rough `ux1` at `eps_cu2` and a guessed equilibrium
/// force) — `ArcLength` does the rest; the actual landing turns out at
/// `lambda ~ -26500`, over 10x the rough guess, which is fine since these
/// are only *scales*, not predictions.
#[test]
fn over_reinforced_strut_lands_on_concrete_crushing() {
    const EPS_CU2: f64 = -0.0035;
    let rebar = Material::csfm_rebar(
        200_000.0, 500.0, 540.0, 0.05, 2.9, 32000.0, 16.0, 0.02, 0.67, true,
    )
    .unwrap();
    let concrete =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    let (domain, concrete_id, bar_id) = one_dof_tie_or_strut(concrete, rebar, 50.0);

    let mut arc = ArcLength::adaptive(
        0.1,
        0.01,
        1.0,
        ArcScales::Explicit {
            displacement: 0.35,
            rotation: None,
            load: 1750.0,
        },
    );
    arc.direction = ArcDirection::Decreasing;
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
        .try_build(domain)
        .unwrap_or_else(|e| panic!("invalid model: {e:?}"));

    let (landed_id, kind, lambda) = 'outer: {
        for _ in 0..200 {
            let result = analysis
                .step()
                .unwrap_or_else(|e| panic!("arc-length step failed before landing: {e:?}"));
            let arc_info = result.arc.expect("ArcLength reports continuation info");
            if let Some(stop) = arc_info.stop {
                assert_eq!(stop.reason, StopReason::MaterialLimit);
                assert!(stop.landed_exactly);
                let (element, utilization, kind) = analysis
                    .domain()
                    .controlling_limit_utilization(&LimitStateConfig::default())
                    .unwrap();
                assert!((utilization - 1.0).abs() < 1e-6, "{utilization}");
                break 'outer (element, kind, result.load_factor);
            }
        }
        panic!("did not land on the material limit within the step budget");
    };
    assert_eq!(kind, LimitKind::ConcreteCrushing);
    assert_eq!(
        landed_id, concrete_id,
        "the concrete must be the controlling element"
    );
    assert!(
        lambda < 0.0,
        "a compressive load must be negative: {lambda}"
    );

    let strain = analysis
        .domain()
        .element_gauss_responses(concrete_id)
        .unwrap()[0]
        .strain;
    let eps2 = strain[1].min(strain[0]);
    assert!((eps2 - EPS_CU2).abs() < 1e-6 * EPS_CU2.abs(), "eps2 {eps2}");

    let (just_below, _) =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap()
            .csfm_concrete_limit_utilization(
                &carapace_core::model::PlaneVector::new(0.0, EPS_CU2 * 0.999, 0.0),
                &LimitStateConfig::default(),
            );
    let (just_above, _) =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap()
            .csfm_concrete_limit_utilization(
                &carapace_core::model::PlaneVector::new(0.0, EPS_CU2 * 1.001, 0.0),
                &LimitStateConfig::default(),
            );
    assert!(just_below < 1.0, "{just_below}");
    assert!(just_above > 1.0, "{just_above}");

    // Over-reinforced: the bar stays far from its own rupture utilization
    // at the strain that crushed the concrete.
    let bar_strain = analysis.domain().element_gauss_responses(bar_id).unwrap()[0].strain[0];
    let carapace_core::model::Element::EmbeddedBar(bar) = analysis.domain().element(bar_id) else {
        panic!()
    };
    let (bar_utilization, bar_kind) = bar
        .material()
        .csfm_rebar_limit_utilization(bar_strain, None);
    assert_eq!(bar_kind, LimitKind::RebarRupture);
    assert!(
        bar_utilization < 1.0,
        "the bar must not also be at its limit: {bar_utilization}"
    );
}
