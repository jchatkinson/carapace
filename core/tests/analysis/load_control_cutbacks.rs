//! `Integrator::LoadControlCutbacks` (`docs/csfm-plan.md` §2.5 point 1): adaptive substepping on
//! `Analysis::step`/`step_prescribed`, ported from CalcsApp's `solveIncrementWithSession`/
//! `solveLoadPath` (`../calcs.app/lib/fem/nonlinear.ts`).
//!
//! Fixture: a `Truss` (k_t = 50) and a `ZeroLength` + `ElasticPP` spring (k_epp = 100 while
//! `|strain| < eyp = 0.01`) in parallel between two nodes, loaded by a reference force of 5.0 — the
//! same system `newton.rs` uses to show Newton resolves a regime crossing. Here
//! `ConvergenceTest::NormUnbalance { max_iter: 1 }` makes that resolution itself the failure mode: a
//! single Newton iteration is *exact* within one material regime (the tangent used to form it never
//! goes stale), so any load factor interval that stays within one regime converges in exactly one
//! iteration, while an interval that crosses the yield point needs a second iteration (to re-form
//! the tangent in the new regime) and so fails with `max_iter: 1`.
//!
//! The two springs share one displacement `u` (same two nodes, same DOF), so while elastic the
//! total external force is `F = u * (k_t + k_epp)`, and the spring alone reaches `fy = u * k_epp`
//! at `u = eyp`. The yield point in *load factor* terms is therefore
//! `lambda_yield = eyp * (k_t + k_epp) / FORCE`, not `fy / FORCE` — the parallel truss carries part
//! of the load too. This gives a deterministic, closed-form-checkable boundary for cutbacks to find
//! without any near-singular or genuinely divergent physics.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, AnalysisError, ConstraintHandler, ConvergenceTest,
    Integrator, TangentStrategy,
};
use carapace_core::model::{Domain, Element, Material, Node, NodeId, Truss, ZeroLength};

const K_TRUSS: f64 = 50.0;
const K_EPP: f64 = 100.0;
const EYP: f64 = 0.01;
const FORCE: f64 = 5.0;
// The pseudo-time (load factor) at which the spring yields: u = EYP while elastic needs
// F = EYP * (K_TRUSS + K_EPP), so lambda_yield = F / FORCE.
const YIELD_LAMBDA: f64 = (EYP * (K_TRUSS + K_EPP)) / FORCE;

fn build() -> (Domain, NodeId) {
    let length = 100.0;
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, FORCE);

    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        1.0,
        Material::Elastic {
            e: K_TRUSS * length,
        },
    )));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::elastic_pp(K_EPP, EYP)),
    ));
    (domain, node_j)
}

fn build_analysis(increment: f64, min_increment: f64) -> (Analysis, NodeId) {
    let (domain, node_j) = build();
    let analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControlCutbacks {
            increment,
            min_increment,
        })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 1,
        })
        .build(domain);
    (analysis, node_j)
}

/// A nominal increment that spans the yield point is halved until it fits inside one regime, and
/// that reduced size — not the nominal one — is what the call accepts. `YIELD_LAMBDA = 0.3`, so
/// `1.0 -> 0.5 -> 0.25` lands the first elastic-only attempt after two retries.
#[test]
fn a_nominal_increment_crossing_a_regime_boundary_is_halved_until_it_fits() {
    let (mut analysis, node_j) = build_analysis(1.0, 0.1);

    let result = analysis.step().expect("0.25 stays below the yield lambda");
    let cutbacks = result
        .cutbacks
        .expect("LoadControlCutbacks reports cutbacks");

    assert_eq!(cutbacks.nominal_increment, 1.0);
    assert_eq!(cutbacks.accepted_increment, 0.25);
    assert_eq!(cutbacks.retries, 2);
    assert!((result.load_factor - 0.25).abs() < 1e-12);

    // Elastic regime only: both springs act, so u = F / (k_t + k_epp).
    let expected_u = (FORCE * 0.25) / (K_TRUSS + K_EPP);
    let u = analysis.domain().node(node_j).displacement[0];
    assert!(
        (u - expected_u).abs() < 1e-10,
        "u = {u}, expected {expected_u}"
    );
}

/// The uncovered remainder from a cutback is requested in full (not pre-halved) by the next call —
/// CalcsApp's `pending.unshift(remainder)`. Starting from the 0.25 accepted above, the remainder is
/// 0.75: `0.75 -> 0.375 -> 0.1875` all still push the total load factor past `YIELD_LAMBDA` (0.3),
/// and the next halving (0.09375) falls below the shared 0.1 floor before it can even be tried, so
/// cutbacks exhaust — and the error reports the load factor actually committed by the first call
/// (0.25), not zero or a failed attempt's target.
#[test]
fn exhausted_cutbacks_report_the_last_committed_load_factor_not_the_failed_attempt() {
    let (mut analysis, _node_j) = build_analysis(1.0, 0.1);
    let first = analysis.step().expect("first call accepts 0.25");
    assert!((first.load_factor - 0.25).abs() < 1e-12);

    let err = analysis
        .step()
        .expect_err("every attempt down to the floor still crosses the yield lambda");
    match err {
        AnalysisError::LoadControlExhausted {
            accepted_load_factor,
            ..
        } => {
            assert!(
                (accepted_load_factor - 0.25).abs() < 1e-12,
                "accepted_load_factor = {accepted_load_factor}"
            );
        }
        other => panic!("expected LoadControlExhausted, got {other:?}"),
    }
}

/// A plain (non-cutback) `LoadControl` on the same model and the same tight `max_iter: 1` fails
/// outright on the first attempt — confirms the fixture's regime-crossing failure mode is real and
/// not an artifact of the cutback plumbing.
#[test]
fn without_cutbacks_the_same_regime_crossing_increment_fails_to_converge() {
    let (domain, _node_j) = build();
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 1,
        })
        .build(domain);

    let err = analysis
        .step()
        .expect_err("one iteration cannot resolve the regime crossing");
    assert!(matches!(err, AnalysisError::FailedToConverge { .. }));
}

// `YIELD_LAMBDA` sanity, checked at compile time: it sits strictly between `1.0 / 4` and `1.0 / 2`,
// which is what makes the two-halving trajectory in the tests above deterministic.
const _: () = assert!(YIELD_LAMBDA > 0.25 && YIELD_LAMBDA < 0.5);
