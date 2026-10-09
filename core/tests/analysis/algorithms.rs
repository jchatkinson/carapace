use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, LineSearch,
    TangentStrategy,
};
use carapace_core::model::{Domain, Element, Material, Node, Truss, ZeroLength};

/// Same `Truss` (k_t=50, elastic) + `ZeroLength`+`ElasticPP` (k_epp=100,
/// fy=1.0) parallel system `materials/state.rs`/`analysis/newton.rs` already use
/// — a small, well-understood nonlinear (past-yield) case, reused here to
/// isolate `docs/algorithms.md` Part I's *iteration-strategy* behavior from
/// element/material correctness (already covered elsewhere).
fn epp_truss_domain() -> (Domain, carapace_core::model::NodeId) {
    let (k_t, length) = (50.0, 100.0);
    let (e_area, e_modulus) = (1.0, k_t * length);
    let (e_epp, eyp) = (100.0, 0.01);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, 1.0);

    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        e_area,
        Material::Elastic { e: e_modulus },
    )));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::elastic_pp(e_epp, eyp)),
    ));
    (domain, node_j)
}

const TOL: f64 = 1e-9;
// Closed form (same derivation as `materials/state.rs`): load 0->5 past
// yield (fy=1.0) gives `u1 = (5-1)/50 = 0.08`.
const EXPECTED_U: f64 = 0.08;

/// `TangentStrategy::Current` re-forms and factors the tangent every iteration. On a step that
/// crosses yield (elastic-perfectly-plastic spring parallel to a truss) the answer matches the
/// closed-form displacement, the factorization count equals the iteration count, and more than one
/// iteration is needed.
#[test]
fn tangent_strategy_current_factors_once_per_iteration() {
    let (domain, node_j) = epp_truss_domain();
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 5.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: TOL,
            max_iter: 20,
        })
        .build(domain);

    let result = analysis.step().expect("past-yield step should converge");
    assert!((analysis.domain().node(node_j).displacement[0] - EXPECTED_U).abs() < TOL);
    assert_eq!(
        result.factorizations, result.iterations,
        "TangentStrategy::Current must factor exactly once per iteration"
    );
    assert!(
        result.iterations > 1,
        "a regime change within the step should need more than one iteration"
    );
}

/// `TangentStrategy::Initial`: forms the tangent once, from the model's
/// very first (undeformed, elastic) state, and never re-factors again for
/// the rest of the `Analysis`'s lifetime — proven here across *two*
/// separate `step()` calls, not just within one, since `docs/algorithms.md`
/// §3.1 specifically calls out "once ever," not "once per step"
/// (`TangentStrategy::ReuseAtStepStart` is the once-per-step variant).
/// Converges to the *same* answer as `Current` (linear convergence rate,
/// more iterations, but the same equilibrium point), while never paying
/// for a second factorization.
#[test]
fn tangent_strategy_initial_never_refactors_across_multiple_steps() {
    let (domain, node_j) = epp_truss_domain();
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 2.5 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Initial,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: TOL,
            max_iter: 100,
        })
        .build(domain);

    let result1 = analysis.step().expect("step 1 should converge");
    assert_eq!(
        result1.factorizations, 1,
        "the very first iteration must factor once"
    );

    let result2 = analysis.step().expect("step 2 should converge");
    assert_eq!(
        result2.factorizations, 0,
        "a later step must reuse the very first factorization, forming none of its own"
    );

    assert!(
        (analysis.domain().node(node_j).displacement[0] - EXPECTED_U).abs() < TOL,
        "Initial-tangent iteration must still converge to the same equilibrium as Current"
    );
}

/// `TangentStrategy::ReuseAtStepStart`: factors once per `step()` call
/// (unlike `Initial`'s once-ever), so a *second* step gets its own fresh
/// factorization at the state that step started from — distinguishing this
/// from `tangent_strategy_initial_never_refactors_across_multiple_steps`
/// above is the whole point of these being two different enum variants.
#[test]
fn tangent_strategy_reuse_at_step_start_refactors_once_per_step() {
    let (domain, node_j) = epp_truss_domain();
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 2.5 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::ReuseAtStepStart,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: TOL,
            max_iter: 100,
        })
        .build(domain);

    let result1 = analysis.step().expect("step 1 should converge");
    assert_eq!(result1.factorizations, 1);
    let result2 = analysis.step().expect("step 2 should converge");
    assert_eq!(
        result2.factorizations, 1,
        "each step gets its own fresh step-start factorization"
    );

    assert!((analysis.domain().node(node_j).displacement[0] - EXPECTED_U).abs() < TOL);
}

/// `LineSearch` must not change the converged answer versus plain Newton —
/// only how it gets there. Both variants checked against the same
/// `Current`-tangent, no-line-search baseline.
#[test]
fn line_search_converges_to_the_same_answer_as_plain_newton() {
    for line_search in [
        LineSearch::Bisection {
            tol: 1e-10,
            max_iter: 30,
            max_eta: 16.0,
        },
        LineSearch::RegulaFalsi {
            tol: 1e-10,
            max_iter: 30,
            max_eta: 16.0,
        },
    ] {
        let (domain, node_j) = epp_truss_domain();
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 5.0 })
            .algorithm(Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: Some(line_search),
            })
            .test(ConvergenceTest::NormUnbalance {
                tol: TOL,
                max_iter: 20,
            })
            .build(domain);

        analysis
            .step()
            .expect("line-search-assisted step should converge");
        let u = analysis.domain().node(node_j).displacement[0];
        assert!(
            (u - EXPECTED_U).abs() < TOL,
            "{line_search:?} gave u={u}, expected {EXPECTED_U}"
        );
    }
}

/// `Algorithm::KrylovNewton` (stale `Initial` tangent + acceleration)
/// converges to the same equilibrium as full `Current`-tangent Newton,
/// while factoring far less often than once per iteration — the whole
/// point of accelerating a stale tangent instead of re-forming it.
#[test]
fn krylov_newton_converges_with_far_fewer_factorizations_than_full_newton() {
    let (domain, node_j) = epp_truss_domain();
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 5.0 })
        .algorithm(Algorithm::KrylovNewton {
            tangent: TangentStrategy::Initial,
            max_dimension: 3,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: TOL,
            max_iter: 100,
        })
        .build(domain);

    let result = analysis
        .step()
        .expect("Krylov-accelerated step should converge");
    let u = analysis.domain().node(node_j).displacement[0];
    assert!(
        (u - EXPECTED_U).abs() < TOL,
        "got u={u}, expected {EXPECTED_U}"
    );
    assert!(
        result.factorizations < result.iterations,
        "acceleration should let most iterations skip factoring entirely (factorizations={}, iterations={})",
        result.factorizations,
        result.iterations
    );
}
