use carapace_core::analysis::{Algorithm, AnalysisBuilder, AnalysisError, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, Material, Node, Truss, ZeroLength};

/// Pre-M7 foundational work: `Material`'s trial/commit split (see its doc
/// comment) and `Domain`'s snapshot-and-restore-on-failure, landed ahead of
/// M7's fiber sections and full material catalog since those need real
/// path-dependent state to be correct.
///
/// Same `Truss` (k_t=50, elastic) + `ZeroLength`+`ElasticPP` (k_epp=100,
/// fy=1.0) parallel system as `m4_analysis.rs`, but driven through *two*
/// separate `Analysis::step()` calls — loading past yield, then partially
/// unloading — rather than one. `Algorithm::NewtonRaphson`'s Newton loop
/// only ever evaluates materials via `Material::trial_stress_tangent`
/// (never mutates), and `Domain::commit` runs once after each step
/// converges: step 1 leaves the EPP spring's plastic strain committed at
/// `ep = u1 - fy/e_epp`, and step 2's unload must react to *that*
/// committed history, not recompute from scratch.
///
/// Closed form: after step 1 (load 0->5, past yield), `u1 = (5-1)/50 =
/// 0.08`, `ep = 0.08 - 1.0/100 = 0.07`. Step 2 unloads to load_factor=3
/// (increment -2.0); solving `F = k_t*u + epp_force(u)` with the EPP now
/// *elastic* about `ep=0.07` (since the trial stress `e_epp*(u-ep)` stays
/// under `fy` for the unload direction): `F = 50*u + 100*(u-0.07) =
/// 150*u - 7 = 3` => `u2 = 10/150 = 1/15`. A model *without* real
/// plasticity (recomputing from strain alone, ignoring `ep`) would instead
/// stay on the plastic plateau and predict `u2 = (3-1)/50 = 0.04` — a
/// materially different, wrong answer.
#[test]
fn multi_step_analysis_shows_real_permanent_set_on_partial_unload() {
    let (k_t, length) = (50.0, 100.0);
    let (e_area, e_modulus) = (1.0, k_t * length);
    let (e_epp, eyp, fy) = (100.0, 0.01, 1.0);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, 1.0);

    domain.add_element(Element::Truss(Truss::new(node_i, node_j, e_area, Material::Elastic { e: e_modulus })));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::elastic_pp(e_epp, eyp)),
    ));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 5.0 })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 20,
        })
        .build(domain);

    // Step 1: load 0 -> 5, past yield.
    analysis.step().expect("step 1 (past yield) should converge");
    let u1 = analysis.domain().node(node_j).displacement[0];
    assert!((u1 - 0.08).abs() < 1e-9, "expected u1=0.08, got {u1}");

    // Step 2: unload to a total force of 3.0. `load_factor` lives on
    // `Analysis`, not `Domain`, so this needs a second, fresh `Analysis`
    // (its own `load_factor` starts at 0) built from a *clone* of the
    // first one's domain — the point is that the clone carries step 1's
    // *committed* material state (`Domain: Clone` is exactly the same
    // derive the snapshot/restore mechanism uses), so step 2 must react to
    // that history, not recompute from scratch. Since the reference load
    // pattern (`node_j.load = 1.0`) is unchanged, `load_factor` directly
    // equals the total applied force — `increment: 3.0` targets F=3.0
    // outright, not a delta from step 1's F=5.0.
    let mut analysis2 = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 3.0 })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 20,
        })
        .build(analysis.domain().clone());

    analysis2.step().expect("step 2 (partial unload) should converge");
    let u2 = analysis2.domain().node(node_j).displacement[0];
    let expected = 1.0 / 15.0;
    assert!((u2 - expected).abs() < 1e-9, "expected u2={expected} (real plasticity), got {u2}");

    // The wrong (path-independent) answer a stateless model would give,
    // named here only to make the contrast explicit — not asserted against.
    let wrong_stateless_answer = (3.0 - fy) / k_t;
    assert!((u2 - wrong_stateless_answer).abs() > 1e-6);
}

/// `Domain`'s snapshot-and-restore: a step that fails to converge must
/// leave the domain exactly as it was before the attempt, not partway
/// through a discarded Newton iteration. A `ConvergenceTest` with
/// `max_iter: 0` can never converge (the loop body never runs), forcing
/// `FailedToConverge` deterministically without needing a contrived
/// numerically-hard system.
#[test]
fn failed_step_leaves_domain_unchanged() {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([100.0, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, 50.0);
    domain.add_element(Element::Truss(Truss::new(node_i, node_j, 2.0, Material::Elastic { e: 30000.0 })));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 0 })
        .build(domain);

    let before = analysis.domain().node(node_j).displacement[0];
    let result = analysis.step();

    assert!(matches!(result, Err(AnalysisError::FailedToConverge { step: 1 })));
    let after = analysis.domain().node(node_j).displacement[0];
    assert_eq!(before, after, "a failed step must not leave displacement partially applied");
    assert_eq!(before, 0.0);
}
