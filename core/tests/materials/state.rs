use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, AnalysisError, ConstraintHandler, ConvergenceTest, Integrator,
    TangentStrategy,
};
use carapace_core::model::{Domain, Element, Material, Node, Truss, ZeroLength};

/// `Material`'s trial/commit split (see its doc
/// comment) and `Domain`'s snapshot-and-restore-on-failure, needed by fiber sections and the full material catalog, which depend on real path-dependent state.
///
/// Same `Truss` (k_t=50, elastic) + `ZeroLength`+`ElasticPP` (k_epp=100,
/// fy=1.0) parallel system as `analysis/newton.rs`, but driven through *two*
/// separate `Analysis::step()` calls — loading past yield, then partially
/// unloading — rather than one. `Algorithm::Newton`'s Newton loop
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

    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        e_area,
        Material::Elastic { e: e_modulus },
    )));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::elastic_pp(e_epp, eyp)),
    ));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 5.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 20,
        })
        .build(domain);

    // Step 1: load 0 -> 5, past yield.
    analysis
        .step()
        .expect("step 1 (past yield) should converge");
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
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 20,
        })
        .build(analysis.domain().clone());

    analysis2
        .step()
        .expect("step 2 (partial unload) should converge");
    let u2 = analysis2.domain().node(node_j).displacement[0];
    let expected = 1.0 / 15.0;
    assert!(
        (u2 - expected).abs() < 1e-9,
        "expected u2={expected} (real plasticity), got {u2}"
    );

    // The wrong (path-independent) answer a stateless model would give,
    // named here only to make the contrast explicit — not asserted against.
    let wrong_stateless_answer = (3.0 - fy) / k_t;
    assert!((u2 - wrong_stateless_answer).abs() > 1e-6);
}

/// `Domain`'s snapshot-and-restore: a step that fails to converge must
/// leave the domain exactly as it was before the attempt, not partway
/// through a discarded Newton iteration. `NormDispIncr` with a single
/// permitted iteration can never converge here (the one correction is the
/// full elastic displacement, far above `tol`), forcing `FailedToConverge`
/// deterministically — after a correction *was* applied to the domain —
/// without needing a contrived numerically-hard system. (A zero iteration
/// budget is rejected outright as `AnalysisError::InvalidOption`.)
#[test]
fn failed_step_leaves_domain_unchanged() {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([100.0, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, 50.0);
    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        2.0,
        Material::Elastic { e: 30000.0 },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-9,
            max_iter: 1,
        })
        .build(domain);

    let before = analysis.domain().node(node_j).displacement[0];
    let result = analysis.step();

    assert!(matches!(
        result,
        Err(AnalysisError::FailedToConverge { step: 1 })
    ));
    let after = analysis.domain().node(node_j).displacement[0];
    assert_eq!(
        before, after,
        "a failed step must not leave displacement partially applied"
    );
    assert_eq!(before, 0.0);
}

/// A step that fails to converge on a model that has already yielded leaves no trace: not in the nodal
/// displacements, and not in any material's committed history.
///
/// `Analysis` restores only nodal displacements when a step fails (materials are never mutated during
/// iteration, see `Material`'s trial/commit split), so this checks that claim on a yielding model. Same
/// `Truss` + `ElasticPP` parallel system as above: step 1 loads past yield and commits `ep = 0.07`. The
/// unload step to F = 3 is then attempted with a single permitted iteration under `NormDispIncr`, which cannot
/// converge (`FailedToConverge`). A domain taken from the failed analysis and a clone taken before the failed attempt
/// must then give bitwise-identical results for the same unload step with a normal iteration budget, and both
/// must land on the closed-form answer `u = 1/15` that depends on the committed plastic strain.
#[test]
fn failed_step_after_yielding_leaves_no_material_trace() {
    let length = 100.0;
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, 1.0);
    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        1.0,
        Material::Elastic { e: 50.0 * length },
    )));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::elastic_pp(100.0, 0.01)),
    ));
    let analysis_with = |domain: Domain, increment: f64, test: ConvergenceTest| {
        AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment })
            .algorithm(Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            })
            .test(test)
            .build(domain)
    };
    let converging = ConvergenceTest::NormUnbalance {
        tol: 1e-9,
        max_iter: 20,
    };
    // The single correction is the full displacement change, far above `tol`, so one iteration cannot pass.
    let never_in_one_iteration = ConvergenceTest::NormDispIncr {
        tol: 1e-9,
        max_iter: 1,
    };

    let mut yielded = analysis_with(domain, 5.0, converging);
    yielded.step().expect("loading past yield converges");
    let before_attempt = yielded.domain().clone();

    let mut failing = analysis_with(before_attempt.clone(), 3.0, never_in_one_iteration);
    assert!(matches!(
        failing.step(),
        Err(AnalysisError::FailedToConverge { step: 1 })
    ));
    let after_failure = failing.into_domain();
    assert_eq!(
        before_attempt.node(node_j).displacement,
        after_failure.node(node_j).displacement,
        "the failed attempt must not leave displacement partially applied"
    );

    let mut control = analysis_with(before_attempt, 3.0, converging);
    let mut retried = analysis_with(after_failure, 3.0, converging);
    control.step().expect("unload step converges");
    retried
        .step()
        .expect("unload step converges after a failed attempt");
    let (u_control, u_retried) = (
        control.domain().node(node_j).displacement[0],
        retried.domain().node(node_j).displacement[0],
    );
    assert_eq!(
        u_control, u_retried,
        "a failed attempt must not change the result of the next step"
    );
    assert!(
        (u_retried - 1.0 / 15.0).abs() < 1e-9,
        "expected u = 1/15, got {u_retried}"
    );
}
