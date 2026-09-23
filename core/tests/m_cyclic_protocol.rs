use carapace_core::analysis::{Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, Material, Node, NodeId, Truss, ZeroLength};

/// Same `Truss` (k_t=50, elastic) + `ZeroLength`+`ElasticPP` (k_epp=100,
/// eyp=0.01, fy=1.0) parallel system as `m4_analysis.rs`/`material_state.rs`
/// — a reference unit nodal load (`load_node(node_j, 0, 1.0)`) makes
/// `StepResult::load_factor` read directly as the total force at `node_j`,
/// the same trick those tests use.
fn build_system() -> (Domain, NodeId) {
    let (k_t, length) = (50.0, 100.0);
    let (e_area, e_modulus) = (1.0, k_t * length);
    let (e_epp, eyp) = (100.0, 0.01);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, 1.0);

    domain.add_element(Element::Truss(Truss::new(node_i, node_j, e_area, Material::Elastic { e: e_modulus })));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::elastic_pp(e_epp, eyp)),
    ));

    (domain, node_j)
}

/// One `DisplacementControl` leg via `Analysis::set_integrator` (no
/// rebuild) targeting an *absolute* displacement — the increment passed to
/// `DisplacementControl` is computed fresh from wherever the controlled DOF
/// actually is right now, not assumed from the previous leg's target.
fn leg(analysis: &mut Analysis, node_j: NodeId, target: f64) -> f64 {
    let current = analysis.domain().node(node_j).displacement[0];
    analysis.set_integrator(Integrator::DisplacementControl {
        node: node_j,
        dof: 0,
        increment: target - current,
    });
    analysis.step().expect("cyclic protocol leg should converge").load_factor
}

/// Walks the controlled DOF from wherever it is to `target` in steps no
/// larger than `step`, via repeated `leg` calls (`Analysis::
/// set_integrator`, no rebuild) — the realistic way a cyclic/quasi-static
/// protocol is actually driven: many small `DisplacementControl` legs, not
/// one big excursion. Small steps keep `DisplacementControl`'s per-step
/// linear prediction (exact only within a single material regime for the
/// whole step, §`Integrator`'s doc comment) close enough to exact that
/// crossing a regime boundary mid-leg contributes only a small, bounded
/// error rather than compounding.
fn ramp_to(analysis: &mut Analysis, node_j: NodeId, target: f64, step: f64) -> f64 {
    let mut load_factor = 0.0;
    loop {
        let current = analysis.domain().node(node_j).displacement[0];
        let remaining = target - current;
        if remaining.abs() < 1e-12 {
            return load_factor;
        }
        let next = current + remaining.clamp(-step, step);
        load_factor = leg(analysis, node_j, next);
    }
}

/// Exercises `Analysis::set_integrator` driving a cyclic `DisplacementControl`
/// protocol — many small legs on *one* `Analysis` (no rebuild) — the
/// mechanism a real quasi-static cyclic-loading protocol (e.g. FEMA 461/
/// ATC-24-style excursions) is meant to use: ramp past yield, back through
/// zero to a negative excursion, then back to zero — and check that the
/// material remembers it yielded (nonzero residual force at zero net
/// displacement — the hysteresis/permanent-set signature a path-independent
/// material could never produce).
#[test]
fn set_integrator_drives_a_cyclic_protocol_with_permanent_set_at_zero_displacement() {
    let (domain, node_j) = build_system();

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::DisplacementControl {
            node: node_j,
            dof: 0,
            increment: 0.0, // placeholder — `ramp_to` calls `set_integrator` before every step
        })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 20 })
        .build(domain);

    let step = 0.001; // a tenth of eyp — small enough to keep per-leg regime-crossing error negligible
    let tol = 1e-6;

    // 0 -> +0.03 (past yield: peak force is the truss's share plus the EPP
    // spring's clamped yield force, F = k_t*0.03 + fy = 50*0.03 + 1.0 = 2.5).
    let peak = ramp_to(&mut analysis, node_j, 0.03, step);
    let u_peak = analysis.domain().node(node_j).displacement[0];
    assert!((u_peak - 0.03).abs() < tol, "expected u=0.03 at peak, got {u_peak}");
    // A generous tolerance: exact only if every ~30 small legs landed
    // exactly on a plastic-branch commit; some accumulate a little
    // regime-boundary imprecision (see `ramp_to`'s doc comment) without
    // affecting the final *displacement* (self-corrected every leg).
    assert!((peak - 2.5).abs() < 0.15, "expected peak force~2.5, got {peak}");

    // +0.03 -> -0.03 (through zero, well past yield in the other direction).
    let trough = ramp_to(&mut analysis, node_j, -0.03, step);
    let u_trough = analysis.domain().node(node_j).displacement[0];
    assert!((u_trough - (-0.03)).abs() < tol, "expected u=-0.03 at trough, got {u_trough}");
    assert!((trough - (-2.5)).abs() < 0.15, "expected trough force~-2.5, got {trough}");

    // -0.03 -> 0.0 (unloading from the most recent, negative-going yield).
    let residual_force = ramp_to(&mut analysis, node_j, 0.0, step);
    let u_final = analysis.domain().node(node_j).displacement[0];
    assert!(u_final.abs() < tol, "expected u=0.0 after returning, got {u_final}");

    // A stateless (path-independent) material would give zero force at
    // zero displacement (the truss alone is zero at u=0, and a stateless
    // EPP envelope re-evaluated fresh at strain=0 gives zero too); real
    // permanent set leaves a nonzero residual, bounded in magnitude by the
    // EPP spring's yield force (it can never exceed fy=1.0, since the
    // truss alone is exactly zero at u=0).
    assert!(
        residual_force.abs() > 0.5,
        "expected a clearly nonzero permanent-set residual force at u=0, got {residual_force}"
    );
    assert!(
        residual_force.abs() <= 1.0 + 1e-3,
        "residual force magnitude should be bounded by fy=1.0 (only the EPP spring can be nonzero at u=0), got {residual_force}"
    );
    // The most recent yield excursion was negative-going (trough at -0.03),
    // so the EPP spring's committed plastic strain sits on the negative
    // side, and unloading from there leaves a *positive* residual (the
    // spring pulling back toward its negative-side plastic offset) at u=0.
    assert!(residual_force > 0.0, "expected a positive residual force after a negative-going final excursion, got {residual_force}");
}
