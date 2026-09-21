use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, Material, Node, Truss, ZeroLength};

/// Builds a `Truss` (elastic, k_t = E*A/L = 50) and a `ZeroLength`+`ElasticPP`
/// (k_epp = 100 while `|strain| < eyp = 0.01`, i.e. fy = 1.0, then zero
/// stiffness / constant force) in parallel between the same two nodes, both
/// acting on node_j's x DOF. This is a genuinely path-dependent-*shaped*
/// nonlinear system within a single step (the EPP spring changes regime
/// partway through the applied load): `Algorithm::Linear`'s one-shot solve
/// (exact only within a single regime, per M1-M3) cannot resolve it, but
/// `Algorithm::NewtonRaphson` — new at M4 — can, by re-forming the tangent
/// each iteration as displacement crosses the yield point.
///
/// Closed form for total applied force F > fy: once yielded, the EPP
/// spring's force is pinned at fy and all further load is carried by the
/// truss alone, so `u = (F - fy) / k_t`.
fn build_elastic_plastic_parallel_system(force: f64) -> (Domain, carapace_core::model::NodeId, carapace_core::model::NodeId) {
    let (k_t, length) = (50.0, 100.0);
    let (e_area, e_modulus) = (1.0, k_t * length); // E*A/L = k_t
    let (e_epp, eyp) = (100.0, 0.01);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2).with_load(0, force));

    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        e_area,
        Material::Elastic { e: e_modulus },
    )));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::elastic_pp(e_epp, eyp)),
    ));

    (domain, node_i, node_j)
}

#[test]
fn newton_raphson_resolves_elastic_perfectly_plastic_regime_crossing() {
    let force = 5.0;
    let (domain, _node_i, node_j) = build_elastic_plastic_parallel_system(force);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 20,
        })
        .build(domain);

    analysis.step().expect("should converge across the EPP yield point");

    let (k_t, fy) = (50.0, 1.0);
    let expected = (force - fy) / k_t;
    let u = analysis.domain().node(node_j).displacement[0];
    assert!((u - expected).abs() < 1e-9, "expected u={expected}, got {u}");
}

#[test]
fn norm_disp_incr_and_energy_incr_also_converge_to_the_same_result() {
    let (k_t, fy, force) = (50.0, 1.0, 5.0);
    let expected = (force - fy) / k_t;

    for test in [
        ConvergenceTest::NormDispIncr { tol: 1e-9, max_iter: 20 },
        ConvergenceTest::EnergyIncr { tol: 1e-9, max_iter: 20 },
    ] {
        let (domain, _node_i, node_j) = build_elastic_plastic_parallel_system(force);

        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::NewtonRaphson)
            .test(test)
            .build(domain);

        analysis.step().expect("should converge under every ConvergenceTest variant");
        let u = analysis.domain().node(node_j).displacement[0];
        assert!((u - expected).abs() < 1e-6, "expected u={expected}, got {u} (test={test:?})");
    }
}

/// `Integrator::DisplacementControl` on the M1 truss case: prescribing the
/// same target displacement (rather than the equivalent load) should drive
/// the analysis to the same load factor M1 verified directly (P=50).
#[test]
fn displacement_control_recovers_the_load_that_produces_the_target_displacement() {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    // Reference load pattern is a unit load — DisplacementControl scales it
    // by whatever load factor produces the target displacement.
    let node_j = domain.add_node(Node::new([100.0, 0.0]).fix(1).fix(2).with_load(0, 1.0));

    domain.add_element(Element::Truss(Truss::new(node_i, node_j, 2.0, Material::Elastic { e: 30000.0 })));

    let target = 0.0833333333;
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::DisplacementControl {
            node: node_j,
            dof: 0,
            increment: target,
        })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    let result = analysis.step().expect("displacement-controlled truss should solve");

    let dx = analysis.domain().node(node_j).displacement[0];
    assert!((dx - target).abs() < 1e-9, "expected dx={target}, got {dx}");
    assert!(
        (result.load_factor - 50.0).abs() < 1e-6,
        "expected the controlling load factor to recover P=50, got {}",
        result.load_factor
    );
}
