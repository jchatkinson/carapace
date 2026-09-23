use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, Material, Node, Truss};

/// A ground-anchored truss oriented along x (like `m1_truss.rs`'s case), so
/// that a load/displacement on the free node's x DOF actually exercises the
/// truss's axial stiffness `k = E*A/L`. `y` offsets each independent
/// ground/free node pair so multiple calls don't collide in space (though
/// only x-topology, via `NodeId`, actually matters to the solve).
fn truss_to_ground(domain: &mut Domain, y: f64, k: f64, load: f64) -> carapace_core::model::NodeId {
    let (area, length) = (1.0, 100.0);
    let ground = domain.add_node(Node::new([0.0, y]).fix(0).fix(1).fix(2));
    let free = domain.add_node(Node::new([length, y]).fix(1).fix(2));
    domain.load_node(free, 0, load);
    domain.add_element(Element::Truss(Truss::new(ground, free, area, Material::Elastic { e: k * length / area })));
    free
}

/// `equal_dof` tying two otherwise-unconnected nodes' x DOF should make two
/// separate ground-trusses act exactly like two springs in parallel on a
/// single node (same closed form `m4_analysis.rs`'s parallel-EPP/Truss case
/// uses for a single shared node, but here the "sharing" is via constraint
/// aliasing across two distinct `NodeId`s instead of two elements on one
/// node) — the real thing under test is that `Domain::number_dofs` gives
/// the constrained node's x DOF the *same* equation number as the retained
/// node's, so both trusses' stiffness contributions land in the same row/
/// column of the assembled system.
#[test]
fn equal_dof_couples_two_nodes_x_translation_like_parallel_springs() {
    let (k1, k2, force) = (50.0, 30.0, 20.0);
    let mut domain = Domain::new();

    let retained = truss_to_ground(&mut domain, 0.0, k1, force);
    let constrained = truss_to_ground(&mut domain, 200.0, k2, 0.0);
    domain.equal_dof(retained, constrained, &[0]);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);

    analysis.step().expect("linear elastic parallel-spring system should solve");

    let expected = force / (k1 + k2);
    let u_retained = analysis.domain().node(retained).displacement[0];
    let u_constrained = analysis.domain().node(constrained).displacement[0];
    assert!((u_retained - expected).abs() < 1e-9, "expected u={expected}, got {u_retained}");
    assert_eq!(u_retained, u_constrained, "constrained node must move exactly with the retained node");
}

/// `rigid_diaphragm` ties every constrained node's x DOF to the retained
/// node's — same aliasing mechanism as `equal_dof`, exercised across more
/// than one constrained node at once (the 2D-frame simplification: only
/// x-translation is tied, not rotation/lever-arm kinematics — see
/// `Domain::rigid_diaphragm`'s doc comment).
#[test]
fn rigid_diaphragm_couples_every_constrained_node_to_the_retained_node() {
    let (k0, k1, k2, force) = (50.0, 30.0, 10.0, 27.0);
    let mut domain = Domain::new();

    let retained = truss_to_ground(&mut domain, 0.0, k0, force);
    let c1 = truss_to_ground(&mut domain, 200.0, k1, 0.0);
    let c2 = truss_to_ground(&mut domain, 400.0, k2, 0.0);
    domain.rigid_diaphragm(retained, &[c1, c2]);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);

    analysis.step().expect("linear elastic rigid-diaphragm system should solve");

    let expected = force / (k0 + k1 + k2);
    let u0 = analysis.domain().node(retained).displacement[0];
    let u1 = analysis.domain().node(c1).displacement[0];
    let u2 = analysis.domain().node(c2).displacement[0];
    assert!((u0 - expected).abs() < 1e-9, "expected u={expected}, got {u0}");
    assert_eq!(u0, u1);
    assert_eq!(u0, u2);
}

/// `ConstraintHandler::Plain` can't resolve multi-point constraints — using
/// it with a domain that has one is a model-construction error caught at
/// `build`, not a silently-wrong solve.
#[test]
#[should_panic(expected = "Plain")]
fn plain_handler_rejects_a_domain_with_multi_point_constraints() {
    let mut domain = Domain::new();
    let retained = truss_to_ground(&mut domain, 0.0, 50.0, 0.0);
    let constrained = truss_to_ground(&mut domain, 200.0, 30.0, 0.0);
    domain.equal_dof(retained, constrained, &[0]);

    let _ = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
}
