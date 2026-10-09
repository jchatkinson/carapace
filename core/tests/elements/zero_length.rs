use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{Domain, Element, Material, Node, ZeroLength};

/// A `ZeroLength` + `Ent` ("no tension") element under a compressive load.
/// The whole step stays within `Ent`'s single linear (engaged) regime, so
/// `Algorithm::Linear` (no Newton iteration) is exact.
#[test]
fn zero_length_ent_matches_hand_calc() {
    let e = 100.0;
    let load = -10.0; // compressive: engages Ent's elastic branch

    // dof 2 (rotation) is fixed at both nodes: no material is assigned to
    // it here, so an unconnected rotation DOF would leave the global
    // system singular (rotations are fixed so the system is not singular).
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, load);

    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::Ent { e }),
    ));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    analysis
        .step()
        .expect("compressive Ent response should solve");

    let dx = analysis.domain().node(node_j).displacement[0];
    // Equilibrium: internal force (e * dx) balances applied load -> dx = load / e.
    assert!(
        (dx - load / e).abs() < 1e-9,
        "expected dx={}, got {dx}",
        load / e
    );
}
