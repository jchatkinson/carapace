use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, Material, Node, Truss};

/// M1 acceptance criterion (implementation-plan §7): same 2-node truss case
/// as M0 (L=100, A=2, E=30000, P=50 -> dx=0.0833333333), computed through
/// the real Domain/Node/Element::Truss/Material::Elastic/Analysis
/// architecture instead of the closed-form placeholder.
#[test]
fn matches_xara_poc_truss_case_via_real_architecture() {
    let mut domain = Domain::new();

    let node_i = domain.add_node(Node::new([0.0, 100.0]).fix(0).fix(1));
    let node_j = domain.add_node(Node::new([100.0, 100.0]).fix(1).with_load(0, 50.0));

    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        2.0,
        Material::Elastic { e: 30000.0 },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    let result = analysis.step().expect("linear elastic truss should solve");
    assert_eq!(result.step, 1);
    assert_eq!(result.load_factor, 1.0);

    let dx = analysis.domain().node(node_j).displacement[0];
    assert!(
        (dx - 0.0833333333).abs() < 1e-9,
        "expected dx=0.0833333333 (PL/AE), got {dx}"
    );
}
