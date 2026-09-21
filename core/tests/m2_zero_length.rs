use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, Material, Node, ZeroLength};

/// M2 acceptance (implementation-plan §7): prove `Element`/`Material` enum
/// dispatch generalizes beyond `Truss`/`Elastic` by running a `ZeroLength` +
/// `Ent` ("no tension") element through the same real `Domain`/`Analysis`
/// pipeline M1 verified. The applied load stays compressive, so the whole
/// step is within `Ent`'s single linear (engaged) regime — `Algorithm::Linear`
/// (no Newton iteration until M4) is exact for that case, same as it was for
/// `Truss`/`Elastic` in M1.
#[test]
fn zero_length_ent_matches_hand_calc_via_real_architecture() {
    let e = 100.0;
    let load = -10.0; // compressive: engages Ent's elastic branch

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1));
    let node_j = domain.add_node(Node::new([0.0, 0.0]).fix(1).with_load(0, load));

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

    analysis.step().expect("compressive Ent response should solve");

    let dx = analysis.domain().node(node_j).displacement[0];
    // Equilibrium: internal force (e * dx) balances applied load -> dx = load / e.
    assert!((dx - load / e).abs() < 1e-9, "expected dx={}, got {dx}", load / e);
}
