use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, AnalysisError, ConstraintHandler, ConvergenceTest,
    Integrator, TangentStrategy,
};
use carapace_core::model::{Domain, Element, Material, Node, NodeId, Truss, ZeroLength};

fn build(domain: Domain) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 25,
        })
        .build(domain)
}

/// The OpenSees uniaxial-material probe: a zero-length element between a
/// grounded node and a node whose only DOF is prescribed — no free DOFs at all.
fn probe(material: Material) -> (Analysis, NodeId, NodeId) {
    let mut domain = Domain::new();
    let a = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let b = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(a, b).with_material(0, material),
    ));
    (build(domain), a, b)
}

fn force(analysis: &Analysis, node: NodeId) -> f64 {
    // Material resisting force on the prescribed node's DOF.
    analysis.domain().reaction(node, 0, 0.0)
}

#[test]
fn elastic_probe_follows_prescribed_strain() {
    let (mut analysis, _a, b) = probe(Material::Elastic { e: 200_000.0 });
    for (target, expected) in [(0.0, 0.0), (0.001, 200.0), (-0.001, -200.0), (0.0, 0.0)] {
        analysis.step_prescribed(&[(b, 0, target)]).unwrap();
        assert_eq!(analysis.domain().node(b).displacement[0], target);
        assert!((force(&analysis, b) - expected).abs() < 1e-9);
    }
}

/// `Ent` carries no tension: tangent and force are both zero at positive
/// strain, which is exactly where `DisplacementControl` has nothing to divide by.
#[test]
fn zero_tangent_branch_reports_zero_force() {
    let (mut analysis, _a, b) = probe(Material::Ent { e: 1000.0 });
    analysis.step_prescribed(&[(b, 0, 0.002)]).unwrap();
    assert_eq!(force(&analysis, b), 0.0);
    analysis.step_prescribed(&[(b, 0, -0.002)]).unwrap();
    assert!((force(&analysis, b) + 2.0).abs() < 1e-12);
}

#[test]
fn committed_history_persists_across_targets() {
    // fy = 100. Load past +yield, reverse past -yield (plastic strain -0.001), then
    // come back to -0.0005: the elastic branch gives +50, a fresh envelope would give -50.
    let (mut analysis, _a, b) = probe(Material::elastic_pp(100_000.0, 0.001));
    let mut forces = vec![];
    for target in [0.002, -0.002, -0.0005] {
        analysis.step_prescribed(&[(b, 0, target)]).unwrap();
        forces.push(force(&analysis, b));
    }
    for (got, want) in forces.iter().zip([100.0, -100.0, 50.0]) {
        assert!((got - want).abs() < 1e-6, "{forces:?}");
    }
}

#[test]
fn free_dofs_respond_to_a_support_displacement() {
    // Support settlement of a truss: node j free in x, node i displaced by 0.5.
    let mut domain = Domain::new();
    let i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let j = domain.add_node(Node::new([100.0, 0.0]).fix(1).fix(2));
    let k = domain.add_node(Node::new([200.0, 0.0]).fix(0).fix(1).fix(2));
    domain.add_element(Element::Truss(Truss::new(
        i,
        j,
        1.0,
        Material::Elastic { e: 100.0 },
    )));
    domain.add_element(Element::Truss(Truss::new(
        j,
        k,
        1.0,
        Material::Elastic { e: 100.0 },
    )));
    let mut analysis = build(domain);
    analysis.step_prescribed(&[(i, 0, 0.5)]).unwrap();
    // Equal stiffness in series: the free node moves half the support's displacement.
    assert!((analysis.domain().node(j).displacement[0] - 0.25).abs() < 1e-9);
}

#[test]
fn invalid_targets_change_nothing() {
    let (mut analysis, _a, b) = probe(Material::Elastic { e: 1.0 });
    analysis.step_prescribed(&[(b, 0, 0.5)]).unwrap();
    for bad in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            analysis.step_prescribed(&[(b, 0, bad)]).unwrap_err(),
            AnalysisError::InvalidConstraint
        );
    }
    assert_eq!(analysis.domain().node(b).displacement[0], 0.5);

    // A free DOF can't be prescribed.
    let mut domain = Domain::new();
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    let mut analysis = build(domain);
    assert_eq!(
        analysis.step_prescribed(&[(free, 0, 1.0)]).unwrap_err(),
        AnalysisError::InvalidConstraint
    );
}
