use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{
    BeamIntegration, DispBeamColumn, Domain, Element, Fiber, ForceBeamColumn, Material, Node,
};

/// `DispBeamColumn::fiber_responses` after a converged pure-axial-load
/// step must reproduce, fiber by fiber, the exact closed-form strain/stress
/// a two-symmetric-elastic-fiber section gives — cross-checked against
/// `fiber_section.rs`'s own `two_symmetric_elastic_fibers_reproduce_ea_ei_
/// exactly`, here read back through the element rather than the bare
/// `FiberSection`: a symmetric section under pure axial extension has zero
/// curvature, so every fiber (regardless of its `y`) should show identical
/// strain `== eps0 == axial_disp/length` and stress `== e*eps0`.
#[test]
fn disp_beam_column_fiber_responses_match_hand_computed_strain_and_stress() {
    let (e, area, iz, length, axial_load): (f64, f64, f64, f64, f64) = (30_000.0, 2.0, 1000.0, 100.0, 60.0);
    let h = (iz / area).sqrt();
    let fibers = vec![
        Fiber::new(h, area / 2.0, Material::Elastic { e }),
        Fiber::new(-h, area / 2.0, Material::Elastic { e }),
    ];

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, axial_load);
    let element_id = domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(
        node_i,
        node_j,
        fibers,
        BeamIntegration::Legendre { points: 2 },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    analysis.step().expect("elastic axial DispBeamColumn should solve");

    let expected_strain = axial_load / (e * area);
    let expected_stress = e * expected_strain;

    let responses = analysis
        .domain()
        .element_fiber_responses(element_id)
        .expect("DispBeamColumn should have fiber responses");
    for point in &responses {
        assert_eq!(point.len(), 2, "two fibers per section");
        for &(strain, stress) in point {
            assert!((strain - expected_strain).abs() < 1e-9, "expected strain {expected_strain}, got {strain}");
            assert!((stress - expected_stress).abs() < 1e-6, "expected stress {expected_stress}, got {stress}");
        }
    }
}

/// `ForceBeamColumn::fiber_responses` after a converged nonlinear step:
/// summing `stress*area` back up across every fiber (the same definition
/// `FiberSection::trial`'s `N` uses) at each integration point must
/// reproduce that point's true axial force (~0 here — pure transverse tip
/// load, no net axial anywhere along a straight cantilever) — a self-
/// consistency check that `fiber_responses` reads the same `e_commit`
/// state `local_force`/`commit` already rely on, not stale or mismatched
/// data.
#[test]
fn force_beam_column_fiber_responses_are_internally_consistent_with_committed_section_state() {
    let (e, area, iz, length, fy): (f64, f64, f64, f64, f64) = (30_000.0, 2.0, 1000.0, 100.0, -10.0);
    let h = (iz / area).sqrt();
    let fibers = vec![
        Fiber::new(h, area / 2.0, Material::elastic_pp(e, 0.01)),
        Fiber::new(-h, area / 2.0, Material::elastic_pp(e, 0.01)),
    ];

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]));
    domain.load_node(node_j, 1, fy);
    let element_id = domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
        node_i,
        node_j,
        fibers,
        BeamIntegration::Lobatto { points: 3 },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-10, max_iter: 30 })
        .build(domain);
    analysis.step().expect("elastic-range force-beam-column cantilever should converge");

    let responses = analysis
        .domain()
        .element_fiber_responses(element_id)
        .expect("ForceBeamColumn should have fiber responses");

    for point in &responses {
        let axial: f64 = point.iter().map(|&(_strain, stress)| stress * (area / 2.0)).sum();
        assert!(axial.abs() < 1e-6, "expected ~0 net axial force, got {axial}");
    }
}
