//! `Truss` / `Truss3` through the full `Domain`/`Analysis` stack: 2D closed-form axial
//! displacement and the symmetric 3D tripod.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    Domain, Domain3, Element, Element3, Material, Node, Node3, SpatialDof, Truss, Truss3,
};

/// A 2-node truss (L=100, A=2, E=30000, P=50) gives the closed-form
/// displacement `P*L/(A*E) = 0.0833333333` through `Domain`/`Analysis`.
#[test]
fn truss_matches_closed_form_axial_displacement() {
    let mut domain = Domain::new();

    // dof 2 (rotation) is fixed at both nodes: a Truss carries no moment,
    // so an unconnected rotation DOF would leave the global system singular
    // (rotations are fixed so the system is not singular).
    let node_i = domain.add_node(Node::new([0.0, 100.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([100.0, 100.0]).fix(1).fix(2));
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

/// Spatial counterpart to `elements/truss.rs`'s acceptance case, run through the
/// real `Domain3`/`Node3`/`Element3::Truss3`/`Analysis3` architecture (the
/// spatial instantiation of the now-generic `Domain`/`Analysis` — see their
/// doc comments) — a
/// symmetric tripod of three truss members meeting at one free apex node,
/// each dropping to a fixed base support spaced 120 degrees apart — the
/// standard space-truss case with a clean closed form, and genuinely
/// exercises all three translational DOFs at once (unlike a single truss,
/// which only stiffens the direction along its own axis and leaves a
/// singular system in the other two — that's not a valid spatial model on
/// its own).
///
/// Base supports lie on a radius-3 circle at z=0; the apex sits at (0,0,4),
/// so every bar is a 3-4-5 triangle (L=5). By the 3-fold symmetry, the
/// assembled stiffness at the apex is diagonal (cross terms cancel), so a
/// purely vertical load produces a purely vertical displacement:
/// `dz = P / (3 * A * E * H^2 / L^3)`.
#[test]
fn symmetric_tripod_matches_analytical_vertical_stiffness_via_analysis3() {
    let mut domain = Domain3::new();

    let radius: f64 = 3.0;
    let height = 4.0;
    let angles = [0.0_f64, 120.0_f64.to_radians(), 240.0_f64.to_radians()];

    let apex = domain.add_node(
        Node3::new([0.0, 0.0, height])
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize),
    );

    let area = 2.0;
    let e = 30_000.0;
    for angle in angles {
        let support = domain.add_node(
            Node3::new([radius * angle.cos(), radius * angle.sin(), 0.0])
                .fix(SpatialDof::Ux as usize)
                .fix(SpatialDof::Uy as usize)
                .fix(SpatialDof::Uz as usize)
                .fix(SpatialDof::Rx as usize)
                .fix(SpatialDof::Ry as usize)
                .fix(SpatialDof::Rz as usize),
        );
        domain.add_element(Element3::Truss3(Truss3::new(
            apex,
            support,
            area,
            Material::Elastic { e },
        )));
    }

    let load = 100.0;
    domain.load_node(apex, SpatialDof::Uz as usize, load);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    let result = analysis
        .step()
        .expect("linear elastic space truss should solve");
    assert_eq!(result.step, 1);
    assert_eq!(result.load_factor, 1.0);

    let length: f64 = (radius * radius + height * height).sqrt();
    let k_zz = 3.0 * area * e * height * height / length.powi(3);
    let expected_dz = load / k_zz;

    let node = analysis.domain().node(apex);
    assert!((node.displacement[SpatialDof::Uz as usize] - expected_dz).abs() < 1e-9);
    assert!(node.displacement[SpatialDof::Ux as usize].abs() < 1e-9);
    assert!(node.displacement[SpatialDof::Uy as usize].abs() < 1e-9);
}
