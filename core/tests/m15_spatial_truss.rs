use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain3, Element3, Material, Node3, SpatialDof, Truss3};

/// Spatial counterpart to `m1_truss.rs`'s acceptance case, run through the
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
        domain.add_element(Element3::Truss3(Truss3::new(apex, support, area, Material::Elastic { e })));
    }

    let load = 100.0;
    domain.load_node(apex, SpatialDof::Uz as usize, load);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);

    let result = analysis.step().expect("linear elastic space truss should solve");
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
