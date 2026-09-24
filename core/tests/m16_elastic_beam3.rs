use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain3, Element3, ElasticBeamColumn3, GeomTransf3, Node3, SpatialDof};

/// M16 acceptance case: a spatial `ElasticBeamColumn3` cantilever, run
/// through the real `Domain3`/`Analysis3` stack (not just the element-level
/// checks in `elastic_beam_column.rs`), matching biaxial bending, axial,
/// and torsion deflections against textbook closed forms simultaneously —
/// the "axial, torsional, and both bending deflections of rotated
/// cantilevers" acceptance criterion from `docs/spatial-architecture.md`'s
/// M16 entry (member here is axis-aligned rather than rotated, since
/// `elastic_beam_column.rs`'s own tests already isolate the `vec_xz`
/// transform on a skew member; this test's job is the `Domain3` assembly
/// path, not the transform).
#[test]
fn cantilever_elastic_beam3_matches_closed_form_under_combined_loading() {
    let mut domain = Domain3::new();

    let fixed = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(SpatialDof::Ux as usize)
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize),
    );

    let length = 100.0;
    let free = domain.add_node(Node3::new([length, 0.0, 0.0]));

    let (e, g, a, j, iy, iz) = (30_000.0, 12_000.0, 20.0, 5.0, 400.0, 800.0);
    domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
        fixed,
        free,
        e,
        g,
        a,
        j,
        iy,
        iz,
        GeomTransf3::linear([0.0, 0.0, 1.0]),
    )));

    let (axial_force, fy, fz, torque) = (100.0, 50.0, 30.0, 20.0);
    domain.load_node(free, SpatialDof::Ux as usize, axial_force);
    domain.load_node(free, SpatialDof::Uy as usize, fy);
    domain.load_node(free, SpatialDof::Uz as usize, fz);
    domain.load_node(free, SpatialDof::Rx as usize, torque);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);

    let result = analysis.step().expect("linear elastic spatial beam should solve");
    assert_eq!(result.load_factor, 1.0);

    let node = analysis.domain().node(free);
    let tol = 1e-9;
    assert!((node.displacement[SpatialDof::Ux as usize] - axial_force * length / (e * a)).abs() < tol);
    assert!((node.displacement[SpatialDof::Uy as usize] - fy * length.powi(3) / (3.0 * e * iz)).abs() < tol);
    assert!((node.displacement[SpatialDof::Rz as usize] - fy * length.powi(2) / (2.0 * e * iz)).abs() < tol);
    assert!((node.displacement[SpatialDof::Uz as usize] - fz * length.powi(3) / (3.0 * e * iy)).abs() < tol);
    assert!((node.displacement[SpatialDof::Ry as usize] - (-fz * length.powi(2) / (2.0 * e * iy))).abs() < tol);
    assert!((node.displacement[SpatialDof::Rx as usize] - torque * length / (g * j)).abs() < tol);
}
