use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain3, Element3, ElementLoad3, ElasticBeamColumn3, GeomTransf3, Node3, SpatialDof};

/// Spatial counterpart to `m3_beam.rs`'s `simply_supported_beam_udl_
/// matches_closed_form_end_rotation` — `ElasticBeamColumn3`'s new local
/// `wy`/`wz` uniform transverse load (`ElementLoad3::UniformTransverse`,
/// the spatial generalization of planar `ElementLoad::UniformTransverse`;
/// see `spatial-architecture.md`'s "Spatial beam loads" milestone),
/// exercised through the full `Domain3`/`Analysis3` stack with a member
/// aligned along global `x` (`vec_xz = [0,0,1]` makes local `[x,y,z]` equal
/// global `[x,y,z]` exactly, isolating the load formula itself from
/// `GeomTransf3`'s rotation — same isolation strategy as
/// `elastic_beam_column.rs`'s `spatial_tests`).
///
/// Both bending planes are loaded simultaneously (`wy` and `wz` both
/// nonzero) to confirm they're independent (no cross term) and each
/// matches the same closed form as the planar case,
/// `theta = w*L^3/(24*E*I)`, with `Iz` governing the `wy`-driven `rz`
/// rotation and `Iy` governing the `wz`-driven `ry` rotation.
///
/// Boundary conditions: node_i is a pin (`Ux,Uy,Uz` fixed) plus torsion
/// fixed (`Rx`, since a doubly-free torsion DOF pair is a zero-energy rigid
/// rotation about the member axis — see `GeomTransf3`'s discussion of why
/// spatial models need every DOF meaningfully restrained); node_j is a
/// roller (`Uy,Uz` fixed, axial and every rotation free) — the direct
/// spatial generalization of the planar test's pin/roller pair.
#[test]
fn axis_aligned_biaxial_udl_matches_closed_form_end_rotations_in_both_planes() {
    let (e, g, area, j, iy, iz) = (30_000.0, 12_000.0, 100.0, 500.0, 1500.0, 1000.0);
    let length = 100.0;
    let (wy, wz) = (-1.0, 0.7);

    let mut domain = Domain3::new();
    let node_i = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(SpatialDof::Ux as usize)
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize),
    );
    let node_j = domain.add_node(
        Node3::new([length, 0.0, 0.0])
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize),
    );

    let beam = domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
        node_i,
        node_j,
        e,
        g,
        area,
        j,
        iy,
        iz,
        GeomTransf3::linear([0.0, 0.0, 1.0]),
    )));
    let pattern = domain.default_pattern();
    domain.add_element_load(pattern, beam, ElementLoad3::UniformTransverse { wy, wz });

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);

    analysis.step().expect("simply supported spatial beam under biaxial UDL should solve");

    let rz_i = analysis.domain().node(node_i).displacement[SpatialDof::Rz as usize];
    let rz_j = analysis.domain().node(node_j).displacement[SpatialDof::Rz as usize];
    let ry_i = analysis.domain().node(node_i).displacement[SpatialDof::Ry as usize];
    let ry_j = analysis.domain().node(node_j).displacement[SpatialDof::Ry as usize];

    let expected_z = (wy.abs() * length.powi(3)) / (24.0 * e * iz);
    let expected_y = (wz.abs() * length.powi(3)) / (24.0 * e * iy);

    assert!(
        (rz_i.abs() - expected_z).abs() < 1e-9,
        "expected |rz_i|={expected_z}, got {}",
        rz_i.abs()
    );
    assert!(
        (rz_j.abs() - expected_z).abs() < 1e-9,
        "expected |rz_j|={expected_z}, got {}",
        rz_j.abs()
    );
    assert!((rz_i + rz_j).abs() < 1e-9, "z-bending end rotations should be equal and opposite by symmetry");

    assert!(
        (ry_i.abs() - expected_y).abs() < 1e-9,
        "expected |ry_i|={expected_y}, got {}",
        ry_i.abs()
    );
    assert!(
        (ry_j.abs() - expected_y).abs() < 1e-9,
        "expected |ry_j|={expected_y}, got {}",
        ry_j.abs()
    );
    assert!((ry_i + ry_j).abs() < 1e-9, "y-bending end rotations should be equal and opposite by symmetry");
}

/// Zero load must still assemble (the `wy == 0.0 && wz == 0.0` fast path in
/// `ElasticBeamColumn3::form_load_vector`), and a pattern with no element
/// load attached at all must behave identically to one whose load is
/// exactly zero — the spatial counterpart of the implicit "no load"
/// coverage every other element-load test gets for free by simply not
/// calling `add_element_load`.
#[test]
fn zero_uniform_transverse_load_matches_no_load_at_all() {
    let (e, g, area, j, iy, iz) = (30_000.0, 12_000.0, 100.0, 500.0, 1500.0, 1000.0);
    let length = 100.0;

    let build = |with_zero_load: bool| {
        let mut domain = Domain3::new();
        let node_i = domain.add_node(
            Node3::new([0.0, 0.0, 0.0])
                .fix(SpatialDof::Ux as usize)
                .fix(SpatialDof::Uy as usize)
                .fix(SpatialDof::Uz as usize)
                .fix(SpatialDof::Rx as usize)
                .fix(SpatialDof::Ry as usize)
                .fix(SpatialDof::Rz as usize),
        );
        let node_j = domain.add_node(Node3::new([length, 0.0, 0.0]));
        let beam = domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
            node_i,
            node_j,
            e,
            g,
            area,
            j,
            iy,
            iz,
            GeomTransf3::linear([0.0, 0.0, 1.0]),
        )));
        if with_zero_load {
            let pattern = domain.default_pattern();
            domain.add_element_load(pattern, beam, ElementLoad3::UniformTransverse { wy: 0.0, wz: 0.0 });
        }
        domain.load_node(node_j, SpatialDof::Uy as usize, 1.0);

        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
            .build(domain);
        analysis.step().expect("cantilever spatial beam should solve");
        analysis.domain().node(node_j).displacement[SpatialDof::Uy as usize]
    };

    let v_with_zero_load = build(true);
    let v_without_load = build(false);
    assert!(
        (v_with_zero_load - v_without_load).abs() < 1e-12,
        "an explicit zero UniformTransverse load must match having no element load at all"
    );
}
