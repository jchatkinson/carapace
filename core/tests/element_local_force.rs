use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{
    BeamIntegration, Domain, Domain3, ElasticBeamColumn, Element, Element3, Fiber, ForceBeamColumn, GeomTransf, Material, Node, Node3,
    SpatialDof,
};

/// Every check here is derived from plain nodal-equilibrium statics (sum of
/// forces/moments = 0), independent of any element's internal formulation —
/// not by re-deriving the element's own formula and checking it agrees with
/// itself. A single-element cantilever is statically determinate, so this
/// pins down `local_force` exactly, not just approximately: at a free,
/// directly-loaded dof, the element's local nodal force must equal the
/// applied load (that's what "solved" means); at a fixed dof, it must equal
/// the reaction the rest of the free body's equilibrium requires.
///
/// Each member below is aligned with global `x`, so its local axes coincide
/// with global — isolating the per-element-type `local_force` logic itself
/// from the separate rotation logic a skew member would also exercise (see
/// `elastic_beam_column_local_force_decomposes_purely_axially_...` below).

#[test]
fn truss_local_force_matches_axial_nodal_equilibrium() {
    let (e, area, length, axial) = (30_000.0, 2.0, 100.0, 500.0);
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, axial);
    let truss = domain.add_element(Element::Truss(carapace_core::model::Truss::new(node_i, node_j, area, Material::Elastic { e })));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    analysis.step().expect("axial truss should solve");

    let f = analysis.domain().element_local_force(truss);
    let expected = [-axial, 0.0, 0.0, axial, 0.0, 0.0];
    for i in 0..6 {
        assert!((f[i] - expected[i]).abs() < 1e-9, "component {i}: got {}, expected {}", f[i], expected[i]);
    }
}

#[test]
fn zero_length_local_force_matches_spring_nodal_equilibrium() {
    let (k, load) = (100.0, 50.0);
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, load);
    let spring = domain.add_element(Element::ZeroLength(
        carapace_core::model::ZeroLength::new(node_i, node_j).with_material(0, Material::Elastic { e: k }),
    ));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    analysis.step().expect("spring should solve");

    let f = analysis.domain().element_local_force(spring);
    let expected = [-load, 0.0, 0.0, load, 0.0, 0.0];
    for i in 0..6 {
        assert!((f[i] - expected[i]).abs() < 1e-9, "component {i}: got {}, expected {}", f[i], expected[i]);
    }
}

/// The expected local force here — `[0, -Fy, -L*Fy, 0, Fy, 0]` — is a pure
/// statics result (see the module doc comment) that must hold regardless of
/// whether the member is `ElasticBeamColumn`, `DispBeamColumn`, or
/// `ForceBeamColumn`: a single-element cantilever is determinate, so its
/// support reactions don't depend on the member's internal stiffness
/// formulation at all.
fn cantilever_local_force_equilibrium(fy: f64, length: f64) -> [f64; 6] {
    [0.0, -fy, -length * fy, 0.0, fy, 0.0]
}

#[test]
fn elastic_beam_column_local_force_matches_cantilever_equilibrium_for_horizontal_member() {
    let (e, area, iz, length, fy): (f64, f64, f64, f64, f64) = (30_000.0, 10.0, 1000.0, 100.0, -10.0);
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]));
    domain.load_node(node_j, 1, fy);
    let beam = domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        node_i,
        node_j,
        e,
        area,
        iz,
        GeomTransf::Linear,
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    analysis.step().expect("cantilever should solve");

    let f = analysis.domain().element_local_force(beam);
    let expected = cantilever_local_force_equilibrium(fy, length);
    for i in 0..6 {
        assert!((f[i] - expected[i]).abs() < 1e-6, "component {i}: got {}, expected {}", f[i], expected[i]);
    }
}

#[test]
fn disp_beam_column_local_force_matches_cantilever_equilibrium_for_horizontal_member() {
    let (e, area, iz, length, fy): (f64, f64, f64, f64, f64) = (30_000.0, 10.0, 1000.0, 100.0, -10.0);
    let h = (iz / area).sqrt();
    let fibers = vec![
        Fiber::new(h, area / 2.0, Material::Elastic { e }),
        Fiber::new(-h, area / 2.0, Material::Elastic { e }),
    ];

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]));
    domain.load_node(node_j, 1, fy);
    let beam = domain.add_element(Element::DispBeamColumn(carapace_core::model::DispBeamColumn::new(
        node_i,
        node_j,
        fibers,
        BeamIntegration::Legendre { points: 3 },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    analysis.step().expect("cantilever should solve");

    let f = analysis.domain().element_local_force(beam);
    let expected = cantilever_local_force_equilibrium(fy, length);
    for i in 0..6 {
        assert!((f[i] - expected[i]).abs() < 1e-6, "component {i}: got {}, expected {}", f[i], expected[i]);
    }
}

/// The highest-risk case: `ForceBeamColumn::local_force` reads a value
/// cached at `commit` time (to avoid re-running the Newton state
/// determination), not recomputed fresh like every other element here — a
/// caching bug (stale value, wrong sign in the cached expression) wouldn't
/// be caught by anything that doesn't actually call `commit`.
#[test]
fn force_beam_column_local_force_matches_cantilever_equilibrium_for_horizontal_member() {
    let (e, area, iz, length, fy): (f64, f64, f64, f64, f64) = (30_000.0, 10.0, 1000.0, 100.0, -10.0);
    let h = (iz / area).sqrt();
    let fibers = vec![
        Fiber::new(h, area / 2.0, Material::Elastic { e }),
        Fiber::new(-h, area / 2.0, Material::Elastic { e }),
    ];

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]));
    domain.load_node(node_j, 1, fy);
    let beam = domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
        node_i,
        node_j,
        fibers,
        BeamIntegration::Lobatto { points: 4 },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    analysis.step().expect("cantilever should solve");

    let f = analysis.domain().element_local_force(beam);
    let expected = cantilever_local_force_equilibrium(fy, length);
    for i in 0..6 {
        assert!((f[i] - expected[i]).abs() < 1e-6, "component {i}: got {}, expected {}", f[i], expected[i]);
    }
}

/// Isolates the local-axis rotation itself (as opposed to the tests above,
/// which use a horizontal member where local coincides with global and so
/// never actually exercise the rotation). A global load applied exactly
/// along a skew member's own axis must decompose to *pure* axial local
/// force — `ElasticBeamColumn`'s local stiffness has zero coupling between
/// its axial row/column and its bending rows/columns, so any nonzero local
/// shear/moment here means the rotation itself (not the element formula)
/// has a sign or axis error.
#[test]
fn elastic_beam_column_local_force_decomposes_purely_axially_for_a_load_along_a_skew_members_own_axis() {
    let (e, area, iz) = (30_000.0, 10.0, 1000.0);
    let (dx, dy, length) = (3.0, 4.0, 5.0);
    let (cx, cy) = (dx / length, dy / length);
    let axial = 100.0;
    let (fx, fy) = (axial * cx, axial * cy);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([dx, dy]));
    domain.load_node(node_j, 0, fx);
    domain.load_node(node_j, 1, fy);
    let beam = domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        node_i,
        node_j,
        e,
        area,
        iz,
        GeomTransf::Linear,
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    analysis.step().expect("skew cantilever should solve");

    let f = analysis.domain().element_local_force(beam);
    assert!((f[0] - -axial).abs() < 1e-6, "local axial reaction: got {}, expected {}", f[0], -axial);
    assert!((f[3] - axial).abs() < 1e-6, "local axial at tip: got {}, expected {}", f[3], axial);
    for i in [1, 2, 4, 5] {
        assert!(f[i].abs() < 1e-6, "component {i} should be ~0 (pure axial), got {}", f[i]);
    }
}

/// Spatial counterpart, exercising `ForceBeamColumn3::local_force`'s cached
/// value (including its decoupled torsion term) against full 3D nodal
/// equilibrium — same technique as the planar cantilever checks above,
/// extended to all six components at each end. Reuses the load case from
/// `m17_force_beam_column3.rs`'s exactness test (axial, biaxial shear,
/// torque all simultaneously) so this also implicitly cross-checks that
/// test's already-verified tip displacements are consistent with these
/// reactions, not just internally self-consistent.
#[test]
fn force_beam_column3_local_force_matches_cantilever_equilibrium_for_horizontal_member() {
    let (e, g, area, iy, iz, j, length): (f64, f64, f64, f64, f64, f64, f64) = (30_000.0, 12_000.0, 4.0, 500.0, 2000.0, 50.0, 100.0);
    let (fy, fz, torque, axial) = (-10.0, -6.0, 20.0, 500.0);

    let hz = (iy / area).sqrt();
    let hy = (iz / area).sqrt();
    let a4 = area / 4.0;
    let fibers = vec![
        carapace_core::model::Fiber3::new(hy, hz, a4, Material::Elastic { e }),
        carapace_core::model::Fiber3::new(hy, -hz, a4, Material::Elastic { e }),
        carapace_core::model::Fiber3::new(-hy, hz, a4, Material::Elastic { e }),
        carapace_core::model::Fiber3::new(-hy, -hz, a4, Material::Elastic { e }),
    ];

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
    domain.load_node(node_j, SpatialDof::Ux as usize, axial);
    domain.load_node(node_j, SpatialDof::Uy as usize, fy);
    domain.load_node(node_j, SpatialDof::Uz as usize, fz);
    domain.load_node(node_j, SpatialDof::Rx as usize, torque);
    let beam = domain.add_element(Element3::ForceBeamColumn3(carapace_core::model::ForceBeamColumn3::new(
        node_i,
        node_j,
        g,
        j,
        [0.0, 0.0, 1.0],
        fibers,
        BeamIntegration::Lobatto { points: 3 },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    analysis.step().expect("spatial cantilever should solve");

    let f = analysis.domain().element_local_force(beam);
    #[rustfmt::skip]
    let expected = [
        -axial, -fy, -fz, -torque, length * fz, -length * fy,
         axial,  fy,  fz,  torque, 0.0,          0.0,
    ];
    for i in 0..12 {
        assert!((f[i] - expected[i]).abs() < 1e-6, "component {i}: got {}, expected {}", f[i], expected[i]);
    }
}
