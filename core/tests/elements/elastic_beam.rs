//! `ElasticBeamColumn` / `ElasticBeamColumn3` through the full stack: closed-form UDL
//! and combined-load deflections, `PDelta` against `Linear`, and axial/transverse element loads.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    Domain, Domain3, ElasticBeamColumn, ElasticBeamColumn3, Element, Element3, ElementLoad,
    ElementLoad3, GeomTransf, GeomTransf3, Node, Node3, SpatialDof,
};

/// A single `ElasticBeamColumn`
/// spanning a simply-supported beam, under a uniform transverse element
/// load, matches the classical closed-form end rotation
/// `theta = w*L^3 / (24*E*I)`.
///
/// A single 2-node beam element with consistent (virtual-work) equivalent
/// nodal loads exactly reproduces the nodal DOFs of the continuous beam
/// solution for element loads with no interior point loads — the same
/// identity underlying the classical slope-deflection / moment-distribution
/// fixed-end-moment method — so this is an exact check, not an
/// approximation, despite using only one element.
#[test]
fn simply_supported_beam_udl_matches_closed_form_end_rotation() {
    let (e, iz, area) = (30000.0, 1000.0, 100.0);
    let length = 100.0;
    let w = -1.0; // local +y transverse load

    let mut domain = Domain::new();
    // Pin at node_i (ux, uy fixed), roller at node_j (uy fixed) — a simply
    // supported beam. Rotations free at both.
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1));

    let beam = domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        node_i,
        node_j,
        e,
        area,
        iz,
        GeomTransf::Linear,
    )));
    let pattern = domain.default_pattern();
    domain.add_element_load(pattern, beam, ElementLoad::uniform(0.0, w));

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
        .expect("simply supported beam under UDL should solve");

    let theta_i = analysis.domain().node(node_i).displacement[2];
    let theta_j = analysis.domain().node(node_j).displacement[2];

    let expected = (w.abs() * length.powi(3)) / (24.0 * e * iz);
    assert!(
        (theta_i.abs() - expected).abs() < 1e-9,
        "expected |theta_i|={expected}, got {}",
        theta_i.abs()
    );
    assert!(
        (theta_j.abs() - expected).abs() < 1e-9,
        "expected |theta_j|={expected}, got {}",
        theta_j.abs()
    );
    assert!(
        (theta_i + theta_j).abs() < 1e-9,
        "end rotations should be equal and opposite by symmetry"
    );
}

/// `GeomTransf::PDelta` sanity check: with zero axial force it must reduce
/// exactly to `Linear` (the geometric-stiffness correction is a function of
/// axial force, so it vanishes at N=0), and under a compressive axial force
/// it must soften the bending stiffness relative to `Linear` (the intended
/// P-Delta direction — see beam.rs's `geometric_stiffness` doc comment).
///
/// Each case runs as *two* equal `LoadControl` increments, not one: within
/// a single `Algorithm::Linear` solve from a zero initial state, the axial
/// force used for the geometric-stiffness correction is read from the
/// *pre-step* displacement (always zero on a first step from rest), so a
/// one-shot solve can never actually exercise it (see beam.rs's doc
/// comment on why path-following P-Delta needs Newton iteration).
/// Splitting into two increments lets the first establish a nonzero axial
/// force that the second increment's tangent stiffness can then react to —
/// a legitimate (if approximate) incremental-linear P-Delta technique.
#[test]
fn pdelta_matches_linear_at_zero_axial_force_and_softens_under_compression() {
    let (e, iz, area) = (30000.0, 1000.0, 100.0);
    let length = 100.0;

    let build = |transform: GeomTransf, axial_load: f64| {
        let mut domain = Domain::new();
        let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let node_j = domain.add_node(Node::new([length, 0.0]));
        domain.load_node(node_j, 0, axial_load);
        domain.load_node(node_j, 1, -1.0);
        domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            node_i, node_j, e, area, iz, transform,
        )));
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 0.5 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            })
            .build(domain);
        analysis
            .step()
            .expect("cantilever beam should solve (increment 1/2)");
        analysis
            .step()
            .expect("cantilever beam should solve (increment 2/2)");
        analysis.domain().node(node_j).displacement[1]
    };

    // Zero axial force: PDelta's geometric-stiffness correction vanishes,
    // so the transverse-tip-load response must match Linear exactly.
    let v_linear_no_axial = build(GeomTransf::Linear, 0.0);
    let v_pdelta_no_axial = build(GeomTransf::PDelta, 0.0);
    assert!(
        (v_linear_no_axial - v_pdelta_no_axial).abs() < 1e-9,
        "PDelta with zero axial force should match Linear exactly"
    );

    // Compressive axial force: PDelta should soften the beam, i.e. produce
    // a larger-magnitude transverse tip deflection than Linear for the
    // same transverse load.
    let v_linear_compression = build(GeomTransf::Linear, -1000.0);
    let v_pdelta_compression = build(GeomTransf::PDelta, -1000.0);
    assert!(
        v_pdelta_compression.abs() > v_linear_compression.abs(),
        "compressive PDelta should soften the beam relative to Linear: got {v_pdelta_compression} vs {v_linear_compression}"
    );
}

/// A 3D `ElasticBeamColumn3` cantilever, run
/// through the real `Domain3`/`Analysis3` stack (not just the element-level
/// checks in `elastic_beam_column.rs`), matching biaxial bending, axial,
/// and torsion deflections against textbook closed forms simultaneously —
/// the "axial, torsional, and both bending deflections of rotated
/// cantilevers" acceptance criterion from the 3D design notes (member here is axis-aligned rather than rotated, since
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
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    let result = analysis
        .step()
        .expect("linear elastic spatial beam should solve");
    assert_eq!(result.load_factor, 1.0);

    let node = analysis.domain().node(free);
    let tol = 1e-9;
    assert!(
        (node.displacement[SpatialDof::Ux as usize] - axial_force * length / (e * a)).abs() < tol
    );
    assert!(
        (node.displacement[SpatialDof::Uy as usize] - fy * length.powi(3) / (3.0 * e * iz)).abs()
            < tol
    );
    assert!(
        (node.displacement[SpatialDof::Rz as usize] - fy * length.powi(2) / (2.0 * e * iz)).abs()
            < tol
    );
    assert!(
        (node.displacement[SpatialDof::Uz as usize] - fz * length.powi(3) / (3.0 * e * iy)).abs()
            < tol
    );
    assert!(
        (node.displacement[SpatialDof::Ry as usize] - (-fz * length.powi(2) / (2.0 * e * iy)))
            .abs()
            < tol
    );
    assert!((node.displacement[SpatialDof::Rx as usize] - torque * length / (g * j)).abs() < tol);
}

/// `GeomTransf3::PDelta3` sanity check, the spatial analogue of
/// `elements/elastic_beam.rs`'s planar `pdelta_matches_linear_at_zero_axial_force_and_
/// softens_under_compression` (see its doc comment for why each case runs
/// as two `LoadControl` increments rather than one): at zero axial force
/// `PDelta3` must reduce exactly to `Linear3`, and under compression it
/// must soften the bending stiffness relative to `Linear3`.
///
/// Run in *both* bending planes (`Uy`, governed by `Iz`, and `Uz`, governed
/// by `Iy`) — not just one — because `ElasticBeamColumn3::geometric_
/// stiffness`'s `y`-plane block is a hand-derived sign flip of the
/// `z`-plane one (`ry = -dw/dx` vs `rz = +dv/dx`), not copied from a
/// reference; a bug there would only show up by actually exercising the
/// `Iy` plane, not just `Iz`.
#[test]
fn pdelta3_matches_linear3_at_zero_axial_force_and_softens_under_compression_in_both_planes() {
    let (e, g, a, j, iy, iz) = (30_000.0, 12_000.0, 100.0, 100.0, 1000.0, 2000.0);
    let length = 100.0;

    let build = |transform: GeomTransf3, axial_load: f64, transverse_dof: usize| {
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
        let free = domain.add_node(Node3::new([length, 0.0, 0.0]));
        domain.load_node(free, SpatialDof::Ux as usize, axial_load);
        domain.load_node(free, transverse_dof, -1.0);
        domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
            fixed, free, e, g, a, j, iy, iz, transform,
        )));

        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 0.5 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            })
            .build(domain);

        analysis
            .step()
            .expect("cantilever beam should solve (increment 1/2)");
        analysis
            .step()
            .expect("cantilever beam should solve (increment 2/2)");
        analysis.domain().node(free).displacement[transverse_dof]
    };

    for transverse_dof in [SpatialDof::Uy as usize, SpatialDof::Uz as usize] {
        let v_linear_no_axial = build(GeomTransf3::linear([0.0, 0.0, 1.0]), 0.0, transverse_dof);
        let v_pdelta_no_axial = build(GeomTransf3::p_delta([0.0, 0.0, 1.0]), 0.0, transverse_dof);
        assert!(
            (v_linear_no_axial - v_pdelta_no_axial).abs() < 1e-9,
            "dof {transverse_dof}: PDelta3 with zero axial force should match Linear3 exactly"
        );

        let v_linear_compression = build(
            GeomTransf3::linear([0.0, 0.0, 1.0]),
            -1000.0,
            transverse_dof,
        );
        let v_pdelta_compression = build(
            GeomTransf3::p_delta([0.0, 0.0, 1.0]),
            -1000.0,
            transverse_dof,
        );
        assert!(
            v_pdelta_compression.abs() > v_linear_compression.abs(),
            "dof {transverse_dof}: compressive PDelta3 should soften the beam relative to Linear3: got {v_pdelta_compression} vs {v_linear_compression}"
        );
    }
}

/// Spatial counterpart to `elements/elastic_beam.rs`'s `simply_supported_beam_udl_
/// matches_closed_form_end_rotation` — `ElasticBeamColumn3`'s new local
/// `wy`/`wz` uniform transverse load (`ElementLoad3::uniform`,
/// the spatial generalization of planar `ElementLoad::Uniform`;
/// see the 3D beam-load design),
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
    domain.add_element_load(pattern, beam, ElementLoad3::uniform(0.0, wy, wz));

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
        .expect("simply supported spatial beam under biaxial UDL should solve");

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
    assert!(
        (rz_i + rz_j).abs() < 1e-9,
        "z-bending end rotations should be equal and opposite by symmetry"
    );

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
    assert!(
        (ry_i + ry_j).abs() < 1e-9,
        "y-bending end rotations should be equal and opposite by symmetry"
    );
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
            domain.add_element_load(pattern, beam, ElementLoad3::uniform(0.0, 0.0, 0.0));
        }
        domain.load_node(node_j, SpatialDof::Uy as usize, 1.0);

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
            .expect("cantilever spatial beam should solve");
        analysis.domain().node(node_j).displacement[SpatialDof::Uy as usize]
    };

    let v_with_zero_load = build(true);
    let v_without_load = build(false);
    assert!(
        (v_with_zero_load - v_without_load).abs() < 1e-12,
        "an explicit zero UniformTransverse load must match having no element load at all"
    );
}

const E: f64 = 30_000.0;
const A: f64 = 10.0;
const IZ: f64 = 1000.0;
const L: f64 = 100.0;

fn solve_planar_cantilever(
    loads: &[(f64, f64)],
) -> (
    Domain,
    carapace_core::model::ElementId,
    carapace_core::model::NodeId,
) {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([L, 0.0]));
    let beam = domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        node_i,
        node_j,
        E,
        A,
        IZ,
        GeomTransf::Linear,
    )));
    let pattern = domain.default_pattern();
    for &(wx, wy) in loads {
        domain.add_element_load(pattern, beam, ElementLoad::uniform(wx, wy));
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    analysis.step().expect("cantilever should solve");
    (analysis.domain().clone(), beam, node_j)
}

/// Axial + transverse UDL on a cantilever: tip `ux = wx L² / 2EA`, tip
/// `uy = wy L⁴ / 8EI`, and the member end forces equal the statics result —
/// the free end carries nothing, the fixed end carries the whole load.
#[test]
fn planar_cantilever_axial_and_transverse_udl_match_closed_form() {
    let (wx, wy) = (2.0, -1.0);
    let (domain, beam, node_j) = solve_planar_cantilever(&[(wx, wy)]);

    let ux_tip = domain.node(node_j).displacement[0];
    assert!((ux_tip - wx * L * L / (2.0 * E * A)).abs() < 1e-12);
    let uy_tip = domain.node(node_j).displacement[1];
    assert!((uy_tip - wy * L.powi(4) / (8.0 * E * IZ)).abs() < 1e-9);

    let f = domain.element_end_force(beam, 1.0);
    let expected = [-wx * L, -wy * L, -wy * L * L / 2.0, 0.0, 0.0, 0.0];
    for i in 0..6 {
        assert!(
            (f[i] - expected[i]).abs() < 1e-6,
            "component {i}: got {}, expected {}",
            f[i],
            expected[i]
        );
    }
}

/// Several loads on one element in one pattern sum rather than the last one
/// winning (OpenSees semantics).
#[test]
fn loads_on_the_same_element_accumulate() {
    let (d, _, n) = solve_planar_cantilever(&[(1.0, 0.0), (1.0, 0.0)]);
    let ux_split = d.node(n).displacement[0];
    let (d, _, n) = solve_planar_cantilever(&[(2.0, 0.0)]);
    let ux_whole = d.node(n).displacement[0];
    assert!(ux_whole.abs() > 0.0);
    assert!((ux_split - ux_whole).abs() < 1e-15);
}

/// Spatial counterpart: a member along global x with `vec_xz = [0,0,1]` has
/// local axes equal to global, so `wx` is plain axial.
#[test]
fn spatial_cantilever_axial_udl_matches_closed_form() {
    let (g, j, iy, iz) = (12_000.0, 500.0, 1500.0, 1000.0);
    let wx = 3.0;
    let mut domain = Domain3::new();
    let mut fixed = Node3::new([0.0, 0.0, 0.0]);
    for dof in [
        SpatialDof::Ux,
        SpatialDof::Uy,
        SpatialDof::Uz,
        SpatialDof::Rx,
        SpatialDof::Ry,
        SpatialDof::Rz,
    ] {
        fixed = fixed.fix(dof as usize);
    }
    let node_i = domain.add_node(fixed);
    let node_j = domain.add_node(Node3::new([L, 0.0, 0.0]));
    let beam = domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
        node_i,
        node_j,
        E,
        g,
        A,
        j,
        iy,
        iz,
        GeomTransf3::linear([0.0, 0.0, 1.0]),
    )));
    let pattern = domain.default_pattern();
    domain.add_element_load(pattern, beam, ElementLoad3::uniform(wx, 0.0, 0.0));
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    analysis.step().expect("spatial cantilever should solve");

    let ux = analysis.domain().node(node_j).displacement[SpatialDof::Ux as usize];
    assert!((ux - wx * L * L / (2.0 * E * A)).abs() < 1e-12);
    let f = analysis.domain().element_end_force(beam, 1.0);
    assert!((f[0] + wx * L).abs() < 1e-6, "fixed-end axial {}", f[0]);
    assert!(f[6].abs() < 1e-6, "free-end axial {}", f[6]);
}
