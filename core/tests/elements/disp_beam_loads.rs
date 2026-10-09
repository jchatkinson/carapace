//! `DispBeamColumn` / `DispBeamColumn3` uniform element loads: elastic
//! equivalence with the closed-form beams, and a yielding cantilever benchmarked
//! against OpenSeesPy.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    BeamIntegration, DispBeamColumn, DispBeamColumn3, Domain, Domain3, ElasticBeamColumn,
    ElasticBeamColumn3, Element, Element3, ElementLoad, ElementLoad3, Fiber, Fiber3, GeomTransf,
    GeomTransf3, Material, Node, Node3, SpatialDof,
};

const E: f64 = 30_000.0;
const L: f64 = 100.0;

fn solve_linear(domain: Domain) -> Domain {
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    analysis.step().expect("linear solve");
    analysis.domain().clone()
}

/// A fiber disp beam with two symmetric elastic fibers reproduces `E·A`/`E·I`
/// exactly, so under the same uniform load it must match `ElasticBeamColumn`:
/// same nodal displacements and same member end forces (the latter only if the
/// load's fixed-end effect is subtracted identically).
#[test]
fn planar_disp_beam_matches_elastic_beam_under_axial_and_transverse_udl() {
    let (area, iz): (f64, f64) = (2.0, 1000.0);
    let (wx, wy) = (2.0, -1.5);
    let h = (iz / area).sqrt();

    let build = |disp: bool| {
        let mut domain = Domain::new();
        let i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let j = domain.add_node(Node::new([L, 0.0]));
        let beam = if disp {
            let fibers = vec![
                Fiber::new(h, area / 2.0, Material::Elastic { e: E }),
                Fiber::new(-h, area / 2.0, Material::Elastic { e: E }),
            ];
            domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(
                i,
                j,
                fibers,
                BeamIntegration::Lobatto { points: 5 },
            )))
        } else {
            domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
                i,
                j,
                E,
                area,
                iz,
                GeomTransf::Linear,
            )))
        };
        let pattern = domain.default_pattern();
        domain.add_element_load(pattern, beam, ElementLoad::uniform(wx, wy));
        (solve_linear(domain), beam, j)
    };

    let (d_disp, b_disp, j_disp) = build(true);
    let (d_elas, b_elas, j_elas) = build(false);
    for dof in 0..3 {
        let (u1, u2) = (
            d_disp.node(j_disp).displacement[dof],
            d_elas.node(j_elas).displacement[dof],
        );
        assert!((u1 - u2).abs() < 1e-9, "dof {dof}: {u1} vs {u2}");
    }
    let (f1, f2) = (
        d_disp.element_end_force(b_disp, 1.0),
        d_elas.element_end_force(b_elas, 1.0),
    );
    for c in 0..6 {
        assert!(
            (f1[c] - f2[c]).abs() < 1e-6,
            "end force {c}: {} vs {}",
            f1[c],
            f2[c]
        );
    }
    // Statics: free end carries nothing, fixed end carries the whole load.
    let expected = [-wx * L, -wy * L, -wy * L * L / 2.0, 0.0, 0.0, 0.0];
    for c in 0..6 {
        assert!((f1[c] - expected[c]).abs() < 1e-6, "statics {c}: {}", f1[c]);
    }
}

/// Spatial counterpart: four fibers at `(±h, ±h)` give `Iy = Iz`; `wx`, `wy`
/// and `wz` all nonzero on a member along global x.
#[test]
fn spatial_disp_beam_matches_elastic_beam_under_triaxial_udl() {
    let (area, i_bend, g, j): (f64, f64, f64, f64) = (4.0, 1000.0, 12_000.0, 500.0);
    let (wx, wy, wz) = (1.0, -2.0, 0.75);
    let h = (i_bend / area).sqrt();

    let fixed = || {
        let mut node = Node3::new([0.0, 0.0, 0.0]);
        for dof in [
            SpatialDof::Ux,
            SpatialDof::Uy,
            SpatialDof::Uz,
            SpatialDof::Rx,
            SpatialDof::Ry,
            SpatialDof::Rz,
        ] {
            node = node.fix(dof as usize);
        }
        node
    };
    let solve = |domain: Domain3| {
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            })
            .build(domain);
        analysis.step().expect("linear solve");
        analysis.domain().clone()
    };
    let build = |disp: bool| {
        let mut domain = Domain3::new();
        let i = domain.add_node(fixed());
        let jn = domain.add_node(Node3::new([L, 0.0, 0.0]));
        let vec_xz = [0.0, 0.0, 1.0];
        let beam = if disp {
            let fibers = [(h, h), (h, -h), (-h, h), (-h, -h)]
                .iter()
                .map(|&(y, z)| Fiber3::new(y, z, area / 4.0, Material::Elastic { e: E }))
                .collect();
            domain.add_element(Element3::DispBeamColumn3(DispBeamColumn3::new(
                i,
                jn,
                g,
                j,
                vec_xz,
                fibers,
                BeamIntegration::Lobatto { points: 5 },
            )))
        } else {
            domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
                i,
                jn,
                E,
                g,
                area,
                j,
                i_bend,
                i_bend,
                GeomTransf3::linear(vec_xz),
            )))
        };
        let pattern = domain.default_pattern();
        domain.add_element_load(pattern, beam, ElementLoad3::Uniform { wx, wy, wz });
        (solve(domain), beam, jn)
    };

    let (d_disp, b_disp, j_disp) = build(true);
    let (d_elas, b_elas, j_elas) = build(false);
    for dof in 0..6 {
        let (u1, u2) = (
            d_disp.node(j_disp).displacement[dof],
            d_elas.node(j_elas).displacement[dof],
        );
        assert!((u1 - u2).abs() < 1e-9, "dof {dof}: {u1} vs {u2}");
    }
    let (f1, f2) = (
        d_disp.element_end_force(b_disp, 1.0),
        d_elas.element_end_force(b_elas, 1.0),
    );
    for c in 0..12 {
        assert!(
            (f1[c] - f2[c]).abs() < 1e-6,
            "end force {c}: {} vs {}",
            f1[c],
            f2[c]
        );
    }
}

/// Yielding cantilever (Steel01, `fy = 30`, `E = 30000`, `b = 0.01`; 10 × 20
/// rectangle, 20 fiber layers) under a ramped UDL `wy = -5.5` (fixed-end
/// moment 27500 vs first yield 20000), `n` equal `DispBeamColumn` elements,
/// Lobatto `points`, 10 load steps. Reference values are OpenSeesPy 3.8
/// `dispBeamColumn` with the identical model (`eleLoad -beamUniform`).
fn yielding_cantilever(n: usize, points: usize) -> (f64, f64) {
    let (b, h) = (10.0, 20.0);
    let layers = 20;
    let fibers = || -> Vec<Fiber> {
        (0..layers)
            .map(|k| {
                let y = -h / 2.0 + (k as f64 + 0.5) * h / layers as f64;
                Fiber::new(
                    y,
                    b * h / layers as f64,
                    Material::steel01(30.0, E, 0.01, 0.0, 1.0, 0.0, 1.0),
                )
            })
            .collect()
    };
    let mut domain = Domain::new();
    let nodes: Vec<_> = (0..=n)
        .map(|k| {
            let node = Node::new([L * k as f64 / n as f64, 0.0]);
            domain.add_node(if k == 0 {
                node.fix(0).fix(1).fix(2)
            } else {
                node
            })
        })
        .collect();
    let pattern = domain.default_pattern();
    for k in 0..n {
        let beam = domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(
            nodes[k],
            nodes[k + 1],
            fibers(),
            BeamIntegration::Lobatto { points },
        )));
        domain.add_element_load(pattern, beam, ElementLoad::uniform(0.0, -5.5));
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.1 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-8,
            max_iter: 50,
        })
        .build(domain);
    for _ in 0..10 {
        analysis
            .step()
            .expect("yielding cantilever should converge");
    }
    let tip = analysis.domain().node(nodes[n]);
    (tip.displacement[1], tip.displacement[2])
}

#[test]
fn yielding_disp_cantilever_under_udl_matches_opensees() {
    // (elements, Lobatto points, OpenSees tip uy, OpenSees tip rz)
    let cases = [
        (1, 3, -0.34985327048913684, -0.00464723780096489),
        (1, 5, -0.3459756573110966, -0.0046084616691844825),
        (8, 3, -0.362380673289648, -0.004779128518230232),
        (8, 5, -0.3628563331496834, -0.004783698676346214),
    ];
    for (n, points, uy_ref, rz_ref) in cases {
        let (uy, rz) = yielding_cantilever(n, points);
        assert!(
            (uy - uy_ref).abs() < 1e-6 * uy_ref.abs() && (rz - rz_ref).abs() < 1e-6 * rz_ref.abs(),
            "{n} elements, {points} points: got ({uy}, {rz}), OpenSees ({uy_ref}, {rz_ref})"
        );
    }
}
