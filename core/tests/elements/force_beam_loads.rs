//! `ForceBeamColumn` / `ForceBeamColumn3` uniform element loads: elastic
//! equivalence with the closed-form beams (displacements *and* member end
//! forces), and yielding members benchmarked against OpenSeesPy.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    BeamIntegration, Domain, Domain3, ElasticBeamColumn, ElasticBeamColumn3, Element, Element3,
    ElementLoad, ElementLoad3, Fiber, Fiber3, ForceBeamColumn, ForceBeamColumn3, GeomTransf,
    GeomTransf3, LoadSeries, Material, Node, Node3, SpatialDof,
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

/// A fiber force beam with two symmetric elastic fibers reproduces `E·A`/`E·I`
/// exactly, so under the same uniform load it must match `ElasticBeamColumn`:
/// same nodal displacements and same member end forces (the latter only if the
/// load's fixed-end effect is subtracted identically).
#[test]
fn planar_force_beam_matches_elastic_beam_under_axial_and_transverse_udl() {
    let (area, iz): (f64, f64) = (2.0, 1000.0);
    let (wx, wy) = (2.0, -1.5);
    let h = (iz / area).sqrt();

    let build = |force: bool| {
        let mut domain = Domain::new();
        let i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let j = domain.add_node(Node::new([L, 0.0]));
        let beam = if force {
            let fibers = vec![
                Fiber::new(h, area / 2.0, Material::Elastic { e: E }),
                Fiber::new(-h, area / 2.0, Material::Elastic { e: E }),
            ];
            domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
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

    let (d_force, b_force, j_force) = build(true);
    let (d_elas, b_elas, j_elas) = build(false);
    for dof in 0..3 {
        let (u1, u2) = (
            d_force.node(j_force).displacement[dof],
            d_elas.node(j_elas).displacement[dof],
        );
        assert!((u1 - u2).abs() < 1e-9, "dof {dof}: {u1} vs {u2}");
    }
    let (f1, f2) = (
        d_force.element_end_force(b_force, 1.0),
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
fn spatial_force_beam_matches_elastic_beam_under_triaxial_udl() {
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
    let build = |force: bool| {
        let mut domain = Domain3::new();
        let i = domain.add_node(fixed());
        let jn = domain.add_node(Node3::new([L, 0.0, 0.0]));
        let vec_xz = [0.0, 0.0, 1.0];
        let beam = if force {
            let fibers = [(h, h), (h, -h), (-h, h), (-h, -h)]
                .iter()
                .map(|&(y, z)| Fiber3::new(y, z, area / 4.0, Material::Elastic { e: E }))
                .collect();
            domain.add_element(Element3::ForceBeamColumn3(ForceBeamColumn3::new(
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
        domain.add_element_load(pattern, beam, ElementLoad3::uniform(wx, wy, wz));
        (solve(domain), beam, jn)
    };

    let (d_force, b_force, j_force) = build(true);
    let (d_elas, b_elas, j_elas) = build(false);
    for dof in 0..6 {
        let (u1, u2) = (
            d_force.node(j_force).displacement[dof],
            d_elas.node(j_elas).displacement[dof],
        );
        assert!((u1 - u2).abs() < 1e-9, "dof {dof}: {u1} vs {u2}");
    }
    let (f1, f2) = (
        d_force.element_end_force(b_force, 1.0),
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

/// `n` equal `ForceBeamColumn`s of the yielding model (Steel01 `fy = 30`,
/// `E = 30000`, `b = 0.01`; 10 x 20 rectangle in 20 fiber layers) under a ramped
/// UDL (`wx`, `wy`), 10 load steps, Lobatto `points`. `fixed_both` clamps both
/// ends; returns the displacement `[ux, uy, rz]` of the free end (cantilever) or
/// the midspan node (clamped-clamped).
fn yielding_member(n: usize, points: usize, wx: f64, wy: f64, fixed_both: bool) -> [f64; 3] {
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
            let clamp = k == 0 || (fixed_both && k == n);
            domain.add_node(if clamp {
                node.fix(0).fix(1).fix(2)
            } else {
                node
            })
        })
        .collect();
    let pattern = domain.default_pattern();
    for k in 0..n {
        let beam = domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
            nodes[k],
            nodes[k + 1],
            fibers(),
            BeamIntegration::Lobatto { points },
        )));
        domain.add_element_load(pattern, beam, ElementLoad::uniform(wx, wy));
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
        analysis.step().expect("yielding member should converge");
    }
    let reported = if fixed_both { n / 2 } else { n };
    let d = &analysis.domain().node(nodes[reported]).displacement;
    [d[0], d[1], d[2]]
}

/// References are OpenSeesPy 3.8 `forceBeamColumn` (Lobatto, `eleLoad
/// -beamUniform`) with the identical model.
#[test]
fn yielding_force_members_under_udl_match_opensees() {
    // (elements, points, wx, wy, clamped both ends, OpenSees [ux, uy, rz])
    let cases: [(usize, usize, f64, f64, bool, [f64; 3]); 6] = [
        (
            1,
            3,
            0.0,
            -5.5,
            false,
            [0.0, -0.4435178662762192, -0.0055838837588357065],
        ),
        (
            1,
            5,
            0.0,
            -5.5,
            false,
            [0.0, -0.37428343005830444, -0.004891539396656561],
        ),
        (
            4,
            5,
            0.0,
            -5.5,
            false,
            [0.0, -0.36382010866362074, -0.004792302190234896],
        ),
        (
            1,
            5,
            1.0,
            -5.5,
            false,
            [
                0.0009150165016501648,
                -0.37428343005830456,
                -0.004891539396656563,
            ],
        ),
        (2, 5, 0.0, -30.0, true, [0.0, -0.04009379574907932, 0.0]),
        (4, 5, 0.0, -30.0, true, [0.0, -0.03964791993503361, 0.0]),
    ];
    for (n, points, wx, wy, fixed_both, reference) in cases {
        let got = yielding_member(n, points, wx, wy, fixed_both);
        for c in 0..3 {
            let tolerance = 1e-6 * reference.iter().map(|v| v.abs()).fold(0.0, f64::max);
            assert!(
                (got[c] - reference[c]).abs() < tolerance,
                "{n} elements, {points} points, wx {wx}, wy {wy}, clamped {fixed_both}, dof {c}: got {}, OpenSees {}",
                got[c],
                reference[c]
            );
        }
    }
}

/// Gravity UDL ramped on a yielding cantilever, frozen with
/// `hold_pattern_constant`, then a tip load pushed under displacement control —
/// the two-phase shape pysees models take. The frozen element load must stay in
/// the elements' section forces through the second phase's fresh `Analysis`
/// (whose pseudo-time restarts at 0) and the committed state carried across.
/// References: OpenSeesPy `forceBeamColumn`, `loadConst`, `DisplacementControl`.
#[test]
fn frozen_udl_then_displacement_controlled_tip_load_matches_opensees() {
    let n = 4;
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
    let gravity = domain.default_pattern();
    for k in 0..n {
        let beam = domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
            nodes[k],
            nodes[k + 1],
            fibers(),
            BeamIntegration::Lobatto { points: 5 },
        )));
        domain.add_element_load(gravity, beam, ElementLoad::uniform(0.0, -4.5));
    }
    let newton = || Algorithm::Newton {
        tangent: TangentStrategy::Current,
        line_search: None,
    };
    let test = || ConvergenceTest::NormUnbalance {
        tol: 1e-8,
        max_iter: 50,
    };

    let mut phase1 = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.2 })
        .algorithm(newton())
        .test(test())
        .build(domain);
    let mut factor = 0.0;
    for _ in 0..5 {
        factor = phase1.step().expect("gravity step").load_factor;
    }
    let tip = *nodes.last().unwrap();
    let after_gravity = phase1.domain().node(tip).displacement;
    assert!((after_gravity[1] + 0.2823002392678316).abs() < 1e-6);
    assert!((after_gravity[2] + 0.0037628520167384666).abs() < 1e-6);

    phase1.domain_mut().hold_pattern_constant(gravity, factor);
    let mut domain = phase1.into_domain();
    let lateral = domain.add_load_pattern(LoadSeries::Linear { slope: 1.0 });
    domain.add_nodal_load(lateral, tip, 1, -1.0);

    let mut phase2 = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::DisplacementControl {
            node: tip,
            dof: 1,
            increment: -0.05,
        })
        .algorithm(newton())
        .test(test())
        .build(domain);
    let mut load_factor = 0.0;
    for _ in 0..6 {
        load_factor = phase2.step().expect("push step").load_factor;
    }
    let d = phase2.domain().node(tip).displacement;
    let close = |got: f64, want: f64| (got - want).abs() < 1e-6 * want.abs();
    assert!(close(d[1], -0.5823002392678316), "uy {}", d[1]);
    assert!(close(d[2], -0.007520147627162879), "rz {}", d[2]);
    assert!(
        close(load_factor, 84.82012353714734),
        "load factor {load_factor}"
    );
}
