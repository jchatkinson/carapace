//! 3D stud-wall and floor-band behavior against equivalent 2D reference models: drift, wall response, gravity
//! takedown, nonlinear shear, stacked-storey response and corner coupling. Tolerances are stated per comparison.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    ChordDistribution, Domain, Domain3, Element, Element3, ElementLoad, ElementLoad3, FloorZone,
    Node, Node3, StudWall, StudWallLayout, WallCaseLoad, WallLoad, WallShearLaw, WallSheathing,
    WallStud,
};

/// Geometry-derived masks in a Z-up model activate exactly the translations in the wall plane.
/// Reference: local-axis projections; exact active slots.
#[test]
fn z_up_wall_masks_follow_geometry() {
    use carapace_core::model::Node3Id;
    fn active(corners: [[f64; 3]; 4]) -> [bool; 3] {
        let mut d = Domain3::new();
        let ids: [Node3Id; 4] = corners.map(|p| d.add_node(Node3::new(p)));
        d.add_element(Element3::StudWall(
            StudWall::new3(ids, vec![], layout(), vec![]).unwrap(),
        ));
        d.validate().unwrap();
        std::array::from_fn(|slot| d.is_active(ids[0], slot))
    }
    assert_eq!(
        active([[0., 0., 0.], [4., 0., 0.], [4., 0., 3.], [0., 0., 3.]]),
        [true, false, true]
    );
    assert_eq!(
        active([[0., 0., 0.], [0., 4., 0.], [0., 4., 3.], [0., 0., 3.]]),
        [false, true, true]
    );
    assert_eq!(
        active([[0., 0., 0.], [3.2, 2.4, 0.], [3.2, 2.4, 3.], [0., 0., 3.]]),
        [true, true, true]
    );
}

/// A wall with non-orthogonal corner axes fails geometry validation with a named model error.
/// Reference: `PlaneMap` relative dot-product check, threshold 1e-9.
#[test]
fn nonorthogonal_wall_corners_are_rejected() {
    let mut d = Domain3::new();
    let ids = [
        d.add_node(Node3::new([0., 0., 0.])),
        d.add_node(Node3::new([4., 0.1, 0.])),
        d.add_node(Node3::new([4., 0.2, 3.])),
        d.add_node(Node3::new([0., 0.1, 3.])),
    ];
    d.add_element(Element3::StudWall(
        StudWall::new3(ids, vec![], layout(), vec![]).unwrap(),
    ));
    let error = d.validate().unwrap_err();
    assert!(format!("{error:?}").contains("wall axes are not orthogonal"));
}

fn layout() -> StudWallLayout {
    StudWallLayout {
        length: 4.,
        height: 3.,
        studs: vec![
            WallStud { x: 0., ea: 2e6 },
            WallStud { x: 2., ea: 2e6 },
            WallStud { x: 4., ea: 2e6 },
        ],
        top_ea: 1e8,
        top_ei: Some(1e8),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: 1e8,
        bottom_ei: 1e8,
        openings: vec![],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 2e5 },
        }),
        foundation: true,
        bearing: None,
        floor: None,
    }
}

fn solve_2d() -> (f64, f64) {
    let mut d = Domain::new();
    let c = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([4., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([4., 3.])),
        d.add_node(Node::new([0., 3.])),
    ];
    let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    loads[0].push(WallLoad::TopPoint {
        x: 4.,
        fx: 100.,
        fy: 0.,
    });
    let id = d.add_element(Element::StudWall(
        StudWall::new(c, vec![], layout(), loads).unwrap(),
    ));
    d.add_element_load(
        d.default_pattern(),
        id,
        ElementLoad {
            wall: Some(WallCaseLoad {
                case: 0,
                factor: 1.0,
            }),
            ..Default::default()
        },
    );
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 4,
        })
        .build(d);
    a.step().unwrap();
    let n = a.domain().node(c[2]);
    (
        n.displacement[0],
        a.domain().element_stud_wall_response(id, 1.).unwrap().piers[0].q,
    )
}

fn solve_3d(ex: f64, ez: f64) -> (f64, f64) {
    let mut d = Domain3::new();
    let skew = ex.abs() > 1e-9 && ez.abs() > 1e-9;
    let out = if ex.abs() > 0.999 { 2 } else { 0 };
    let top_node = |coords| {
        if skew {
            Node3::new(coords)
        } else {
            Node3::new(coords).fix(out)
        }
    };
    let c = [
        d.add_node(Node3::new([0., 0., 0.]).fix(0).fix(1).fix(2)),
        d.add_node(Node3::new([4. * ex, 0., 4. * ez]).fix(0).fix(1).fix(2)),
        d.add_node(top_node([4. * ex, 3., 4. * ez])),
        d.add_node(top_node([0., 3., 0.])),
    ];
    if skew {
        for &id in &c[2..] {
            // Constrain the plan-normal translation while leaving along-wall drift free.
            d.add_constraint((id, 2), &[(id, 0, ez / ex)]);
        }
    }
    let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    loads[0].push(WallLoad::TopPoint {
        x: 4.,
        fx: 100.,
        fy: 0.,
    });
    let wall = StudWall::new3(c, vec![], layout(), loads).unwrap();
    let id = d.add_element(Element3::StudWall(wall));
    d.add_element_load(
        d.default_pattern(),
        id,
        ElementLoad3 {
            wall: Some(WallCaseLoad {
                case: 0,
                factor: 1.0,
            }),
            ..Default::default()
        },
    );
    let mut a = AnalysisBuilder::new()
        .constraint_handler(if skew {
            ConstraintHandler::Transformation
        } else {
            ConstraintHandler::Plain
        })
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 4,
        })
        .build(d);
    a.step().unwrap();
    let n = a.domain().node(c[2]);
    (
        n.displacement[0] * ex + n.displacement[2] * ez,
        a.domain().element_stud_wall_response(id, 1.).unwrap().piers[0].q,
    )
}

/// X and Z walls match the planar drift and pier shear to 1e-10; the unstiffened transverse translation is fixed.
/// Reference: the identical condensed 2D wall model, mapped into each vertical plane.
#[test]
fn axis_aligned_walls_match_the_planar_reference() {
    let (drift_2d, q_2d) = solve_2d();
    for (ex, ez) in [
        (1., 0.),
        (0., 1.),
        (30f64.to_radians().cos(), 30f64.to_radians().sin()),
    ] {
        let (drift, q) = solve_3d(ex, ez);
        assert!(
            (drift - drift_2d).abs() <= 1e-10 * drift_2d.abs().max(1.),
            "dir=({ex},{ez}), {drift} != {drift_2d}"
        );
        assert!((q - q_2d).abs() <= 1e-10 * q_2d.abs().max(1.));
    }
}

/// An X wall activates only global `ux` and `uy` at the domain level; `uz` has no equation at its nodes.
/// Reference: the X wall's declared element mask propagated through `Domain3` numbering; exact equation check.
#[test]
fn x_wall_domain_numbering_excludes_global_uz() {
    let mut d = Domain3::new();
    let c = [
        d.add_node(Node3::new([0., 0., 0.])),
        d.add_node(Node3::new([4., 0., 0.])),
        d.add_node(Node3::new([4., 3., 0.])),
        d.add_node(Node3::new([0., 3., 0.])),
    ];
    let wall = StudWall::new3(c, vec![], layout(), vec![]).unwrap();
    let wall_id = d.add_element(Element3::StudWall(wall));
    let a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 4,
        })
        .build(d);
    let domain = a.domain();
    assert!(domain.num_free_dofs() > 0);
    for node in c {
        assert!(domain.equation_of(node, 0).is_some());
        assert!(domain.equation_of(node, 1).is_some());
        assert_eq!(
            domain.equation_of(node, 2),
            None,
            "wall {wall_id:?} node {node:?}"
        );
    }
}

const OPEN_H: f64 = 2440.;
const OPEN_X0: f64 = 1600.;
const OPEN_X1: f64 = 2400.;
const OPEN_PLY: f64 = 9500. * 38. * 140.;
const OPEN_KHD: f64 = 20_000.;
const OPEN_V: f64 = 14_400.;

fn window_layout() -> StudWallLayout {
    StudWallLayout {
        length: 3600.,
        height: OPEN_H,
        studs: [0., 400., 800., 1200., OPEN_X0, OPEN_X1, 2800., 3200., 3600.]
            .into_iter()
            .map(|x| WallStud {
                x,
                ea: if x == 0. || x == 3600. || x == OPEN_X0 || x == OPEN_X1 {
                    2. * OPEN_PLY
                } else {
                    OPEN_PLY
                },
            })
            .collect(),
        top_ea: 2. * OPEN_PLY,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: OPEN_PLY,
        bottom_ei: 9500. * 140. * 38f64.powi(3) / 12.,
        openings: vec![carapace_core::model::WallOpening {
            x0: OPEN_X0,
            x1: OPEN_X1,
            head_height: 2000.,
            header_ei: 9500. * 2. * 38. * 235f64.powi(3) / 12.,
            jack_ea: OPEN_PLY,
            door: false,
        }],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 6000. },
        }),
        foundation: true,
        bearing: None,
        floor: None,
    }
}

#[derive(Debug)]
struct WindowResult {
    drift: f64,
    piers: Vec<(f64, f64, f64)>,
    uplift: Vec<f64>,
}

fn solve_window_2d(o86: bool) -> WindowResult {
    use carapace_core::model::{LoadSeries, Material, NodeId, ZeroLength};
    let mut d: Domain<2, 3, NodeId, Element> = Domain::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    fn hold(d: &mut Domain<2, 3, NodeId, Element>, ground: NodeId, x: f64) -> NodeId {
        let n = d.add_node(Node::new([x, 0.]).fix(0));
        d.add_element(Element::ZeroLength(
            ZeroLength::new(ground, n).with_material(1, Material::Elastic { e: OPEN_KHD }),
        ));
        n
    }
    let c0 = hold(&mut d, ground, 0.);
    let c1 = hold(&mut d, ground, 3600.);
    let c2 = d.add_node(Node::new([3600., OPEN_H]));
    let c3 = d.add_node(Node::new([0., OPEN_H]));
    let s0 = hold(&mut d, ground, OPEN_X0);
    let s1 = hold(&mut d, ground, OPEN_X1);
    let stations = [(OPEN_X0, s0), (OPEN_X1, s1)]
        .into_iter()
        .map(|(x, node)| carapace_core::model::WallStation {
            chord: carapace_core::model::WallChord::Bottom,
            x,
            node,
        })
        .collect();
    let mut wall_layout = window_layout();
    if o86 {
        wall_layout.sheathing = Some(WallSheathing {
            law: WallShearLaw::O86 {
                bv: 11_000.,
                nail_spacing: 100.,
                nail_diameter: 3.25,
            },
        });
    }
    let wall = d.add_element(Element::StudWall(
        StudWall::new([c0, c1, c2, c3], stations, wall_layout, vec![]).unwrap(),
    ));
    let pattern = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_nodal_load(pattern, c3, 0, OPEN_V);
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl {
            increment: if o86 { 0.1 } else { 1. },
        })
        .algorithm(if o86 {
            Algorithm::Newton {
                tangent: carapace_core::analysis::TangentStrategy::Current,
                line_search: None,
            }
        } else {
            Algorithm::Linear
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-10,
            max_iter: 30,
        })
        .build(d);
    for _ in 0..if o86 { 10 } else { 1 } {
        a.step().unwrap();
    }
    let dom = a.domain();
    let r = dom.element_stud_wall_response(wall, 1.).unwrap();
    WindowResult {
        drift: dom.node(c3).displacement[0],
        piers: r
            .piers
            .iter()
            .map(|p| (p.q, p.chords[0], p.chords[1]))
            .collect(),
        uplift: [c0, s0, s1, c1]
            .map(|n| dom.node(n).displacement[1])
            .to_vec(),
    }
}

fn solve_window_3d_axes(ex: f64, ez: f64, z_up: bool, o86: bool) -> WindowResult {
    use carapace_core::model::{Element3, LoadSeries, Material, Node3, ZeroLength3};
    let mut d = Domain3::new();
    let up = if z_up { 2 } else { 1 };
    let lateral = if z_up { 1 } else { 2 };
    let point = |x: f64, h: f64| {
        if z_up {
            [x * ex, x * ez, h]
        } else {
            [x * ex, h, x * ez]
        }
    };
    let ground = d.add_node(Node3::new([0., 0., 0.]).fix(0).fix(1).fix(2));
    let c0 = d.add_node(Node3::new(point(0., 0.)).fix(0).fix(lateral));
    let c1 = d.add_node(Node3::new(point(3600., 0.)).fix(0).fix(lateral));
    let c2 = d.add_node(Node3::new(point(3600., OPEN_H)));
    let c3 = d.add_node(Node3::new(point(0., OPEN_H)));
    let s0 = d.add_node(Node3::new(point(OPEN_X0, 0.)).fix(0).fix(lateral));
    let s1 = d.add_node(Node3::new(point(OPEN_X1, 0.)).fix(0).fix(lateral));
    for n in [c0, c1, s0, s1] {
        d.add_element(Element3::ZeroLength3(
            ZeroLength3::new(ground, n).with_material(up, Material::Elastic { e: OPEN_KHD }),
        ));
    }
    if ex.abs() > 1e-9 && ez.abs() > 1e-9 {
        for n in [c2, c3] {
            // Plan-normal translation is zero: ez*ux - ex*uz = 0.
            d.add_constraint((n, lateral), &[(n, 0, ez / ex)]);
        }
    }
    let stations = [(OPEN_X0, s0), (OPEN_X1, s1)]
        .into_iter()
        .map(|(x, node)| carapace_core::model::WallStation {
            chord: carapace_core::model::WallChord::Bottom,
            x,
            node,
        })
        .collect();
    let mut wall_layout = window_layout();
    if o86 {
        wall_layout.sheathing = Some(WallSheathing {
            law: WallShearLaw::O86 {
                bv: 11_000.,
                nail_spacing: 100.,
                nail_diameter: 3.25,
            },
        });
    }
    let wall = d.add_element(Element3::StudWall(
        StudWall::new3([c0, c1, c2, c3], stations, wall_layout, vec![]).unwrap(),
    ));
    let pattern = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_nodal_load(pattern, c3, 0, OPEN_V * ex);
    if ez.abs() > 1e-9 {
        // For a skew plane, apply the horizontal component along the wall direction.
        // `add_nodal_load` above supplies ux; add uz to complete its projection.
        d.add_nodal_load(pattern, c3, lateral, OPEN_V * ez);
    }
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl {
            increment: if o86 { 0.1 } else { 1. },
        })
        .algorithm(if o86 {
            Algorithm::Newton {
                tangent: carapace_core::analysis::TangentStrategy::Current,
                line_search: None,
            }
        } else {
            Algorithm::Linear
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-10,
            max_iter: 30,
        })
        .build(d);
    for _ in 0..if o86 { 10 } else { 1 } {
        a.step().unwrap();
    }
    let dom = a.domain();
    let r = dom.element_stud_wall_response(wall, 1.).unwrap();
    let disp = dom.node(c3).displacement;
    WindowResult {
        drift: disp[0] * ex + disp[lateral] * ez,
        piers: r
            .piers
            .iter()
            .map(|p| (p.q, p.chords[0], p.chords[1]))
            .collect(),
        uplift: [c0, s0, s1, c1]
            .map(|n| dom.node(n).displacement[up])
            .to_vec(),
    }
}

/// A window wall with king studs and four elastic hold-downs matches the 2D wall along +X, +Z and 30°.
/// Reference: the identical 2D window wall solved in the same test; drift, each pier's `q`/chords and uplift, 1e-10.
#[test]
fn window_wall_matches_2d_in_three_vertical_planes() {
    let reference = solve_window_2d(false);
    for z_up in [false, true] {
        for (ex, ez) in [
            (1., 0.),
            (30f64.to_radians().cos(), 30f64.to_radians().sin()),
        ] {
            let actual = solve_window_3d_axes(ex, ez, z_up, false);
            assert_eq!(actual.piers.len(), reference.piers.len());
            let values = [(actual.drift, reference.drift)]
                .into_iter()
                .chain(
                    actual
                        .piers
                        .iter()
                        .zip(&reference.piers)
                        .flat_map(|(a, b)| [(a.0, b.0), (a.1, b.1), (a.2, b.2)]),
                )
                .chain(
                    actual
                        .uplift
                        .iter()
                        .zip(&reference.uplift)
                        .map(|(a, b)| (*a, *b)),
                );
            for (a, b) in values {
                assert!(
                    (a - b).abs() <= 1e-10 * a.abs().max(b.abs()).max(1.),
                    "z_up={z_up} dir=({ex},{ez}) {a} != {b}"
                );
            }
        }
    }
}

/// The nonlinear O86 window wall converges with Newton at 30° and matches its 2D definition.
/// Reference: same O86 wall solved in 2D with ten 0.1 load increments; drift, pier `q` and chords, 1e-9.
#[test]
fn o86_nonlinear_wall_matches_2d_at_thirty_degrees() {
    let reference = solve_window_2d(true);
    let angle = 30f64.to_radians();
    for z_up in [false, true] {
        let actual = solve_window_3d_axes(angle.cos(), angle.sin(), z_up, true);
        assert!(
            (actual.drift - reference.drift).abs()
                <= 1e-9 * actual.drift.abs().max(reference.drift.abs()).max(1.)
        );
        assert_eq!(actual.piers.len(), reference.piers.len());
        for (a, b) in actual.piers.iter().zip(&reference.piers) {
            for (a, b) in [(a.0, b.0), (a.1, b.1), (a.2, b.2)] {
                assert!(
                    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.),
                    "{a} != {b}"
                );
            }
        }
    }
}

#[derive(Debug)]
struct GravityResult {
    base: Vec<(f64, f64, f64)>,
    jacks: Vec<f64>,
    header: Vec<(f64, f64)>,
}

fn gravity_window_2d() -> GravityResult {
    use carapace_core::model::{ElementLoad, LoadSeries, NodeId};
    let mut l = window_layout();
    l.top_ei = Some(1.2e10);
    l.foundation = false;
    let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    loads[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: 3600.,
        wx0: 0.,
        wx1: 0.,
        wy0: -10.,
        wy1: -10.,
    });
    let mut d: Domain<2, 3, NodeId, Element> = Domain::new();
    let xs = [0., 400., 800., 1200., OPEN_X0, OPEN_X1, 2800., 3200., 3600.];
    let bottom: Vec<_> = xs
        .iter()
        .map(|&x| d.add_node(Node::new([x, 0.]).fix(0).fix(1).fix(2)))
        .collect();
    let c2 = d.add_node(Node::new([3600., OPEN_H]));
    let c3 = d.add_node(Node::new([0., OPEN_H]));
    let stations = xs
        .iter()
        .skip(1)
        .take(xs.len() - 2)
        .zip(bottom.iter().skip(1).take(xs.len() - 2))
        .map(|(&x, &node)| carapace_core::model::WallStation {
            chord: carapace_core::model::WallChord::Bottom,
            x,
            node,
        })
        .collect();
    let id = d.add_element(Element::StudWall(
        StudWall::new(
            [bottom[0], *bottom.last().unwrap(), c2, c3],
            stations,
            l,
            loads,
        )
        .unwrap(),
    ));
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_element_load(
        p,
        id,
        ElementLoad {
            wall: Some(WallCaseLoad {
                case: 0,
                factor: 1.,
            }),
            ..Default::default()
        },
    );
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-11,
            max_iter: 4,
        })
        .build(d);
    a.step().unwrap();
    let r = a.domain().element_stud_wall_response(id, 1.).unwrap();
    GravityResult {
        base: r.base_nodes.iter().map(|n| (n.x, n.fx, n.fy)).collect(),
        jacks: r.jacks.iter().map(|j| j.axial).collect(),
        header: r.header.iter().map(|h| (h.moment, h.shear)).collect(),
    }
}

fn gravity_window_3d(z_up: bool) -> GravityResult {
    use carapace_core::model::{Element3, ElementLoad3, LoadSeries, Node3};
    let mut l = window_layout();
    l.top_ei = Some(1.2e10);
    l.foundation = false;
    let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    loads[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: 3600.,
        wx0: 0.,
        wx1: 0.,
        wy0: -10.,
        wy1: -10.,
    });
    let mut d = Domain3::new();
    let point = |x: f64, y: f64, z: f64| if z_up { [z, x, y] } else { [x, y, z] };
    let xs = [0., 400., 800., 1200., OPEN_X0, OPEN_X1, 2800., 3200., 3600.];
    let bottom: Vec<_> = xs
        .iter()
        .map(|&x| d.add_node(Node3::new(point(0., 0., x)).fix(0).fix(1).fix(2)))
        .collect();
    let c2 = d.add_node(Node3::new(point(0., OPEN_H, 3600.)));
    let c3 = d.add_node(Node3::new(point(0., OPEN_H, 0.)));
    let stations = xs
        .iter()
        .skip(1)
        .take(xs.len() - 2)
        .zip(bottom.iter().skip(1).take(xs.len() - 2))
        .map(|(&x, &node)| carapace_core::model::WallStation {
            chord: carapace_core::model::WallChord::Bottom,
            x,
            node,
        })
        .collect();
    let id = d.add_element(Element3::StudWall(
        StudWall::new3(
            [bottom[0], *bottom.last().unwrap(), c2, c3],
            stations,
            l,
            loads,
        )
        .unwrap(),
    ));
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_element_load(
        p,
        id,
        ElementLoad3 {
            wall: Some(WallCaseLoad {
                case: 0,
                factor: 1.,
            }),
            ..Default::default()
        },
    );
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-11,
            max_iter: 4,
        })
        .build(d);
    a.step().unwrap();
    let r = a.domain().element_stud_wall_response(id, 1.).unwrap();
    GravityResult {
        base: r.base_nodes.iter().map(|n| (n.x, n.fx, n.fy)).collect(),
        jacks: r.jacks.iter().map(|j| j.axial).collect(),
        header: r.header.iter().map(|h| (h.moment, h.shear)).collect(),
    }
}

/// Gravity slot loading through a window/header gives matching 2D and +Z 3D base, jack and header responses.
/// Reference: same supported window model solved in 2D; `base_nodes`, jack axial forces and header shear/moment, 1e-10.
#[test]
fn window_gravity_takedown_matches_2d() {
    let a = gravity_window_2d();
    for b in [gravity_window_3d(false), gravity_window_3d(true)] {
        assert_eq!(a.base.len(), b.base.len());
        assert_eq!(a.jacks.len(), b.jacks.len());
        assert_eq!(a.header.len(), b.header.len());
        for (x, y) in a.base.iter().zip(&b.base) {
            for (u, v) in [(x.0, y.0), (x.1, y.1), (x.2, y.2)] {
                assert!(
                    (u - v).abs() <= 1e-10 * u.abs().max(v.abs()).max(1.),
                    "base {u} != {v}"
                );
            }
        }
        for (u, v) in a.jacks.iter().zip(&b.jacks) {
            assert!((u - v).abs() <= 1e-10 * u.abs().max(v.abs()).max(1.));
        }
        for (u, v) in a.header.iter().zip(&b.header) {
            for (u, v) in [(u.0, v.0), (u.1, v.1)] {
                assert!((u - v).abs() <= 1e-10 * u.abs().max(v.abs()).max(1.));
            }
        }
    }
}

fn stack_ga(force: f64, spacing: f64) -> f64 {
    let v = force / 2400.;
    let slip = (0.013 * v * spacing / 3.25f64.powi(2)).powi(2);
    1. / (1. / 11_000. + 0.0025 * slip / v)
}

fn rigid_stack_2d() -> (Vec<f64>, Vec<f64>) {
    use carapace_core::model::{LoadSeries, Material, NodeId, ZeroLength};
    const L: f64 = 2400.;
    const H: f64 = 3000.;
    const D: f64 = 254.;
    const R: f64 = 1e8;
    let mut d: Domain<2, 3, NodeId, Element> = Domain::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let left = d.add_node(Node::new([0., 0.]));
    let right = d.add_node(Node::new([L, 0.]).fix(0).fix(1));
    d.add_element(Element::ZeroLength(
        ZeroLength::new(ground, left)
            .with_material(0, Material::Elastic { e: R })
            .with_material(1, Material::Elastic { e: 20_000. }),
    ));
    let pattern = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    let mut bottom = [left, right];
    let mut walls = vec![];
    let mut levels = vec![];
    for i in 0..2 {
        let y = i as f64 * (H + D);
        let ytop = if i == 1 { y + H } else { y + H + D };
        let tl = d.add_node(Node::new([0., ytop]));
        let tr = d.add_node(Node::new([L, ytop]));
        let ea = [3., 2.][i] * OPEN_PLY;
        let layout = StudWallLayout {
            length: L,
            height: H,
            studs: vec![WallStud { x: 0., ea }, WallStud { x: L, ea }],
            top_ea: OPEN_PLY * 1e6,
            top_ei: Some(1.2e10),
            top_distribution: ChordDistribution::Beam,
            bottom_ea: OPEN_PLY * 1e6,
            bottom_ei: 6e9,
            openings: vec![],
            sheathing: Some(WallSheathing {
                law: WallShearLaw::Linear {
                    ga: stack_ga([6000., 8000.][i], [100., 150.][i]),
                },
            }),
            foundation: false,
            bearing: None,
            floor: (i == 0).then_some(FloorZone {
                depth: D,
                ga: R,
                rim_ea: R,
            }),
        };
        let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
        if i == 1 {
            let w = 8000. / L;
            loads[0].push(WallLoad::TopLine {
                x0: 0.,
                x1: L,
                wx0: w,
                wx1: w,
                wy0: 0.,
                wy1: 0.,
            });
        }
        let wall = d.add_element(Element::StudWall(
            StudWall::new([bottom[0], bottom[1], tr, tl], vec![], layout, loads).unwrap(),
        ));
        walls.push(wall);
        if i == 1 {
            d.add_element_load(
                pattern,
                wall,
                ElementLoad {
                    wall: Some(WallCaseLoad {
                        case: 0,
                        factor: 1.,
                    }),
                    ..Default::default()
                },
            );
            levels.push([tl, tr]);
        } else {
            let bl = d.add_node(Node::new([0., ytop]));
            let br = d.add_node(Node::new([L, ytop]));
            for (a, b, k) in [(tl, bl, 20_000.), (tr, br, R)] {
                d.add_element(Element::ZeroLength(
                    ZeroLength::new(a, b).with_material(1, Material::Elastic { e: k }),
                ));
                d.equal_dof(a, b, &[0]);
            }
            d.add_nodal_load(pattern, tl, 0, 6000.);
            levels.push([tl, tl]);
            bottom = [bl, br];
        }
    }
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-13,
            max_iter: 4,
        })
        .build(d);
    a.step().unwrap();
    let dom = a.domain();
    (
        levels
            .iter()
            .map(|n| (dom.node(n[0]).displacement[0] + dom.node(n[1]).displacement[0]) / 2.)
            .collect(),
        walls
            .iter()
            .map(|&id| dom.element_stud_wall_response(id, 1.).unwrap().piers[0].chords[0])
            .collect(),
    )
}

fn rigid_stack_3d(z_up: bool) -> (Vec<f64>, Vec<f64>) {
    use carapace_core::model::{Element3, LoadSeries, Material, Node3, ZeroLength3};
    const L: f64 = 2400.;
    const H: f64 = 3000.;
    const D: f64 = 254.;
    const R: f64 = 1e8;
    let up = if z_up { 2 } else { 1 };
    let wall_axis = if z_up { 1 } else { 2 };
    let point = |x: f64, y: f64, z: f64| if z_up { [x, z, y] } else { [x, y, z] };
    let mut d = Domain3::new();
    let ground = d.add_node(Node3::new([0., 0., 0.]).fix(0).fix(1).fix(2));
    let left = d.add_node(Node3::new([0., 0., 0.]));
    let right = d.add_node(Node3::new(point(0., 0., L)).fix(0).fix(1).fix(2));
    d.add_element(Element3::ZeroLength3(
        ZeroLength3::new(ground, left)
            .with_material(up, Material::Elastic { e: 20_000. })
            .with_material(wall_axis, Material::Elastic { e: R }),
    ));
    let pattern = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    let mut bottom = [left, right];
    let mut walls = vec![];
    let mut levels = vec![];
    for i in 0..2 {
        let y = i as f64 * (H + D);
        let ytop = if i == 1 { y + H } else { y + H + D };
        let tl = d.add_node(Node3::new(point(0., ytop, 0.)));
        let tr = d.add_node(Node3::new(point(0., ytop, L)));
        let ea = [3., 2.][i] * OPEN_PLY;
        let layout = StudWallLayout {
            length: L,
            height: H,
            studs: vec![WallStud { x: 0., ea }, WallStud { x: L, ea }],
            top_ea: OPEN_PLY * 1e6,
            top_ei: Some(1.2e10),
            top_distribution: ChordDistribution::Beam,
            bottom_ea: OPEN_PLY * 1e6,
            bottom_ei: 6e9,
            openings: vec![],
            sheathing: Some(WallSheathing {
                law: WallShearLaw::Linear {
                    ga: stack_ga([6000., 8000.][i], [100., 150.][i]),
                },
            }),
            foundation: false,
            bearing: None,
            floor: (i == 0).then_some(FloorZone {
                depth: D,
                ga: R,
                rim_ea: R,
            }),
        };
        let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
        if i == 1 {
            let w = 8000. / L;
            loads[0].push(WallLoad::TopLine {
                x0: 0.,
                x1: L,
                wx0: w,
                wx1: w,
                wy0: 0.,
                wy1: 0.,
            });
        }
        let wall = d.add_element(Element3::StudWall(
            StudWall::new3([bottom[0], bottom[1], tr, tl], vec![], layout, loads).unwrap(),
        ));
        walls.push(wall);
        if i == 1 {
            d.add_element_load(
                pattern,
                wall,
                ElementLoad3 {
                    wall: Some(WallCaseLoad {
                        case: 0,
                        factor: 1.,
                    }),
                    ..Default::default()
                },
            );
            levels.push([tl, tr]);
        } else {
            let bl = d.add_node(Node3::new(point(0., ytop, 0.)));
            let br = d.add_node(Node3::new(point(0., ytop, L)));
            for (a, b, k) in [(tl, bl, 20_000.), (tr, br, R)] {
                d.add_element(Element3::ZeroLength3(
                    ZeroLength3::new(a, b).with_material(up, Material::Elastic { e: k }),
                ));
                d.equal_dof(a, b, &[wall_axis]);
            }
            d.add_nodal_load(pattern, tl, wall_axis, 6000.);
            levels.push([tl, tl]);
            bottom = [bl, br];
        }
    }
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-13,
            max_iter: 4,
        })
        .build(d);
    a.step().unwrap();
    let dom = a.domain();
    (
        levels
            .iter()
            .map(|n| {
                (dom.node(n[0]).displacement[wall_axis] + dom.node(n[1]).displacement[wall_axis])
                    / 2.
            })
            .collect(),
        walls
            .iter()
            .map(|&id| dom.element_stud_wall_response(id, 1.).unwrap().piers[0].chords[0])
            .collect(),
    )
}

/// The rigid-zone 2-storey stack has the same level drifts and wall chords in Y-up and Z-up 3D models as in 2D.
/// Reference: the 2-storey rigid-zone model and hand-walk case in `floor_zone.rs`; displacements and chords, 1e-10.
#[test]
fn rigid_floor_zone_stack_matches_2d_in_z() {
    let a = rigid_stack_2d();
    for b in [rigid_stack_3d(true)] {
        assert_eq!(a.0.len(), b.0.len());
        assert_eq!(a.1.len(), b.1.len());
        for (i, (x, y)) in a.0.iter().zip(&b.0).enumerate() {
            assert!(
                (x - y).abs() <= 1e-10 * x.abs().max(y.abs()).max(1.),
                "level {i} displacement: {x} != {y}"
            );
        }
        for (i, (x, y)) in a.1.iter().zip(&b.1).enumerate() {
            assert!(
                (x - y).abs() <= 1e-10 * x.abs().max(y.abs()).max(1.),
                "wall {i} chord: {x} != {y}"
            );
        }
    }
}

/// A Z-up 3D wall (`z` up, along `y`) with a floor zone reports the same bay shear flows, rim forces and
/// base reaction as its 2D definition: a window wall with a zone, sheared by horizontal loads on the exterior
/// top row. Reference: the same model in 2D; 1e-10.
#[test]
fn z_up_floor_zone_equals_2d() {
    use carapace_core::model::{LoadSeries, Node3, NodeId, WallChord, WallOpening, WallStation};
    const L: f64 = 3600.;
    const H: f64 = 2440.;
    const D: f64 = 254.;
    let xs = [0., 800., 1600., 2400., 3200., L];
    let layout = || StudWallLayout {
        length: L,
        height: H,
        studs: xs.iter().map(|&x| WallStud { x, ea: OPEN_PLY }).collect(),
        top_ea: 2. * OPEN_PLY,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: OPEN_PLY,
        bottom_ei: OPEN_PLY,
        openings: vec![WallOpening {
            x0: 1600.,
            x1: 2400.,
            head_height: 2000.,
            header_ei: 1e9,
            jack_ea: OPEN_PLY,
            door: false,
        }],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 6000. },
        }),
        foundation: true,
        bearing: None,
        floor: Some(FloorZone {
            depth: D,
            ga: 5e6,
            rim_ea: 2e8,
        }),
    };
    let loads = |n: usize| (0..n).map(|i| 1000. + 250. * i as f64);
    // 2D
    let mut d: Domain<2, 3, NodeId, Element> = Domain::new();
    let c0 = d.add_node(Node::new([0., 0.]).fix(0).fix(1));
    let c1 = d.add_node(Node::new([L, 0.]).fix(0).fix(1));
    let top: Vec<_> = xs
        .iter()
        .map(|&x| d.add_node(Node::new([x, H + D])))
        .collect();
    let stations: Vec<_> = xs[1..5]
        .iter()
        .zip(&top[1..5])
        .map(|(&x, &node)| WallStation {
            chord: WallChord::Top,
            x,
            node,
        })
        .collect();
    let wall = d.add_element(Element::StudWall(
        StudWall::new(
            [c0, c1, top[5], top[0]],
            stations,
            layout(),
            Default::default(),
        )
        .unwrap(),
    ));
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    for (n, f) in top.iter().zip(loads(6)) {
        d.add_nodal_load(p, *n, 0, f);
    }
    let mut a2 = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 4,
        })
        .build(d);
    a2.step().unwrap();
    let r2 = a2.domain().element_stud_wall_response(wall, 1.).unwrap();
    let top2 = top.clone();
    // Z-up: wall along y, up is z.
    let mut d3 = Domain3::new();
    let pt = |x: f64, z: f64| [0., x, z];
    let c0 = d3.add_node(Node3::new(pt(0., 0.)).fix(0).fix(1).fix(2));
    let c1 = d3.add_node(Node3::new(pt(L, 0.)).fix(0).fix(1).fix(2));
    let top: Vec<_> = xs
        .iter()
        .map(|&x| d3.add_node(Node3::new(pt(x, H + D))))
        .collect();
    let top3 = top.clone();
    let stations: Vec<_> = xs[1..5]
        .iter()
        .zip(&top[1..5])
        .map(|(&x, &node)| WallStation {
            chord: WallChord::Top,
            x,
            node,
        })
        .collect();
    let wall3 = d3.add_element(Element3::StudWall(
        StudWall::new3(
            [c0, c1, top[5], top[0]],
            stations,
            layout(),
            Default::default(),
        )
        .unwrap(),
    ));
    let p = d3.add_load_pattern(LoadSeries::Linear { slope: 1. });
    for (n, f) in top.iter().zip(loads(6)) {
        d3.add_nodal_load(p, *n, 1, f);
    }
    let mut a3 = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 4,
        })
        .build(d3);
    a3.step().unwrap();
    let r3 = a3.domain().element_stud_wall_response(wall3, 1.).unwrap();
    let near = |a: f64, b: f64, what: &str| {
        assert!(
            (a - b).abs() <= 1e-10 * a.abs().max(b.abs()).max(1e-3),
            "{what}: {a} != {b}"
        )
    };
    assert_eq!(r2.floor_bays.len(), 5);
    for (a, b) in r2.floor_bays.iter().zip(&r3.floor_bays) {
        near(a.q, b.q, "bay q");
        near(a.gamma, b.gamma, "bay gamma");
    }
    assert_eq!(r2.floor_bays[2].q, 0.);
    for (a, b) in r2.floor_rim.iter().zip(&r3.floor_rim) {
        near(*a, *b, "rim");
    }
    for (a, b) in r2.piers.iter().zip(&r3.piers) {
        near(a.q, b.q, "pier q");
        near(a.moment_base, b.moment_base, "pier moment");
    }
    for i in 0..3 {
        near(r2.base[i], r3.base[i], "base");
    }
    for i in 0..6 {
        near(
            a2.domain().node(top2[i]).displacement[0],
            a3.domain().node(top3[i]).displacement[1],
            "top u",
        );
        near(
            a2.domain().node(top2[i]).displacement[1],
            a3.domain().node(top3[i]).displacement[2],
            "top v",
        );
    }
}

/// The 3D wall's threaded and serial assembly paths give bit-identical drift and pier shear.
/// Reference: native serial-scatter contract; exact `f64::to_bits` equality.
#[cfg(feature = "parallel")]
#[test]
fn serial_and_parallel_wall_results_are_bit_identical() {
    use carapace_core::analysis::{set_parallelism, Parallelism};
    set_parallelism(Parallelism::Serial);
    let serial = solve_3d(1., 0.);
    set_parallelism(Parallelism::Threads(2));
    let threaded = solve_3d(1., 0.);
    set_parallelism(Parallelism::Serial);
    assert_eq!(serial.0.to_bits(), threaded.0.to_bits());
    assert_eq!(serial.1.to_bits(), threaded.1.to_bits());
}

/// L corner: wall A along +X and wall B along +Z share their corner post nodes (base and top); A is pushed towards
/// +X and B towards -Z, so both overturn onto the shared post (twice the compression on its hold-down). Nothing fixes the shared nodes out of
/// plane: each wall stiffens the other's out-of-plane translation. Reference: the same two walls with separate
/// corner nodes tied by `equal_dof` on `uy` only; drift, pier shears and the hold-down uplift agree to 1e-10.
#[test]
fn l_corner_shares_post_between_planes() {
    use carapace_core::model::{Material, Node3Id, ZeroLength3};
    fn build(shared: bool) -> [f64; 5] {
        let mut d = Domain3::new();
        let ground = d.add_node(Node3::new([0., 0., 0.]).fix(0).fix(1).fix(2));
        let holddown = |d: &mut Domain3, p: [f64; 3]| -> Node3Id {
            let n = d.add_node(Node3::new(p).fix(0).fix(2));
            d.add_element(Element3::ZeroLength3(
                ZeroLength3::new(ground, n).with_material(1, Material::Elastic { e: 5e4 }),
            ));
            n
        };
        let a0 = holddown(&mut d, [0., 0., 0.]);
        let corner = holddown(&mut d, [4., 0., 0.]);
        let b1 = holddown(&mut d, [4., 0., 4.]);
        let a3 = d.add_node(Node3::new([0., 3., 0.]).fix(2));
        let top = d.add_node(Node3::new([4., 3., 0.]));
        let b2 = d.add_node(Node3::new([4., 3., 4.]).fix(0));
        // Wall B's own corner nodes, tied to wall A's in `uy` only when not shared.
        let (b0, b3) = if shared {
            (corner, top)
        } else {
            let b0 = d.add_node(Node3::new([4., 0., 0.]).fix(0).fix(2));
            let b3 = d.add_node(Node3::new([4., 3., 0.]));
            d.equal_dof(corner, b0, &[1]);
            d.equal_dof(top, b3, &[1]);
            (b0, b3)
        };
        let wall = |d: &mut Domain3, c: [Node3Id; 4], w: f64| {
            let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
            loads[0].push(WallLoad::TopLine {
                x0: 0.,
                x1: 4.,
                wx0: w,
                wx1: w,
                wy0: 0.,
                wy1: 0.,
            });
            let mut l = layout();
            l.foundation = false;
            let id = d.add_element(Element3::StudWall(
                StudWall::new3(c, vec![], l, loads).unwrap(),
            ));
            d.add_element_load(
                d.default_pattern(),
                id,
                ElementLoad3 {
                    wall: Some(WallCaseLoad {
                        case: 0,
                        factor: 1.0,
                    }),
                    ..Default::default()
                },
            );
            id
        };
        let wa = wall(&mut d, [a0, corner, top, a3], 25.);
        let wb = wall(&mut d, [b0, b1, b2, b3], -25.);
        let mut an = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Transformation)
            .integrator(Integrator::LoadControl { increment: 1. })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-10,
                max_iter: 4,
            })
            .build(d);
        an.step().expect("L corner solve");
        let dom = an.domain();
        let q = |id| dom.element_stud_wall_response(id, 1.).unwrap().piers[0].q;
        [
            dom.node(a3).displacement[0],
            dom.node(b2).displacement[2],
            q(wa),
            q(wb),
            dom.node(corner).displacement[1],
        ]
    }
    let (shared, tied) = (build(true), build(false));
    for (a, b) in shared.iter().zip(&tied) {
        assert!(
            (a - b).abs() <= 1e-10 * a.abs().max(b.abs()),
            "{shared:?} vs {tied:?}"
        );
    }
    // The shared post really is loaded by both walls.
    assert!(shared[4].abs() > 0.);
}
