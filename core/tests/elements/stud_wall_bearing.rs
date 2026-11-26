//! Bearing base (`StudWallLayout::bearing`, `docs/stud-wall-bearing.md`): a wall standing on base nodes through
//! internal `kv` / `kh` springs at every bottom-plate position and pier-end hold-downs, interpolated linearly
//! between base positions.
//!
//! - On two fixed ground nodes it equals `foundation: true` with external pier-end links (rigid-support stand-in
//!   `kv = 1e12`): 1e-6.
//! - On the floor-zone top row of the wall below it equals the station + `ZeroLength` link model with the same
//!   springs: 1e-8.
//! - Offset studs: each base node's reaction is the lever-rule share of the spring forces, coincident base
//!   nodes share equally, and the reactions sum to the load: 1e-9.
//! - O86-crush hold-downs cycled between bearing and uplift (lateral load reversed under constant gravity) follow
//!   the same link model step by step: 1e-6. The element's local solve crosses the hold-down kinks without cycling.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    AnchorageBearing, ChordDistribution, Domain, Element, ElementId, ElementLoad, FloorZone,
    LoadSeries, Material, Node, NodeId, StudWall, StudWallLayout, WallBearing, WallCaseLoad,
    WallChord, WallLoad, WallOpening, WallShearLaw, WallSheathing, WallStation, WallStud,
    ZeroLength,
};

type D2 = Domain<2, 3, NodeId, Element>;
const EA: f64 = 9500. * 38. * 140.;
const H: f64 = 2440.;
const DEPTH: f64 = 254.;
const KH: f64 = 1e9;
const KV: f64 = 1e6;
/// Stand-in for the foundation's rigid internal supports.
const RIGID: f64 = 1e12;

fn near(a: f64, b: f64, tol: f64, scale: f64, what: &str) {
    assert!(
        (a - b).abs() <= tol * scale.abs().max(1e-30),
        "{what}: {a} vs {b} (scale {scale})"
    );
}
fn anchor() -> Material {
    Material::anchorage(15_000., 0., AnchorageBearing::Linear { k: 1e6 }).unwrap()
}
fn xs(length: f64) -> Vec<f64> {
    let mut v: Vec<f64> = (0..)
        .map(|i| i as f64 * 400.)
        .take_while(|x| *x < length)
        .collect();
    v.push(length);
    v
}
fn layout(length: f64, window: bool, floor: bool) -> StudWallLayout {
    StudWallLayout {
        length,
        height: H,
        studs: xs(length)
            .into_iter()
            .map(|x| WallStud {
                x,
                ea: if x == 0. || x == length { 3. * EA } else { EA },
            })
            .collect(),
        top_ea: 2. * EA,
        top_ei: Some(2. * 9500. * 140. * 38f64.powi(3) / 12.),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: EA,
        bottom_ei: 9500. * 140. * 38f64.powi(3) / 12.,
        openings: if window {
            vec![WallOpening {
                x0: 1600.,
                x1: 2400.,
                head_height: 2000.,
                header_ei: 1e11,
                jack_ea: 2. * EA,
                door: false,
            }]
        } else {
            vec![]
        },
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 6000. },
        }),
        foundation: false,
        bearing: None,
        floor: floor.then_some(FloorZone {
            depth: DEPTH,
            ga: 1e7,
            rim_ea: 9500. * 38. * 235.,
        }),
    }
}
/// Gravity `wy` and lateral `wx` along the top, as one case.
fn top_line(length: f64, wx: f64, wy: f64) -> Vec<Vec<WallLoad>> {
    vec![vec![WallLoad::TopLine {
        x0: 0.,
        x1: length,
        wx0: wx,
        wx1: wx,
        wy0: wy,
        wy1: wy,
    }]]
}
fn load_wall(d: &mut D2, wall: ElementId) {
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_element_load(
        p,
        wall,
        ElementLoad {
            wall: Some(WallCaseLoad {
                case: 0,
                factor: 1.,
            }),
            ..Default::default()
        },
    );
}
fn solve(d: D2) -> carapace_core::analysis::Analysis {
    solve_steps(d, 0.25, 4)
}
fn solve_steps(d: D2, increment: f64, steps: usize) -> carapace_core::analysis::Analysis {
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-10,
            max_iter: 100,
        })
        .build(d);
    for _ in 0..steps {
        a.step().expect("bearing wall step");
    }
    a
}

/// A 3.6 m window wall (pier ends 0, 1600, 2400, 3600) under gravity and a lateral load that lifts the left
/// end: bearing on two fixed ground nodes against the foundation with external links at the same four pier ends.
#[test]
fn bearing_on_ground_matches_foundation_with_links() {
    let length = 3600.;
    let loads = || top_line(length, 6., -4.);
    let ends = [0., 1600., 2400., length];
    let run = |bearing: bool| {
        let mut d = D2::new();
        let c: Vec<NodeId> = [(0., 0.), (length, 0.), (length, H), (0., H)]
            .iter()
            .map(|&(x, y)| d.add_node(Node::new([x, y])))
            .collect();
        let mut l = layout(length, true, false);
        let mut stations = Vec::new();
        if bearing {
            l.bearing = Some(WallBearing {
                kv: RIGID,
                kh: RIGID,
                hold_down: Some(anchor()),
            });
            for x in [0., length] {
                let n = d.add_node(Node::new([x, 0.]).fix(0).fix(1).fix(2));
                stations.push(WallStation {
                    chord: WallChord::Base,
                    x,
                    node: n,
                });
            }
        } else {
            l.foundation = true;
            let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
            for &x in &ends {
                let n = if x == 0. {
                    c[0]
                } else if x == length {
                    c[1]
                } else {
                    let n = d.add_node(Node::new([x, 0.]));
                    stations.push(WallStation {
                        chord: WallChord::Bottom,
                        x,
                        node: n,
                    });
                    n
                };
                d.add_element(Element::ZeroLength(
                    ZeroLength::new(ground, n)
                        .with_material(0, Material::Elastic { e: RIGID })
                        .with_material(1, anchor()),
                ));
            }
        }
        let wall = d.add_element(Element::StudWall(
            StudWall::new([c[0], c[1], c[2], c[3]], stations, l, loads()).unwrap(),
        ));
        load_wall(&mut d, wall);
        let a = solve(d);
        let r = a.domain().element_stud_wall_response(wall, 1.).unwrap();
        let top = [c[2], c[3]].map(|n| a.domain().node(n).displacement);
        (r, top)
    };
    let (b, tb) = run(true);
    let (f, tf) = run(false);
    let scale = tf[0][0].abs();
    for (p, q) in tb.iter().zip(&tf) {
        near(p[0], q[0], 1e-6, scale, "top ux");
        near(p[1], q[1], 1e-6, scale, "top uy");
    }
    for (p, q) in b.piers.iter().zip(&f.piers) {
        near(p.q, q.q, 1e-6, q.q, "pier q");
    }
    let fmax = f.base_nodes.iter().map(|n| n.fy.abs()).fold(0., f64::max);
    assert_eq!(b.base_nodes.len(), f.base_nodes.len());
    for (p, q) in b.base_nodes.iter().zip(&f.base_nodes) {
        near(p.x, q.x, 1e-12, 1., "base x");
        near(p.fy, q.fy, 1e-6, fmax, "base fy");
        near(p.fx, q.fx, 1e-6, fmax, "base fx");
    }
    // The hold-downs are the pier-end entries; the left end is in tension.
    assert_eq!(b.hold_downs.iter().map(|h| h.x).collect::<Vec<_>>(), ends);
    assert!(b.hold_downs[0].force > 0. && b.hold_downs[0].deformation > 0.);
    near(
        b.hold_downs[0].force,
        15_000. * b.hold_downs[0].deformation,
        1e-9,
        b.hold_downs[0].force,
        "hold-down law",
    );
    near(b.base[1], 4. * length, 1e-9, 4. * length, "gravity sum");
}

/// The upper wall of a two-storey stack on the floor-zone top row of the wall below, studs aligned: bearing on
/// the lower top nodes against bottom stations at every stud joined by `ZeroLength`s with the same springs.
#[test]
fn aligned_stack_matches_station_links() {
    let length = 3600.;
    let (b, tb) = stack(true, anchor(), StackLoad::Proportional);
    let (s, ts) = stack(false, anchor(), StackLoad::Proportional);
    let scale = ts[0][0].abs();
    for (p, q) in tb.iter().zip(&ts) {
        near(p[0], q[0], 1e-8, scale, "top ux");
        near(p[1], q[1], 1e-8, scale, "top uy");
    }
    for (p, q) in b.piers.iter().zip(&s.piers) {
        near(p.q, q.q, 1e-8, q.q, "pier q");
        near(p.chords[0], q.chords[0], 1e-8, q.chords[0], "pier chord");
    }
    for (p, q) in b.studs.iter().zip(&s.studs) {
        near(p.bottom, q.bottom, 1e-8, 6. * length, "stud force");
    }
    near(b.base[1], 6. * length, 1e-9, 6. * length, "gravity sum");
}

/// A 2x6 S-P-F window wall on the wall below with O86-crush hold-downs (tension 15 kN/mm, end posts Qr = 84.3 kN,
/// studs at the crush curve's initial 47.4 kN/mm), constant gravity and a lateral load cycled +1, -1, +1: every
/// pier end crosses between bearing and uplift. The bearing wall follows the equivalent link model at every step
/// (1e-6 of the top drift). Before the local solve's line search, full Newton steps cycled across the hold-down
/// kinks and the step failed.
#[test]
fn crush_hold_downs_cycle_between_bearing_and_uplift() {
    let hold =
        Material::anchorage(15_000., 0., AnchorageBearing::O86Crush { qr: 84_323. }).unwrap();
    let (b, tb) = stack(true, hold.clone(), StackLoad::Cyclic);
    let (s, ts) = stack(false, hold, StackLoad::Cyclic);
    let scale = ts.iter().map(|p| p[0].abs()).fold(0., f64::max);
    for (p, q) in tb.iter().zip(&ts) {
        near(p[0], q[0], 1e-6, scale, "top ux");
        near(p[1], q[1], 1e-6, scale, "top uy");
    }
    for (p, q) in b.piers.iter().zip(&s.piers) {
        near(p.chords[0], q.chords[0], 1e-6, q.chords[0], "pier chord");
    }
}

#[derive(Clone, Copy, PartialEq)]
enum StackLoad {
    /// Gravity and lateral together, in four steps.
    Proportional,
    /// Gravity constant; lateral (12 N/mm) along a path 0, 1, -1, 1 in 24 steps.
    Cyclic,
}

/// The aligned stack: a 3.6 m wall with a floor zone on the ground (foundation), and a window wall on its top row,
/// either bearing on it (`bearing`) or on bottom stations linked to it at every stud (`hold` at pier ends, `KV`
/// elsewhere). Returns the upper wall's response and its top corners' displacements at the end.
fn stack(
    bearing: bool,
    hold: Material,
    load: StackLoad,
) -> (carapace_core::model::StudWallResponse, [[f64; 3]; 2]) {
    let length = 3600.;
    let anchor = || hold.clone();
    let run = |bearing: bool| {
        let mut d = D2::new();
        let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
        let yt = H + DEPTH;
        let lc: Vec<NodeId> = [(0., 0.), (length, 0.)]
            .iter()
            .map(|&(x, y)| d.add_node(Node::new([x, y])))
            .collect();
        for &n in &lc {
            d.add_element(Element::ZeroLength(
                ZeroLength::new(ground, n)
                    .with_material(0, Material::Elastic { e: 1e10 })
                    .with_material(1, Material::Elastic { e: 1e10 }),
            ));
        }
        let lower_top: Vec<(f64, NodeId)> = xs(length)
            .into_iter()
            .map(|x| (x, d.add_node(Node::new([x, yt]))))
            .collect();
        let mut ll = layout(length, false, true);
        ll.foundation = true;
        let inner = |v: &[(f64, NodeId)], chord| {
            v.iter()
                .filter(|(x, _)| *x > 0. && *x < length)
                .map(|&(x, node)| WallStation { chord, x, node })
                .collect::<Vec<_>>()
        };
        d.add_element(Element::StudWall(
            StudWall::new(
                [lc[0], lc[1], lower_top.last().unwrap().1, lower_top[0].1],
                inner(&lower_top, WallChord::Top),
                ll,
                vec![],
            )
            .unwrap(),
        ));
        let uc: Vec<NodeId> = [(0., yt), (length, yt), (length, yt + H), (0., yt + H)]
            .iter()
            .map(|&(x, y)| d.add_node(Node::new([x, y])))
            .collect();
        let mut ul = layout(length, true, false);
        let mut stations = Vec::new();
        if bearing {
            ul.bearing = Some(WallBearing {
                kv: KV,
                kh: KH,
                hold_down: Some(anchor()),
            });
            stations.extend(lower_top.iter().map(|&(x, node)| WallStation {
                chord: WallChord::Base,
                x,
                node,
            }));
        } else {
            // Bottom stations at every bearing position (studs and jacks), each linked to the node below.
            let ends = ul.pier_ends();
            for &(x, below) in lower_top
                .iter()
                .filter(|(x, _)| !(*x > 1600. && *x < 2400.))
            {
                let above = if x == 0. {
                    uc[0]
                } else if x == length {
                    uc[1]
                } else {
                    let n = d.add_node(Node::new([x, yt]));
                    stations.push(WallStation {
                        chord: WallChord::Bottom,
                        x,
                        node: n,
                    });
                    n
                };
                let v = if ends.contains(&x) {
                    anchor()
                } else {
                    Material::Elastic { e: KV }
                };
                d.add_element(Element::ZeroLength(
                    ZeroLength::new(below, above)
                        .with_material(0, Material::Elastic { e: KH })
                        .with_material(1, v),
                ));
            }
        }
        let cases = match load {
            StackLoad::Proportional => top_line(length, 5., -6.),
            StackLoad::Cyclic => [top_line(length, 0., -6.), top_line(length, 12., 0.)].concat(),
        };
        let wall = d.add_element(Element::StudWall(
            StudWall::new([uc[0], uc[1], uc[2], uc[3]], stations, ul, cases).unwrap(),
        ));
        let a = match load {
            StackLoad::Proportional => {
                load_wall(&mut d, wall);
                solve(d)
            }
            StackLoad::Cyclic => {
                for (case, series) in [
                    (0, LoadSeries::Constant),
                    (
                        1,
                        LoadSeries::Path {
                            times: vec![0., 1., 3., 5.],
                            factors: vec![0., 1., -1., 1.],
                        },
                    ),
                ] {
                    let p = d.add_load_pattern(series);
                    d.add_element_load(
                        p,
                        wall,
                        ElementLoad {
                            wall: Some(WallCaseLoad { case, factor: 1. }),
                            ..Default::default()
                        },
                    );
                }
                solve_steps(d, 5. / 24., 24)
            }
        };
        let r = a.domain().element_stud_wall_response(wall, 1.).unwrap();
        (r, [uc[2], uc[3]].map(|n| a.domain().node(n).displacement))
    };
    run(bearing)
}

/// Studs at 400 on base nodes offset by 150 (and two nodes at one position): each base node's reaction is the
/// lever-rule share of the spring forces above it, and the reactions sum to the load.
#[test]
fn offset_studs_split_by_lever_rule() {
    let length = 2400.;
    let mut d = D2::new();
    let c: Vec<NodeId> = [(0., 0.), (length, 0.), (length, H), (0., H)]
        .iter()
        .enumerate()
        .map(|(i, &(x, y))| {
            // An unsheathed racking mode does not arise (linear sheathing), but fix one top ux against drift anyway.
            d.add_node(if i == 3 {
                Node::new([x, y]).fix(0)
            } else {
                Node::new([x, y])
            })
        })
        .collect();
    let mut bx: Vec<f64> = (0..8).map(|i| -150. + 400. * i as f64).collect();
    bx.push(1050.); // a second node at 1050
    let base: Vec<(f64, NodeId)> = bx
        .iter()
        .map(|&x| (x, d.add_node(Node::new([x, 0.]).fix(0).fix(1).fix(2))))
        .collect();
    let mut l = layout(length, false, false);
    l.top_distribution = ChordDistribution::Direct;
    l.top_ei = None;
    l.bearing = Some(WallBearing {
        kv: KV,
        kh: KH,
        hold_down: None,
    });
    let stations = base
        .iter()
        .map(|&(x, node)| WallStation {
            chord: WallChord::Base,
            x,
            node,
        })
        .collect();
    let wall = d.add_element(Element::StudWall(
        StudWall::new(
            [c[0], c[1], c[2], c[3]],
            stations,
            l,
            top_line(length, 0., -5.),
        )
        .unwrap(),
    ));
    load_wall(&mut d, wall);
    let a = solve(d);
    let dom = a.domain();
    let r = dom.element_stud_wall_response(wall, 1.).unwrap();
    assert_eq!(r.base_nodes.len(), 7);
    let total = 5. * length;
    near(r.base[1], total, 1e-9, total, "spring sum");
    // Expected reaction of each base position: every spring's force split by the lever rule.
    let pos: Vec<f64> = (0..8).map(|i| -150. + 400. * i as f64).collect();
    let mut expect = vec![0.; pos.len()];
    for n in &r.base_nodes {
        let i = pos
            .windows(2)
            .position(|w| w[0] <= n.x && n.x <= w[1])
            .unwrap();
        let t = (n.x - pos[i]) / (pos[i + 1] - pos[i]);
        expect[i] += n.fy * (1. - t);
        expect[i + 1] += n.fy * t;
    }
    let reaction = |node: NodeId| dom.reaction(node, 1, 1.);
    let mut sum = 0.;
    for (i, &x) in pos.iter().enumerate() {
        let nodes: Vec<NodeId> = base
            .iter()
            .filter(|(bx, _)| *bx == x)
            .map(|b| b.1)
            .collect();
        let got: f64 = nodes.iter().map(|&n| reaction(n)).sum();
        near(got, expect[i], 1e-9, total, "lever-rule reaction");
        if nodes.len() == 2 {
            near(
                reaction(nodes[0]),
                reaction(nodes[1]),
                1e-12,
                total,
                "coincident nodes share",
            );
        }
        sum += got;
    }
    near(sum, total, 1e-9, total, "reactions sum to the load");
}

#[test]
fn bearing_errors_are_named() {
    let mut d = D2::new();
    let ids: Vec<NodeId> = (0..6)
        .map(|i| d.add_node(Node::new([i as f64, 0.])))
        .collect();
    let n = |i: usize| ids[i];
    let base = |x: f64, i: usize| WallStation {
        chord: WallChord::Base,
        x,
        node: n(i),
    };
    let corners = [n(0), n(1), n(2), n(3)];
    let mut l = layout(2400., false, false);
    l.bearing = Some(WallBearing {
        kv: KV,
        kh: KH,
        hold_down: None,
    });
    let e = StudWall::new(
        corners,
        vec![base(100., 4), base(2400., 5)],
        l.clone(),
        vec![],
    )
    .unwrap_err();
    assert!(e.0.contains("bearing position x = 0 lies outside"), "{e}");
    let mut f = l.clone();
    f.foundation = true;
    let e = StudWall::new(corners, vec![base(0., 4), base(2400., 5)], f, vec![]).unwrap_err();
    assert!(e.0.contains("exclusive"), "{e}");
    let mut g = l.clone();
    g.bearing = None;
    let e = StudWall::new(corners, vec![base(0., 4)], g, vec![]).unwrap_err();
    assert!(e.0.contains("base station without bearing"), "{e}");
}

/// 3D: the wall stiffens only its own plane, so a bearing wall ties its bottom corners to the base along its normal (`kh`), and
/// only those. A wall along Y on three fixed base nodes (y = -1000, 2000, 5000) pushed out of plane by a prescribed `ux` linear in
/// y: its bottom corners follow the interpolated base exactly (nothing else resists), its top corners' `ux` is no equation.
#[test]
fn bearing_ties_bottom_corners_out_of_plane() {
    use carapace_core::model::{Domain3, Element3, Node3};
    let mut d = Domain3::new();
    let fixed = |p: [f64; 3]| (0..6).fold(Node3::new(p), |n, k| n.fix(k));
    let base: Vec<_> = [-1000., 2000., 5000.]
        .map(|y| d.add_node(fixed([0., y, 0.])))
        .to_vec();
    let c = [[0., 0., 0.], [0., 4000., 0.], [0., 4000., H], [0., 0., H]]
        .map(|p| d.add_node(Node3::new(p)));
    let stations = base
        .iter()
        .zip([-1000., 2000., 5000.])
        .map(|(&node, x)| WallStation {
            chord: WallChord::Base,
            x,
            node,
        })
        .collect();
    let l = StudWallLayout {
        bearing: Some(WallBearing {
            kv: KV,
            kh: KH,
            hold_down: None,
        }),
        ..layout(4000., false, false)
    };
    d.add_element(Element3::StudWall(
        StudWall::new3(c, stations, l, vec![]).unwrap(),
    ));
    d.validate().unwrap();
    assert!(
        d.is_active(c[0], 0) && d.is_active(c[1], 0),
        "bottom corners are tied out of plane"
    );
    assert!(
        !d.is_active(c[2], 0) && !d.is_active(c[3], 0),
        "top corners are not"
    );
    for (&b, y) in base.iter().zip([-1000., 2000., 5000.]) {
        assert!(d.prescribe_displacement(b, 0, (y + 1000.) / 1000.));
    }
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
    near(
        a.domain().node(c[0]).displacement[0],
        1.,
        1e-12,
        1.,
        "c0 ux",
    );
    near(
        a.domain().node(c[1]).displacement[0],
        5.,
        1e-12,
        1.,
        "c1 ux",
    );
    near(
        a.domain().node(c[2]).displacement[1],
        0.,
        1e-12,
        1.,
        "no in-plane movement",
    );
}

/// 3D: the tie only joins an out-of-plane equation something else makes. A wall along X on a wall along X with nothing holding
/// either out of plane (no diaphragm) solves under gravity as before (no mechanism): the lower top corner's `uy` stays no equation
/// and the upper bottom corner, tied to it, stays at zero.
#[test]
fn bearing_tie_makes_no_mechanism_out_of_plane() {
    use carapace_core::model::{Domain3, Element3, ElementLoad3, Node3};
    let mut d = Domain3::new();
    let fixed = |p: [f64; 3]| (0..6).fold(Node3::new(p), |n, k| n.fix(k));
    let ground = [0., 4000.].map(|x| d.add_node(fixed([x, 0., 0.])));
    let lo = [[0., 0., 0.], [4000., 0., 0.], [4000., 0., H], [0., 0., H]]
        .map(|p| d.add_node(Node3::new(p)));
    let up = [
        [0., 0., H],
        [4000., 0., H],
        [4000., 0., 2. * H],
        [0., 0., 2. * H],
    ]
    .map(|p| d.add_node(Node3::new(p)));
    let bearing = Some(WallBearing {
        kv: KV,
        kh: KH,
        hold_down: None,
    });
    let base = |nodes: [carapace_core::model::Node3Id; 2]| {
        nodes
            .iter()
            .zip([0., 4000.])
            .map(|(&node, x)| WallStation {
                chord: WallChord::Base,
                x,
                node,
            })
            .collect::<Vec<_>>()
    };
    let lower = StudWall::new3(
        lo,
        base(ground),
        StudWallLayout {
            bearing: bearing.clone(),
            ..layout(4000., false, false)
        },
        vec![],
    )
    .unwrap();
    let upper = StudWall::new3(
        up,
        base([lo[3], lo[2]]),
        StudWallLayout {
            bearing,
            ..layout(4000., false, false)
        },
        top_line(4000., 0., -5.),
    )
    .unwrap();
    d.add_element(Element3::StudWall(lower));
    let wall = d.add_element(Element3::StudWall(upper));
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_element_load(
        p,
        wall,
        ElementLoad3 {
            wall: Some(WallCaseLoad {
                case: 0,
                factor: 1.,
            }),
            ..Default::default()
        },
    );
    d.validate().unwrap();
    assert!(
        !d.is_active(lo[3], 1),
        "nothing holds the lower top corner out of plane"
    );
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-6,
            max_iter: 4,
        })
        .build(d);
    a.step().expect("no out-of-plane mechanism");
    assert_eq!(a.domain().node(up[0]).displacement[1], 0.);
    assert!(
        a.domain().node(up[2]).displacement[2] < 0.,
        "the load goes down"
    );
}
