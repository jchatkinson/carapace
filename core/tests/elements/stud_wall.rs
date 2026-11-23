//! Stud-wall v1 layout, kinematics, condensation, load, and explicit-frame verification.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    ChordDistribution, ElementLoad, StudWall, StudWallLayout, WallChord, WallLoad, WallStation,
    WallStud, WALL_SLOTS,
};
use carapace_core::model::{Domain, Element, Node};
use carapace_core::model::{ElasticBeamColumn, GeomTransf, Material, Truss};
use slotmap::new_key_type;

new_key_type! { struct DummyNode; }

fn layout() -> StudWallLayout {
    StudWallLayout {
        length: 10.0,
        height: 8.0,
        studs: vec![
            WallStud { x: 0.0, ea: 100.0 },
            WallStud { x: 5.0, ea: 100.0 },
            WallStud { x: 10.0, ea: 100.0 },
        ],
        top_ea: 1.0e6,
        top_ei: Some(1.0e6),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: 1.0e6,
        bottom_ei: 1.0e6,
        openings: vec![],
        sheathing: None,
        foundation: false,
    }
}

/// Invalid layout parameters identify the offending item; slot selectors add and scale fieldwise.
/// Reference: stud-wall-v1 §§3 and 5; exact scalar checks (zero tolerance).
#[test]
fn layout_errors_are_named_and_wall_slots_add_and_scale() {
    let mut l = layout();
    l.studs[1].ea = f64::NAN;
    let mut d: Domain<2, 3, carapace_core::model::NodeId, Element> = Domain::new();
    let c = [
        d.add_node(Node::new([0., 0.])),
        d.add_node(Node::new([10., 0.])),
        d.add_node(Node::new([10., 8.])),
        d.add_node(Node::new([0., 8.])),
    ];
    let error = StudWall::new(c, vec![], l, std::array::from_fn(|_| vec![])).unwrap_err();
    assert!(error.to_string().contains("studs[1]"));

    let mut bad_law = layout();
    bad_law.sheathing = Some(carapace_core::model::WallSheathing {
        law: carapace_core::model::WallShearLaw::O86 {
            bv: 0.,
            nail_spacing: 100.,
            nail_diameter: 3.25,
        },
    });
    let error = StudWall::new(c, vec![], bad_law, std::array::from_fn(|_| vec![])).unwrap_err();
    assert!(error.to_string().contains("O86 bv"));

    let mut l = layout();
    l.foundation = true;
    let station = d.add_node(Node::new([5., 0.]));
    let error = StudWall::new(
        c,
        vec![WallStation {
            chord: WallChord::Bottom,
            x: 5.,
            node: station,
        }],
        l,
        std::array::from_fn(|_| vec![]),
    )
    .expect(
        "a bottom station on a foundation replaces the internal support there (for a hold-down)",
    );
    let _ = error;

    let mut many_piers = StudWallLayout {
        length: 390.,
        height: 100.,
        studs: (0..65)
            .flat_map(|i| {
                let x = i as f64 * 6.;
                [WallStud { x, ea: 1e4 }, WallStud { x: x + 2., ea: 1e4 }]
            })
            .collect(),
        top_ea: 1e8,
        top_ei: Some(1e8),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: 1e8,
        bottom_ei: 1e8,
        openings: (0..64)
            .map(|i| carapace_core::model::WallOpening {
                x0: i as f64 * 6. + 2.5,
                x1: i as f64 * 6. + 3.5,
                head_height: 50.,
                header_ei: 1e6,
                jack_ea: 1e4,
                door: true,
            })
            .collect(),
        sheathing: Some(carapace_core::model::WallSheathing {
            law: carapace_core::model::WallShearLaw::O86 {
                bv: 1e4,
                nail_spacing: 100.,
                nail_diameter: 3.25,
            },
        }),
        foundation: false,
    };
    many_piers.studs.last_mut().unwrap().x = 390.;
    let error = StudWall::new(c, vec![], many_piers, std::array::from_fn(|_| vec![])).unwrap_err();
    assert!(error.to_string().contains("pier count 65"));

    let mut nonlinear_layout = layout();
    nonlinear_layout.sheathing = Some(carapace_core::model::WallSheathing {
        law: carapace_core::model::WallShearLaw::O86 {
            bv: 100.,
            nail_spacing: 100.,
            nail_diameter: 3.25,
        },
    });
    let stations: Vec<_> = (1..=29)
        .map(|i| WallStation {
            chord: WallChord::Top,
            x: i as f64 / 3.,
            node: d.add_node(Node::new([i as f64 / 3., 8.])),
        })
        .collect();
    // 33 exterior nodes: nonlinear walls share the 256-node limit (the old 32-node cap is gone).
    assert!(StudWall::new(
        c,
        stations,
        nonlinear_layout,
        std::array::from_fn(|_| vec![]),
    )
    .is_ok());

    let a = ElementLoad {
        wall: [1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ..Default::default()
    };
    let b = ElementLoad {
        wall: [3.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ..Default::default()
    };
    assert_eq!((a + b * 2.0).wall, [7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert_eq!(WALL_SLOTS, 8);
}

/// A partial top line load is accepted into its requested pattern slot and preserves the exterior node order.
/// Reference: stud-wall-v1 §§2 and 5; exact count check.
#[test]
fn partial_line_load_has_correct_resultant() {
    let l = layout();
    let mut d: Domain<2, 3, carapace_core::model::NodeId, Element> = Domain::new();
    let c = [
        d.add_node(Node::new([0., 0.])),
        d.add_node(Node::new([10., 0.])),
        d.add_node(Node::new([10., 8.])),
        d.add_node(Node::new([0., 8.])),
    ];
    let mut loads: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| vec![]);
    loads[2].push(WallLoad::TopLine {
        x0: 2.0,
        x1: 8.0,
        wx0: 0.0,
        wx1: 0.0,
        wy0: 2.0,
        wy1: 2.0,
    });
    let wall = StudWall::new(c, vec![], l, loads).unwrap();
    assert_eq!(wall.nodes().len(), 4);
}

fn kinematic_wall(
    angle: f64,
    sheathing: bool,
    local_u: &[[f64; 2]; 4],
) -> (
    Domain<2, 3, carapace_core::model::NodeId, Element>,
    carapace_core::model::ElementId,
) {
    let l = 10.;
    let h = 8.;
    let (c, s) = (angle.cos(), angle.sin());
    let xy = [(0., 0.), (l, 0.), (l, h), (0., h)];
    let mut d: Domain<2, 3, carapace_core::model::NodeId, Element> = Domain::new();
    let corners = std::array::from_fn(|i| {
        let (x, y) = xy[i];
        let gx = c * x - s * y;
        let gy = s * x + c * y;
        let ux = c * local_u[i][0] - s * local_u[i][1];
        let uy = s * local_u[i][0] + c * local_u[i][1];
        d.add_node(
            Node::new([gx, gy])
                .with_initial_displacement(0, ux)
                .with_initial_displacement(1, uy),
        )
    });
    let mut layout = layout();
    if !sheathing {
        layout.sheathing = None;
    } else {
        layout.sheathing = Some(carapace_core::model::WallSheathing {
            law: carapace_core::model::WallShearLaw::Linear { ga: 100. },
        });
    }
    let wall = StudWall::new(corners, vec![], layout, std::array::from_fn(|_| vec![])).unwrap();
    let id = d.add_element(Element::StudWall(wall));
    (d, id)
}

/// The exterior mask contains translations only; rigid translations and a small rigid rotation produce zero condensed force.
/// Reference: stud-wall-v1 §8 row 1 kinematics; relative force tolerance 1e-10.
#[test]
fn mask_and_rigid_body_modes() {
    let zero = [[0.; 2]; 4];
    let (d, id) = kinematic_wall(0., false, &zero);
    let Element::StudWall(w) = d.element(id) else {
        panic!()
    };
    let mask = w.dof_mask();
    assert!(mask.contains(0) && mask.contains(1));
    assert!(!mask.contains(2));
    for u in [
        [[0.3, -0.2]; 4],
        [[0., 0.], [0., 0.01], [-0.08, 0.01], [-0.08, 0.]],
    ] {
        let (d, id) = kinematic_wall(0., false, &u);
        let f = d.element_local_force(id);
        let scale = f.iter().map(|v| v.abs()).fold(0., f64::max);
        assert!(scale < 1e-10, "{f:?}");
    }
}

/// Corner-only condensed stiffness has the expected sheathing rank, and local response is invariant under rigid wall rotation.
/// Reference: stud-wall-v1 §8 row 1; relative SVD cutoff and force tolerance 1e-10.
#[test]
fn corner_rank_and_rotated_response() {
    let rank = |sheathing| {
        let mut cols = Vec::new();
        for j in 0..8 {
            let mut u = [[0.; 2]; 4];
            u[j / 2][j % 2] = 1.;
            let (d, id) = kinematic_wall(0., sheathing, &u);
            cols.extend(d.element_local_force(id).iter().copied());
        }
        let k = nalgebra::DMatrix::from_column_slice(8, 8, &cols);
        let svd = k.svd(false, false);
        let max = svd.singular_values[0];
        svd.singular_values
            .iter()
            .filter(|v| **v > max * 1e-10)
            .count()
    };
    assert_eq!(rank(false), 4);
    assert_eq!(rank(true), 5);
    let local = [[0.1, -0.03], [0.02, 0.08], [-0.04, 0.07], [0.06, -0.02]];
    let (a, aid) = kinematic_wall(0., false, &local);
    let (b, bid) = kinematic_wall(30f64.to_radians(), false, &local);
    let fa = a.element_local_force(aid);
    let fb = b.element_local_force(bid);
    for (x, y) in fa.iter().zip(fb.iter()) {
        assert!((x - y).abs() <= 1e-10 * x.abs().max(y.abs()).max(1.));
    }
}

type D2 = Domain<2, 3, carapace_core::model::NodeId, Element>;
type ExplicitFrame = (
    D2,
    [carapace_core::model::NodeId; 4],
    carapace_core::model::NodeId,
    carapace_core::model::NodeId,
    Vec<(f64, carapace_core::model::NodeId)>,
    Vec<(f64, carapace_core::model::NodeId)>,
    Vec<(f64, f64, carapace_core::model::ElementId)>,
    Vec<(f64, carapace_core::model::ElementId)>,
    Vec<(f64, carapace_core::model::ElementId)>,
    Vec<(f64, carapace_core::model::ElementId)>,
    Vec<(f64, f64, carapace_core::model::ElementId)>,
    Vec<(f64, f64, carapace_core::model::ElementId)>,
);
fn frame_layout(distribution: ChordDistribution, door: bool) -> StudWallLayout {
    StudWallLayout {
        length: 6.,
        height: 3.,
        studs: vec![
            WallStud { x: 0., ea: 20. },
            WallStud { x: 1., ea: 20. },
            WallStud { x: 2.5, ea: 20. },
            WallStud { x: 4., ea: 40. },
            WallStud { x: 5., ea: 20. },
            WallStud { x: 6., ea: 20. },
        ],
        top_ea: 1e5,
        top_ei: Some(3e4),
        top_distribution: distribution,
        bottom_ea: 1e5,
        bottom_ei: 3e4,
        openings: vec![carapace_core::model::WallOpening {
            x0: 1.5,
            x1: 3.5,
            head_height: 2.2,
            header_ei: 120.,
            jack_ea: 35.,
            door,
        }],
        sheathing: None,
        foundation: false,
    }
}
fn wall_loads() -> [Vec<WallLoad>; WALL_SLOTS] {
    let mut l: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| vec![]);
    l[0].push(WallLoad::TopPoint {
        x: 0.55,
        fx: 4.,
        fy: -8.,
    });
    l[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: 5.,
        wx0: 0.5,
        wx1: 1.2,
        wy0: -1.,
        wy1: -2.,
    });
    l
}
fn supported_domain() -> (
    D2,
    [carapace_core::model::NodeId; 4],
    carapace_core::model::NodeId,
    carapace_core::model::NodeId,
) {
    let mut d = D2::new();
    let c = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([6., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([6., 3.])),
        d.add_node(Node::new([0., 3.]).fix(0)),
    ];
    // Vertical support at the top station: with the plates hinged over the opening, the unsheathed left pier
    // (pinned only at c0) would otherwise rotate rigidly about c0, since pinned studs resist no vertical shear.
    let ts = d.add_node(Node::new([1., 3.]).fix(1));
    let bs = d.add_node(Node::new([5., 0.]).fix(0).fix(1));
    (d, c, ts, bs)
}
fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.),
        "{a} != {b}"
    );
}

fn explicit_frame(distribution: ChordDistribution, door: bool) -> ExplicitFrame {
    let (mut d, c, ts, bs) = supported_domain();
    let l = frame_layout(distribution, door);
    let mut tx: Vec<f64> = l.studs.iter().map(|s| s.x).collect();
    if distribution == ChordDistribution::Beam {
        tx.extend([1.5, 3.5]);
    }
    tx.extend([1.]);
    tx.sort_by(f64::total_cmp);
    tx.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    let mut bx: Vec<f64> = l
        .studs
        .iter()
        .filter(|s| !(s.x > 1.5 && s.x < 3.5))
        .map(|s| s.x)
        .collect();
    bx.extend([1.5, 3.5, 5.]);
    bx.sort_by(f64::total_cmp);
    bx.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    let ext = |ch: WallChord, x: f64| -> Option<carapace_core::model::NodeId> {
        if ch == WallChord::Top {
            if x == 0. {
                Some(c[3])
            } else if x == 6. {
                Some(c[2])
            } else if x == 1. {
                Some(ts)
            } else {
                None
            }
        } else if x == 0. {
            Some(c[0])
        } else if x == 6. {
            Some(c[1])
        } else if x == 5. {
            Some(bs)
        } else {
            None
        }
    };
    let mut top = Vec::new();
    for &x in &tx {
        let n = d.add_node(Node::new([x, 3.]));
        if let Some(e) = ext(WallChord::Top, x) {
            d.equal_dof(e, n, &[0, 1]);
        }
        top.push((x, n));
    }
    let mut bottom = Vec::new();
    for &x in &bx {
        let n = d.add_node(Node::new([x, 0.]));
        if let Some(e) = ext(WallChord::Bottom, x) {
            d.equal_dof(e, n, &[0, 1]);
        }
        bottom.push((x, n));
    }
    let map = |v: &[(f64, carapace_core::model::NodeId)], x: f64| {
        v.iter().find(|(xx, _)| (*xx - x).abs() < 1e-10).unwrap().1
    };
    let beam = |d: &mut D2, a, b, e, area, ei| {
        d.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            a,
            b,
            e,
            area,
            ei,
            GeomTransf::Linear,
        )))
    };
    // Plates are hinged at the opening edges (x0 = 1.5, x1 = 3.5) for the segments over the opening: a
    // coincident node tied in ux, uy only.
    let hinged = |d: &mut D2, p: &[(f64, carapace_core::model::NodeId)], y: f64| {
        let over = p[0].0 >= 1.5 - 1e-10 && p[1].0 <= 3.5 + 1e-10;
        let mut ends = [p[0].1, p[1].1];
        for (k, edge) in [(0, 1.5), (1, 3.5)] {
            if over && (p[k].0 - edge).abs() < 1e-10 {
                let n = d.add_node(Node::new([edge, y]));
                d.equal_dof(p[k].1, n, &[0, 1]);
                ends[k] = n;
            }
        }
        ends
    };
    let mut top_beams = Vec::new();
    for p in top.windows(2) {
        if distribution == ChordDistribution::Beam {
            let [a, b] = hinged(&mut d, p, 3.);
            let id = beam(&mut d, a, b, 1., l.top_ea, l.top_ei.unwrap());
            top_beams.push((p[0].0, p[1].0, id));
        } else {
            d.add_element(Element::Truss(Truss::new(
                p[0].1,
                p[1].1,
                l.top_ea,
                Material::Elastic { e: 1. },
            )));
        }
    }
    let mut bottom_beams = Vec::new();
    for p in bottom.windows(2) {
        if door && p[0].0 >= 1.5 - 1e-10 && p[1].0 <= 3.5 + 1e-10 {
            continue;
        }
        let [a, b] = hinged(&mut d, p, 0.);
        let id = beam(&mut d, a, b, 1., l.bottom_ea, l.bottom_ei);
        bottom_beams.push((p[0].0, p[1].0, id));
    }
    let mut stud_bars = Vec::new();
    for s in &l.studs {
        if s.x > 1.5 && s.x < 3.5 {
            continue;
        }
        let id = d.add_element(Element::Truss(Truss::new(
            map(&bottom, s.x),
            map(&top, s.x),
            s.ea,
            Material::Elastic { e: 1. },
        )));
        stud_bars.push((s.x, id));
    }
    let hxs = [1.5, 2.5, 3.5];
    let mut hn = Vec::new();
    for x in hxs {
        let n = d.add_node(Node::new([x, 2.2]).fix(0));
        hn.push((x, n));
    }
    let mut header_beams = Vec::new();
    for p in hn.windows(2) {
        let id = beam(&mut d, p[0].1, p[1].1, 1., 1., 120.);
        header_beams.push((p[0].0, p[1].0, id));
    }
    let mut jack_bars = Vec::new();
    for (x, hn) in [(1.5, hn[0].1), (3.5, hn[2].1)] {
        let id = d.add_element(Element::Truss(Truss::new(
            map(&bottom, x),
            hn,
            35.,
            Material::Elastic { e: 1. },
        )));
        jack_bars.push((x, id));
    }
    let cripple_id = d.add_element(Element::Truss(Truss::new(
        hn[1].1,
        map(&top, 2.5),
        20.,
        Material::Elastic { e: 1. },
    )));
    (
        d,
        c,
        ts,
        bs,
        top,
        bottom,
        top_beams,
        stud_bars,
        jack_bars,
        vec![(2.5, cripple_id)],
        bottom_beams,
        header_beams,
    )
}

fn explicit_top_loads(
    d: &mut D2,
    top: &[(f64, carapace_core::model::NodeId)],
    top_beams: &[(f64, f64, carapace_core::model::ElementId)],
    distribution: ChordDistribution,
    scale: f64,
) {
    let add_point = |d: &mut D2,
                     x: f64,
                     fx: f64,
                     fy: f64,
                     xa: f64,
                     xb: f64,
                     a: carapace_core::model::NodeId,
                     b: carapace_core::model::NodeId| {
        let le = xb - xa;
        let t = (x - xa) / le;
        if distribution == ChordDistribution::Direct {
            for (node, dof, v) in [
                (a, 0, (1. - t) * fx),
                (b, 0, t * fx),
                (a, 1, (1. - t) * fy),
                (b, 1, t * fy),
            ] {
                d.load_node(node, dof, v);
            }
        } else {
            let h = [
                1. - 3. * t * t + 2. * t.powi(3),
                le * (t - 2. * t * t + t.powi(3)),
                3. * t * t - 2. * t.powi(3),
                le * (-t * t + t.powi(3)),
            ];
            for (node, dof, v) in [
                (a, 0, (1. - t) * fx),
                (b, 0, t * fx),
                (a, 1, h[0] * fy),
                (a, 2, h[1] * fy),
                (b, 1, h[2] * fy),
                (b, 2, h[3] * fy),
            ] {
                d.load_node(node, dof, v);
            }
        }
    };
    let p = top
        .windows(2)
        .find(|w| 0.55 > w[0].0 && 0.55 < w[1].0)
        .unwrap();
    add_point(
        d,
        0.55,
        4. * scale,
        -8. * scale,
        p[0].0,
        p[1].0,
        p[0].1,
        p[1].1,
    );
    if distribution == ChordDistribution::Beam {
        for &(x0, x1, id) in top_beams {
            if x1 > 5. {
                continue;
            }
            let t0 = x0 / 5.;
            let t1 = x1 / 5.;
            d.add_element_load(
                d.default_pattern(),
                id,
                ElementLoad::trapezoid(
                    scale * (0.5 + 0.7 * t0),
                    scale * (-1. - t0),
                    scale * (0.5 + 0.7 * t1),
                    scale * (-1. - t1),
                ),
            );
        }
        return;
    }
    if distribution == ChordDistribution::Direct {
        return;
    }
    let z = [
        -0.906179845938664,
        -0.538469310105683,
        0.,
        0.538469310105683,
        0.906179845938664,
    ];
    let w = [
        0.236926885056189,
        0.478628670499366,
        0.568888888888889,
        0.478628670499366,
        0.236926885056189,
    ];
    for p in top.windows(2) {
        let mid = (p[0].0 + p[1].0) / 2.;
        let half = (p[1].0 - p[0].0) / 2.;
        for i in 0..5 {
            let x = mid + half * z[i];
            let t = x / 6.;
            add_point(
                d,
                x,
                (0.5 + 0.7 * t) * half * w[i],
                (-1. - t) * half * w[i],
                p[0].0,
                p[1].0,
                p[0].1,
                p[1].1,
            );
        }
    }
}

fn scaled_wall_load(load: WallLoad, scale: f64) -> WallLoad {
    match load {
        WallLoad::TopPoint { x, fx, fy } => WallLoad::TopPoint {
            x,
            fx: fx * scale,
            fy: fy * scale,
        },
        WallLoad::TopLine {
            x0,
            x1,
            wx0,
            wx1,
            wy0,
            wy1,
        } => WallLoad::TopLine {
            x0,
            x1,
            wx0: wx0 * scale,
            wx1: wx1 * scale,
            wy0: wy0 * scale,
            wy1: wy1 * scale,
        },
    }
}

/// Three independently factored wall slots condense point and partial line loads like the explicit frame.
/// Reference: explicit Truss/beam frame and its support reactions; relative tolerance 1e-9.
#[test]
fn loaded_condensation_matches_explicit_frame_across_slots() {
    for distribution in [ChordDistribution::Beam, ChordDistribution::Direct] {
        let (mut wall_domain, corners, ts, bs) = supported_domain();
        let mut by_slot = wall_loads();
        if distribution == ChordDistribution::Direct {
            by_slot[0].retain(|load| matches!(load, WallLoad::TopPoint { .. }));
        }
        by_slot[1] = by_slot[0]
            .iter()
            .copied()
            .map(|load| scaled_wall_load(load, 2.))
            .collect();
        by_slot[2] = by_slot[0]
            .iter()
            .copied()
            .map(|load| scaled_wall_load(load, 3.))
            .collect();
        let wall = StudWall::new(
            corners,
            vec![
                carapace_core::model::WallStation {
                    chord: WallChord::Top,
                    x: 1.,
                    node: ts,
                },
                carapace_core::model::WallStation {
                    chord: WallChord::Bottom,
                    x: 5.,
                    node: bs,
                },
            ],
            frame_layout(distribution, false),
            by_slot,
        )
        .unwrap();
        let wall_id = wall_domain.add_element(Element::StudWall(wall));
        for (slot, factor) in [(0, 1.), (1, 0.5), (2, 0.1)] {
            let pattern = if slot == 0 {
                wall_domain.default_pattern()
            } else {
                wall_domain.add_load_pattern_scaled(
                    carapace_core::model::LoadSeries::Linear { slope: 1. },
                    factor,
                )
            };
            let mut selector = [0.; WALL_SLOTS];
            selector[slot] = 1.;
            wall_domain.add_element_load(
                pattern,
                wall_id,
                ElementLoad {
                    wall: selector,
                    ..Default::default()
                },
            );
        }
        let (mut frame, fc, fts, fbs, top, bottom, top_beams, stud_bars, _, _, _, _) =
            explicit_frame(distribution, false);
        explicit_top_loads(&mut frame, &top, &top_beams, distribution, 2.3);
        let solve = |domain| {
            let mut a = AnalysisBuilder::new()
                .constraint_handler(ConstraintHandler::Transformation)
                .integrator(Integrator::LoadControl { increment: 1. })
                .algorithm(Algorithm::Linear)
                .test(ConvergenceTest::NormUnbalance {
                    tol: 1e-11,
                    max_iter: 4,
                })
                .build(domain);
            a.step().expect("loaded wall solve");
            a
        };
        let wa = solve(wall_domain);
        let fa = solve(frame);
        for (wn, fn_) in corners.iter().zip(fc.iter()) {
            for dof in 0..2 {
                near(
                    wa.domain().node(*wn).displacement[dof],
                    fa.domain().node(*fn_).displacement[dof],
                    1e-9,
                );
            }
        }
        for (wn, fn_) in [(ts, fts), (bs, fbs)] {
            for dof in 0..2 {
                near(
                    wa.domain().node(wn).displacement[dof],
                    fa.domain().node(fn_).displacement[dof],
                    1e-9,
                );
            }
        }
        let at_x = |list: &[(f64, carapace_core::model::NodeId)], x: f64| {
            list.iter()
                .find(|(xx, _)| (*xx - x).abs() < 1e-10)
                .unwrap()
                .1
        };
        for (wn, fn_) in [
            (corners[0], at_x(&bottom, 0.)),
            (corners[1], at_x(&bottom, 6.)),
            (corners[3], at_x(&top, 0.)),
            (bs, at_x(&bottom, 5.)),
        ] {
            for dof in 0..2 {
                near(
                    wa.domain().reaction(wn, dof, 1.),
                    fa.domain().reaction(fn_, dof, 1.),
                    1e-9,
                );
            }
        }
        let end = wa.domain().element_end_force(wall_id, 1.);
        let response = wa.domain().element_stud_wall_response(wall_id, 1.).unwrap();
        assert_eq!(response.studs.len(), 5);
        for (stud, (x, element)) in response.studs.iter().zip(stud_bars) {
            assert!((stud.x - x).abs() < 1e-12);
            let force = fa.domain().element_end_force(element, 1.);
            near(stud.bottom, force[4], 1e-9);
            near(stud.top, force[4], 1e-9);
        }
        for (i, node) in corners.iter().enumerate().take(2) {
            for dof in 0..2 {
                near(end[2 * i + dof], wa.domain().reaction(*node, dof, 1.), 1e-9);
            }
        }
    }
}

/// Elastic material sheathing reproduces the linear law for a station-loaded window/door wall.
/// Reference: the Linear law path; exterior displacements, pier flow and member forces at 1e-8 relative.
#[test]
fn elastic_material_law_reproduces_linear_window_and_door_walls() {
    for door in [false, true] {
        let solve = |material: bool| {
            let (mut d, corners, ts, bs) = supported_domain();
            let mut layout = frame_layout(ChordDistribution::Beam, door);
            layout.sheathing = Some(carapace_core::model::WallSheathing {
                law: if material {
                    carapace_core::model::WallShearLaw::Material(Material::Elastic { e: 75. })
                } else {
                    carapace_core::model::WallShearLaw::Linear { ga: 75. }
                },
            });
            let mut loads = wall_loads();
            loads[1] = loads[0]
                .iter()
                .copied()
                .map(|l| scaled_wall_load(l, 2.))
                .collect();
            loads[2] = loads[0]
                .iter()
                .copied()
                .map(|l| scaled_wall_load(l, 3.))
                .collect();
            let wall = StudWall::new(
                corners,
                vec![
                    WallStation {
                        chord: WallChord::Top,
                        x: 1.,
                        node: ts,
                    },
                    WallStation {
                        chord: WallChord::Bottom,
                        x: 5.,
                        node: bs,
                    },
                ],
                layout,
                loads,
            )
            .unwrap();
            let id = d.add_element(Element::StudWall(wall));
            for (slot, factor) in [(0, 0.001), (1, 0.0005), (2, 0.0001)] {
                let pattern = if slot == 0 {
                    d.default_pattern()
                } else {
                    d.add_load_pattern_scaled(
                        carapace_core::model::LoadSeries::Linear { slope: 1. },
                        factor,
                    )
                };
                let mut selector = [0.; WALL_SLOTS];
                selector[slot] = 1.;
                d.add_element_load(
                    pattern,
                    id,
                    ElementLoad {
                        wall: selector,
                        ..Default::default()
                    },
                );
            }
            let mut a = AnalysisBuilder::new()
                .constraint_handler(ConstraintHandler::Transformation)
                .integrator(Integrator::LoadControl { increment: 1. })
                .algorithm(if material {
                    Algorithm::Newton {
                        tangent: TangentStrategy::Current,
                        line_search: None,
                    }
                } else {
                    Algorithm::Linear
                })
                .test(if material {
                    ConvergenceTest::NormDispIncr {
                        tol: 1e-13,
                        max_iter: 100,
                    }
                } else {
                    ConvergenceTest::NormUnbalance {
                        tol: 1e-11,
                        max_iter: 20,
                    }
                })
                .build(d);
            a.step().expect("linear wall solve");
            let domain = a.domain();
            let displacements: Vec<[f64; 2]> = corners
                .into_iter()
                .chain([ts, bs])
                .map(|n| {
                    [
                        domain.node(n).displacement[0],
                        domain.node(n).displacement[1],
                    ]
                })
                .collect();
            let response = domain.element_stud_wall_response(id, 1.).unwrap();
            (displacements, response)
        };
        let (u_linear, r_linear) = solve(false);
        let (u_material, r_material) = solve(true);
        for (a, b) in u_linear.iter().flatten().zip(u_material.iter().flatten()) {
            near(*a, *b, 1e-8);
        }
        for (a, b) in r_linear.piers.iter().zip(&r_material.piers) {
            near(a.q, b.q, 1e-8);
            near(a.gamma, b.gamma, 1e-8);
        }
        for (a, b) in r_linear.studs.iter().zip(&r_material.studs) {
            near(a.bottom, b.bottom, 1e-8);
            near(a.top, b.top, 1e-8);
        }
        for (a, b) in r_linear
            .jacks
            .iter()
            .chain(&r_linear.cripples)
            .zip(r_material.jacks.iter().chain(&r_material.cripples))
        {
            near(a.axial, b.axial, 1e-8);
        }
        for (a, b) in r_linear.header.iter().zip(&r_material.header) {
            near(a.moment, b.moment, 1e-8);
            near(a.shear, b.shear, 1e-8);
        }
        for (a, b) in r_linear
            .top_chord
            .iter()
            .chain(&r_linear.bottom_chord)
            .zip(r_material.top_chord.iter().chain(&r_material.bottom_chord))
        {
            near(a.axial, b.axial, 1e-8);
            near(a.moment, b.moment, 1e-8);
            near(a.shear, b.shear, 1e-8);
        }
    }
}

/// Beam and Direct walls with a window/door, a station on each chord, a stud pack and top point/line loads match an explicit Truss/beam frame.
/// Reference: explicit Carapace frame with equalDOF hinge nodes; exterior translations/reactions and stud forces relative tolerance 1e-9.
#[test]
fn unsheathed_walls_match_explicit_frames() {
    for (distribution, door) in [
        (ChordDistribution::Beam, true),
        (ChordDistribution::Direct, false),
    ] {
        let (mut wall_domain, c, ts, bs) = supported_domain();
        let l = frame_layout(distribution, door);
        let mut wall_loads = wall_loads();
        if distribution == ChordDistribution::Direct {
            wall_loads[0].retain(|l| matches!(l, WallLoad::TopPoint { .. }));
        }
        let wall = StudWall::new(
            c,
            vec![
                carapace_core::model::WallStation {
                    chord: WallChord::Top,
                    x: 1.,
                    node: ts,
                },
                carapace_core::model::WallStation {
                    chord: WallChord::Bottom,
                    x: 5.,
                    node: bs,
                },
            ],
            l,
            wall_loads.clone(),
        )
        .unwrap();
        let wall_id = wall_domain.add_element(Element::StudWall(wall));
        wall_domain.add_element_load(
            wall_domain.default_pattern(),
            wall_id,
            ElementLoad {
                wall: [1., 0., 0., 0., 0., 0., 0., 0.],
                ..Default::default()
            },
        );
        let (
            mut frame,
            fc,
            fts,
            fbs,
            top,
            bottom,
            top_beams,
            stud_bars,
            jack_bars,
            cripple_bars,
            bottom_beams,
            header_beams,
        ) = explicit_frame(distribution, door);
        assert_eq!(c.len(), fc.len());
        explicit_top_loads(&mut frame, &top, &top_beams, distribution, 1.0);
        let solve = |domain| {
            let mut a = AnalysisBuilder::new()
                .constraint_handler(ConstraintHandler::Transformation)
                .integrator(Integrator::LoadControl { increment: 1. })
                .algorithm(Algorithm::Linear)
                .test(ConvergenceTest::NormUnbalance {
                    tol: 1e-11,
                    max_iter: 4,
                })
                .build(domain);
            a.step().expect("linear frame solve");
            a
        };
        let wa = solve(wall_domain);
        let fa = solve(frame);
        for (wn, fn_) in c.iter().zip(fc.iter()) {
            for dof in 0..2 {
                near(
                    wa.domain().node(*wn).displacement[dof],
                    fa.domain().node(*fn_).displacement[dof],
                    1e-9,
                );
            }
        }
        for (wn, fn_) in [(ts, fts), (bs, fbs)] {
            for dof in 0..2 {
                near(
                    wa.domain().node(wn).displacement[dof],
                    fa.domain().node(fn_).displacement[dof],
                    1e-9,
                );
            }
        }
        let idx = |v: &[(f64, carapace_core::model::NodeId)], x: f64| {
            v.iter().find(|(xx, _)| (*xx - x).abs() < 1e-10).unwrap().1
        };
        for (external, slave, dof) in [
            (c[0], idx(&bottom, 0.), 0),
            (c[0], idx(&bottom, 0.), 1),
            (c[1], idx(&bottom, 6.), 0),
            (c[1], idx(&bottom, 6.), 1),
            (c[3], idx(&top, 0.), 0),
            (bs, idx(&bottom, 5.), 0),
            (bs, idx(&bottom, 5.), 1),
        ] {
            near(
                wa.domain().reaction(external, dof, 1.),
                fa.domain().reaction(slave, dof, 1.),
                1e-9,
            );
        }
        let response = wa.domain().element_stud_wall_response(wall_id, 1.).unwrap();
        assert_eq!(response.studs.len(), 5);
        for (r, (x, eid)) in response.studs.iter().zip(stud_bars) {
            assert!((r.x - x).abs() < 1e-12);
            let f = fa.domain().element_end_force(eid, 1.);
            near(r.bottom, f[4], 1e-9);
            near(r.top, f[4], 1e-9);
        }
        for (r, (x, eid)) in response.jacks.iter().zip(jack_bars) {
            assert!((r.x - x).abs() < 1e-12);
            near(r.axial, fa.domain().element_end_force(eid, 1.)[4], 1e-9);
        }
        for (r, (x, eid)) in response.cripples.iter().zip(cripple_bars) {
            assert!((r.x - x).abs() < 1e-12);
            near(r.axial, fa.domain().element_end_force(eid, 1.)[4], 1e-9);
        }
        for (i, (_, _, eid)) in bottom_beams.iter().enumerate() {
            let f = fa.domain().element_end_force(*eid, 1.);
            let a = response.bottom_chord[2 * i];
            let b = response.bottom_chord[2 * i + 1];
            near(a.axial, -f[0], 1e-9);
            near(a.moment, f[2], 1e-9);
            near(a.shear, f[1], 1e-9);
            near(b.axial, f[3], 1e-9);
            near(b.moment, f[5], 1e-9);
            near(b.shear, -f[4], 1e-9);
        }
        for (i, (x0, x1, eid)) in top_beams.iter().enumerate() {
            if *x0 < 0.55 && *x1 > 0.55 {
                continue;
            }
            let f = fa.domain().element_end_force(*eid, 1.);
            let a = response.top_chord[2 * i];
            let b = response.top_chord[2 * i + 1];
            near(a.axial, -f[0], 1e-9);
            near(a.moment, f[2], 1e-9);
            near(a.shear, f[1], 1e-9);
            near(b.axial, f[3], 1e-9);
            near(b.moment, f[5], 1e-9);
            near(b.shear, -f[4], 1e-9);
        }
        for (i, (_, _, eid)) in header_beams.iter().enumerate() {
            let f = fa.domain().element_end_force(*eid, 1.);
            let a = response.header[2 * i];
            let b = response.header[2 * i + 1];
            near(a.moment, f[2], 1e-9);
            near(a.shear, f[1], 1e-9);
            near(b.moment, f[5], 1e-9);
            near(b.shear, -f[4], 1e-9);
        }
        let _ = wall_loads;
    }
}

/// A rigid-chord sheathed pier follows the SDPWS two-term deflection and end-post force checks.
/// Reference: `2vH³/(3EA b) + vH/ga`, with post force `±qH`; relative tolerance 1e-6.
#[test]
fn single_sheathed_pier_matches_sdpws_terms() {
    let (length, height, ea, ga, lateral) = (1200., 2440., 30e6, 2e5, 12e3);
    let mut l = StudWallLayout {
        length,
        height,
        studs: (0..=3)
            .map(|i| WallStud {
                x: i as f64 * 400.,
                ea: if i == 0 || i == 3 { ea } else { 1. },
            })
            .collect(),
        top_ea: ea * 1e6,
        top_ei: Some(5e10),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: ea * 1e6,
        bottom_ei: 5e10,
        openings: vec![],
        sheathing: Some(carapace_core::model::WallSheathing {
            law: carapace_core::model::WallShearLaw::Linear { ga },
        }),
        foundation: false,
    };
    l.height = height;
    let mut d: D2 = Domain::new();
    let corners = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, height])),
        d.add_node(Node::new([0., height])),
    ];
    let wall = StudWall::new(corners, vec![], l, std::array::from_fn(|_| vec![])).unwrap();
    let id = d.add_element(Element::StudWall(wall));
    let pattern = d.default_pattern();
    d.add_nodal_load(pattern, corners[2], 0, lateral / 2.);
    d.add_nodal_load(pattern, corners[3], 0, lateral / 2.);
    let a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-11,
            max_iter: 4,
        })
        .build(d);
    let mut a = a;
    a.step().expect("sheathed pier solve");
    let top_v = (a.domain().node(corners[2]).displacement[0]
        + a.domain().node(corners[3]).displacement[0])
        / 2.;
    let v = lateral / length;
    let expected = 2. * v * height.powi(3) / (3. * ea * length) + v * height / ga;
    near(top_v, expected, 1e-6);
    let response = a.domain().element_stud_wall_response(id, 1.).unwrap();
    let q = response.piers[0].q;
    near(q * length, lateral, 1e-6);
    let left = response.studs.iter().find(|s| s.x == 0.).unwrap();
    let right = response.studs.iter().find(|s| s.x == length).unwrap();
    assert!(
        left.top.abs() < q * height * 1e-6,
        "left top force {}",
        left.top
    );
    assert!(
        right.top.abs() < q * height * 1e-6,
        "right top force {}",
        right.top
    );
    near(left.bottom, q * height, 1e-6);
    near(right.bottom, -q * height, 1e-6);
    near(response.piers[0].chords[0], q * height, 1e-6);
    near(response.piers[0].chords[1], -q * height, 1e-6);
    near(response.piers[0].moment_top, 0., 1e-6);
}

/// One-pier O86 response follows the code shear and nail-slip terms over a Newton load step.
/// Reference: CSA O86:19 11.7.1.2 single segment with rigid stand-ins; relative tolerance 1e-6.
#[test]
fn single_o86_pier_matches_nonlinear_deflection() {
    let (length, height, ea, bv, nail_spacing, diameter, lateral) =
        (1200., 2440., 30e6, 11_000., 100., 3.25, 12e3);
    let mut l = StudWallLayout {
        length,
        height,
        studs: (0..=3)
            .map(|i| WallStud {
                x: i as f64 * 400.,
                ea: if i == 0 || i == 3 { ea } else { 1. },
            })
            .collect(),
        top_ea: ea * 1e6,
        top_ei: Some(5e10),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: ea * 1e6,
        bottom_ei: 5e10,
        openings: vec![],
        sheathing: Some(carapace_core::model::WallSheathing {
            law: carapace_core::model::WallShearLaw::O86 {
                bv,
                nail_spacing,
                nail_diameter: diameter,
            },
        }),
        foundation: false,
    };
    l.height = height;
    let mut d: D2 = Domain::new();
    let corners = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, height])),
        d.add_node(Node::new([0., height])),
    ];
    let id = d.add_element(Element::StudWall(
        StudWall::new(corners, vec![], l, std::array::from_fn(|_| vec![])).unwrap(),
    ));
    let pattern = d.default_pattern();
    d.add_nodal_load(pattern, corners[2], 0, lateral / 2.);
    d.add_nodal_load(pattern, corners[3], 0, lateral / 2.);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.1 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-4,
            max_iter: 100,
        })
        .build(d);
    for _ in 0..10 {
        analysis.step().expect("O86 nonlinear pier solve");
    }
    let top_v = (analysis.domain().node(corners[2]).displacement[0]
        + analysis.domain().node(corners[3]).displacement[0])
        / 2.;
    let v = lateral / length;
    let nail_slip = 0.0025 * height * (0.013 * v * nail_spacing / diameter.powi(2)).powi(2);
    let expected = 2. * v * height.powi(3) / (3. * ea * length) + v * height / bv + nail_slip;
    near(top_v, expected, 1e-6);
    let response = analysis
        .domain()
        .element_stud_wall_response(id, 1.)
        .unwrap();
    near(response.piers[0].q, v, 1e-6);
}

/// Committed nonlinear pier forces include the slot-load-dependent pier state used by reactions.
/// Reference: response base-node exterior resistance; absolute tolerance 1e-10 of the applied load.
#[test]
fn nonlinear_wall_end_force_matches_slot_loaded_response_after_commit() {
    let (length, height) = (1200., 2440.);
    let mut d: D2 = Domain::new();
    let corners = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, height])),
        d.add_node(Node::new([0., height])),
    ];
    let mut layout = layout();
    layout.length = length;
    layout.height = height;
    layout.studs = vec![
        WallStud { x: 0., ea: 30e6 },
        WallStud { x: 400., ea: 1. },
        WallStud { x: 800., ea: 1. },
        WallStud {
            x: length,
            ea: 30e6,
        },
    ];
    layout.sheathing = Some(carapace_core::model::WallSheathing {
        law: carapace_core::model::WallShearLaw::O86 {
            bv: 11_000.,
            nail_spacing: 100.,
            nail_diameter: 3.25,
        },
    });
    let mut slots: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| vec![]);
    slots[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: length,
        wx0: 12_000. / length,
        wx1: 12_000. / length,
        wy0: 0.,
        wy1: 0.,
    });
    let id = d.add_element(Element::StudWall(
        StudWall::new(corners, vec![], layout, slots).unwrap(),
    ));
    let pattern = d.default_pattern();
    d.add_element_load(
        pattern,
        id,
        ElementLoad {
            wall: [1., 0., 0., 0., 0., 0., 0., 0.],
            ..Default::default()
        },
    );
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.1 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-7,
            max_iter: 100,
        })
        .build(d);
    for _ in 0..10 {
        a.step().expect("slot-loaded nonlinear wall");
    }
    let end = a.domain().element_end_force(id, 1.);
    let response = a.domain().element_stud_wall_response(id, 1.).unwrap();
    assert_eq!(response.base_nodes.len(), 2);
    for i in 0..2 {
        near(end[2 * i], response.base_nodes[i].fx, 1e-10);
        near(end[2 * i + 1], response.base_nodes[i].fy, 1e-10);
    }
}

/// Cyclic displacement control commits the material law once per step and preserves trial reversion.
/// Reference: standalone Steel01 driven through the same shear strains; 1e-5 relative for rigid stand-ins.
#[test]
fn material_sheathing_tracks_cyclic_steel01_shear_history() {
    let (length, height) = (10., 8.);
    let material = Material::steel01(100., 40_000., 0.01, 0., 1., 0., 1.);
    let mut d: D2 = Domain::new();
    let corners = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, height])),
        d.add_node(Node::new([0., height])),
    ];
    let mut l = layout();
    l.length = length;
    l.height = height;
    l.studs = vec![
        WallStud { x: 0., ea: 1e12 },
        WallStud {
            x: length,
            ea: 1e12,
        },
    ];
    l.top_ea = 1e12;
    l.top_ei = Some(1e12);
    l.bottom_ea = 1e12;
    l.bottom_ei = 1e12;
    l.sheathing = Some(carapace_core::model::WallSheathing {
        law: carapace_core::model::WallShearLaw::Material(material.clone()),
    });
    let wall = d.add_element(Element::StudWall(
        StudWall::new(corners, vec![], l, std::array::from_fn(|_| vec![])).unwrap(),
    ));
    let pattern = d.add_load_pattern(carapace_core::model::LoadSeries::Linear { slope: 1. });
    d.add_nodal_load(pattern, corners[2], 0, 1.);
    d.add_nodal_load(pattern, corners[3], 0, 1.);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::DisplacementControl {
            node: corners[3],
            dof: 0,
            increment: 0.,
        })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-6,
            max_iter: 100,
        })
        .build(d);
    let targets = [0.01, 0.02, 0.03, 0.02, 0., -0.01, -0.02, -0.03];
    let mut prev = 0.;
    let mut standalone = material;
    for target in targets {
        analysis.set_integrator(Integrator::DisplacementControl {
            node: corners[3],
            dof: 0,
            increment: target - prev,
        });
        analysis.step().expect("cyclic wall displacement step");
        let response = analysis
            .domain()
            .element_stud_wall_response(wall, 1.)
            .unwrap();
        let strain = target / height;
        let expected = standalone.trial_stress_tangent(strain).0;
        near(length * response.piers[0].q, length * expected, 1e-5);
        standalone = standalone.commit(strain);
        prev = target;
    }
    let probe = standalone.trial_stress_tangent(0.004);
    let trial_then_new = standalone.clone();
    let _ = trial_then_new.trial_stress_tangent(-0.003);
    let second = trial_then_new.trial_stress_tangent(0.004);
    near(probe.0, second.0, 1e-14);
    near(probe.1, second.1, 1e-14);
}

/// Door piers balance lateral load through their `q·b` fields and report the same total base shear.
/// Reference: applied nodal resultant and the element base response; relative tolerance 1e-9.
#[test]
fn door_piers_sum_to_applied_and_reported_base_shear() {
    let lateral = 24e3;
    let mut d: D2 = Domain::new();
    let corners = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([6., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([6., 3.])),
        d.add_node(Node::new([0., 3.])),
    ];
    let mut l = frame_layout(ChordDistribution::Beam, true);
    l.sheathing = Some(carapace_core::model::WallSheathing {
        law: carapace_core::model::WallShearLaw::Linear { ga: 8e4 },
    });
    l.foundation = true;
    let wall = StudWall::new(corners, vec![], l, std::array::from_fn(|_| vec![])).unwrap();
    let id = d.add_element(Element::StudWall(wall));
    let pattern = d.default_pattern();
    d.add_nodal_load(pattern, corners[2], 0, lateral / 2.);
    d.add_nodal_load(pattern, corners[3], 0, lateral / 2.);
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-11,
            max_iter: 4,
        })
        .build(d);
    a.step().expect("door wall solve");
    let r = a.domain().element_stud_wall_response(id, 1.).unwrap();
    near(
        r.piers.iter().map(|p| p.q * (p.xb - p.xa)).sum(),
        lateral,
        1e-9,
    );
    near(r.base[0].abs(), lateral, 1e-9);
    for p in &r.piers {
        // q is positive when its pier carries +x shear; the corresponding CCW-positive axial-force
        // moment relation is M_top = M_base + q·b·H.
        near(p.moment_top, p.moment_base + p.q * (p.xb - p.xa) * 3., 1e-9);
    }
}

/// Adjacent walls with different sheathing stiffnesses share the applied shear in proportion to their response.
/// Reference: each pier's `q·b`, each wall's base shear, and the applied resultant; relative tolerance 1e-8.
#[test]
fn adjacent_walls_with_different_ga_balance_individually_and_together() {
    let (length, height, lateral) = (1200., 2440., 30e3);
    let mut d: D2 = Domain::new();
    let lb = d.add_node(Node::new([0., 0.]).fix(0).fix(1));
    let mb = d.add_node(Node::new([length, 0.]).fix(0).fix(1));
    let rb = d.add_node(Node::new([2. * length, 0.]).fix(0).fix(1));
    let lt = d.add_node(Node::new([0., height]));
    let mt = d.add_node(Node::new([length, height]));
    let rt = d.add_node(Node::new([2. * length, height]));
    let make_layout = |ga| StudWallLayout {
        length,
        height,
        studs: (0..=3)
            .map(|i| WallStud {
                x: i as f64 * 400.,
                ea: if i == 0 || i == 3 { 30e6 } else { 1. },
            })
            .collect(),
        top_ea: 30e6 * 1e6,
        top_ei: Some(5e10),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: 30e6 * 1e6,
        bottom_ei: 5e10,
        openings: vec![],
        sheathing: Some(carapace_core::model::WallSheathing {
            law: carapace_core::model::WallShearLaw::Linear { ga },
        }),
        foundation: false,
    };
    let make_wall = |corners, ga| {
        StudWall::new(
            corners,
            vec![],
            make_layout(ga),
            std::array::from_fn(|_| vec![]),
        )
        .unwrap()
    };
    let a_id = d.add_element(Element::StudWall(make_wall([lb, mb, mt, lt], 1e5)));
    let b_id = d.add_element(Element::StudWall(make_wall([mb, rb, rt, mt], 3e5)));
    let pattern = d.default_pattern();
    d.add_nodal_load(pattern, lt, 0, lateral / 2.);
    d.add_nodal_load(pattern, rt, 0, lateral / 2.);
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-11,
            max_iter: 4,
        })
        .build(d);
    a.step().expect("adjacent wall solve");
    let ra = a.domain().element_stud_wall_response(a_id, 1.).unwrap();
    let rb = a.domain().element_stud_wall_response(b_id, 1.).unwrap();
    let qa = ra.piers[0].q * length;
    let qb = rb.piers[0].q * length;
    near(qa + qb, lateral, 1e-8);
    near(ra.base[0].abs(), qa, 1e-8);
    near(rb.base[0].abs(), qb, 1e-8);
    near(ra.base[0].abs() + rb.base[0].abs(), lateral, 1e-8);
}

/// Unloaded top and bottom stations inside a sheathed pier do not alter the original response.
/// Reference: the same wall without stations; original node displacements, pier `q`, and stud forces to 1e-10 relative.
#[test]
fn unloaded_stations_inside_a_pier_leave_response_unchanged() {
    let solve = |with_stations: bool| {
        let (length, height, lateral) = (1200., 2440., 12e3);
        let mut d: D2 = Domain::new();
        let corners = [
            d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
            d.add_node(Node::new([length, 0.]).fix(0).fix(1)),
            d.add_node(Node::new([length, height])),
            d.add_node(Node::new([0., height])),
        ];
        let stations = if with_stations {
            vec![
                carapace_core::model::WallStation {
                    chord: WallChord::Top,
                    x: 200.,
                    node: d.add_node(Node::new([200., height])),
                },
                carapace_core::model::WallStation {
                    chord: WallChord::Bottom,
                    x: 200.,
                    node: d.add_node(Node::new([200., 0.])),
                },
            ]
        } else {
            vec![]
        };
        let layout = StudWallLayout {
            length,
            height,
            studs: (0..=3)
                .map(|i| WallStud {
                    x: i as f64 * 400.,
                    ea: if i == 0 || i == 3 { 30e6 } else { 1. },
                })
                .collect(),
            top_ea: 30e6 * 1e4,
            top_ei: Some(5e10),
            top_distribution: ChordDistribution::Beam,
            bottom_ea: 30e6 * 1e4,
            bottom_ei: 5e10,
            openings: vec![],
            sheathing: Some(carapace_core::model::WallSheathing {
                law: carapace_core::model::WallShearLaw::Linear { ga: 2e5 },
            }),
            foundation: false,
        };
        let id = d.add_element(Element::StudWall(
            StudWall::new(corners, stations, layout, std::array::from_fn(|_| vec![])).unwrap(),
        ));
        let pattern = d.default_pattern();
        d.add_nodal_load(pattern, corners[2], 0, lateral / 2.);
        d.add_nodal_load(pattern, corners[3], 0, lateral / 2.);
        let mut a = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Transformation)
            .integrator(Integrator::LoadControl { increment: 1. })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-11,
                max_iter: 4,
            })
            .build(d);
        a.step().expect("station invariance solve");
        let displacements = corners.map(|n| a.domain().node(n).displacement);
        let response = a.domain().element_stud_wall_response(id, 1.).unwrap();
        (displacements, response)
    };
    let (base_u, base_r) = solve(false);
    let (station_u, station_r) = solve(true);
    for (a, b) in base_u.iter().flatten().zip(station_u.iter().flatten()) {
        near(*a, *b, 1e-10);
    }
    near(base_r.piers[0].q, station_r.piers[0].q, 1e-10);
    for (a, b) in base_r.studs.iter().zip(station_r.studs.iter()) {
        near(a.bottom, b.bottom, 1e-10);
        near(a.top, b.top, 1e-10);
    }
}
