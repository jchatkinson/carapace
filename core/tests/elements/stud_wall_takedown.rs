//! Stud-wall gravity takedowns against hand calculations.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    ChordDistribution, ElasticBeamColumn, Element, ElementLoad, GeomTransf, LoadSeries, Node,
    StudWall, StudWallLayout, WallCaseLoad, WallChord, WallLoad, WallOpening, WallStation,
    WallStud,
};
use carapace_core::model::{Domain, ElementId, NodeId};

type D2 = Domain<2, 3, NodeId, Element>;

struct BuiltWall {
    id: ElementId,
    corners: [NodeId; 4],
    bottom: Vec<(f64, NodeId)>,
}

fn layout(
    length: f64,
    studs: &[(f64, f64)],
    distribution: ChordDistribution,
    top_ei: Option<f64>,
    openings: Vec<WallOpening>,
) -> StudWallLayout {
    StudWallLayout {
        length,
        height: 2440.,
        studs: studs.iter().map(|&(x, ea)| WallStud { x, ea }).collect(),
        top_ea: 100e6,
        top_ei,
        top_distribution: distribution,
        bottom_ea: 50e6,
        bottom_ei: 6.0e9,
        openings,
        sheathing: None,
        foundation: false,
        bearing: None,
        floor: None,
    }
}

fn add_wall(
    d: &mut D2,
    corners: [NodeId; 4],
    mut l: StudWallLayout,
    loads: Vec<Vec<WallLoad>>,
    bottom: Vec<(f64, NodeId)>,
    top: Vec<(f64, NodeId)>,
) -> BuiltWall {
    let mut stations = Vec::new();
    for (chord, nodes) in [(WallChord::Bottom, &bottom), (WallChord::Top, &top)] {
        for &(x, node) in nodes {
            stations.push(WallStation { chord, x, node });
        }
    }
    // This helper is deliberately explicit about the shared chord constants.
    l.height = 2440.;
    let wall = StudWall::new(corners, stations, l, loads).unwrap();
    let id = d.add_element(Element::StudWall(wall));
    BuiltWall {
        id,
        corners,
        bottom,
    }
}

fn supported_wall(
    length: f64,
    studs: &[(f64, f64)],
    distribution: ChordDistribution,
    top_ei: Option<f64>,
    openings: Vec<WallOpening>,
    loads: Vec<Vec<WallLoad>>,
) -> (D2, BuiltWall) {
    supported_wall_mode(length, studs, distribution, top_ei, openings, loads, false)
}

fn supported_wall_mode(
    length: f64,
    studs: &[(f64, f64)],
    distribution: ChordDistribution,
    top_ei: Option<f64>,
    openings: Vec<WallOpening>,
    loads: Vec<Vec<WallLoad>>,
    foundation: bool,
) -> (D2, BuiltWall) {
    let mut d = D2::new();
    let h = 2440.;
    let corners = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, h])),
        d.add_node(Node::new([0., h]).fix(0)),
    ];
    let mut xs: Vec<f64> = studs
        .iter()
        .map(|&(x, _)| x)
        .filter(|&x| x > 0. && x < length && !openings.iter().any(|o| x > o.x0 && x < o.x1))
        .collect();
    for o in &openings {
        xs.extend([o.x0, o.x1]);
    }
    xs.sort_by(f64::total_cmp);
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-10);
    let bottom: Vec<_> = if foundation {
        vec![]
    } else {
        xs.into_iter()
            .map(|x| (x, d.add_node(Node::new([x, 0.]).fix(0).fix(1))))
            .collect()
    };
    let l = layout(length, studs, distribution, top_ei, openings);
    let mut l = l;
    l.foundation = foundation;
    let wall = add_wall(&mut d, corners, l, loads, bottom, vec![]);
    (d, wall)
}

fn empty_loads() -> Vec<Vec<WallLoad>> {
    vec![Vec::new(); 8]
}

fn add_wall_selector(
    d: &mut D2,
    id: ElementId,
    slot: usize,
    pattern: carapace_core::model::LoadPatternId,
) {
    d.add_element_load(
        pattern,
        id,
        ElementLoad {
            wall: Some(WallCaseLoad {
                case: slot as u32,
                factor: 1.,
            }),
            ..Default::default()
        },
    );
}

fn solve(d: D2) -> AnalysisBuilderResult {
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-11,
            max_iter: 4,
        })
        .build(d);
    a.step().expect("linear static takedown");
    a
}

type AnalysisBuilderResult = carapace_core::analysis::Analysis<2, 3, NodeId, Element>;

fn near(actual: f64, expected: f64, tol: f64, label: &str) {
    let error = (actual - expected).abs();
    let limit = if expected == 0. {
        1e-6
    } else {
        tol * actual.abs().max(expected.abs())
    };
    assert!(error <= limit, "{label}: computed {actual:.12e}, expected {expected:.12e}, abs error {error:.4e}, relative error {:.4e}", error / actual.abs().max(expected.abs()).max(1.));
}

fn stud_forces(a: &AnalysisBuilderResult, wall: &BuiltWall, time: f64) -> Vec<(f64, f64)> {
    a.domain()
        .element_stud_wall_response(wall.id, time)
        .unwrap()
        .studs
        .iter()
        .map(|r| (r.x, r.top))
        .collect()
}

fn at(values: &[(f64, f64)], x: f64) -> f64 {
    values
        .iter()
        .find(|(xx, _)| (*xx - x).abs() < 1e-8)
        .unwrap()
        .1
}

/// Uniform `Direct` top load follows stud tributary widths; vertical support reactions balance the 36 kN resultant.
/// Hand calculation: tributary width times −10 N/mm; tolerance 1e-9 relative (1e-6 N absolute for zeros).
#[test]
fn uniform_direct_blocking_tributaries() {
    let studs: Vec<_> = (0..=9).map(|i| (i as f64 * 400., 50e6)).collect();
    let mut loads = empty_loads();
    loads[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: 3600.,
        wx0: 0.,
        wx1: 0.,
        wy0: -10.,
        wy1: -10.,
    });
    let (mut d, wall) = supported_wall(
        3600.,
        &studs,
        ChordDistribution::Direct,
        None,
        vec![],
        loads,
    );
    let pattern = d.default_pattern();
    add_wall_selector(&mut d, wall.id, 0, pattern);
    let a = solve(d);
    let forces = stud_forces(&a, &wall, 1.);
    for &(x, expected) in &[
        (0., -2000.),
        (400., -4000.),
        (800., -4000.),
        (1200., -4000.),
        (1600., -4000.),
        (2000., -4000.),
        (2400., -4000.),
        (2800., -4000.),
        (3200., -4000.),
        (3600., -2000.),
    ] {
        near(at(&forces, x), expected, 1e-9, "stud axial force");
        let n = if x == 0. {
            wall.corners[0]
        } else if x == 3600. {
            wall.corners[1]
        } else {
            wall.bottom
                .iter()
                .find(|(xx, _)| (*xx - x).abs() < 1e-8)
                .unwrap()
                .1
        };
        near(
            a.domain().reaction(n, 1, 1.),
            -expected,
            1e-9,
            "bottom vertical reaction",
        );
    }
    let reactions: f64 = wall
        .bottom
        .iter()
        .map(|&(_, n)| a.domain().reaction(n, 1, 1.))
        .sum::<f64>()
        + a.domain().reaction(wall.corners[0], 1, 1.)
        + a.domain().reaction(wall.corners[1], 1, 1.);
    near(reactions, 36000., 1e-9, "total vertical reactions");
    near(
        a.domain().reaction(wall.corners[3], 0, 1.),
        0.,
        0.,
        "top-left horizontal reaction",
    );
}

/// A continuous top beam on effectively rigid studs gives the classic three-span reactions and hogging moments, distinct from tributary blocking.
/// Hand calculation: three-span beam coefficients (0.4ws,1.1ws) and (−0.1ws²); tolerance 1e-6 relative. Moment convention: recovered positive moment is sagging, so hogging is negative.
#[test]
fn uniform_beam_chord_continuous_reactions() {
    let studs = [(0., 50e12), (400., 50e12), (800., 50e12), (1200., 50e12)];
    let mut loads = empty_loads();
    loads[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: 1200.,
        wx0: 0.,
        wx1: 0.,
        wy0: -10.,
        wy1: -10.,
    });
    let (mut d, wall) = supported_wall(
        1200.,
        &studs,
        ChordDistribution::Beam,
        Some(1.216e10),
        vec![],
        loads,
    );
    let pattern = d.default_pattern();
    add_wall_selector(&mut d, wall.id, 0, pattern);
    let a = solve(d);
    let r = a.domain().element_stud_wall_response(wall.id, 1.).unwrap();
    for (x, expected) in [
        (0., -1600.),
        (400., -4400.),
        (800., -4400.),
        (1200., -1600.),
    ] {
        near(
            at(&stud_forces(&a, &wall, 1.), x),
            expected,
            1e-6,
            "beam-supported stud force",
        );
    }
    for x in [400., 800.] {
        let m = r
            .top_chord
            .iter()
            .find(|v| (v.x - x).abs() < 1e-8)
            .unwrap()
            .moment;
        near(
            m,
            -160000.,
            1e-6,
            "top chord interior moment; positive is sagging, negative hogging",
        );
    }
}

/// A point load at the middle stud is shared by beam action and three stud springs, while Direct blocking sends it to that stud alone.
/// Hand calculation: symmetric two-span beam reactions (P/(3+ks³/(3EI))); Direct tributary result; tolerance 1e-9 relative.
#[test]
fn point_on_stud_beam_and_direct() {
    let studs = [(0., 50e6), (400., 50e6), (800., 50e6)];
    let mut loads = empty_loads();
    loads[0].push(WallLoad::TopPoint {
        x: 400.,
        fx: 0.,
        fy: -10000.,
    });
    let (mut d, wall) = supported_wall(
        800.,
        &studs,
        ChordDistribution::Beam,
        Some(4.03e11),
        vec![],
        loads,
    );
    let pattern = d.default_pattern();
    add_wall_selector(&mut d, wall.id, 0, pattern);
    let a = solve(d);
    let s: f64 = 400.;
    let k = 50e6 / 2440.;
    let p = -10000.;
    let ei = 4.03e11;
    let rend = p / (3. + k * s.powi(3) / (3. * ei));
    let rm = p - 2. * rend;
    let f = stud_forces(&a, &wall, 1.);
    near(at(&f, 0.), rend, 1e-9, "left end stud");
    near(at(&f, 400.), rm, 1e-9, "middle stud");
    near(at(&f, 800.), rend, 1e-9, "right end stud");

    let mut loads = empty_loads();
    loads[0].push(WallLoad::TopPoint {
        x: 400.,
        fx: 0.,
        fy: -10000.,
    });
    let (mut d, wall) =
        supported_wall(800., &studs, ChordDistribution::Direct, None, vec![], loads);
    let pattern = d.default_pattern();
    add_wall_selector(&mut d, wall.id, 0, pattern);
    let a = solve(d);
    let f = stud_forces(&a, &wall, 1.);
    near(at(&f, 0.), 0., 1e-9, "left direct stud");
    near(at(&f, 400.), -10000., 1e-9, "middle direct stud");
    near(at(&f, 800.), 0., 1e-9, "right direct stud");
}

/// Direct point loads between studs follow the lever rule, and a load over a three-ply stud pack is carried by the pack.
/// Hand calculation: linear lever split between adjacent stud lines; tolerance 1e-9 relative (1e-6 N absolute for zeros).
#[test]
fn point_between_studs_and_pack() {
    let studs = [(0., 50e6), (400., 50e6), (800., 150e6), (1200., 50e6)];
    for (x, expected) in [
        (
            500.,
            vec![(0., 0.), (400., -6750.), (800., -2250.), (1200., 0.)],
        ),
        (
            800.,
            vec![(0., 0.), (400., 0.), (800., -9000.), (1200., 0.)],
        ),
    ] {
        let mut loads = empty_loads();
        loads[0].push(WallLoad::TopPoint {
            x,
            fx: 0.,
            fy: -9000.,
        });
        let (mut d, wall) = supported_wall(
            1200.,
            &studs,
            ChordDistribution::Direct,
            None,
            vec![],
            loads,
        );
        let pattern = d.default_pattern();
        add_wall_selector(&mut d, wall.id, 0, pattern);
        let a = solve(d);
        let f = stud_forces(&a, &wall, 1.);
        for (x, e) in expected {
            near(at(&f, x), e, 1e-9, "point-load stud force");
        }
    }
}

fn window_studs() -> Vec<(f64, f64)> {
    [
        0., 400., 800., 1100., 1200., 1600., 2000., 2400., 2500., 2800., 3200., 3600.,
    ]
    .into_iter()
    .map(|x| (x, 50e6))
    .collect()
}

/// A window header and cripples transfer full-height tributary load through the opening; repeating as a door leaves gravity forces unchanged.
/// Hand calculation: lever rule on the hinged chord plus simply supported header reactions and span moments; tolerance 1e-9 relative. Header moment convention: positive is sagging.
#[test]
fn window_and_door_header_takedown() {
    for door in [false, true] {
        let opening = WallOpening {
            x0: 1140.,
            x1: 2460.,
            head_height: 2000.,
            header_ei: 7.8e11,
            jack_ea: 50e6,
            door,
        };
        let mut loads = empty_loads();
        loads[0].push(WallLoad::TopLine {
            x0: 0.,
            x1: 3600.,
            wx0: 0.,
            wx1: 0.,
            wy0: -10.,
            wy1: -10.,
        });
        let (mut d, wall) = supported_wall(
            3600.,
            &window_studs(),
            ChordDistribution::Direct,
            None,
            vec![opening],
            loads,
        );
        let pattern = d.default_pattern();
        add_wall_selector(&mut d, wall.id, 0, pattern);
        let a = solve(d);
        let r = a.domain().element_stud_wall_response(wall.id, 1.).unwrap();
        for (x, e) in [
            (0., -2000.),
            (400., -4000.),
            (800., -3500.),
            (1100., -2000.),
            (2500., -2000.),
            (2800., -3500.),
            (3200., -4000.),
            (3600., -2000.),
        ] {
            near(
                at(&stud_forces(&a, &wall, 1.), x),
                e,
                1e-9,
                "full-height stud force",
            );
        }
        for (x, e) in [
            (1200., -2500.),
            (1600., -4000.),
            (2000., -4000.),
            (2400., -2500.),
        ] {
            near(
                r.cripples
                    .iter()
                    .find(|v| (v.x - x).abs() < 1e-8)
                    .unwrap()
                    .axial,
                e,
                1e-9,
                "cripple force",
            );
        }
        for j in &r.jacks {
            near(j.axial, -6500., 1e-9, "jack force");
        }
        for (x, e) in [
            (1140., 0.),
            (1200., 390000.),
            (1600., 1990000.),
            (2000., 1990000.),
            (2400., 390000.),
            (2460., 0.),
        ] {
            let m = r
                .header
                .iter()
                .find(|v| (v.x - x).abs() < 1e-8)
                .unwrap()
                .moment;
            near(m, e, 1e-9, "header moment; positive is sagging");
        }
        let sum: f64 = wall
            .bottom
            .iter()
            .map(|&(_, n)| a.domain().reaction(n, 1, 1.))
            .sum::<f64>()
            + a.domain().reaction(wall.corners[0], 1, 1.)
            + a.domain().reaction(wall.corners[1], 1, 1.);
        near(sum, 36000., 1e-9, "window/door vertical reactions");
    }
}

/// Internal foundation supports reproduce fixed bottom-station takedowns for solid, window, and door walls.
/// Reference: paired Carapace station model; member forces, top-corner displacements, and per-x reactions agree to 1e-10 relative.
#[test]
fn foundation_matches_fixed_bottom_stations() {
    for opening in [None, Some(false), Some(true)] {
        let openings = opening
            .map(|door| {
                vec![WallOpening {
                    x0: 1140.,
                    x1: 2460.,
                    head_height: 2000.,
                    header_ei: 7.8e11,
                    jack_ea: 50e6,
                    door,
                }]
            })
            .unwrap_or_default();
        let mut loads = empty_loads();
        loads[0].push(WallLoad::TopLine {
            x0: 0.,
            x1: 3600.,
            wx0: 0.,
            wx1: 0.,
            wy0: -10.,
            wy1: -10.,
        });
        let studs = if opening.is_some() {
            window_studs()
        } else {
            (0..=9).map(|i| (i as f64 * 400., 50e6)).collect()
        };
        let solve_mode = |foundation| {
            let (mut d, wall) = supported_wall_mode(
                3600.,
                &studs,
                ChordDistribution::Direct,
                None,
                openings.clone(),
                loads.clone(),
                foundation,
            );
            let pattern = d.default_pattern();
            add_wall_selector(&mut d, wall.id, 0, pattern);
            let a = solve(d);
            let r = a.domain().element_stud_wall_response(wall.id, 1.).unwrap();
            let corners = wall.corners.map(|n| a.domain().node(n).displacement);
            (r, corners)
        };
        let (stations, u_station) = solve_mode(false);
        let (foundation, u_foundation) = solve_mode(true);
        for (a, b) in u_station
            .iter()
            .flatten()
            .zip(u_foundation.iter().flatten())
        {
            near(*a, *b, 1e-10, "corner displacement station/foundation");
        }
        assert_eq!(stations.studs.len(), foundation.studs.len(), "stud count");
        assert_eq!(stations.cripples.len(), foundation.cripples.len());
        assert_eq!(stations.jacks.len(), foundation.jacks.len());
        assert_eq!(stations.header.len(), foundation.header.len());
        assert_eq!(stations.top_chord.len(), foundation.top_chord.len());
        assert_eq!(stations.bottom_chord.len(), foundation.bottom_chord.len());
        for (a, b) in stations.studs.iter().zip(&foundation.studs) {
            near(a.bottom, b.bottom, 1e-10, "stud bottom station/foundation");
            near(a.top, b.top, 1e-10, "stud top station/foundation");
        }
        for (a, b) in stations.cripples.iter().zip(&foundation.cripples) {
            near(a.axial, b.axial, 1e-10, "cripple station/foundation");
        }
        for (a, b) in stations.jacks.iter().zip(&foundation.jacks) {
            near(a.axial, b.axial, 1e-10, "jack station/foundation");
        }
        for (a, b) in stations.header.iter().zip(&foundation.header) {
            near(a.moment, b.moment, 1e-10, "header station/foundation");
            near(a.shear, b.shear, 1e-10, "header shear station/foundation");
        }
        assert_eq!(stations.base_nodes.len(), foundation.base_nodes.len());
        for a in &stations.base_nodes {
            let b = foundation
                .base_nodes
                .iter()
                .find(|b| (a.x - b.x).abs() < 1e-8)
                .unwrap();
            near(a.fx, b.fx, 1e-10, "base Fx station/foundation");
            near(a.fy, b.fy, 1e-10, "base Fy station/foundation");
        }
        for (a, b) in stations.top_chord.iter().zip(&foundation.top_chord) {
            near(a.axial, b.axial, 1e-10, "header station/foundation");
        }
        for (a, b) in stations.bottom_chord.iter().zip(&foundation.bottom_chord) {
            near(a.axial, b.axial, 1e-10, "bottom chord station/foundation");
        }
    }
}

/// Three independently factored wall-load slots combine dead, live, and linearly varying snow gravity at each stud.
/// Hand calculation: tributaries for dead/live and segment trapezoid reactions for snow; tolerance 1e-9 relative.
#[test]
fn factored_combination_across_slots() {
    let studs: Vec<_> = (0..=6)
        .map(|i| (i as f64 * 400., if i == 3 { 150e6 } else { 50e6 }))
        .collect();
    let mut loads = empty_loads();
    loads[0].extend([
        WallLoad::TopLine {
            x0: 0.,
            x1: 2400.,
            wx0: 0.,
            wx1: 0.,
            wy0: -3.,
            wy1: -3.,
        },
        WallLoad::TopPoint {
            x: 1200.,
            fx: 0.,
            fy: -8000.,
        },
    ]);
    loads[1].extend([
        WallLoad::TopLine {
            x0: 0.,
            x1: 1200.,
            wx0: 0.,
            wx1: 0.,
            wy0: -5.,
            wy1: -5.,
        },
        WallLoad::TopPoint {
            x: 1200.,
            fx: 0.,
            fy: -12000.,
        },
    ]);
    loads[2].push(WallLoad::TopLine {
        x0: 0.,
        x1: 2400.,
        wx0: 0.,
        wx1: 0.,
        wy0: 0.,
        wy1: -6.,
    });
    let (mut d, wall) = supported_wall(
        2400.,
        &studs,
        ChordDistribution::Direct,
        None,
        vec![],
        loads,
    );
    let p0 = d.add_load_pattern_scaled(LoadSeries::Linear { slope: 1. }, 1.2);
    let p1 = d.add_load_pattern_scaled(LoadSeries::Linear { slope: 1. }, 1.6);
    let p2 = d.add_load_pattern_scaled(LoadSeries::Linear { slope: 1. }, 0.5);
    add_wall_selector(&mut d, wall.id, 0, p0);
    add_wall_selector(&mut d, wall.id, 1, p1);
    add_wall_selector(&mut d, wall.id, 2, p2);
    let a = solve(d);
    let f = stud_forces(&a, &wall, 1.);
    let snow_w: Vec<f64> = (0..=6).map(|i| -6. * (i as f64 / 6.)).collect();
    let mut snow = [0.; 7];
    for i in 0..6 {
        let ell = 400.;
        let wa = snow_w[i];
        let wb = snow_w[i + 1];
        snow[i] += ell * (2. * wa + wb) / 6.;
        snow[i + 1] += ell * (wa + 2. * wb) / 6.;
    }
    let dl = [-2320., -4640., -4640., -31840., -1440., -1440., -720.];
    for i in 0..7 {
        near(
            at(&f, i as f64 * 400.),
            dl[i] + 0.5 * snow[i],
            1e-9,
            "factored stud force",
        );
    }
    near(
        f.iter().map(|(_, v)| *v).sum(),
        -50640.,
        1e-9,
        "factored total stud force",
    );
}

/// Two directly stacked walls transfer upper-storey gravity through shared floor nodes into the lower-storey
/// studs, with the lower end studs at half the interior `EA` so every lower stud shortens equally.
/// Hand calculation: each wall's uniform-load tributaries, the lower wall carrying its own plus the upper storey's
/// gravity. With uniform settlement the upper sole plate does not bend, so the hand calc is exact: tolerance 1e-9.
#[test]
fn stacked_storeys_takedown_uniform_settlement() {
    stacked_storeys(25e6, 1e-9);
}

/// As above with equal `EA` in every lower stud, as built. The lower end studs carry half the load, so they
/// shorten half as much; the upper wall's sole plate, continuous over the shared floor-line nodes, bends and
/// sheds load toward the ends. Tributary hand calc within 1% (the redistribution is about 0.7% here); upper
/// studs stay exact (1e-9) because the hinged top chord makes them statically determinate.
#[test]
fn stacked_storeys_takedown_sole_plate_redistribution() {
    stacked_storeys(50e6, 1e-2);
}

fn stacked_storeys(lower_end_ea: f64, lower_tol: f64) {
    let studs: Vec<_> = (0..=6).map(|i| (i as f64 * 400., 50e6)).collect();
    let mut lower_studs = studs.clone();
    lower_studs[0].1 = lower_end_ea;
    lower_studs[6].1 = lower_end_ea;
    let h = 2440.;
    let length = 2400.;
    let mut d = D2::new();
    let lb = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, h])),
        d.add_node(Node::new([0., h]).fix(0)),
    ];
    let ub = [
        lb[3],
        lb[2],
        d.add_node(Node::new([length, 2. * h])),
        d.add_node(Node::new([0., 2. * h]).fix(0)),
    ];
    let bnodes: Vec<_> = (1..6)
        .map(|i| {
            (
                i as f64 * 400.,
                d.add_node(Node::new([i as f64 * 400., 0.]).fix(0).fix(1)),
            )
        })
        .collect();
    let shared_top: Vec<_> = (1..6)
        .map(|i| (i as f64 * 400., d.add_node(Node::new([i as f64 * 400., h]))))
        .collect();
    let mut ll = empty_loads();
    ll[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: length,
        wx0: 0.,
        wx1: 0.,
        wy0: -8.,
        wy1: -8.,
    });
    let lower = add_wall(
        &mut d,
        lb,
        layout(
            length,
            &lower_studs,
            ChordDistribution::Direct,
            None,
            vec![],
        ),
        ll,
        bnodes.clone(),
        shared_top.clone(),
    );
    let mut ul = empty_loads();
    ul[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: length,
        wx0: 0.,
        wx1: 0.,
        wy0: -5.,
        wy1: -5.,
    });
    let upper = add_wall(
        &mut d,
        ub,
        layout(length, &studs, ChordDistribution::Direct, None, vec![]),
        ul,
        shared_top,
        vec![],
    );
    let p = d.default_pattern();
    add_wall_selector(&mut d, lower.id, 0, p);
    add_wall_selector(&mut d, upper.id, 0, p);
    let a = solve(d);
    for (x, e) in [
        (0., -1000.),
        (400., -2000.),
        (800., -2000.),
        (1200., -2000.),
        (1600., -2000.),
        (2000., -2000.),
        (2400., -1000.),
    ] {
        near(
            at(&stud_forces(&a, &upper, 1.), x),
            e,
            1e-9,
            "upper stud force",
        );
    }
    for (x, e) in [
        (0., -2600.),
        (400., -5200.),
        (800., -5200.),
        (1200., -5200.),
        (1600., -5200.),
        (2000., -5200.),
        (2400., -2600.),
    ] {
        near(
            at(&stud_forces(&a, &lower, 1.), x),
            e,
            lower_tol,
            "lower stud force",
        );
    }
}

/// A Direct wall bears on a simply supported beam; stud tributaries remain determinate and beam action shares midspan moment with the sole plate.
/// Hand calculation: stud tributaries, 12 kN beam reactions, and parallel flexural-stiffness moment split; split tolerance 1e-9, ignored-plate approximation within 1e-3.
#[test]
fn wall_on_support_beam() {
    let studs: Vec<_> = (0..=6).map(|i| (i as f64 * 400., 50e6)).collect();
    let length = 2400.;
    let h = 2440.;
    let mut d = D2::new();
    let c0 = d.add_node(Node::new([0., 0.]).fix(0).fix(1));
    let c1 = d.add_node(Node::new([length, 0.]).fix(1));
    let mut bottom = vec![(0., c0), (length, c1)];
    for i in 1..6 {
        let x = i as f64 * 400.;
        bottom.push((x, d.add_node(Node::new([x, 0.]))));
    }
    bottom.sort_by(|a, b| a.0.total_cmp(&b.0));
    let corners = [
        c0,
        c1,
        d.add_node(Node::new([length, h])),
        d.add_node(Node::new([0., h]).fix(0)),
    ];
    let mut loads = empty_loads();
    loads[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: length,
        wx0: 0.,
        wx1: 0.,
        wy0: -10.,
        wy1: -10.,
    });
    let wall_stations: Vec<_> = bottom
        .iter()
        .copied()
        .filter(|(x, _)| *x > 0. && *x < length)
        .collect();
    let wall = add_wall(
        &mut d,
        corners,
        layout(length, &studs, ChordDistribution::Direct, None, vec![]),
        loads,
        wall_stations,
        vec![],
    );
    for w in bottom.windows(2) {
        d.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            w[0].1,
            w[1].1,
            1.,
            1e9,
            2e13,
            GeomTransf::Linear,
        )));
    }
    // Beam supports use ux at the left and uy at both ends; the right wall node is otherwise supported by the beam.
    // The wall's bottom nodes are shared with the beam mesh.
    let pattern = d.default_pattern();
    add_wall_selector(&mut d, wall.id, 0, pattern);
    let a = solve(d);
    let f = stud_forces(&a, &wall, 1.);
    for (x, e) in [
        (0., -2000.),
        (400., -4000.),
        (800., -4000.),
        (1200., -4000.),
        (1600., -4000.),
        (2000., -4000.),
        (2400., -2000.),
    ] {
        near(at(&f, x), e, 1e-9, "wall stud force");
    }
    near(
        a.domain().reaction(c0, 1, 1.),
        12000.,
        1e-9,
        "left beam reaction",
    );
    near(
        a.domain().reaction(c1, 1, 1.),
        12000.,
        1e-9,
        "right beam reaction",
    );
    let total = 7.2e6;
    let expected = total * 2e13 / (2e13 + 6e9);
    let beam_elems: Vec<_> = a
        .domain()
        .elements()
        .filter_map(|(id, e)| matches!(e, Element::ElasticBeamColumn(_)).then_some(id))
        .collect();
    let mid = beam_elems[2];
    let end = a.domain().element_end_force(mid, 1.);
    // At the j end of a beam element, the counter-clockwise end moment acting on the element equals the
    // sagging-positive section moment there.
    let m = end[5];
    near(
        m,
        expected,
        1e-9,
        "beam midspan sagging moment stiffness split (j-end moment)",
    );
    assert!(
        (m / total - 1.).abs() < 1e-3,
        "sole-plate-ignored hand estimate: M_beam/M_total={}, expected near 1",
        m / total
    );
}
