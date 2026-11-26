//! Native `BuildingSpec` generation against hand statics: storey shear in the only resisting line, the
//! symmetric two-storey box with floor band (shears and the `H + d` overturning levers), the rigid-diaphragm
//! stiffness + torsion split, flexible tributary shares, a factored gravity takedown through a window header,
//! and input errors. Tolerances per test (1e-8 to 1e-12).
use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_input::building::*;

pub fn layout(window: bool) -> WallLayoutSpec {
    WallLayoutSpec {
        length: 4000.,
        height: 2440.,
        studs: if window {
            // Studs at 400 mm: king studs at 1200 and 2800, cripples over the window at 1600, 2000, 2400.
            (0..=10)
                .map(|i| WallStudSpec {
                    x: 400. * i as f64,
                    ea: 1e8,
                })
                .collect()
        } else {
            vec![
                WallStudSpec { x: 0., ea: 1e8 },
                WallStudSpec { x: 2000., ea: 1e8 },
                WallStudSpec { x: 4000., ea: 1e8 },
            ]
        },
        top_ea: 1e12,
        top_ei: Some(1e12),
        top_distribution: ChordDistributionSpec::Beam,
        bottom_ea: 1e12,
        bottom_ei: 1e12,
        openings: if window {
            vec![WallOpeningSpec {
                x0: 1200.,
                x1: 2800.,
                head_height: 2000.,
                header_ei: 1e12,
                jack_ea: 1e8,
                door: false,
            }]
        } else {
            vec![]
        },
        sheathing: Some(WallSheathingSpec {
            law: WallShearLawSpec::Linear { ga: 1e7 },
        }),
        foundation: false,
        bearing: None,
        floor: None,
    }
}
pub fn wall(name: &str, start: [f64; 2], direction: [f64; 2], window: bool) -> WallSpec {
    WallSpec {
        name: name.into(),
        start,
        direction,
        layout: layout(window),
        loads: vec![],
    }
}
pub fn box_walls(window: bool) -> Vec<WallSpec> {
    vec![
        wall("front", [0., 0.], [1., 0.], window),
        wall("back", [0., 4000.], [1., 0.], false),
        wall("left", [0., 0.], [0., 1.], false),
        wall("right", [4000., 0.], [0., 1.], false),
    ]
}
pub fn anchorage() -> AnchorageSpec {
    AnchorageSpec {
        tension: 15000.,
        slack: 0.,
        bearing: AnchorageBearingSpec::Linear { k: 1e9 },
    }
}
pub fn storey(walls: Vec<WallSpec>, depth: f64) -> StoreySpec {
    StoreySpec {
        wall_height: 2440.,
        floor_depth: depth,
        band_ga: 1e9,
        rim_ea: 1e12,
        anchorage: anchorage(),
        bearing: 1e9,
        walls,
        floors: vec![],
        lateral: vec![],
    }
}
pub fn spec(storeys: Vec<StoreySpec>, diaphragm: DiaphragmKind) -> BuildingSpec {
    BuildingSpec {
        storeys,
        foundation_anchorage: anchorage(),
        diaphragm,
        combination: [1.; 8],
    }
}
fn build_error(s: &BuildingSpec) -> BuildError {
    match build(s) {
        Ok(_) => panic!("expected building error"),
        Err(e) => e,
    }
}
fn solve(b: Building) -> carapace_core::model::Domain3 {
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-8,
            max_iter: 8,
        })
        .build(b.domain);
    a.step().expect("linear building solve");
    a.into_domain()
}

/// One storey: one X wall and two Z return walls under a rigid diaphragm, 10 kN in X: the X wall, the only line
/// resisting X, carries the whole storey shear (statics, 1e-9). Guards the generated foundation anchorage (a free
/// ground DOF makes the walls rock).
#[test]
fn one_storey_wall_and_returns() {
    let mut s = storey(
        vec![
            wall("x", [0., 0.], [1., 0.], false),
            wall("z0", [0., 0.], [0., 1.], false),
            wall("z4", [4000., 0.], [0., 1.], false),
        ],
        0.,
    );
    s.lateral.push(LateralLoad {
        slot: 0,
        fx: 10000.,
        fz: 0.,
        at: [2000., 4000. / 3.],
    });
    let mut b = build(&spec(vec![s], DiaphragmKind::Rigid)).unwrap();
    let loadfactors = [1.; 8];
    b.apply_combination(loadfactors).unwrap();
    let id = b.walls[0]["x"].element;
    let d = solve(b);
    let r = d.element_stud_wall_response(id, 1.).unwrap();
    let shear = r.base[0];
    assert!(
        (shear.abs() - 10000.).abs() <= 10000. * 1e-9,
        "base shear {shear}"
    );
}

/// A symmetric two-storey box (254 mm floor band, anchorage at every pier end, rigid diaphragms, 10 kN at each
/// level): each X wall takes half of each storey's shear, and the lower front wall's base moment is the statics
/// value with the true levers, `5000·(2H + d) + 5000·(H + d)` (O86 Fig. 11.7). 1e-8.
#[test]
fn two_storey_box_band_and_symmetric_shear() {
    let mut a = storey(box_walls(false), 254.);
    a.lateral.push(LateralLoad {
        slot: 0,
        fx: 10000.,
        fz: 0.,
        at: [2000., 2000.],
    });
    let mut b = storey(box_walls(false), 0.);
    b.lateral.push(LateralLoad {
        slot: 0,
        fx: 10000.,
        fz: 0.,
        at: [2000., 2000.],
    });
    let built = build(&spec(vec![a, b], DiaphragmKind::Rigid)).unwrap();
    assert!(!built.links[0].is_empty());
    assert_eq!(built.walls[1].len(), 4);
    let front = built.walls[1]["front"].element;
    let back = built.walls[1]["back"].element;
    let lower = built.walls[0]["front"].element;
    let d = solve(built);
    let m = d.element_stud_wall_response(lower, 1.).unwrap().piers[0]
        .moment_base
        .abs();
    let statics = 5000. * (2. * 2440. + 254.) + 5000. * (2440. + 254.);
    assert!(
        (m - statics).abs() <= statics * 1e-8,
        "lower base moment {m} vs {statics}"
    );
    let shear = |id| d.element_stud_wall_response(id, 1.).unwrap().base[0].abs();
    assert!(
        (shear(front) - 5000.).abs() <= 5000. * 1e-8,
        "front {} expected 5000",
        shear(front)
    );
    assert!(
        (shear(back) - 5000.).abs() <= 5000. * 1e-8,
        "back {} expected 5000",
        shear(back)
    );
}

/// Eccentric X force (10 kN at z = 3500) on a rigid diaphragm over four identical, separate walls (X lines at
/// z = 0 and 4000, Z lines at x = -500 and 4500, so no corner posts are shared): the master rotates positively about
/// +Y (`M = fx·dz`), and the wall shears equal the hand rigid-body split, which for equal line stiffnesses is
/// `V = F/2 ∓ M d/J`, `J = 2·2000² + 2·2500²` (k cancels). 1e-8.
#[test]
fn rigid_torsion_load() {
    let mut s = storey(
        vec![
            wall("front", [0., 0.], [1., 0.], false),
            wall("back", [0., 4000.], [1., 0.], false),
            wall("left", [-500., 0.], [0., 1.], false),
            wall("right", [4500., 0.], [0., 1.], false),
        ],
        0.,
    );
    s.lateral.push(LateralLoad {
        slot: 0,
        fx: 10000.,
        fz: 0.,
        at: [2000., 3500.],
    });
    let b = build(&spec(vec![s], DiaphragmKind::Rigid)).unwrap();
    let ids = ["front", "back", "left", "right"].map(|n| b.walls[0][n].element);
    let master = b.levels[0].master.unwrap();
    let d = solve(b);
    assert!(d.node(master).displacement[4] > 0., "rotation sign");
    let (f, moment, j) = (
        10000.,
        10000. * 1500.,
        2. * 2000f64.powi(2) + 2. * 2500f64.powi(2),
    );
    let want = [
        f / 2. - moment * 2000. / j,
        f / 2. + moment * 2000. / j,
        moment * 2500. / j,
        moment * 2500. / j,
    ];
    for (id, w) in ids.into_iter().zip(want) {
        let g = d.element_stud_wall_response(id, 1.).unwrap().base[0].abs();
        assert!((g - w).abs() <= f * 1e-8, "{g} vs {w}");
    }
}

/// Three flexible wall lines receive the summed uniform-area force in their tributary widths, as the hand strip reactions predict (1e-9).
#[test]
fn flexible_tributary_loads_sum() {
    let walls = vec![
        wall("a", [0., 0.], [1., 0.], false),
        wall("b", [0., 2000.], [1., 0.], false),
        wall("c", [0., 4000.], [1., 0.], false),
    ];
    let mut s = storey(walls, 0.);
    s.lateral.push(LateralLoad {
        slot: 0,
        fx: 12000.,
        fz: 0.,
        at: [0., 0.],
    });
    let b = build(&spec(vec![s], DiaphragmKind::Flexible)).unwrap();
    assert_eq!(b.levels[0].master, None);
    let ids = [
        b.walls[0]["a"].element,
        b.walls[0]["b"].element,
        b.walls[0]["c"].element,
    ];
    let d = solve(b);
    let got: Vec<f64> = ids
        .iter()
        .map(|id| d.element_stud_wall_response(*id, 1.).unwrap().base[0].abs())
        .collect();
    for (g, w) in got.into_iter().zip([3000., 6000., 3000.]) {
        assert!((g - w).abs() <= 12000. * 1e-9, "{g} != {w}");
    }
}

/// A floor spanning X between the two Z walls (D = 2, L = 1 per area) at 1.25D + 1.5L: the reported reactions
/// integrate to the applied load (1e-12), the left wall's base carries 32e6 N, and its window header gets the
/// cripple tributaries, 4.8e6 N per end (hand takedown, 1e-9).
#[test]
fn gravity_reactions_to_window_wall() {
    let mut s = storey(
        vec![
            wall("left", [0., 0.], [0., 1.], true),
            wall("right", [4000., 0.], [0., 1.], false),
            wall("front", [0., 0.], [1., 0.], false),
            wall("back", [0., 4000.], [1., 0.], false),
        ],
        0.,
    );
    s.floors.push(FloorSpec {
        bounds: [0., 4000., 0., 4000.],
        supports: vec!["left".into(), "right".into()],
        span_axis: SpanAxis::X,
        loads: vec![
            FloorLoadSpec::Area {
                slot: 0,
                w: 2.,
                region: None,
            },
            FloorLoadSpec::Area {
                slot: 1,
                w: 1.,
                region: None,
            },
        ],
        continuity: FloorContinuitySpec::Simple,
    });
    // Hand takedown: the opening's tributary goes to the header (top plate spans stud to stud: `Direct`).
    s.walls[0].layout.top_distribution = ChordDistributionSpec::Direct;
    s.walls[0].layout.top_ei = None;
    let mut b = build(&spec(vec![s], DiaphragmKind::Rigid)).unwrap();
    let r = &b.floor_reactions[0];
    let sum: f64 = r
        .per_support
        .iter()
        .flat_map(|s| s.lines.iter())
        .map(|x| x.w * (x.t1 - x.t0))
        .sum();
    assert!((sum - 48_000_000.).abs() <= 48_000_000. * 1e-12);
    let mut factors = [0.; 8];
    factors[0] = 1.25;
    factors[1] = 1.5;
    b.apply_combination(factors).unwrap();
    let id = b.walls[0]["left"].element;
    let d = solve(b);
    let r = d.element_stud_wall_response(id, 1.).unwrap();
    // Hand takedown: 1.25·(2·2000) + 1.5·(1·2000) = 8000 N/mm on the left wall's top plate; the cripples at 1600,
    // 2000 and 2400 (Direct tributaries 1400..2600) put 8000·1200 on the simply supported header: 4.8e6 per end.
    assert!(
        (r.base[1] - 32e6).abs() <= 32e6 * 1e-9,
        "base {}",
        r.base[1]
    );
    for end in &r.header_reactions[0..2] {
        assert!(
            (end[1].abs() - 4.8e6).abs() <= 4.8e6 * 1e-9,
            "header end {end:?} expected 4.8e6"
        );
    }
}

/// Misaligned upper walls, unknown floor supports, non-perpendicular supports, and non-unit directions all return named build errors.
#[test]
fn invalid_geometry_errors() {
    let low = storey(box_walls(false), 0.);
    let mut high = storey(box_walls(false), 0.);
    high.walls[0].start[0] = 1.;
    let e = build_error(&spec(vec![low, high], DiaphragmKind::Rigid));
    assert!(e.to_string().contains("front"));
    let mut s = storey(box_walls(false), 0.);
    s.floors.push(FloorSpec {
        bounds: [0., 4., 0., 4.],
        supports: vec!["missing".into()],
        span_axis: SpanAxis::X,
        loads: vec![],
        continuity: FloorContinuitySpec::Simple,
    });
    let e = build_error(&spec(vec![s], DiaphragmKind::Rigid));
    assert!(e.to_string().contains("missing"));
    let mut s = storey(vec![wall("bad", [0., 0.], [1., 0.], false)], 0.);
    s.floors.push(FloorSpec {
        bounds: [0., 4., 0., 4.],
        supports: vec!["bad".into()],
        span_axis: SpanAxis::X,
        loads: vec![],
        continuity: FloorContinuitySpec::Simple,
    });
    let e = build_error(&spec(vec![s], DiaphragmKind::Rigid));
    assert!(e.to_string().contains("bad"));
    let s = storey(vec![wall("bad", [0., 0.], [2., 0.], false)], 0.);
    assert!(build_error(&spec(vec![s], DiaphragmKind::Rigid))
        .to_string()
        .contains("direction"));
}

/// Regression for two limits the PySees template hit. (1) O86 nonlinear walls with more than 32 exterior nodes:
/// 8 m walls with studs at 400 mm get a band station at every stud, 40+ exterior nodes in the upper storey.
/// (2) A roof band: its links must hold the band's top row vertically. A symmetric two-storey box with 10 kN per
/// level in X, solved by Newton: each X wall takes half of each storey shear, and the upper front wall's base moment
/// is `5000 (H + d)` (statics, 1e-8).
#[test]
fn nonlinear_long_walls_with_roof_band() {
    let long = |name: &str, start: [f64; 2], direction: [f64; 2], length: f64| {
        let mut w = wall(name, start, direction, false);
        w.layout.length = length;
        w.layout.studs = (0..=(length / 400.) as usize)
            .map(|i| WallStudSpec {
                x: 400. * i as f64,
                ea: 1e8,
            })
            .collect();
        w.layout.sheathing = Some(WallSheathingSpec {
            law: WallShearLawSpec::O86 {
                bv: 11000.,
                nail_spacing: 150.,
                nail_diameter: 3.25,
            },
        });
        w
    };
    let walls = || {
        vec![
            long("front", [0., 0.], [1., 0.], 8000.),
            long("back", [0., 4000.], [1., 0.], 8000.),
            long("left", [0., 0.], [0., 1.], 4000.),
            long("right", [8000., 0.], [0., 1.], 4000.),
        ]
    };
    let level = |depth| {
        let mut s = storey(walls(), depth);
        s.lateral.push(LateralLoad {
            slot: 0,
            fx: 10000.,
            fz: 0.,
            at: [4000., 2000.],
        });
        s
    };
    let mut built = build(&spec(vec![level(254.), level(254.)], DiaphragmKind::Rigid)).unwrap();
    built.apply_combination([1.; 8]).unwrap();
    let ids = [0, 1].map(|s| {
        [
            built.walls[s]["front"].element,
            built.walls[s]["back"].element,
        ]
    });
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.25 })
        .algorithm(Algorithm::Newton {
            tangent: carapace_core::analysis::TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-7,
            max_iter: 100,
        })
        .build(built.domain);
    for _ in 0..4 {
        a.step().expect("nonlinear building step");
    }
    let d = a.into_domain();
    for (s, v) in [(0, 10000.), (1, 5000.)] {
        for id in ids[s] {
            let shear = d.element_stud_wall_response(id, 1.).unwrap().base[0].abs();
            assert!((shear - v).abs() <= v * 1e-8, "storey {s}: {shear} vs {v}");
        }
    }
    let m = d.element_stud_wall_response(ids[1][0], 1.).unwrap().piers[0]
        .moment_base
        .abs();
    let statics = 5000. * (2440. + 254.);
    assert!(
        (m - statics).abs() <= statics * 1e-8,
        "upper base moment {m} vs {statics}"
    );
}

/// A floor spanning Z on a non-square plan (8 m in X, 5 m in Z; bounds `[x0, x1, z0, z1]`): each X wall gets
/// `w L_z / 2` per unit length over its full 8 m, and the base of each carries `w L_x L_z / 2` (statics, 1e-9).
/// Regression: Z spans once read the X bounds as the span.
#[test]
fn z_span_floor_on_rectangular_plan() {
    let long = |name: &str, start: [f64; 2], direction: [f64; 2], length: f64| {
        let mut w = wall(name, start, direction, false);
        w.layout.length = length;
        w.layout.studs = vec![
            WallStudSpec { x: 0., ea: 1e8 },
            WallStudSpec {
                x: length / 2.,
                ea: 1e8,
            },
            WallStudSpec { x: length, ea: 1e8 },
        ];
        w
    };
    let mut s = storey(
        vec![
            long("front", [0., 0.], [1., 0.], 8000.),
            long("back", [0., 5000.], [1., 0.], 8000.),
            long("left", [0., 0.], [0., 1.], 5000.),
            long("right", [8000., 0.], [0., 1.], 5000.),
        ],
        0.,
    );
    s.floors.push(FloorSpec {
        bounds: [0., 8000., 0., 5000.],
        supports: vec!["front".into(), "back".into()],
        span_axis: SpanAxis::Z,
        loads: vec![FloorLoadSpec::Area {
            slot: 0,
            w: 0.002,
            region: None,
        }],
        continuity: FloorContinuitySpec::Simple,
    });
    let mut b = build(&spec(vec![s], DiaphragmKind::Rigid)).unwrap();
    for support in &b.floor_reactions[0].per_support {
        let (w, len): (f64, f64) = (
            support.lines[0].w,
            support.lines.iter().map(|l| l.t1 - l.t0).sum(),
        );
        assert!(
            (w - 5.).abs() <= 5e-9 && (len - 8000.).abs() <= 8e-6,
            "line {w} over {len}"
        );
    }
    b.apply_combination([1.; 8]).unwrap();
    let ids = ["front", "back"].map(|n| b.walls[0][n].element);
    let d = solve(b);
    for id in ids {
        let v = d.element_stud_wall_response(id, 1.).unwrap().base[1];
        assert!((v - 40_000.).abs() <= 40_000. * 1e-9, "base {v}");
    }
}
