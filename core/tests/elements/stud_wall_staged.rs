//! Staged sheathing (`StudWall::stage_sheathing`, `Domain::install_stud_wall_sheathing`; `docs/stud-wall-panel.md` §6.1):
//! gravity applied before installation goes down the framing alone, and the installed sheathing carries only what comes
//! after.
//!
//! A 4.8 m wall with a window (1200–2400) and a floor zone stands on a rigid foundation at every stud. Gravity is a
//! line load plus a 10 kN point over the header, enough to pull load through the sheathing when it is installed.
//! - Before installation its stud, jack and cripple forces equal the same wall without sheathing and with no floor-zone
//!   shear (a stick takedown, braced at the top): 1e-5 (the uninstalled fields keep 1e-6 of their stiffness).
//! - Installing changes nothing at the committed state but the uninstalled fields' own tiny shear, dropped: 1e-6.
//! - A lateral load after installation gives the same increments (top drift, stud forces, pier shear) as the always-
//!   sheathed wall under the lateral load alone (linear, so superposition is exact): 1e-8.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    ChordDistribution, Domain, Element, ElementId, ElementLoad, FloorZone, LoadSeries, Node,
    NodeId, StudWall, StudWallLayout, StudWallResponse, WallCaseLoad, WallLoad, WallOpening,
    WallShearLaw, WallSheathing, WallStation, WallStud,
};

type D2 = Domain<2, 3, NodeId, Element>;
const L: f64 = 4800.;
const H: f64 = 2440.;
const DEPTH: f64 = 254.;
const EA: f64 = 9500. * 38. * 140.;

fn layout(sheathing: bool, floor_ga: f64) -> StudWallLayout {
    let xs: Vec<f64> = (0..=12).map(|i| i as f64 * 400.).collect();
    StudWallLayout {
        length: L,
        height: H,
        studs: xs
            .iter()
            .map(|&x| WallStud {
                x,
                ea: if x == 0. || x == L { 3. * EA } else { EA },
            })
            .collect(),
        top_ea: 2. * EA,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: EA,
        bottom_ei: 9500. * 140. * 38f64.powi(3) / 12.,
        openings: vec![WallOpening {
            x0: 1200.,
            x1: 2400.,
            head_height: 2100.,
            header_ei: 2. * 9500. * 38. * 235f64.powi(3) / 12.,
            jack_ea: EA,
            door: false,
        }],
        sheathing: sheathing.then_some(WallSheathing {
            law: WallShearLaw::Linear { ga: 11_000. },
        }),
        foundation: true,
        bearing: None,
        floor: Some(FloorZone {
            depth: DEPTH,
            ga: floor_ga,
            rim_ea: 9500. * 38. * DEPTH,
        }),
    }
}

/// Gravity (case 0): 4 N/mm down the top and 10 kN at x = 1500 over the header. Lateral (case 1): 3 N/mm along the top.
fn cases() -> Vec<Vec<WallLoad>> {
    vec![
        vec![
            WallLoad::TopLine {
                x0: 0.,
                x1: L,
                wx0: 0.,
                wx1: 0.,
                wy0: -4.,
                wy1: -4.,
            },
            WallLoad::TopPoint {
                x: 1500.,
                fx: 0.,
                fy: -10_000.,
            },
        ],
        vec![WallLoad::TopLine {
            x0: 0.,
            x1: L,
            wx0: 3.,
            wx1: 3.,
            wy0: 0.,
            wy1: 0.,
        }],
    ]
}

/// The wall on its foundation (bottom corners pinned), its top-left corner held along the wall when `braced`. Returns the
/// domain, the wall and the top corners.
fn build(wall: impl FnOnce([NodeId; 4]) -> StudWall, braced: bool) -> (D2, ElementId, [NodeId; 2]) {
    let mut d = D2::new();
    let top = H + DEPTH;
    let tl = Node::new([0., top]);
    let c = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([L, 0.]).fix(0).fix(1)),
        d.add_node(Node::new([L, top])),
        d.add_node(if braced { tl.fix(0) } else { tl }),
    ];
    let id = d.add_element(Element::StudWall(wall(c)));
    (d, id, [c[2], c[3]])
}

fn new_wall(c: [NodeId; 4], l: StudWallLayout) -> StudWall {
    StudWall::new(c, Vec::<WallStation>::new(), l, cases()).unwrap()
}

fn select(
    d: &mut D2,
    wall: ElementId,
    case: u32,
    series: LoadSeries,
) -> carapace_core::model::LoadPatternId {
    let p = d.add_load_pattern(series);
    d.add_element_load(
        p,
        wall,
        ElementLoad {
            wall: Some(WallCaseLoad { case, factor: 1. }),
            ..Default::default()
        },
    );
    p
}

fn step(d: D2) -> D2 {
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 4,
        })
        .build(d);
    a.step().expect("linear step");
    a.into_domain()
}

fn members(r: &StudWallResponse) -> Vec<(f64, f64)> {
    let mut v: Vec<(f64, f64)> = r.studs.iter().map(|s| (s.x, s.bottom)).collect();
    v.extend(r.jacks.iter().map(|j| (j.x + 0.1, j.axial)));
    v.extend(r.cripples.iter().map(|c| (c.x + 0.2, c.axial)));
    v.sort_by(|a, b| a.0.total_cmp(&b.0));
    v
}

fn close(a: &[(f64, f64)], b: &[(f64, f64)], tol: f64, scale: f64, what: &str) {
    assert_eq!(a.len(), b.len(), "{what}: member count");
    for (p, q) in a.iter().zip(b) {
        assert!(
            (p.0 - q.0).abs() < 1e-9 && (p.1 - q.1).abs() <= tol * scale,
            "{what} at x = {}: {} vs {} (scale {scale})",
            p.0,
            p.1,
            q.1
        );
    }
}

#[test]
fn staged_sheathing_carries_no_gravity_and_all_later_lateral() {
    // Stick reference: no sheathing, (practically) no floor-zone shear, braced at the top.
    let (mut stick, sw, _) = build(|c| new_wall(c, layout(false, 1e-9)), true);
    select(&mut stick, sw, 0, LoadSeries::Linear { slope: 1. });
    let stick = step(stick);
    let reference = members(&stick.element_stud_wall_response(sw, 1.).unwrap());

    // The always-sheathed wall shares gravity: a stud force differs from the stick takedown by kN.
    let (mut full, fw, _) = build(|c| new_wall(c, layout(true, 10_000.)), false);
    select(&mut full, fw, 0, LoadSeries::Linear { slope: 1. });
    let full = step(full);
    let shared = members(&full.element_stud_wall_response(fw, 1.).unwrap());
    let moved = shared
        .iter()
        .zip(&reference)
        .map(|(a, b)| (a.1 - b.1).abs())
        .fold(0., f64::max);
    assert!(
        moved > 1000.,
        "the sheathed wall should share gravity (largest change {moved} N)"
    );

    // Staged: gravity before installation is the stick takedown.
    let (mut d, w, top) = build(
        |c| {
            new_wall(c, layout(true, 10_000.))
                .stage_sheathing()
                .unwrap()
        },
        false,
    );
    let g = select(&mut d, w, 0, LoadSeries::Linear { slope: 1. });
    let lat = select(&mut d, w, 1, LoadSeries::Linear { slope: 1. });
    // The lateral pattern stays at zero for the gravity step: a separate analysis below ramps it.
    d.hold_pattern_constant(lat, 0.);
    let mut d = step(d);
    let before = d.element_stud_wall_response(w, 1.).unwrap();
    let scale = 10_000.;
    close(
        &members(&before),
        &reference,
        1e-5,
        scale,
        "staged gravity vs stick",
    );

    // Installing changes nothing at the committed state.
    d.hold_pattern_constant(g, 1.);
    assert_eq!(d.install_stud_wall_sheathing(0.), 1);
    let after = d.element_stud_wall_response(w, 0.).unwrap();
    close(
        &members(&after),
        &members(&before),
        1e-6,
        scale,
        "install continuity",
    );
    assert!(
        after.piers.iter().all(|p| p.q.abs() < 1e-9),
        "installed sheathing starts unstressed"
    );
    let drift0 = top.map(|n| d.node(n).displacement[0]);

    // Lateral after installation: the increments of the always-sheathed wall under lateral alone.
    select(&mut d, w, 1, LoadSeries::Linear { slope: 1. });
    let d = step(d);
    let staged = d.element_stud_wall_response(w, 1.).unwrap();
    let (mut ref_lat, rw, rtop) = build(|c| new_wall(c, layout(true, 10_000.)), false);
    select(&mut ref_lat, rw, 1, LoadSeries::Linear { slope: 1. });
    let ref_lat = step(ref_lat);
    let lateral = ref_lat.element_stud_wall_response(rw, 1.).unwrap();
    let drift_ref = rtop.map(|n| ref_lat.node(n).displacement[0]);
    for k in 0..2 {
        let inc = d.node(top[k]).displacement[0] - drift0[k];
        assert!(
            (inc - drift_ref[k]).abs() <= 1e-8 * drift_ref[k].abs(),
            "top drift increment {inc} vs {}",
            drift_ref[k]
        );
    }
    let inc: Vec<(f64, f64)> = members(&staged)
        .iter()
        .zip(members(&after))
        .map(|(a, b)| (a.0, a.1 - b.1))
        .collect();
    close(&inc, &members(&lateral), 1e-8, 3. * L, "lateral increments");
    for (p, q) in staged.piers.iter().zip(&lateral.piers) {
        assert!(
            (p.q - q.q).abs() <= 1e-8 * q.q.abs(),
            "pier shear {} vs {}",
            p.q,
            q.q
        );
    }
}

/// The nonlinear path: the window wall bearing on two ground nodes through O86-crush hold-downs at its pier ends, with
/// O86 nail-slip sheathing. Staged gravity equals the stick takedown (no sheathing, braced): 1e-5. After installation a
/// wind load that lifts the left pier end converges (Newton), the sheathing carries the whole storey shear, and the
/// hold-down there goes into tension.
#[test]
fn staged_o86_bearing_wall_takes_gravity_by_framing_and_wind_by_sheathing() {
    use carapace_core::model::{AnchorageBearing, Material, WallBearing};
    let hold =
        Material::anchorage(15_000., 0., AnchorageBearing::O86Crush { qr: 84_323. }).unwrap();
    let make = |sheathing: bool, braced: bool, staged: bool| {
        let mut d = D2::new();
        let top = H + DEPTH;
        let tl = Node::new([0., top]);
        let c = [
            d.add_node(Node::new([0., 0.])),
            d.add_node(Node::new([L, 0.])),
            d.add_node(Node::new([L, top])),
            d.add_node(if braced { tl.fix(0) } else { tl }),
        ];
        let mut l = layout(false, if sheathing { 10_000. } else { 1e-9 });
        l.foundation = false;
        l.sheathing = sheathing.then_some(WallSheathing {
            law: WallShearLaw::O86 {
                bv: 11_000.,
                nail_spacing: 150.,
                nail_diameter: 3.25,
            },
        });
        l.bearing = Some(WallBearing {
            kv: 47_366.,
            kh: 47_366.,
            hold_down: Some(hold.clone()),
        });
        let stations: Vec<WallStation> = [0., L]
            .iter()
            .map(|&x| WallStation {
                chord: carapace_core::model::WallChord::Base,
                x,
                node: d.add_node(Node::new([x, 0.]).fix(0).fix(1)),
            })
            .collect();
        let wall = StudWall::new(c, stations, l, cases()).unwrap();
        let w = d.add_element(Element::StudWall(if staged {
            wall.stage_sheathing().unwrap()
        } else {
            wall
        }));
        (d, w, [c[2], c[3]])
    };
    let newton = |d: D2, steps: usize| {
        let mut a = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl {
                increment: 1. / steps as f64,
            })
            .algorithm(Algorithm::Newton {
                tangent: carapace_core::analysis::TangentStrategy::Current,
                line_search: None,
            })
            .test(ConvergenceTest::NormDispIncr {
                tol: 1e-10,
                max_iter: 50,
            })
            .build(d);
        for _ in 0..steps {
            a.step().expect("Newton step");
        }
        a.into_domain()
    };
    let (mut stick, sw, _) = make(false, true, false);
    select(&mut stick, sw, 0, LoadSeries::Linear { slope: 1. });
    let stick = newton(stick, 4);
    let reference = members(&stick.element_stud_wall_response(sw, 1.).unwrap());

    let (mut d, w, _) = make(true, false, true);
    let g = select(&mut d, w, 0, LoadSeries::Linear { slope: 1. });
    let mut d = newton(d, 4);
    close(
        &members(&d.element_stud_wall_response(w, 1.).unwrap()),
        &reference,
        1e-5,
        10_000.,
        "staged O86 gravity vs stick",
    );

    d.hold_pattern_constant(g, 1.);
    assert_eq!(d.install_stud_wall_sheathing(0.), 1);
    select(&mut d, w, 1, LoadSeries::Linear { slope: 1. });
    let d = newton(d, 8);
    let r = d.element_stud_wall_response(w, 1.).unwrap();
    let shear: f64 = r.piers.iter().map(|p| p.q * (p.xb - p.xa)).sum();
    assert!(
        (shear - 3. * L).abs() <= 1e-3 * 3. * L,
        "sheathing carries the storey shear: {shear} vs {}",
        3. * L
    );
    assert!(
        r.hold_downs[0].force > 0.,
        "the left pier end lifts: {}",
        r.hold_downs[0].force
    );
}
