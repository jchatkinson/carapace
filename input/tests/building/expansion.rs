//! `building::expand`: the explicit model it returns, run through the wire decoder, against the native builder
//! (`build` = `expand` + decode, with its own analysis), the Y-up and Z-up expansions of one building, and the
//! index map. Tolerances: 1e-12 for the wire round trip and the map (the same model solved twice; this needs
//! serde_json's `float_roundtrip`, whose default parser can move a coordinate by an ulp, which the 1e12 link
//! stiffnesses amplify to 1e-11), 1e-10 for Y-up against Z-up (a permuted DOF numbering, so another elimination order).
use crate::builder::{anchorage, box_walls, spec, storey, wall};
use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_input::building::*;
use carapace_input::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
    StudWallBaseComponent, StudWallQuantity,
};
use carapace_input::input_v1::tables::ElementKind;
use carapace_input::input_v1::{decode, CarapaceInputV1, DecodeError};

/// Serializes `expansion.input` to JSON and back (the wire a client sees), appends one linear static stage with
/// `recorders`, decodes, runs, and returns each recorder's one sample.
fn run_wire(expansion: &Expansion, recorders: Vec<RecorderSpec>) -> Vec<f64> {
    let text = serde_json::to_string(expansion).unwrap();
    let mut input = serde_json::from_str::<Expansion>(&text).unwrap().input;
    let count = recorders.len();
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "expanded".into(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1. },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-8,
                max_iter: 8,
            }),
            hold_patterns_after: vec![],
        }],
        recorders,
    };
    let mut session = decode(input).unwrap();
    let out = session.advance(10);
    assert!(out.done && out.error.is_none(), "{out:?}");
    (0..count)
        .map(|i| {
            out.recorder_batches
                .iter()
                .find(|b| b.recorder_index == i)
                .unwrap()
                .samples[0]
                .1
        })
        .collect()
}

fn pier(element: u32, quantity: StudWallQuantity) -> RecorderSpec {
    RecorderSpec::StudWallPier {
        element_kind: ElementKind::StudWall,
        element_index: element,
        pier: 0,
        quantity,
    }
}
fn base(element: u32, component: StudWallBaseComponent) -> RecorderSpec {
    RecorderSpec::StudWallBase {
        element_kind: ElementKind::StudWall,
        element_index: element,
        component,
    }
}
/// Master-node displacement relative to the level below, in model DOF `dof`.
fn drift(map: &ExpansionMap, storey: usize, dof: u8) -> RecorderSpec {
    RecorderSpec::NodeDifference {
        node: map.storeys[storey].master.unwrap(),
        other_node: storey
            .checked_sub(1)
            .map(|s| map.storeys[s].master.unwrap()),
        dof,
    }
}

/// Per wall (storey order, spec order): pier q, base moment of the pier, and the three base reactions.
fn wall_recorders(map: &ExpansionMap) -> Vec<RecorderSpec> {
    map.storeys
        .iter()
        .flat_map(|s| &s.walls)
        .flat_map(|w| {
            [
                pier(w.element, StudWallQuantity::Q),
                pier(w.element, StudWallQuantity::MomentBase),
                base(w.element, StudWallBaseComponent::Shear),
                base(w.element, StudWallBaseComponent::Axial),
                base(w.element, StudWallBaseComponent::Moment),
            ]
        })
        .collect()
}

fn near(got: f64, want: f64, scale: f64, tol: f64, label: &str) {
    assert!(
        (got - want).abs() <= tol * scale,
        "{label}: {got} vs {want} (scale {scale:e}, tol {tol:e})"
    );
}

/// The native builder's linear solution of `spec` (its own `AnalysisBuilder`, not the wire session).
fn solve_native(b: &Building) -> carapace_core::model::Domain3 {
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-8,
            max_iter: 8,
        })
        .build(b.domain.clone());
    a.step().expect("linear building solve");
    a.into_domain()
}

/// A two-storey rigid box with a floor zone on the lower storey, an eccentric X force on each level (so the
/// diaphragms rotate), a Z-span gravity floor on the lower storey (slot 1) and combination factors `[1, 1.25, ..]`.
fn eccentric_box() -> BuildingSpec {
    let mut lower = storey(box_walls(false), 254.);
    lower.lateral.push(LateralLoad {
        slot: 0,
        fx: 10000.,
        fz: 2000.,
        at: [2500., 3000.],
    });
    lower.floors.push(FloorSpec {
        bounds: [0., 4000., 0., 4000.],
        supports: vec!["front".into(), "back".into()],
        span_axis: SpanAxis::Z,
        loads: vec![FloorLoadSpec::Area {
            slot: 1,
            w: 0.002,
            region: None,
        }],
        continuity: FloorContinuitySpec::Simple,
    });
    let mut upper = storey(box_walls(false), 0.);
    upper.lateral.push(LateralLoad {
        slot: 0,
        fx: 6000.,
        fz: -3000.,
        at: [1500., 2200.],
    });
    let mut s = spec(vec![lower, upper], DiaphragmKind::Rigid);
    s.combination = [1., 1.25, 1., 1., 1., 1., 1., 1.];
    s
}

/// Round trip, rigid: `expand` -> JSON -> `Expansion` -> decode -> solve gives the native builder's wall pier q,
/// pier base moment, base reactions and storey drift (X, Z, rotation) to 1e-12 of the largest value of that kind.
#[test]
fn round_trip_rigid_matches_native_builder() {
    let s = eccentric_box();
    let ex = expand(&s, &ExpandOptions::default()).unwrap();
    let mut recorders = wall_recorders(&ex.map);
    for storey in 0..2 {
        for dof in [0, 2, 4] {
            recorders.push(drift(&ex.map, storey, dof));
        }
    }
    let got = run_wire(&ex, recorders);
    let b = build(&s).unwrap();
    let d = solve_native(&b);
    let mut want = Vec::new();
    for (si, st) in s.storeys.iter().enumerate() {
        for w in &st.walls {
            let r = b.wall_response(&d, si, &w.name).unwrap();
            want.extend([
                r.piers[0].q,
                r.piers[0].moment_base,
                r.base[0],
                r.base[1],
                r.base[2],
            ]);
        }
    }
    for si in 0..2 {
        let StoreyDrift::Rigid { ux, uz, ry } = b.storey_drift(&d, si).unwrap() else {
            panic!("rigid level")
        };
        want.extend([ux, uz, ry]);
    }
    assert_eq!(got.len(), want.len());
    for kind in 0..5 {
        let scale = (0..8).map(|w| want[5 * w + kind].abs()).fold(0., f64::max);
        for w in 0..8 {
            near(
                got[5 * w + kind],
                want[5 * w + kind],
                scale,
                1e-12,
                "wall quantity",
            );
        }
    }
    let scale = want[40..].iter().fold(0f64, |m, v| m.max(v.abs()));
    for i in 40..want.len() {
        near(got[i], want[i], scale, 1e-12, "drift");
    }
    assert!(
        want[45].abs() > 0. && want[41].abs() > 0.,
        "load is eccentric"
    );
}

/// Round trip, flexible: three X wall lines under a 12 kN X load (tributary shares) plus a Z-span gravity floor
/// on two of them: wall base reactions and a mid-span top node's X displacement equal the native builder's, 1e-12.
#[test]
fn round_trip_flexible_matches_native_builder() {
    let mut st = storey(
        vec![
            wall("a", [0., 0.], [1., 0.], false),
            wall("b", [0., 2000.], [1., 0.], false),
            wall("c", [0., 4000.], [1., 0.], false),
        ],
        0.,
    );
    st.lateral.push(LateralLoad {
        slot: 0,
        fx: 12000.,
        fz: 0.,
        at: [0., 0.],
    });
    st.floors.push(FloorSpec {
        bounds: [0., 4000., 0., 4000.],
        supports: vec!["a".into(), "c".into()],
        span_axis: SpanAxis::Z,
        loads: vec![FloorLoadSpec::Area {
            slot: 1,
            w: 0.003,
            region: None,
        }],
        continuity: FloorContinuitySpec::Simple,
    });
    let s = spec(vec![st], DiaphragmKind::Flexible);
    let ex = expand(&s, &ExpandOptions::default()).unwrap();
    assert_eq!(ex.map.storeys[0].master, None);
    assert!(ex.input.rigid_diaphragms.retained.is_empty());
    let mut recorders = wall_recorders(&ex.map);
    let mid = |w: usize| ex.map.storeys[0].walls[w].top_nodes[1].node;
    recorders.extend((0..3).map(|w| RecorderSpec::NodeDisp {
        node: mid(w),
        dof: 0,
    }));
    let got = run_wire(&ex, recorders);
    let b = build(&s).unwrap();
    let d = solve_native(&b);
    let mut want = Vec::new();
    for w in &s.storeys[0].walls {
        let r = b.wall_response(&d, 0, &w.name).unwrap();
        want.extend([
            r.piers[0].q,
            r.piers[0].moment_base,
            r.base[0],
            r.base[1],
            r.base[2],
        ]);
    }
    for w in 0..3 {
        want.push(
            d.node(b.walls[0][s.storeys[0].walls[w].name.as_str()].top_nodes[1].1)
                .displacement[0],
        );
    }
    for kind in 0..5 {
        let scale = (0..3).map(|w| want[5 * w + kind].abs()).fold(0., f64::max);
        for w in 0..3 {
            near(
                got[5 * w + kind],
                want[5 * w + kind],
                scale,
                1e-12,
                "wall quantity",
            );
        }
    }
    let scale = want[15..].iter().fold(0f64, |m, v| m.max(v.abs()));
    for i in 15..18 {
        near(got[i], want[i], scale, 1e-12, "top-node ux");
    }
    // The lines carry 3/6/3 kN of the lateral load.
    near(got[2].abs(), 3000., 12000., 1e-9, "line a shear");
    near(got[7].abs(), 6000., 12000., 1e-9, "line b shear");
}

/// The torsion case of `builder::rigid_torsion_load` (four separate equal walls, 10 kN in X at z = 3500) and a
/// flexible diaphragm.
fn torsion_spec(diaphragm: DiaphragmKind) -> BuildingSpec {
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
    spec(vec![s], diaphragm)
}

/// Z-up equals Y-up: the same building expanded with `up = Z` (plan `(x, z)` -> model `(x, y)`) gives the same
/// node positions (permuted), wall pier q, base moments and reactions, and storey drift in plan X and Z, to 1e-10
/// of the largest value. Cases: the eccentric two-storey box (floor zone, gravity floor, rigid), the torsion case
/// (rigid) and a flexible diaphragm. The rotation about the vertical has the opposite sign: plan `(x, z)` with `y`
/// up is a left-handed plan frame seen from above, `(x, y)` with `z` up a right-handed one, so `rz = -ry` for the
/// same physical motion (see `docs/building.md`).
#[test]
fn z_up_matches_y_up() {
    let cases = [
        eccentric_box(),
        torsion_spec(DiaphragmKind::Rigid),
        torsion_spec(DiaphragmKind::Flexible),
    ];
    for (case, s) in cases.iter().enumerate() {
        let y = expand(s, &ExpandOptions { up: UpAxis::Y }).unwrap();
        let z = expand(s, &ExpandOptions { up: UpAxis::Z }).unwrap();
        // Same rows, permuted coordinates.
        assert_eq!(y.input.nodes.fixed.len(), z.input.nodes.fixed.len());
        for (a, b) in y
            .input
            .nodes
            .coords
            .chunks(3)
            .zip(z.input.nodes.coords.chunks(3))
        {
            assert_eq!([a[0], a[2], a[1]], [b[0], b[1], b[2]], "case {case}");
        }
        assert_eq!(
            serde_json::to_value(&y.map).unwrap(),
            serde_json::to_value(&z.map).unwrap()
        );
        let rigid = s.diaphragm == DiaphragmKind::Rigid;
        let extra = |ex: &Expansion, dofs: [u8; 3]| {
            let mut r = wall_recorders(&ex.map);
            if rigid {
                for storey in 0..s.storeys.len() {
                    r.extend(dofs.iter().map(|&d| drift(&ex.map, storey, d)));
                }
            }
            r
        };
        let (ry, rz) = (extra(&y, [0, 2, 4]), extra(&z, [0, 1, 5]));
        let (gy, gz) = (run_wire(&y, ry), run_wire(&z, rz));
        let walls = 5 * s.storeys.iter().map(|s| s.walls.len()).sum::<usize>();
        let scale = gy.iter().fold(0f64, |m, v| m.max(v.abs()));
        for i in 0..gy.len() {
            let is_rotation = rigid && i >= walls && (i - walls) % 3 == 2;
            let want = if is_rotation { -gy[i] } else { gy[i] };
            near(gz[i], want, scale, 1e-10, &format!("case {case} value {i}"));
        }
        if rigid {
            // Torsion sign: the Y-up master rotates positively about +Y (`rigid_torsion_load`), the Z-up one
            // negatively about +Z, and it is a real rotation (not noise).
            let r = gy[walls + 2];
            assert!(
                r.abs() > 1e-9 && r.signum() == -gz[walls + 2].signum(),
                "case {case}: ry {r} vs rz {}",
                gz[walls + 2]
            );
            if case == 1 {
                assert!(r > 0. && gz[walls + 2] < 0., "torsion sign");
            }
        }
    }
}

/// The map addresses the right rows: a `studWallPier` recorder at `map.storeys[s].walls[i].element` equals
/// `Building::wall_response`, the pier-link `zeroLengths` rows carry the native hold-down forces, the level
/// masters' differences are `Building::storey_drift`, and the table rows sit where the map says (1e-12).
#[test]
fn map_rows_equal_native_handles() {
    let s = eccentric_box();
    let ex = expand(&s, &ExpandOptions::default()).unwrap();
    let b = build(&s).unwrap();
    let d = solve_native(&b);
    let mut recorders = Vec::new();
    let mut want = Vec::new();
    for (si, sm) in ex.map.storeys.iter().enumerate() {
        for wm in &sm.walls {
            let r = b.wall_response(&d, si, &wm.name).unwrap();
            recorders.push(pier(wm.element, StudWallQuantity::Q));
            want.push(r.piers[0].q);
            recorders.push(base(wm.element, StudWallBaseComponent::Shear));
            want.push(r.base[0]);
            let holddowns = b.holddown_forces(&d, si, &wm.name).unwrap();
            for link in &wm.pier_links {
                let native = holddowns
                    .iter()
                    .find(|h| (h.x - link.x).abs() < 1e-7)
                    .unwrap();
                recorders.push(RecorderSpec::ElementForce {
                    element_kind: ElementKind::ZeroLength,
                    element_index: link.element,
                    component: 1,
                });
                want.push(d.element_local_force(native.link).as_slice()[1]);
            }
        }
        let StoreyDrift::Rigid { ux, uz, ry } = b.storey_drift(&d, si).unwrap() else {
            panic!("rigid level")
        };
        for (dof, v) in [(0, ux), (2, uz), (4, ry)] {
            recorders.push(drift(&ex.map, si, dof));
            want.push(v);
        }
    }
    let got = run_wire(&ex, recorders);
    let scale = want.iter().fold(0f64, |m, v| m.max(v.abs()));
    for (i, (g, w)) in got.iter().zip(&want).enumerate() {
        near(*g, *w, scale, 1e-12, &format!("recorder {i}"));
    }
    // The rows are where the map says: wall corners, anchorage links on the pier-end nodes, level nodes.
    let nodes = |i: u32| {
        let c = &ex.input.nodes.coords;
        [
            c[3 * i as usize],
            c[3 * i as usize + 1],
            c[3 * i as usize + 2],
        ]
    };
    for (si, sm) in ex.map.storeys.iter().enumerate() {
        assert_eq!(sm.walls.len(), s.storeys[si].walls.len());
        assert_eq!(sm.links_below.is_empty(), si == 0);
        for (wm, w) in sm.walls.iter().zip(&s.storeys[si].walls) {
            assert_eq!(wm.name, w.name);
            let row = &ex.input.stud_walls[wm.element as usize];
            assert_eq!(row.corners, wm.corners);
            let c1 = nodes(wm.corners[1]);
            assert_eq!(c1[0], w.start[0] + w.direction[0] * w.layout.length);
            assert_eq!(c1[2], w.start[1] + w.direction[1] * w.layout.length);
            for link in &wm.pier_links {
                let upper = ex.input.zero_lengths.node_j[link.element as usize];
                let bottom = wm.bottom_nodes.iter().find(|n| n.x == link.x).unwrap();
                assert_eq!(upper, bottom.node);
            }
        }
        assert_eq!(sm.nodes.len(), 12);
    }
    assert_eq!(ex.map.slot_patterns[0], Some(0));
    assert_eq!(ex.map.slot_patterns[1], Some(1));
    assert!(ex.map.slot_patterns[2..].iter().all(Option::is_none));
}

/// A Z-span floor (2 kPa) on a 4 m x 4 m box, expanded to a `floorLoads` row on the front and back walls: through
/// the wire decoder the base axial of those walls is the takedown `w Lx Lz / 2 = 16 kN` and the side walls carry
/// none (statics, 1e-9).
#[test]
fn floor_load_rows_give_the_takedown() {
    let mut s = storey(box_walls(false), 0.);
    s.floors.push(FloorSpec {
        bounds: [0., 4000., 0., 4000.],
        supports: vec!["front".into(), "back".into()],
        span_axis: SpanAxis::Z,
        loads: vec![FloorLoadSpec::Area {
            slot: 0,
            w: 0.002,
            region: None,
        }],
        continuity: FloorContinuitySpec::Simple,
    });
    let ex = expand(
        &spec(vec![s], DiaphragmKind::Rigid),
        &ExpandOptions::default(),
    )
    .unwrap();
    assert_eq!(ex.input.floor_loads.len(), 1);
    assert_eq!(ex.input.floor_loads[0].pattern, 0);
    let recorders = (0..4)
        .map(|w| base(w, StudWallBaseComponent::Axial))
        .collect();
    let got = run_wire(&ex, recorders);
    // Rows: front, back, left, right.
    near(got[0], 16_000., 16_000., 1e-9, "front takedown");
    near(got[1], 16_000., 16_000., 1e-9, "back takedown");
    assert!(
        got[2].abs() < 1e-6 && got[3].abs() < 1e-6,
        "side walls {got:?}"
    );
}

/// The expansion is plain tables: combination factors become the pattern scale factors of the used slots only, every
/// load row is at stage 0 of an empty sequence (which the decoder refuses until a stage exists), diaphragm normals
/// follow the up axis, ground nodes are fixed, and bad slots or factors are `BuildError`s.
#[test]
fn expansion_is_explicit_tables() {
    let s = eccentric_box();
    for (up, normal) in [
        (UpAxis::Y, carapace_input::input_v1::tables::Axis3Spec::Y),
        (UpAxis::Z, carapace_input::input_v1::tables::Axis3Spec::Z),
    ] {
        let ex = expand(&s, &ExpandOptions { up }).unwrap();
        let i = &ex.input;
        assert!(i.sequence.stages.is_empty() && i.sequence.recorders.is_empty());
        assert_eq!(i.load_patterns.scale_factor, vec![1., 1.25]);
        assert_eq!(i.stud_walls.len(), 8);
        assert_eq!(i.rigid_diaphragms.normal, vec![normal; 2]);
        assert!(i.nodal_loads.stage.iter().all(|&s| s == 0));
        assert_eq!(i.nodal_loads.node.len(), i.nodal_loads.value.len());
        assert_eq!(i.floor_loads.len(), 1);
        assert!(i.stud_walls[0].layout.foundation && !i.stud_walls[4].layout.foundation);
        assert_eq!(i.stud_walls[0].layout.floor.unwrap().depth, 254.);
        // Eight ground nodes (two pier ends per lower wall), each fixed in all three translations.
        assert_eq!(i.nodes.fixed.iter().filter(|&&f| f == 0b111).count(), 8);
        assert!(matches!(
            decode(i.clone()),
            Err(DecodeError::UnknownStageIndex { .. })
        ));
    }
    let mut bad = eccentric_box();
    bad.storeys[0].lateral[0].slot = 8;
    assert!(expand(&bad, &ExpandOptions::default())
        .unwrap_err()
        .to_string()
        .contains("lateral"));
    let mut bad = eccentric_box();
    bad.combination[3] = f64::NAN;
    assert!(expand(&bad, &ExpandOptions::default())
        .unwrap_err()
        .to_string()
        .contains("combination"));
    let mut bad = eccentric_box();
    bad.foundation_anchorage = AnchorageSpec {
        tension: -1.,
        ..anchorage()
    };
    assert!(expand(&bad, &ExpandOptions::default()).is_err());
}

/// A wall load in a spec slot becomes a `wallLoad` row of that slot's pattern: a 1 N/mm top line on a single wall
/// puts `L` of vertical reaction on its base, scaled by the combination factor (2.5: 10 kN for a 4 m wall; 1e-9).
#[test]
fn wall_loads_become_pattern_rows() {
    let mut w = wall("only", [0., 0.], [1., 0.], false);
    w.loads.push(SlotWallLoad {
        slot: 3,
        load: WallLoadSpec::TopLine {
            x0: 0.,
            x1: 4000.,
            wx0: 0.,
            wx1: 0.,
            wy0: -1.,
            wy1: -1.,
        },
    });
    let mut s = spec(vec![storey(vec![w], 0.)], DiaphragmKind::Flexible);
    s.combination[3] = 2.5;
    let ex = expand(&s, &ExpandOptions::default()).unwrap();
    assert_eq!(ex.map.slot_patterns[3], Some(0));
    assert_eq!(ex.input.element_loads.element_index, vec![0]);
    assert_eq!(ex.input.load_patterns.scale_factor, vec![2.5]);
    let got = run_wire(&ex, vec![base(0, StudWallBaseComponent::Axial)]);
    near(got[0], 10_000., 10_000., 1e-9, "base axial");
}

/// The unused-field guard from the other side: a spec with the legacy `building` input in a model JSON is an
/// unknown field, and `CarapaceInputV1` still parses what `expand` printed.
#[test]
fn expanded_input_parses_and_legacy_field_does_not() {
    let ex = expand(&eccentric_box(), &ExpandOptions::default()).unwrap();
    let mut value = serde_json::to_value(&ex.input).unwrap();
    assert!(serde_json::from_value::<CarapaceInputV1>(value.clone()).is_ok());
    value["building"] = serde_json::json!({});
    assert!(serde_json::from_value::<CarapaceInputV1>(value).is_err());
}
