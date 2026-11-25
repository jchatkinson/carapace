//! Wire tables for stud walls (with their floor zone), diaphragms, the anchorage material and floor loads:
//! each decoded model equals the same model built directly through `carapace-core`
//! (window wall with floor loads, anchorage drift: 1e-12 or 1e-10 as stated per test), specs
//! round-trip through serde, and every malformed row or misuse is a `DecodeError`.
use crate::common::empty_input;
use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    AnchorageBearing, ChordDistribution, Domain, Element, ElementLoad, FloorBounds,
    FloorContinuity, FloorLoad, FloorSupport, LoadSeries, Material, Node, StudWall, StudWallLayout,
    WallCaseLoad, WallLoad, WallShearLaw, WallSheathing, WallStud, ZeroLength,
};
use carapace_input::building::{AnchorageBearingSpec, BuildingSpec, FloorContinuitySpec};
use carapace_input::building::{
    ChordDistributionSpec, FloorZoneSpec, WallLayoutSpec, WallLoadSpec, WallShearLawSpec,
    WallSheathingSpec, WallStudSpec,
};
use carapace_input::input_v1::{
    decode,
    materials::MaterialSpec,
    sequence::{
        AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
        TangentStrategySpec,
    },
    tables::*,
};

fn layout() -> WallLayoutSpec {
    WallLayoutSpec {
        length: 1000.,
        height: 1000.,
        studs: vec![
            WallStudSpec { x: 0., ea: 1e8 },
            WallStudSpec { x: 1000., ea: 1e8 },
        ],
        top_ea: 1e9,
        top_ei: Some(1e9),
        top_distribution: ChordDistributionSpec::Beam,
        bottom_ea: 1e9,
        bottom_ei: 1e9,
        openings: vec![],
        sheathing: Some(WallSheathingSpec {
            law: WallShearLawSpec::Linear { ga: 1e6 },
        }),
        foundation: false,
        floor: None,
    }
}
fn input_wall() -> carapace_input::input_v1::CarapaceInputV1 {
    let mut i = empty_input(2);
    i.nodes = NodeTable {
        coords: vec![0., 0., 1000., 0., 1000., 1000., 0., 1000.],
        fixed: vec![7, 7, 2, 2],
        mass_node_index: vec![],
        mass: vec![],
    };
    i.stud_walls = vec![StudWallSpec {
        corners: [0, 1, 2, 3],
        stations: vec![],
        layout: layout(),
    }];
    i.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1. }],
        scale_factor: vec![1.],
    };
    i.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::StudWall],
        element_index: vec![0],
        load: vec![ElementLoadSpec::WallLoad {
            load: WallLoadSpec::TopPoint {
                x: 500.,
                fx: 100.,
                fy: 0.,
            },
        }],
        stage: vec![0],
    };
    i.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "s".into(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1. },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-10,
                max_iter: 10,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 3, dof: 0 }],
    };
    i
}
fn core_wall_values() -> (f64, f64) {
    let mut d: Domain = Domain::new();
    let n0 = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let n1 = d.add_node(Node::new([1000., 0.]).fix(0).fix(1).fix(2));
    let n2 = d.add_node(Node::new([1000., 1000.]).fix(1));
    let n3 = d.add_node(Node::new([0., 1000.]).fix(1));
    let law = StudWallLayout {
        length: 1000.,
        height: 1000.,
        studs: vec![WallStud { x: 0., ea: 1e8 }, WallStud { x: 1000., ea: 1e8 }],
        top_ea: 1e9,
        top_ei: Some(1e9),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: 1e9,
        bottom_ei: 1e9,
        openings: vec![],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 1e6 },
        }),
        foundation: false,
        floor: None,
    };
    let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    loads[0].push(WallLoad::TopPoint {
        x: 500.,
        fx: 100.,
        fy: 0.,
    });
    let id = d.add_element(Element::StudWall(
        StudWall::new([n0, n1, n2, n3], vec![], law, loads).unwrap(),
    ));
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_element_load(
        p,
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
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 10,
        })
        .build(d);
    a.step().unwrap();
    (
        a.domain().node(n3).displacement[0],
        a.domain().element_stud_wall_response(id, 1.).unwrap().piers[0].q,
    )
}
#[test]
fn stud_wall_wire_solves_like_direct_core_constructor_and_roundtrips() {
    let input = input_wall();
    let value = serde_json::to_value(&input).unwrap();
    let decoded: carapace_input::input_v1::CarapaceInputV1 =
        serde_json::from_value(value.clone()).unwrap();
    assert_eq!(serde_json::to_value(decoded.clone()).unwrap(), value);
    let mut s = decode(decoded).unwrap();
    let o = s.advance(10);
    assert!(o.done && o.error.is_none(), "{o:?}");
    let got = o.recorder_batches[0].samples[0].1;
    assert!(
        (got - core_wall_values().0).abs() < 1e-12,
        "got={got} expected={}",
        core_wall_values().0
    );
}
#[test]
fn stud_wall_pier_recorder_matches_core_and_validates_pier_index() {
    let mut input = input_wall();
    input.sequence.recorders = vec![RecorderSpec::StudWallPier {
        element_kind: ElementKind::StudWall,
        element_index: 0,
        pier: 0,
        quantity: carapace_input::input_v1::sequence::StudWallQuantity::Q,
    }];
    let mut session = decode(input.clone()).unwrap();
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none());
    let value = outcome.recorder_batches[0].samples[0].1;
    assert!((value - core_wall_values().1).abs() < 1e-12);
    input.sequence.recorders = vec![RecorderSpec::StudWallPier {
        element_kind: ElementKind::StudWall,
        element_index: 0,
        pier: 1,
        quantity: carapace_input::input_v1::sequence::StudWallQuantity::Q,
    }];
    assert!(matches!(
        decode(input),
        Err(carapace_input::input_v1::DecodeError::InvalidRecorderPier {
            pier: 1,
            count: 1,
            ..
        })
    ));
}

/// A wall's `floor` zone decodes (exterior top nodes at `H + depth`), serializes back unchanged, solves like the
/// same wall built through `carapace-core`, and a zone with a non-positive depth is a row error. Reference: the
/// direct core construction; top displacement agrees to 1e-12 relative.
#[test]
fn stud_wall_floor_zone_wire_matches_direct_core() {
    let zone = FloorZoneSpec {
        depth: 200.,
        ga: 5e5,
        rim_ea: 1e8,
    };
    let mut input = input_wall();
    input.nodes.coords = vec![0., 0., 1000., 0., 1000., 1200., 0., 1200.];
    input.nodes.fixed = vec![7, 7, 0, 0];
    input.stud_walls[0].layout.floor = Some(zone);
    let value = serde_json::to_value(&input).unwrap();
    assert_eq!(value["studWalls"][0]["layout"]["floor"]["rimEa"], 1e8);
    let decoded: carapace_input::input_v1::CarapaceInputV1 =
        serde_json::from_value(value.clone()).unwrap();
    assert_eq!(serde_json::to_value(decoded.clone()).unwrap(), value);
    let mut s = decode(decoded).unwrap();
    let o = s.advance(10);
    assert!(o.done && o.error.is_none(), "{o:?}");
    let got = o.recorder_batches[0].samples[0].1;

    let mut d: Domain = Domain::new();
    let n0 = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let n1 = d.add_node(Node::new([1000., 0.]).fix(0).fix(1).fix(2));
    let n2 = d.add_node(Node::new([1000., 1200.]));
    let n3 = d.add_node(Node::new([0., 1200.]));
    let mut l = StudWallLayout {
        length: 1000.,
        height: 1000.,
        studs: vec![WallStud { x: 0., ea: 1e8 }, WallStud { x: 1000., ea: 1e8 }],
        top_ea: 1e9,
        top_ei: Some(1e9),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: 1e9,
        bottom_ei: 1e9,
        openings: vec![],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 1e6 },
        }),
        foundation: false,
        floor: None,
    };
    l.floor = Some(carapace_core::model::FloorZone {
        depth: 200.,
        ga: 5e5,
        rim_ea: 1e8,
    });
    let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    loads[0].push(WallLoad::TopPoint {
        x: 500.,
        fx: 100.,
        fy: 0.,
    });
    let id = d.add_element(Element::StudWall(
        StudWall::new([n0, n1, n2, n3], vec![], l, loads).unwrap(),
    ));
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_element_load(
        p,
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
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 10,
        })
        .build(d);
    a.step().unwrap();
    near_relative(got, a.domain().node(n3).displacement[0], 1e-12, "top ux");
    // The load acts on the exterior row: the base reactions' moment is V (H + d).
    let r = a.domain().element_stud_wall_response(id, 1.).unwrap();
    near_relative(r.base[2].abs(), 100. * 1200., 1e-9, "base moment");

    let mut bad = input_wall();
    bad.nodes.coords = vec![0., 0., 1000., 0., 1000., 1200., 0., 1200.];
    bad.stud_walls[0].layout.floor = Some(FloorZoneSpec { depth: 0., ..zone });
    assert!(matches!(
        decode(bad),
        Err(carapace_input::input_v1::DecodeError::InvalidMacroElement {
            table: "stud_walls",
            ..
        })
    ));
}
#[test]
fn macro_table_errors_are_reported() {
    let mut bad = input_wall();
    bad.stud_walls[0].corners[0] = 99;
    assert!(matches!(
        decode(bad),
        Err(carapace_input::input_v1::DecodeError::UnknownNodeIndex { .. })
    ));
    let mut bad = input_wall();
    bad.element_loads.pattern[0] = 3;
    assert!(matches!(
        decode(bad),
        Err(carapace_input::input_v1::DecodeError::UnknownPatternIndex { .. })
    ));
    let mut bad = input_wall();
    bad.stud_walls = vec![StudWallSpec {
        corners: [0; 4],
        stations: vec![],
        layout: layout(),
    }];
    assert!(matches!(
        decode(bad),
        Err(carapace_input::input_v1::DecodeError::InvalidMacroElement {
            table: "stud_walls",
            ..
        })
    ));
}

fn near_relative(got: f64, expected: f64, tolerance: f64, label: &str) {
    let relative = (got - expected).abs() / got.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
    assert!(
        relative <= tolerance,
        "{label}: {got} vs {expected} (relative {relative:e})"
    );
}

fn window_layout(foundation: bool) -> WallLayoutSpec {
    WallLayoutSpec {
        length: 4000.,
        height: 2440.,
        studs: (0..=10)
            .map(|i| WallStudSpec {
                x: 400. * i as f64,
                ea: 1e8,
            })
            .collect(),
        top_ea: 1e12,
        top_ei: None,
        top_distribution: ChordDistributionSpec::Direct,
        bottom_ea: 1e12,
        bottom_ei: 1e12,
        openings: vec![carapace_input::building::WallOpeningSpec {
            x0: 1200.,
            x1: 2800.,
            head_height: 2000.,
            header_ei: 1e12,
            jack_ea: 1e8,
            door: false,
        }],
        sheathing: Some(WallSheathingSpec {
            law: WallShearLawSpec::Linear { ga: 1e7 },
        }),
        foundation,
        floor: None,
    }
}

/// Decodes a two-slot floor transfer onto a window wall with 400 mm studs and Direct top-plate distribution.
/// Reference: the same `transfer_floor`/`wall_loads`/`StudWall::new` core construction; pier q and base response
/// agree to 1e-12 relative. Direct header reactions also match the 4.8 MN hand takedown to 1e-12 relative.
#[test]
fn floor_loads_window_wall_matches_direct_core() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![
            0., 0., 4000., 0., 4000., 2440., 0., 2440., 0., 0., 4000., 0., 4000., 2440., 0., 2440.,
        ],
        fixed: vec![7, 7, 2, 2, 7, 7, 2, 2],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.stud_walls = [0, 1]
        .map(|wall| StudWallSpec {
            corners: [wall * 4, wall * 4 + 1, wall * 4 + 2, wall * 4 + 3],
            stations: vec![],
            layout: window_layout(wall == 0),
        })
        .to_vec();
    input.floor_loads = vec![FloorTransferSpec {
        bounds: [0., 4000., 0., 4000.],
        continuity: FloorContinuitySpec::Simple,
        supports: [0., 4000.]
            .into_iter()
            .enumerate()
            .map(|(i, s)| FloorTransferSupportSpec {
                s,
                t0: 0.,
                t1: 4000.,
                wall: WallTargetSpec::StudWall { index: i as u32 },
                x_at_t0: 0.,
                reversed: false,
            })
            .collect(),
        pattern: 0,
        loads: vec![FloorTransferLoadSpec::Area {
            w: 2.,
            region: None,
        }],
    }];
    // Live load is a second floor row on pattern 1 (the rows share their supports).
    let live = FloorTransferSpec {
        pattern: 1,
        loads: vec![FloorTransferLoadSpec::Area {
            w: 1.,
            region: None,
        }],
        ..input.floor_loads[0].clone()
    };
    input.floor_loads.push(live);
    input.load_patterns = LoadPatternTable {
        series: vec![
            TimeSeriesSpec::Linear { slope: 1. },
            TimeSeriesSpec::Linear { slope: 1. },
        ],
        scale_factor: vec![1.25, 1.5],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "gravity".into(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1. },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-8,
                max_iter: 8,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::StudWallPier {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                pier: 0,
                quantity: carapace_input::input_v1::sequence::StudWallQuantity::Q,
            },
            RecorderSpec::StudWallPier {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                pier: 0,
                quantity: carapace_input::input_v1::sequence::StudWallQuantity::MomentBase,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: 0,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: 1,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: 2,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: 3,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: 4,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: 5,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: 6,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: 7,
            },
        ],
    };
    let mut session = decode(input).unwrap();
    let output = session.advance(10);
    assert!(output.done && output.error.is_none(), "{output:?}");

    let make_direct_wall = |d: &mut Domain, offset: f64, layout: &WallLayoutSpec, loads| {
        let nodes = [
            d.add_node(Node::new([offset, 0.]).fix(0).fix(1).fix(2)),
            d.add_node(Node::new([offset + 4000., 0.]).fix(0).fix(1).fix(2)),
            d.add_node(Node::new([offset + 4000., 2440.]).fix(1)),
            d.add_node(Node::new([offset, 2440.]).fix(1)),
        ];
        let core_layout = StudWallLayout {
            length: layout.length,
            height: layout.height,
            studs: layout
                .studs
                .iter()
                .map(|s| WallStud { x: s.x, ea: s.ea })
                .collect(),
            top_ea: layout.top_ea,
            top_ei: layout.top_ei,
            top_distribution: ChordDistribution::Direct,
            bottom_ea: layout.bottom_ea,
            bottom_ei: layout.bottom_ei,
            openings: layout
                .openings
                .iter()
                .map(|o| carapace_core::model::WallOpening {
                    x0: o.x0,
                    x1: o.x1,
                    head_height: o.head_height,
                    header_ei: o.header_ei,
                    jack_ea: o.jack_ea,
                    door: o.door,
                })
                .collect(),
            sheathing: Some(WallSheathing {
                law: WallShearLaw::Linear { ga: 1e7 },
            }),
            foundation: layout.foundation,
            floor: None,
        };
        d.add_element(Element::StudWall(
            StudWall::new(nodes, vec![], core_layout, loads).unwrap(),
        ))
    };
    let reactions = carapace_core::model::transfer_floor(
        FloorBounds {
            s0: 0.,
            s1: 4000.,
            t0: 0.,
            t1: 4000.,
        },
        &[
            FloorSupport {
                s: 0.,
                t0: 0.,
                t1: 4000.,
            },
            FloorSupport {
                s: 4000.,
                t0: 0.,
                t1: 4000.,
            },
        ],
        &[
            FloorLoad::Area {
                slot: 0,
                w: 2.,
                region: None,
            },
            FloorLoad::Area {
                slot: 1,
                w: 1.,
                region: None,
            },
        ],
        FloorContinuity::Simple,
    )
    .unwrap();
    let mut d = Domain::new();
    let mut loads_by_support: Vec<Vec<Vec<WallLoad>>> = vec![vec![vec![]; 2]; 2];
    for (support_index, support) in reactions.per_support.iter().enumerate() {
        for (slot, loads) in support.wall_loads(0., false) {
            loads_by_support[support_index][slot].extend(loads);
        }
    }
    let id = make_direct_wall(&mut d, 0., &window_layout(true), loads_by_support.remove(0));
    make_direct_wall(
        &mut d,
        10000.,
        &window_layout(false),
        loads_by_support.remove(0),
    );
    for (case, factor) in [(0, 1.25), (1, 1.5)] {
        let p = d.add_load_pattern_scaled(LoadSeries::Linear { slope: 1. }, factor);
        d.add_element_load(
            p,
            id,
            ElementLoad {
                wall: Some(WallCaseLoad { case, factor: 1. }),
                ..Default::default()
            },
        );
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-8,
            max_iter: 8,
        })
        .build(d);
    analysis.step().unwrap();
    let response = analysis
        .domain()
        .element_stud_wall_response(id, 1.)
        .unwrap();
    near_relative(
        output.recorder_batches[0].samples[0].1,
        response.piers[0].q,
        1e-12,
        "pier q",
    );
    near_relative(
        output.recorder_batches[1].samples[0].1,
        response.piers[0].moment_base,
        1e-12,
        "base moment",
    );
    let direct_end_forces = analysis.domain().element_end_force(id, 1.);
    for (component, expected) in direct_end_forces.iter().enumerate() {
        near_relative(
            output.recorder_batches[component + 2].samples[0].1,
            *expected,
            1e-12,
            &format!("wall end force component {component}"),
        );
    }
    for end in &response.header_reactions[0..2] {
        near_relative(end[1].abs(), 4.8e6, 1e-12, "header reaction");
    }
}

/// Decodes the O86Crush anchored single-storey wall with slack at both corners, dead and lateral slots, and Newton steps.
/// Reference: `core/tests/elements/stud_wall_anchorage.rs` direct core construction; top-corner ux agrees to 1e-10 relative.
#[test]
fn anchorage_wire_drift_matches_direct_core() {
    const LENGTH: f64 = 2400.;
    const HEIGHT: f64 = 2440.;
    const POST: f64 = 2. * 9500. * 38. * 140.;
    const GA: f64 = 6000.;
    const TENSION: f64 = 15000.;
    const QR: f64 = 60000.;
    const DEAD: f64 = 3.;
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0., 0., 0., 0., LENGTH, 0., LENGTH, HEIGHT, 0., HEIGHT],
        fixed: vec![7, 0, 0, 0, 0],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![
        MaterialSpec::Elastic { e: 1e10 },
        MaterialSpec::Anchorage {
            tension: TENSION,
            slack: 0.5,
            bearing: AnchorageBearingSpec::O86Crush { qr: QR },
        },
    ];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0, 0],
        node_j: vec![1, 2],
        materials: vec![(0, 0, 0), (0, 1, 1), (1, 0, 0), (1, 1, 1)],
        ..Default::default()
    };
    input.stud_walls = vec![StudWallSpec {
        corners: [1, 2, 3, 4],
        stations: vec![],
        layout: WallLayoutSpec {
            length: LENGTH,
            height: HEIGHT,
            studs: vec![
                WallStudSpec { x: 0., ea: POST },
                WallStudSpec {
                    x: LENGTH,
                    ea: POST,
                },
            ],
            top_ea: POST * 1e6,
            top_ei: Some(1.2e10),
            top_distribution: ChordDistributionSpec::Beam,
            bottom_ea: POST * 1e6,
            bottom_ei: 6e9,
            openings: vec![],
            sheathing: Some(WallSheathingSpec {
                law: WallShearLawSpec::Linear { ga: GA },
            }),
            foundation: false,
            floor: None,
        },
    }];
    // Dead load and lateral shear are two patterns, ramped together.
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1. }; 2],
        scale_factor: vec![1.; 2],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0, 1],
        element_kind: vec![ElementKind::StudWall; 2],
        element_index: vec![0; 2],
        load: vec![
            ElementLoadSpec::WallLoad {
                load: WallLoadSpec::TopLine {
                    x0: 0.,
                    x1: LENGTH,
                    wx0: 0.,
                    wx1: 0.,
                    wy0: -DEAD,
                    wy1: -DEAD,
                },
            },
            ElementLoadSpec::WallLoad {
                load: WallLoadSpec::TopLine {
                    x0: 0.,
                    x1: LENGTH,
                    wx0: 5.,
                    wx1: 5.,
                    wy0: 0.,
                    wy1: 0.,
                },
            },
        ],
        stage: vec![0; 2],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "pushover".into(),
            steps: 10,
            integrator: IntegratorSpec::LoadControl { increment: 0.1 },
            algorithm: AlgorithmSpec::Newton {
                tangent: TangentStrategySpec::Current,
                line_search: None,
            },
            convergence: Some(ConvergenceSpec::NormDispIncr {
                tol: 1e-7,
                max_iter: 500,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 3, dof: 0 }],
    };
    let mut session = decode(input).unwrap();
    let output = session.advance(20);
    assert!(output.done && output.error.is_none(), "{output:?}");

    let mut d: Domain = Domain::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let anchor = Material::anchorage(TENSION, 0.5, AnchorageBearing::O86Crush { qr: QR }).unwrap();
    let base = [0., LENGTH].map(|x| {
        let n = d.add_node(Node::new([x, 0.]));
        d.add_element(Element::ZeroLength(
            ZeroLength::new(ground, n)
                .with_material(0, Material::Elastic { e: 1e10 })
                .with_material(1, anchor.clone()),
        ));
        n
    });
    let c2 = d.add_node(Node::new([LENGTH, HEIGHT]));
    let c3 = d.add_node(Node::new([0., HEIGHT]));
    let mut slots: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    slots[0].push(WallLoad::TopLine {
        x0: 0.,
        x1: LENGTH,
        wx0: 0.,
        wx1: 0.,
        wy0: -DEAD,
        wy1: -DEAD,
    });
    slots[1].push(WallLoad::TopLine {
        x0: 0.,
        x1: LENGTH,
        wx0: 5.,
        wx1: 5.,
        wy0: 0.,
        wy1: 0.,
    });
    let layout = StudWallLayout {
        length: LENGTH,
        height: HEIGHT,
        studs: vec![
            WallStud { x: 0., ea: POST },
            WallStud {
                x: LENGTH,
                ea: POST,
            },
        ],
        top_ea: POST * 1e6,
        top_ei: Some(1.2e10),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: POST * 1e6,
        bottom_ei: 6e9,
        openings: vec![],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: GA },
        }),
        foundation: false,
        floor: None,
    };
    let wall = d.add_element(Element::StudWall(
        StudWall::new([base[0], base[1], c2, c3], vec![], layout, slots).unwrap(),
    ));
    for case in 0..2 {
        let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
        d.add_element_load(
            p,
            wall,
            ElementLoad {
                wall: Some(WallCaseLoad { case, factor: 1. }),
                ..Default::default()
            },
        );
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.1 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-7,
            max_iter: 500,
        })
        .build(d);
    for _ in 0..10 {
        analysis.step().unwrap();
    }
    near_relative(
        output.recorder_batches[0].samples[9].1,
        analysis.domain().node(c2).displacement[0],
        1e-10,
        "top-corner ux",
    );
}

#[test]
fn shared_stud_wall_table_decodes_in_3d_and_rejects_nonorthogonal_geometry() {
    let mut i = empty_input(3);
    i.nodes = NodeTable {
        coords: vec![0., 0., 0., 1000., 0., 0., 1000., 1000., 0., 0., 1000., 0.],
        fixed: vec![63; 4],
        mass_node_index: vec![],
        mass: vec![],
    };
    i.stud_walls = vec![StudWallSpec {
        corners: [0, 1, 2, 3],
        stations: vec![],
        layout: layout(),
    }];
    let json = serde_json::to_value(&i).unwrap();
    let back: carapace_input::input_v1::CarapaceInputV1 =
        serde_json::from_value(json.clone()).unwrap();
    assert_eq!(serde_json::to_value(back.clone()).unwrap(), json);
    assert!(decode(back).is_ok());
    i.nodes.coords[4] = 100.;
    assert!(matches!(
        decode(i),
        Err(carapace_input::input_v1::DecodeError::InvalidMacroElement {
            table: "stud_walls",
            ..
        })
    ));
}

#[test]
fn wall_load_is_rejected_on_a_non_wall_and_wall_constructor_errors_surface() {
    let mut i = empty_input(2);
    i.nodes = NodeTable {
        coords: vec![0., 0., 1., 0.],
        fixed: vec![7, 7],
        mass_node_index: vec![],
        mass: vec![],
    };
    i.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![1.],
        material: vec![0],
        density: vec![0.],
    };
    i.materials = vec![MaterialSpec::Elastic { e: 1. }];
    i.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1. }],
        scale_factor: vec![1.],
    };
    i.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::Truss],
        element_index: vec![0],
        load: vec![ElementLoadSpec::WallLoad {
            load: WallLoadSpec::TopPoint {
                x: 0.,
                fx: 1.,
                fy: 0.,
            },
        }],
        stage: vec![0],
    };
    assert!(matches!(
        decode(i),
        Err(carapace_input::input_v1::DecodeError::UnsupportedElementLoad { .. })
    ));

    let mut bad = input_wall();
    bad.stud_walls[0].layout.length = -1.;
    assert!(matches!(
        decode(bad),
        Err(carapace_input::input_v1::DecodeError::InvalidMacroElement {
            table: "stud_walls",
            ..
        })
    ));
}
#[test]
fn anchorage_material_roundtrips_and_invalid_values_fail() {
    let good = MaterialSpec::Anchorage {
        tension: 1e5,
        slack: 0.,
        bearing: AnchorageBearingSpec::Linear { k: 1e6 },
    };
    let value = serde_json::to_value(&good).unwrap();
    assert_eq!(
        serde_json::to_value(serde_json::from_value::<MaterialSpec>(value.clone()).unwrap())
            .unwrap(),
        value
    );
    let mut i = empty_input(2);
    i.materials = vec![MaterialSpec::Anchorage {
        tension: -1.,
        slack: 0.,
        bearing: AnchorageBearingSpec::Linear { k: 1. },
    }];
    assert!(matches!(
        decode(i),
        Err(carapace_input::input_v1::DecodeError::InvalidAnchorageMaterial { index: 0 })
    ));
}

#[test]
fn floor_transfer_errors_are_not_dropped() {
    let mut i = input_wall();
    i.floor_loads = vec![FloorTransferSpec {
        bounds: [1., 0., 0., 1000.],
        continuity: FloorContinuitySpec::Simple,
        supports: vec![],
        pattern: 0,
        loads: vec![],
    }];
    assert!(matches!(
        decode(i),
        Err(carapace_input::input_v1::DecodeError::InvalidFloorTransfer { row: 0, .. })
    ));
}

/// A one-way floor's two simple-span reactions become `wallLoad`s of the row's pattern and the decoded wall solves.
/// Reference: statics gives 500 N per support for 1 kN total load; tolerance 1e-12 on total reaction.
#[test]
fn floor_transfer_loads_become_wall_loads_of_the_pattern() {
    let mut i = input_wall();
    i.element_loads = ElementLoadTable::default();
    i.floor_loads = vec![FloorTransferSpec {
        bounds: [0., 1000., 0., 1000.],
        continuity: FloorContinuitySpec::Simple,
        supports: [0., 1000.]
            .into_iter()
            .map(|s| FloorTransferSupportSpec {
                s,
                t0: 0.,
                t1: 1000.,
                wall: WallTargetSpec::StudWall { index: 0 },
                x_at_t0: 0.,
                reversed: false,
            })
            .collect(),
        pattern: 0,
        loads: vec![FloorTransferLoadSpec::Area {
            w: 0.001,
            region: None,
        }],
    }];
    let mut session = decode(i).unwrap();
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    let ux = outcome.recorder_batches[0].samples[0].1;
    assert!(ux.is_finite());
    let reactions = carapace_core::model::transfer_floor(
        carapace_core::model::FloorBounds {
            s0: 0.,
            s1: 1000.,
            t0: 0.,
            t1: 1000.,
        },
        &[
            carapace_core::model::FloorSupport {
                s: 0.,
                t0: 0.,
                t1: 1000.,
            },
            carapace_core::model::FloorSupport {
                s: 1000.,
                t0: 0.,
                t1: 1000.,
            },
        ],
        &[carapace_core::model::FloorLoad::Area {
            slot: 0,
            w: 0.001,
            region: None,
        }],
        carapace_core::model::FloorContinuity::Simple,
    )
    .unwrap();
    for support in &reactions.per_support {
        let got: f64 = support
            .lines
            .iter()
            .map(|l| l.w * (l.t1 - l.t0))
            .sum::<f64>()
            + support.points.iter().map(|p| p.p).sum::<f64>();
        assert!((got - 500.).abs() < 1e-12);
    }
}

#[test]
fn diaphragm_table_decodes_and_roundtrips() {
    let mut dia = empty_input(3);
    dia.nodes = NodeTable {
        coords: vec![0., 0., 0., 100., 0., 0., 0., 0., 100., 100., 0., 100.],
        fixed: vec![63; 4],
        mass_node_index: vec![],
        mass: vec![],
    };
    dia.diaphragms = vec![DiaphragmSpec {
        nodes: vec![0, 2, 1, 3],
        layout: DiaphragmLayoutSpec {
            u: vec![0., 100.],
            v: vec![0., 100.],
            level: 0.,
            bays: vec![Some(WallShearLawSpec::Linear { ga: 1000. })],
            ea_u: vec![1000., 1000.],
            ea_v: vec![1000., 1000.],
        },
    }];
    let json = serde_json::to_value(&dia).unwrap();
    let back: carapace_input::input_v1::CarapaceInputV1 =
        serde_json::from_value(json.clone()).unwrap();
    assert_eq!(serde_json::to_value(back.clone()).unwrap(), json);
    assert!(decode(back).is_ok());
    dia.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "dia".into(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1. },
            algorithm: AlgorithmSpec::Linear,
            convergence: None,
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::DiaphragmBay {
            element_index: 0,
            bay: 0,
            quantity: carapace_input::input_v1::sequence::WallBayQuantity::Q,
        }],
    };
    let mut session = decode(dia.clone()).unwrap();
    let outcome = session.advance(5);
    assert!(outcome.done && outcome.error.is_none());
    assert_eq!(outcome.recorder_batches[0].samples[0].1, 0.);
    dia.sequence.recorders = vec![RecorderSpec::DiaphragmBay {
        element_index: 0,
        bay: 1,
        quantity: carapace_input::input_v1::sequence::WallBayQuantity::Q,
    }];
    assert!(matches!(
        decode(dia),
        Err(carapace_input::input_v1::DecodeError::InvalidRecorderBay {
            bay: 1,
            count: 1,
            ..
        })
    ));
}

/// The wire format has no `building` table any more: an input that carries one fails to parse (the top level denies
/// unknown fields) instead of silently dropping it, and `combination` still defaults to ones in the spec.
#[test]
fn building_field_is_rejected_and_combination_defaults_to_one() {
    let value = serde_json::json!({"storeys":[],"foundationAnchorage":{"tension":1.0,"slack":0.0,"bearing":{"kind":"linear","k":1.0}},"diaphragm":"rigid"});
    let spec: BuildingSpec = serde_json::from_value(value.clone()).unwrap();
    assert_eq!(spec.combination, [1.; 8]);
    let mut input = serde_json::to_value(empty_input(3)).unwrap();
    input["building"] = value;
    let error = serde_json::from_value::<carapace_input::input_v1::CarapaceInputV1>(input)
        .unwrap_err()
        .to_string();
    assert!(error.contains("unknown field `building`"), "{error}");
}

fn wall_load_row(
    input: &mut carapace_input::input_v1::CarapaceInputV1,
    pattern: u32,
    fx: f64,
    x: f64,
) {
    let table = &mut input.element_loads;
    table.pattern.push(pattern);
    table.element_kind.push(ElementKind::StudWall);
    table.element_index.push(0);
    table.load.push(ElementLoadSpec::WallLoad {
        load: WallLoadSpec::TopPoint { x, fx, fy: -fx },
    });
    table.stage.push(0);
}

/// Top-corner ux of the decoded model after the (linear, single-step) run.
fn top_ux(input: carapace_input::input_v1::CarapaceInputV1) -> f64 {
    let mut session = decode(input).unwrap();
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    outcome.recorder_batches[0].samples[0].1
}

/// `wallLoad` rows of one pattern on one wall form one load case (the decoder emits one element load per pattern
/// and wall; two different cases from one pattern would fail `Domain::validate`): the rows give the same response
/// in either order, and the same rows split over two patterns give the same response (superposition, 1e-12).
#[test]
fn wall_loads_of_one_pattern_merge_into_one_case() {
    let combined = {
        let mut i = input_wall();
        i.element_loads = ElementLoadTable::default();
        wall_load_row(&mut i, 0, 100., 500.);
        wall_load_row(&mut i, 0, 40., 200.);
        i
    };
    // The same two loads as two `wallLoad` rows of one pattern, in the opposite order.
    let mut reordered = combined.clone();
    for column in [
        &mut reordered.element_loads.pattern,
        &mut reordered.element_loads.element_index,
        &mut reordered.element_loads.stage,
    ] {
        column.reverse();
    }
    reordered.element_loads.load.reverse();
    let (a, b) = (top_ux(combined.clone()), top_ux(reordered));
    assert!(a.abs() > 1e-9);
    near_relative(b, a, 1e-12, "reordered rows");
    // Two patterns (both factor 1) carry the loads as two cases.
    let mut two = combined.clone();
    two.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1. }; 2],
        scale_factor: vec![1.; 2],
    };
    two.element_loads.pattern[1] = 1;
    near_relative(top_ux(two), a, 1e-12, "two patterns");
}

/// Twelve patterns load one wall through the wire (no cap on patterns per wall); each has its own scale factor and
/// the response equals one pattern carrying the factored loads (1e-12).
#[test]
fn twelve_wall_load_patterns_equal_one_combined_pattern() {
    let factors: Vec<f64> = (0..12).map(|k| 0.5 + 0.2 * k as f64).collect();
    let mut many = input_wall();
    many.element_loads = ElementLoadTable::default();
    many.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1. }; 12],
        scale_factor: factors.clone(),
    };
    let mut one = input_wall();
    one.element_loads = ElementLoadTable::default();
    for (k, f) in factors.iter().enumerate() {
        wall_load_row(&mut many, k as u32, 10. + k as f64, 100. + 60. * k as f64);
        wall_load_row(&mut one, 0, f * (10. + k as f64), 100. + 60. * k as f64);
    }
    near_relative(top_ux(many), top_ux(one), 1e-12, "12 patterns");
}

/// A `floorLoads` row loads the pattern it names: its reactions become `wallLoad`s of that pattern. Splitting one
/// floor load into two rows of the same pattern, or putting the floor and a `wallLoad` in different patterns,
/// gives the same response as the merged form (linear superposition, 1e-12).
#[test]
fn floor_load_rows_load_their_pattern() {
    let floor = |pattern: u32, w: f64| FloorTransferSpec {
        bounds: [0., 1000., 0., 1000.],
        continuity: FloorContinuitySpec::Simple,
        supports: [0., 1000.]
            .into_iter()
            .map(|s| FloorTransferSupportSpec {
                s,
                t0: 0.,
                t1: 1000.,
                wall: WallTargetSpec::StudWall { index: 0 },
                x_at_t0: 0.,
                reversed: false,
            })
            .collect(),
        pattern,
        loads: vec![FloorTransferLoadSpec::Area { w, region: None }],
    };
    let base = || {
        let mut i = input_wall();
        i.load_patterns = LoadPatternTable {
            series: vec![TimeSeriesSpec::Linear { slope: 1. }; 2],
            scale_factor: vec![1.; 2],
        };
        i
    };
    // Reference: one floor row, pattern 0, plus the wallLoad of `input_wall` in pattern 0.
    let mut reference = base();
    reference.floor_loads = vec![floor(0, 0.002)];
    let expected = top_ux(reference);
    let mut split = base();
    split.floor_loads = vec![floor(0, 0.001), floor(0, 0.001)];
    near_relative(top_ux(split), expected, 1e-12, "split floor rows");
    let mut apart = base();
    apart.floor_loads = vec![floor(1, 0.002)];
    near_relative(top_ux(apart), expected, 1e-12, "floor in another pattern");
}

/// Misuse is an error, never dropped: a floor row naming a missing pattern, a `wallLoad` on a wall row past the
/// table, and `wallLoad` rows of one pattern and wall in different stages.
#[test]
fn wall_load_and_floor_load_misuse_is_an_error() {
    use carapace_input::input_v1::DecodeError;
    let mut bad = input_wall();
    bad.floor_loads = vec![FloorTransferSpec {
        bounds: [0., 1000., 0., 1000.],
        continuity: FloorContinuitySpec::Simple,
        supports: vec![],
        pattern: 1,
        loads: vec![],
    }];
    assert!(matches!(
        decode(bad),
        Err(DecodeError::UnknownPatternIndex {
            table: "floor_loads",
            row: 1
        })
    ));
    let mut bad = input_wall();
    bad.element_loads.element_index[0] = 4;
    assert!(matches!(
        decode(bad),
        Err(DecodeError::UnknownElementIndex { .. })
    ));
    let mut bad = input_wall();
    wall_load_row(&mut bad, 0, 1., 100.);
    bad.element_loads.stage[1] = 1;
    assert!(matches!(
        decode(bad),
        Err(DecodeError::InvalidRow {
            table: "element_loads",
            ..
        })
    ));
}

/// The `wallLoad` element load and the `pattern` of a `floorLoads` row serialize to the documented JSON and read
/// back unchanged.
#[test]
fn wall_load_and_floor_pattern_round_trip() {
    let load = ElementLoadSpec::WallLoad {
        load: WallLoadSpec::TopLine {
            x0: 0.,
            x1: 100.,
            wx0: 1.,
            wx1: 2.,
            wy0: -3.,
            wy1: -4.,
        },
    };
    let value = serde_json::to_value(load).unwrap();
    assert_eq!(
        value,
        serde_json::json!({"kind": "wallLoad", "load": {"kind": "topLine", "x0": 0.0, "x1": 100.0,
            "wx0": 1.0, "wx1": 2.0, "wy0": -3.0, "wy1": -4.0}})
    );
    assert_eq!(
        serde_json::to_value(serde_json::from_value::<ElementLoadSpec>(value.clone()).unwrap())
            .unwrap(),
        value
    );
    let floor = FloorTransferSpec {
        bounds: [0., 1., 0., 1.],
        continuity: FloorContinuitySpec::Continuous,
        supports: vec![],
        pattern: 3,
        loads: vec![
            FloorTransferLoadSpec::Area {
                w: 1.,
                region: Some([0., 1., 0., 1.]),
            },
            FloorTransferLoadSpec::Point {
                p: 2.,
                s: 0.5,
                t: 0.5,
            },
        ],
    };
    let value = serde_json::to_value(&floor).unwrap();
    assert_eq!(value["pattern"], 3);
    assert_eq!(value["loads"][0]["kind"], "area");
    assert!(value["loads"][0].get("slot").is_none());
    let back: FloorTransferSpec = serde_json::from_value(value.clone()).unwrap();
    assert_eq!(serde_json::to_value(back).unwrap(), value);
    // The retired wire forms are rejected rather than ignored (a `loads` field on a `studWalls` row too).
    let mut row = serde_json::to_value(&input_wall().stud_walls[0]).unwrap();
    row["loads"] = serde_json::json!([]);
    assert!(serde_json::from_value::<StudWallSpec>(row).is_err());
    assert!(serde_json::from_value::<ElementLoadSpec>(
        serde_json::json!({"kind": "wallSlots", "factors": [1, 0, 0, 0, 0, 0, 0, 0]})
    )
    .is_err());
}

/// `memberForce` bulk recorder on the window wall under a top line load: counts and positions follow the layout
/// (8 full-height studs, the 2 jacks at the opening edges, 3 cripples inside it), and the members reaching the
/// bottom plate (full-height studs at their bottom end, jacks) carry the base reaction (statics, 1e-9 relative).
#[test]
fn stud_wall_member_force_bulk_recorder() {
    use carapace_input::input_v1::sequence::{BulkQuantity, BulkTable};
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0., 0., 4000., 0., 4000., 2440., 0., 2440.],
        // Top corners free vertically (rotation only), so the whole top load runs down the framing.
        fixed: vec![7, 7, 4, 4],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.stud_walls = vec![StudWallSpec { corners: [0, 1, 2, 3], stations: vec![], layout: window_layout(false) }];
    input.load_patterns = LoadPatternTable { series: vec![TimeSeriesSpec::Linear { slope: 1. }], scale_factor: vec![1.] };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::StudWall],
        element_index: vec![0],
        load: vec![ElementLoadSpec::WallLoad { load: WallLoadSpec::TopLine { x0: 0., x1: 4000., wx0: 0., wx1: 0., wy0: -3., wy1: -3. } }],
        stage: vec![0],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "g".into(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1. },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-8, max_iter: 8 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::Bulk { table: BulkTable::Elements, quantity: BulkQuantity::MemberForce, element_kind: Some(ElementKind::StudWall), rows: Some(vec![0]) },
            RecorderSpec::StudWallBase {
                element_kind: ElementKind::StudWall,
                element_index: 0,
                component: carapace_input::input_v1::sequence::StudWallBaseComponent::Axial,
            },
        ],
    };
    let mut session = decode(input.clone()).unwrap();
    let output = session.advance(10);
    assert!(output.done && output.error.is_none(), "{output:?}");
    let batch = &output.bulk_batches[0];
    let v = &batch.values;
    assert_eq!((batch.points, batch.components), (1, 3 + 3 * 13));
    assert_eq!([v[0], v[1], v[2]], [8., 2., 3.]);
    let members: Vec<[f64; 3]> = v[3..].chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let sorted = |r: std::ops::Range<usize>| {
        let mut x: Vec<f64> = members[r].iter().map(|m| m[0]).collect();
        x.sort_by(f64::total_cmp);
        x
    };
    assert_eq!(sorted(0..8), [0., 400., 800., 1200., 2800., 3200., 3600., 4000.]);
    assert_eq!(sorted(8..10), [1200., 2800.]);
    assert_eq!(sorted(10..13), [1600., 2000., 2400.]);
    for m in &members[8..] {
        assert_eq!(m[1], m[2]);
    }
    // Direct distribution: each stud's top force is its tributary share of the 3 N/mm (end studs half); the opening's
    // share goes to the cripples. Bottoms differ by the sheathing's shear lag, but with the jacks they carry the base.
    for m in &members[..8] {
        let tributary = if m[0] == 0. || m[0] == 4000. { 600. } else { 1200. };
        assert!((m[2] + tributary).abs() <= 1e-9 * tributary, "stud top at {}: {}", m[0], m[2]);
    }
    let into_base: f64 = members[..10].iter().map(|m| m[1]).sum::<f64>();
    let base = output.recorder_batches.iter().find(|b| b.recorder_index == 1).unwrap().samples[0].1;
    assert!((into_base + base).abs() <= 1e-9 * base.abs(), "{into_base} vs base {base}");
    assert!((base - 12000.).abs() <= 1e-9 * 12000., "base {base}");
    // Anything but a stud-wall elements table is not applicable.
    input.sequence.recorders = vec![RecorderSpec::Bulk { table: BulkTable::Nodes, quantity: BulkQuantity::MemberForce, element_kind: None, rows: None }];
    assert!(decode(input).is_err());
}
