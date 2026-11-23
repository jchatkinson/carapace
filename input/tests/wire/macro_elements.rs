//! Wire tables for stud walls, floor bands, diaphragms, the anchorage material, floor loads and the building:
//! each decoded model equals the same model built directly through `carapace-core` / `carapace_input::building`
//! (window wall with floor loads, anchorage drift, building recorders: 1e-12 or 1e-10 as stated per test), specs
//! round-trip through serde, and every malformed row or misuse is a `DecodeError`.
use crate::common::empty_input;
use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    AnchorageBearing, ChordDistribution, Domain, Element, ElementLoad, FloorBounds,
    FloorContinuity, FloorLoad, FloorSupport, LoadSeries, Material, Node, StudWall, StudWallLayout,
    WallLoad, WallShearLaw, WallSheathing, WallStud, ZeroLength, WALL_SLOTS,
};
use carapace_input::building::{
    AnchorageBearingSpec, AnchorageSpec, BuildingSpec, DiaphragmKind, FloorContinuitySpec,
    FloorLoadSpec, FloorSpec, LateralLoad, SlotWallLoad, SpanAxis, StoreySpec, WallSpec,
};
use carapace_input::building::{
    ChordDistributionSpec, WallLayoutSpec, WallLoadSpec, WallShearLawSpec, WallSheathingSpec,
    WallStudSpec,
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
        loads: vec![SlotWallLoad {
            slot: 0,
            load: WallLoadSpec::TopPoint {
                x: 500.,
                fx: 100.,
                fy: 0.,
            },
        }],
    }];
    i.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1. }],
        scale_factor: vec![1.],
    };
    i.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::StudWall],
        element_index: vec![0],
        load: vec![ElementLoadSpec::WallSlots {
            factors: [1., 0., 0., 0., 0., 0., 0., 0.],
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
    };
    let mut loads: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| vec![]);
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
            wall: [1., 0., 0., 0., 0., 0., 0., 0.],
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

#[test]
fn macro_table_errors_are_reported() {
    let mut bad = input_wall();
    bad.stud_walls[0].corners[0] = 99;
    assert!(matches!(
        decode(bad),
        Err(carapace_input::input_v1::DecodeError::UnknownNodeIndex { .. })
    ));
    let mut bad = input_wall();
    bad.stud_walls[0].loads[0].slot = 8;
    assert!(matches!(
        decode(bad),
        Err(carapace_input::input_v1::DecodeError::InvalidWallSlot { .. })
    ));
    let mut bad = input_wall();
    bad.stud_walls = vec![StudWallSpec {
        corners: [0; 4],
        stations: vec![],
        layout: layout(),
        loads: vec![],
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
            loads: vec![],
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
    }];
    input.load_patterns = LoadPatternTable {
        series: vec![
            TimeSeriesSpec::Linear { slope: 1. },
            TimeSeriesSpec::Linear { slope: 1. },
        ],
        scale_factor: vec![1., 1.],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0, 1],
        element_kind: vec![ElementKind::StudWall; 2],
        element_index: vec![0, 0],
        load: vec![
            ElementLoadSpec::WallSlots {
                factors: [1.25, 0., 0., 0., 0., 0., 0., 0.],
            },
            ElementLoadSpec::WallSlots {
                factors: [0., 1.5, 0., 0., 0., 0., 0., 0.],
            },
        ],
        stage: vec![0, 0],
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
    let mut loads_by_support: Vec<[Vec<WallLoad>; WALL_SLOTS]> =
        (0..2).map(|_| std::array::from_fn(|_| vec![])).collect();
    for (support_index, support) in reactions.per_support.iter().enumerate() {
        let mapped = support.wall_loads(0., false);
        for slot in 0..WALL_SLOTS {
            loads_by_support[support_index][slot].extend(mapped[slot].iter().cloned());
        }
    }
    let id = make_direct_wall(&mut d, 0., &window_layout(true), loads_by_support.remove(0));
    make_direct_wall(
        &mut d,
        10000.,
        &window_layout(false),
        loads_by_support.remove(0),
    );
    for (slot, factor) in [(0, 1.25), (1, 1.5)] {
        let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
        let mut factors = [0.; WALL_SLOTS];
        factors[slot] = factor;
        d.add_element_load(
            p,
            id,
            ElementLoad {
                wall: factors,
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
        },
        loads: vec![
            SlotWallLoad {
                slot: 0,
                load: WallLoadSpec::TopLine {
                    x0: 0.,
                    x1: LENGTH,
                    wx0: 0.,
                    wx1: 0.,
                    wy0: -DEAD,
                    wy1: -DEAD,
                },
            },
            SlotWallLoad {
                slot: 1,
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
    }];
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1. }],
        scale_factor: vec![1.],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::StudWall],
        element_index: vec![0],
        load: vec![ElementLoadSpec::WallSlots {
            factors: [1., 1., 0., 0., 0., 0., 0., 0.],
        }],
        stage: vec![0],
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
    let mut slots: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| vec![]);
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
    };
    let wall = d.add_element(Element::StudWall(
        StudWall::new([base[0], base[1], c2, c3], vec![], layout, slots).unwrap(),
    ));
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_element_load(
        p,
        wall,
        ElementLoad {
            wall: [1., 1., 0., 0., 0., 0., 0., 0.],
            ..Default::default()
        },
    );
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

fn building_wall(name: &str, start: [f64; 2], direction: [f64; 2]) -> WallSpec {
    WallSpec {
        name: name.into(),
        start,
        direction,
        layout: WallLayoutSpec {
            length: 4000.,
            height: 2440.,
            studs: vec![
                WallStudSpec { x: 0., ea: 1e8 },
                WallStudSpec { x: 2000., ea: 1e8 },
                WallStudSpec { x: 4000., ea: 1e8 },
            ],
            top_ea: 1e12,
            top_ei: Some(1e12),
            top_distribution: ChordDistributionSpec::Beam,
            bottom_ea: 1e12,
            bottom_ei: 1e12,
            openings: vec![],
            sheathing: Some(WallSheathingSpec {
                law: WallShearLawSpec::Linear { ga: 1e7 },
            }),
            foundation: false,
        },
        loads: vec![],
    }
}

fn box_storey(floor_depth: f64) -> StoreySpec {
    StoreySpec {
        wall_height: 2440.,
        floor_depth,
        band_ga: 1e9,
        rim_ea: 1e12,
        anchorage: AnchorageSpec {
            tension: 15000.,
            slack: 0.,
            bearing: AnchorageBearingSpec::Linear { k: 1e9 },
        },
        bearing: 1e9,
        walls: vec![
            building_wall("front", [0., 0.], [1., 0.]),
            building_wall("back", [0., 4000.], [1., 0.]),
            building_wall("left", [0., 0.], [0., 1.]),
            building_wall("right", [4000., 0.], [0., 1.]),
        ],
        floors: vec![],
        lateral: vec![LateralLoad {
            slot: 0,
            fx: 10000.,
            fz: 0.,
            at: [2000., 2000.],
        }],
    }
}

/// Decodes a symmetric two-storey rigid-diaphragm box and records the lower front wall and first-storey drift.
/// Reference: `building::build` + `apply_combination` and the builder's direct analysis; q, base moment, and drift
/// agree with `wall_response` and the master-node displacement to 1e-12 relative.
#[test]
fn building_wire_recorders_match_direct_builder() {
    let spec = BuildingSpec {
        storeys: vec![box_storey(254.), box_storey(0.)],
        foundation_anchorage: AnchorageSpec {
            tension: 15000.,
            slack: 0.,
            bearing: AnchorageBearingSpec::Linear { k: 1e9 },
        },
        diaphragm: DiaphragmKind::Rigid,
        combination: [1.; 8],
    };
    let mut input = empty_input(3);
    input.building = Some(spec.clone());
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "building".into(),
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
            RecorderSpec::BuildingWallPier {
                storey: 0,
                wall: "front".into(),
                pier: 0,
                quantity: carapace_input::input_v1::sequence::StudWallQuantity::Q,
            },
            RecorderSpec::BuildingWallPier {
                storey: 0,
                wall: "front".into(),
                pier: 0,
                quantity: carapace_input::input_v1::sequence::StudWallQuantity::MomentBase,
            },
            RecorderSpec::BuildingStoreyDrift {
                storey: 0,
                component: carapace_input::input_v1::sequence::BuildingDriftComponent::X,
            },
            RecorderSpec::BuildingWallBase {
                storey: 0,
                wall: "front".into(),
                component: carapace_input::input_v1::sequence::StudWallBaseComponent::Shear,
            },
            RecorderSpec::BuildingWallBase {
                storey: 0,
                wall: "front".into(),
                component: carapace_input::input_v1::sequence::StudWallBaseComponent::Moment,
            },
        ],
    };
    let mut session = decode(input).unwrap();
    let output = session.advance(10);
    assert!(output.done && output.error.is_none(), "{output:?}");

    let mut built = carapace_input::building::build(&spec).unwrap();
    built.apply_combination(spec.combination).unwrap();
    let master = built.levels[0].master.unwrap();
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-8,
            max_iter: 8,
        })
        .build(built.domain.clone());
    analysis.step().unwrap();
    let response = built.wall_response(analysis.domain(), 0, "front").unwrap();
    near_relative(
        output.recorder_batches[0].samples[0].1,
        response.piers[0].q,
        1e-12,
        "building pier q",
    );
    near_relative(
        output.recorder_batches[1].samples[0].1,
        response.piers[0].moment_base,
        1e-12,
        "building base moment",
    );
    near_relative(
        output.recorder_batches[2].samples[0].1,
        analysis.domain().node(master).displacement[0],
        1e-12,
        "first-storey drift",
    );
    for (i, k) in [(3, 0), (4, 2)] {
        near_relative(
            output.recorder_batches[i].samples[0].1,
            response.base[k],
            1e-12,
            "building wall base",
        );
    }
}

/// `buildingWallBase` axial is the gravity takedown: a Z-span floor (2 kPa) on a 4 m × 4 m box puts
/// `w L_x L_z / 2 = 16 kN` on each of the front and back walls (statics, 1e-9), and none on the side walls.
#[test]
fn building_wall_base_axial_is_the_takedown() {
    let mut storey = box_storey(0.);
    storey.lateral.clear();
    storey.floors.push(FloorSpec {
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
    let mut input = empty_input(3);
    input.building = Some(BuildingSpec {
        storeys: vec![storey],
        foundation_anchorage: AnchorageSpec {
            tension: 15000.,
            slack: 0.,
            bearing: AnchorageBearingSpec::Linear { k: 1e9 },
        },
        diaphragm: DiaphragmKind::Rigid,
        combination: [1.; 8],
    });
    let base = |wall: &str| RecorderSpec::BuildingWallBase {
        storey: 0,
        wall: wall.into(),
        component: carapace_input::input_v1::sequence::StudWallBaseComponent::Axial,
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "gravity".into(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1. },
            algorithm: AlgorithmSpec::Linear,
            convergence: None,
            hold_patterns_after: vec![],
        }],
        recorders: ["front", "back", "left"].map(base).to_vec(),
    };
    let output = decode(input).unwrap().advance(10);
    assert!(output.done && output.error.is_none(), "{output:?}");
    let v = |i: usize| output.recorder_batches[i].samples[0].1;
    near_relative(v(0), 16_000., 1e-9, "front takedown");
    near_relative(v(1), 16_000., 1e-9, "back takedown");
    assert!(v(2).abs() < 1e-6, "side wall takedown {}", v(2));
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
        loads: vec![],
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
fn wall_slot_load_is_rejected_on_a_non_wall_and_wall_constructor_errors_surface() {
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
        load: vec![ElementLoadSpec::WallSlots { factors: [1.; 8] }],
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
        loads: vec![],
    }];
    assert!(matches!(
        decode(i),
        Err(carapace_input::input_v1::DecodeError::InvalidFloorTransfer { row: 0, .. })
    ));
}

/// A one-way floor's two simple-span reactions become ordinary wall slot loads and the decoded wall solves.
/// Reference: statics gives 500 N per support for 1 kN total load; tolerance 1e-12 on total reaction.
#[test]
fn floor_transfer_loads_are_appended_to_wall_slots() {
    let mut i = input_wall();
    i.stud_walls[0].loads.clear();
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
        loads: vec![FloorLoadSpec::Area {
            slot: 0,
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
fn floor_band_and_diaphragm_tables_decode_and_roundtrip() {
    let mut band = empty_input(2);
    band.nodes = NodeTable {
        coords: vec![0., 0., 100., 0., 0., 100., 100., 100.],
        fixed: vec![7; 4],
        mass_node_index: vec![],
        mass: vec![],
    };
    band.floor_bands = vec![FloorBandSpec {
        bottom: vec![0, 1],
        top: vec![2, 3],
        layout: FloorBandLayoutSpec {
            length: 100.,
            depth: 100.,
            xs: vec![0., 100.],
            bays: vec![Some(1e3)],
            rim_ea: 1e4,
        },
    }];
    let json = serde_json::to_value(&band).unwrap();
    let back: carapace_input::input_v1::CarapaceInputV1 =
        serde_json::from_value(json.clone()).unwrap();
    assert_eq!(serde_json::to_value(back.clone()).unwrap(), json);
    assert!(decode(back).is_ok());
    band.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "band".into(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1. },
            algorithm: AlgorithmSpec::Linear,
            convergence: None,
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::FloorBandBay {
            element_index: 0,
            bay: 0,
            quantity: carapace_input::input_v1::sequence::WallBayQuantity::Q,
        }],
    };
    let mut session = decode(band.clone()).unwrap();
    let outcome = session.advance(5);
    assert!(outcome.done && outcome.error.is_none());
    assert_eq!(outcome.recorder_batches[0].samples[0].1, 0.);
    band.sequence.recorders = vec![RecorderSpec::FloorBandBay {
        element_index: 0,
        bay: 1,
        quantity: carapace_input::input_v1::sequence::WallBayQuantity::Q,
    }];
    assert!(matches!(
        decode(band),
        Err(carapace_input::input_v1::DecodeError::InvalidRecorderBay {
            bay: 1,
            count: 1,
            ..
        })
    ));

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

#[test]
fn building_is_exclusive_with_explicit_tables_and_combination_defaults_to_one() {
    let value = serde_json::json!({"storeys":[],"foundationAnchorage":{"tension":1.0,"slack":0.0,"bearing":{"kind":"linear","k":1.0}},"diaphragm":"rigid"});
    let spec: BuildingSpec = serde_json::from_value(value).unwrap();
    assert_eq!(spec.combination, [1.; 8]);
    let mut input = empty_input(3);
    input.building = Some(spec);
    input.nodes = NodeTable {
        coords: vec![0., 0., 0.],
        fixed: vec![63],
        mass_node_index: vec![],
        mass: vec![],
    };
    assert!(matches!(
        decode(input),
        Err(carapace_input::input_v1::DecodeError::BuildingWithTables { table: "nodes" })
    ));
    let mut input = empty_input(2);
    input.building = Some(BuildingSpec {
        storeys: vec![],
        foundation_anchorage: carapace_input::building::AnchorageSpec {
            tension: 1.,
            slack: 0.,
            bearing: AnchorageBearingSpec::Linear { k: 1. },
        },
        diaphragm: DiaphragmKind::Rigid,
        combination: [1.; 8],
    });
    assert!(matches!(
        decode(input),
        Err(carapace_input::input_v1::DecodeError::TableNotInProfile { table: "building" })
    ));
}
