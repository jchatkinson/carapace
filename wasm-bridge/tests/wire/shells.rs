//! The `shellSections` and `shell4s` tables, shell pressure/self-weight element loads and
//! shell Gauss-point recorders, decoded and run end to end. The reference values are from
//! OpenSeesPy (`ShellMITC4` + `ElasticMembranePlateSection`) on the same model; the native
//! engine test (`core/tests/shells`) covers the element itself.

use crate::common::empty_input;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, GaussQuantity, IntegratorSpec, RecorderSpec, SequenceSpec,
    StageSpec,
};
use carapace_wasm::input_v1::tables::*;
use carapace_wasm::input_v1::{decode, CarapaceInputV1, DecodeError, StepOutcome};

fn one_step(recorders: Vec<RecorderSpec>) -> SequenceSpec {
    SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 5,
            }),
            hold_patterns_after: vec![],
        }],
        recorders,
    }
}

fn sample(outcome: &StepOutcome, recorder: usize) -> f64 {
    outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == recorder)
        .and_then(|batch| batch.samples.last())
        .map(|&(_, value)| value)
        .unwrap_or_else(|| panic!("recorder {recorder} recorded nothing"))
}

/// Two distorted quads on the tilted plane `z = 0.3 x + 0.1 y`; nodes 0 and 3 clamped.
fn tilted_patch(self_weight: bool) -> CarapaceInputV1 {
    let xy = [
        (0.0, 0.0),
        (2.0, 0.2),
        (2.3, 1.7),
        (-0.1, 1.4),
        (4.2, 0.1),
        (4.5, 1.6),
    ];
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: xy
            .iter()
            .flat_map(|&(x, y)| [x, y, 0.3 * x + 0.1 * y])
            .collect(),
        fixed: vec![63, 0, 0, 63, 0, 0],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.shell_sections = vec![ShellSectionSpec::ElasticMembranePlate {
        e: 3e4,
        nu: 0.25,
        h: 0.4,
        rho: 2e-3,
    }];
    input.shell4s = Shell4Table {
        node_ids: vec![0, 1, 2, 3, 1, 4, 5, 2],
        section: vec![0, 0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    let loads = [
        (4, [1.0, -2.0, 3.0, 0.5, -0.7, 0.2]),
        (5, [-0.5, 1.5, -4.0, 0.1, 0.3, -0.9]),
    ];
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0; 12],
        node: loads.iter().flat_map(|&(n, _)| [n; 6]).collect(),
        dof: (0..2).flat_map(|_| 0..6).collect(),
        value: loads.iter().flat_map(|&(_, v)| v).collect(),
        stage: vec![0; 12],
    };
    if self_weight {
        input.element_loads = ElementLoadTable {
            pattern: vec![0, 0],
            element_kind: vec![ElementKind::Shell4; 2],
            element_index: vec![0, 1],
            // OpenSees applies the opposite sign for `-selfWeight`.
            load: vec![
                ElementLoadSpec::ShellBody {
                    bx: 0.0,
                    by: 0.0,
                    bz: 9.81
                };
                2
            ],
            stage: vec![0, 0],
        };
    }
    input
}

fn node_disp_recorders(nodes: &[u32]) -> Vec<RecorderSpec> {
    nodes
        .iter()
        .flat_map(|&node| (0..6).map(move |dof| RecorderSpec::NodeDisp { node, dof }))
        .collect()
}

/// The tilted patch decoded from the wire reproduces OpenSees' displacements, with and
/// without self-weight, and the first element's Gauss-point resultants.
#[test]
fn a_shell_patch_runs_end_to_end_and_matches_opensees() {
    for (self_weight, want) in [
        (
            false,
            [
                0.032863282016,
                0.011036245391,
                -0.108117982242,
                -0.054959707678,
                0.05198534551,
                -0.010401219149,
            ],
        ),
        (
            true,
            [
                0.032155942744,
                0.010806864944,
                -0.105747978579,
                -0.054893757487,
                0.05108379486,
                -0.010469481701,
            ],
        ),
    ] {
        let mut input = tilted_patch(self_weight);
        let mut recorders = node_disp_recorders(&[4]);
        for component in 0..8 {
            recorders.push(RecorderSpec::GaussPoint {
                element_kind: ElementKind::Shell4,
                element_index: 0,
                point: 1,
                quantity: GaussQuantity::Stress,
                component,
            });
        }
        input.sequence = one_step(recorders);
        let mut session = decode(input).expect("shell patch decodes");
        let outcome = session.advance(10);
        assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
        for (k, w) in want.iter().enumerate() {
            let got = sample(&outcome, k);
            assert!(
                (got - w).abs() < 1e-8 * 0.11,
                "self_weight {self_weight}: dof {k}: {got} vs {w}"
            );
        }
        if !self_weight {
            let resultants = [
                0.17242928940233437,
                0.1851868064060666,
                -0.37431818354232427,
                2.645743301963054,
                -0.15656032296915945,
                1.3122547152760908,
                3.9950292388576747,
                0.008179010645220794,
            ];
            for (c, w) in resultants.iter().enumerate() {
                assert!(
                    (sample(&outcome, 6 + c) - w).abs() < 1e-8 * 4.0,
                    "resultant {c}"
                );
            }
        }
    }
}

/// A pressure through the wire is the same load as the native engine's: the clamp reactions of
/// a fully clamped quad sum to `p A` along the normal.
#[test]
fn shell_pressure_loads_through_the_wire() {
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 3.0, 0.0, 0.0, 3.0, 0.0],
        fixed: vec![63; 4],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.shell_sections = vec![ShellSectionSpec::ElasticMembranePlate {
        e: 1e4,
        nu: 0.3,
        h: 0.1,
        rho: 0.0,
    }];
    input.shell4s = Shell4Table {
        node_ids: vec![0, 1, 2, 3],
        section: vec![0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::Shell4],
        element_index: vec![0],
        load: vec![ElementLoadSpec::ShellPressure { pressure: 2.0 }],
        stage: vec![0],
    };
    input.sequence = one_step(
        (0..4)
            .map(|node| RecorderSpec::Reaction { node, dof: 2 })
            .collect(),
    );
    let mut session = decode(input).unwrap();
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    let total: f64 = (0..4).map(|k| sample(&outcome, k)).sum();
    assert!((total + 2.0 * 6.0).abs() < 1e-12, "total reaction {total}");
}

fn decode_error(input: CarapaceInputV1) -> DecodeError {
    decode(input).expect_err("expected a decode error")
}

/// The tables serialize with the expected camelCase shapes and round-trip.
#[test]
fn the_shell_tables_round_trip_through_serde() {
    let mut input = tilted_patch(true);
    input.sequence = one_step(vec![RecorderSpec::GaussPoint {
        element_kind: ElementKind::Shell4,
        element_index: 1,
        point: 3,
        quantity: GaussQuantity::Strain,
        component: 7,
    }]);
    let value = serde_json::to_value(&input).unwrap();
    let back: CarapaceInputV1 = serde_json::from_value(value.clone()).unwrap();
    assert_eq!(serde_json::to_value(&back).unwrap(), value);
    assert_eq!(value["shellSections"][0]["kind"], "elasticMembranePlate");
    assert_eq!(value["elementLoads"]["load"][0]["kind"], "shellBody");
    assert_eq!(value["elementLoads"]["elementKind"][0], "shell4");
    assert!(decode(back).is_ok());
}

/// Shells are 3D only; a shell load targets shells only; bad sections, indices and recorder
/// components are structured errors.
#[test]
fn shell_misuse_is_a_decode_error() {
    // A shell table in a 2D model.
    let mut input = empty_input(2);
    input.shell_sections = vec![ShellSectionSpec::ElasticMembranePlate {
        e: 1.0,
        nu: 0.3,
        h: 1.0,
        rho: 0.0,
    }];
    assert!(matches!(
        decode_error(input),
        DecodeError::TableNotInProfile {
            table: "shell_sections"
        }
    ));

    // Bad section constants.
    let mut input = tilted_patch(false);
    input.shell_sections = vec![ShellSectionSpec::ElasticMembranePlate {
        e: 3e4,
        nu: 0.25,
        h: 0.0,
        rho: 0.0,
    }];
    assert!(matches!(
        decode_error(input),
        DecodeError::InvalidShellSection { index: 0, .. }
    ));

    // Unknown section row.
    let mut input = tilted_patch(false);
    input.shell4s.section = vec![0, 3];
    assert!(matches!(
        decode_error(input),
        DecodeError::UnknownMaterialIndex {
            table: "shell4s",
            row: 3
        }
    ));

    // Wrong node stride.
    let mut input = tilted_patch(false);
    input.shell4s.node_ids.pop();
    assert!(matches!(
        decode_error(input),
        DecodeError::InvalidRow {
            table: "shell4s",
            ..
        }
    ));

    // A shell pressure on a non-shell element.
    let mut input = tilted_patch(false);
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![1.0],
        material: vec![0],
        density: vec![0.0],
    };
    input.materials = vec![carapace_wasm::input_v1::materials::MaterialSpec::Elastic { e: 1.0 }];
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::Truss],
        element_index: vec![0],
        load: vec![ElementLoadSpec::ShellPressure { pressure: 1.0 }],
        stage: vec![0],
    };
    input.sequence = one_step(vec![]);
    let error = decode_error(input);
    assert!(
        matches!(error, DecodeError::UnsupportedElementLoad { .. }),
        "{error:?}"
    );

    // A Gauss component past the 8 shell components.
    let mut input = tilted_patch(false);
    input.sequence = one_step(vec![RecorderSpec::GaussPoint {
        element_kind: ElementKind::Shell4,
        element_index: 0,
        point: 0,
        quantity: GaussQuantity::Stress,
        component: 8,
    }]);
    assert!(matches!(
        decode_error(input),
        DecodeError::InvalidRecorderComponent {
            component: 8,
            width: 8,
            ..
        }
    ));
}

/// A triangle patch decoded from the wire format reproduces OpenSees' `ShellDKGT` displacements and the
/// first triangle's resultants, and the table survives serde.
#[test]
fn a_triangle_patch_runs_end_to_end_and_matches_opensees() {
    let xy = [(0.0, 0.0), (2.0, 0.2), (2.3, 1.7), (-0.1, 1.4), (4.2, 0.1)];
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: xy
            .iter()
            .flat_map(|&(x, y)| [x, y, 0.3 * x + 0.1 * y])
            .collect(),
        fixed: vec![63, 0, 0, 63, 0],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.shell_sections = vec![ShellSectionSpec::ElasticMembranePlate {
        e: 3e4,
        nu: 0.25,
        h: 0.4,
        rho: 2e-3,
    }];
    input.shell3s = Shell3Table {
        node_ids: vec![0, 1, 2, 0, 2, 3, 1, 4, 2],
        section: vec![0, 0, 0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    let loads = [
        (4, [1.0, -2.0, 3.0, 0.5, -0.7, 0.2]),
        (2, [-0.5, 1.5, -4.0, 0.1, 0.3, -0.9]),
    ];
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0; 12],
        node: loads.iter().flat_map(|&(n, _)| [n; 6]).collect(),
        dof: (0..2).flat_map(|_| 0..6).collect(),
        value: loads.iter().flat_map(|&(_, v)| v).collect(),
        stage: vec![0; 12],
    };
    let mut recorders = node_disp_recorders(&[4]);
    for component in 0..8 {
        recorders.push(RecorderSpec::GaussPoint {
            element_kind: ElementKind::Shell3,
            element_index: 0,
            point: 3,
            quantity: GaussQuantity::Stress,
            component,
        });
    }
    input.sequence = one_step(recorders);
    let value = serde_json::to_value(&input).unwrap();
    let mut session =
        decode(serde_json::from_value(value).unwrap()).expect("triangle patch decodes");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    let want = [
        -0.06623758349489814,
        -0.028278257772653304,
        0.21662847499794943,
        -0.03716998390490014,
        -0.11139964702985417,
        -0.02700235006241788,
    ];
    for (k, w) in want.iter().enumerate() {
        assert!((sample(&outcome, k) - w).abs() < 1e-8 * 0.22, "dof {k}");
    }
    let resultants = [
        0.6268500586104436,
        1.419419397302228,
        -1.7079804732230155,
        -2.6876534877985727,
        0.188300211150203,
        1.6752693729004784,
        0.0,
        0.0,
    ];
    for (c, w) in resultants.iter().enumerate() {
        assert!(
            (sample(&outcome, 6 + c) - w).abs() < 1e-8 * 3.0,
            "resultant {c}"
        );
    }
}

/// A `shell3s` table in a 2D model is rejected, and its section index is checked.
#[test]
fn triangle_misuse_is_a_decode_error() {
    let mut input = empty_input(2);
    input.shell3s = Shell3Table {
        node_ids: vec![0, 1, 2],
        section: vec![0],
    };
    assert!(matches!(
        decode_error(input),
        DecodeError::TableNotInProfile { table: "shell3s" }
    ));
    let mut input = tilted_patch(false);
    input.shell3s = Shell3Table {
        node_ids: vec![0, 1, 2],
        section: vec![5],
    };
    assert!(matches!(
        decode_error(input),
        DecodeError::UnknownMaterialIndex {
            table: "shell3s",
            row: 5
        }
    ));
    let mut input = tilted_patch(false);
    input.shell3s = Shell3Table {
        node_ids: vec![0, 1],
        section: vec![0],
    };
    assert!(matches!(
        decode_error(input),
        DecodeError::InvalidRow {
            table: "shell3s",
            ..
        }
    ));
}
