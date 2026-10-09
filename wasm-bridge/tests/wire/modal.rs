//! Modal stages through the wire format: eigenvalues, shapes and participation, in 2D and 3D.

use crate::common::{empty_input, last_sample};
use carapace_wasm::input_v1::decode;
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{RecorderSpec, SequenceSpec, StageSpec};
use carapace_wasm::input_v1::tables::*;

/// The golden-ratio closed form (core/tests/dynamics/modal.rs), decoded
/// entirely from wire tables: a `Modal` stage's mode frequencies come back
/// via `ModeShape` recorders, whose "time" component (`StepOutcome::
/// load_factor`'s doc comment) is the mode's own frequency rather than a
/// pseudo-time or load factor.
#[test]
fn decodes_a_modal_stage_and_matches_the_golden_ratio_closed_form() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0, 2.0, 0.0],
        fixed: vec![0b111, 0b110, 0b110],
        mass_node_index: vec![1, 2],
        mass: vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 1.0 }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0, 1],
        node_j: vec![1, 2],
        materials: vec![(0, 0, 0), (1, 0, 0)],
        friction: vec![],
        orient: vec![],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 2,
        }],
        recorders: vec![
            RecorderSpec::ModeShape {
                mode: 0,
                node: 1,
                dof: 0,
            },
            RecorderSpec::ModeShape {
                mode: 1,
                node: 1,
                dof: 0,
            },
        ],
    };

    let mut session = decode(input).expect("well-formed modal input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let (mode0_freq, _) = last_sample(&outcome, 0).expect("mode 0 sample");
    let (mode1_freq, _) = last_sample(&outcome, 1).expect("mode 1 sample");
    assert!(
        (mode0_freq - 1.0 / phi).abs() < 1e-9,
        "expected {}, got {mode0_freq}",
        1.0 / phi
    );
    assert!(
        (mode1_freq - phi).abs() < 1e-9,
        "expected {phi}, got {mode1_freq}"
    );
}

/// `Session::modal_results` returns frequencies, full-node shapes and participation for a
/// finished `Modal` stage — no recorders needed. For the 2-mass chain every mode's effective
/// mass sums to the total x mass (all modes requested), shapes are M-normalized and sign-fixed.
#[test]
fn modal_results_report_shapes_and_participation_without_recorders() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0, 2.0, 0.0],
        fixed: vec![0b111, 0b110, 0b110],
        mass_node_index: vec![1, 2],
        mass: vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 1.0 }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0, 1],
        node_j: vec![1, 2],
        materials: vec![(0, 0, 0), (1, 0, 0)],
        friction: vec![],
        orient: vec![],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 2,
        }],
        recorders: vec![],
    };

    let mut session = decode(input).expect("well-formed modal input should decode");
    assert!(session.modal_results().stages.is_empty());
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");

    let report = session.modal_results();
    assert_eq!(report.stages.len(), 1);
    let stage = &report.stages[0];
    assert_eq!(
        (stage.stage_index, stage.stage_id.as_str(), stage.ndf),
        (0, "modes", 3)
    );
    assert_eq!(stage.modes.len(), 2);
    assert!((stage.total_mass[0] - 2.0).abs() < 1e-12);
    assert_eq!(stage.total_mass[1], 0.0);

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    assert!((stage.modes[0].frequency - 1.0 / phi).abs() < 1e-9);
    let ratio_sum: f64 = stage.modes.iter().map(|m| m.mass_ratio[0]).sum();
    assert!(
        (ratio_sum - 1.0).abs() < 1e-9,
        "x mass ratios sum to {ratio_sum}"
    );
    for mode in &stage.modes {
        assert_eq!(mode.shape.len(), 3 * 3);
        // Ground node fixed -> zero; largest entry positive; M-normalized.
        assert!(mode.shape[..3].iter().all(|v| *v == 0.0));
        let peak = mode
            .shape
            .iter()
            .copied()
            .fold(0.0_f64, |b, v| if v.abs() > b.abs() { v } else { b });
        assert!(peak > 0.0);
        let norm: f64 = mode.shape[3] * mode.shape[3] + mode.shape[6] * mode.shape[6];
        assert!((norm - 1.0).abs() < 1e-9);
        assert_eq!(mode.mass_ratio[1], 0.0);
    }
}

/// An out-of-range `ModeShape.mode` records nothing this call, rather than
/// erroring — `decode` can't validate it upfront (a mode index isn't a
/// table row; it's only meaningful once `modal_analysis` has actually run
/// and produced however many modes it was asked for).
#[test]
fn a_mode_shape_recorder_past_the_computed_mode_count_records_nothing() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![1],
        mass: vec![1.0, 0.0, 0.0],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 1.0 }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)],
        friction: vec![],
        orient: vec![],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 1,
        }],
        recorders: vec![RecorderSpec::ModeShape {
            mode: 5,
            node: 1,
            dof: 0,
        }],
    };

    let mut session = decode(input).expect("well-formed modal input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );
    assert!(
        last_sample(&outcome, 0).is_none(),
        "out-of-range mode should record nothing"
    );
}

/// The spatial SDOF-truss closed form (core/tests/dynamics/free_vibration.rs's
/// `spatial_truss_mass_matches_sdof_closed_form_frequency_via_modal_analysis3`),
/// decoded entirely from 3D wire tables — proves `Modal`/`ModeShape` work
/// through the 3D decode path, not just the 2D one already
/// checked in `wire/elements.rs`.
#[test]
fn decodes_a_modal_stage3_and_matches_the_sdof_truss_closed_form() {
    let (e, area, length, density) = (30_000.0_f64, 2.0, 100.0, 0.5);
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, length, 0.0, 0.0],
        fixed: vec![0b111111, 0b111110], // node 1 free only in ux
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![area],
        material: vec![0],
        density: vec![density],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 1,
        }],
        recorders: vec![RecorderSpec::ModeShape {
            mode: 0,
            node: 1,
            dof: 0,
        }],
    };

    let mut session = decode(input).expect("well-formed spatial modal input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let k = e * area / length;
    let m = density * area * length / 2.0;
    let expected = (k / m).sqrt();
    let (got_freq, _) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!(
        (got_freq - expected).abs() < 1e-9,
        "expected omega={expected}, got {got_freq}"
    );
}

/// A rigid-diaphragm slave has no equation of its own, so its mode-shape
/// components must come from its constraint: reading them through
/// `equation_of` (as the session used to) reported a slave that never moves.
/// Mass and a rotational spring sit on the retained node; the translational
/// springs sit on the slave, so every mode moves it.
#[test]
fn modal_shapes_report_a_rigid_diaphragm_slaves_motion_through_the_constraint() {
    let (dx, dz) = (4.0_f64, -6.0_f64);
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        // node 0: ground; node 1: retained (free ux, uz, ry); node 2: slave (free ux, uz).
        coords: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dx, 0.0, dz],
        fixed: vec![0b111111, 0b101010, 0b111010],
        mass_node_index: vec![1],
        mass: vec![2.0, 0.0, 2.0, 0.0, 5.0, 0.0],
    };
    input.materials = vec![
        MaterialSpec::Elastic { e: 400.0 },
        MaterialSpec::Elastic { e: 250.0 },
        MaterialSpec::Elastic { e: 1000.0 },
    ];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0, 0],
        node_j: vec![2, 1],
        // Row 0: slave springs in ux and uz; row 1: a rotational spring about y on the retained node.
        materials: vec![(0, 0, 0), (0, 2, 1), (1, 4, 2)],
        friction: vec![],
        orient: vec![],
    };
    input.rigid_diaphragms = RigidDiaphragmTable {
        retained: vec![1],
        normal: vec![Axis3Spec::Y],
        constrained: vec![(0, 2)],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 3,
        }],
        recorders: vec![],
    };

    let mut session = decode(input).expect("well-formed spatial modal input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let report = session.modal_results();
    let modal = &report.stages[0];
    assert_eq!(modal.ndf, 6);
    assert_eq!(modal.modes.len(), 3);
    for mode in &modal.modes {
        let node = |n: usize, dof: usize| mode.shape[n * 6 + dof];
        let ry = node(1, 4);
        // Diaphragm about Y: u_c[z] = u_r[z] - ry * (x_c - x_r), u_c[x] = u_r[x] + ry * (z_c - z_r).
        let expected_x = node(1, 0) + ry * dz;
        let expected_z = node(1, 2) - ry * dx;
        assert!(
            (node(2, 0) - expected_x).abs() < 1e-9,
            "slave ux {} vs {expected_x}",
            node(2, 0)
        );
        assert!(
            (node(2, 2) - expected_z).abs() < 1e-9,
            "slave uz {} vs {expected_z}",
            node(2, 2)
        );
        assert!(
            node(2, 0).abs() + node(2, 2).abs() > 1e-9,
            "the slave must move"
        );
    }
}
