//! Arc-length continuation through the wire format and `Session`
//! (`docs/arclength.md` §9 "WASM session / storage ordering"): both
//! profiles decode arc settings, decreasing and repeated load factors keep
//! monotonically increasing sample indices, the step budget doesn't change
//! accepted samples, stop criteria end a stage early, and invalid
//! combinations are rejected at decode.

use crate::common::empty_input;
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmConfigSpec, AlgorithmSpec, ArcLengthSpec, ArcScalesSpec, ArcSeedComponentSpec,
    ArcSeedSpec, ArcStopSpec, ConvergenceSpec, DisplacementTargetSpec, IntegratorSpec,
    LoadFactorTargetSpec, RecorderSpec, SequenceSpec, StageSpec, TangentStrategySpec,
};
use carapace_wasm::input_v1::tables::*;
use carapace_wasm::input_v1::{
    decode, AnalysisErrorDetail, CarapaceInputV1, DecodeError, StopReasonDetail,
};

fn arc_spec(radius: f64) -> ArcLengthSpec {
    ArcLengthSpec {
        initial_radius: radius,
        min_radius: None,
        max_radius: None,
        target_iterations: None,
        max_retries: None,
        direction: Default::default(),
        scales: ArcScalesSpec::Explicit {
            displacement: 0.01,
            rotation: None,
            load: 10.0,
        },
        predictor: Default::default(),
        seed: None,
        arc_tolerance: Some(1e-9),
        correction_tolerance: None,
        backtracking: None,
        stop: None,
    }
}

fn backbone(x: f64) -> f64 {
    if x <= 0.01 {
        1000.0 * x
    } else {
        10.0 - 400.0 * (x - 0.01)
    }
}

/// The §8 series snap-back model on the planar wire format: a degrading
/// `hysteretic` spring (node 0 -> 1) in series with an elastic `380`
/// spring (node 1 -> 2), unit load at node 2, stop exactly at `x = 0.028`.
fn snap_back_input(
    stop: Option<ArcStopSpec>,
    convergence: Option<ConvergenceSpec>,
) -> CarapaceInputV1 {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0; 6],
        fixed: vec![0b111, 0b110, 0b110],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![
        MaterialSpec::Hysteretic {
            mom1p: 10.0,
            rot1p: 0.01,
            mom2p: 6.0,
            rot2p: 0.02,
            mom3p: 2.0,
            rot3p: 0.03,
            mom1n: -10.0,
            rot1n: -0.01,
            mom2n: -6.0,
            rot2n: -0.02,
            mom3n: -2.0,
            rot3n: -0.03,
            pinch_x: 1.0,
            pinch_y: 1.0,
            damfc1: 0.0,
            damfc2: 0.0,
            beta: 0.0,
        },
        MaterialSpec::Elastic { e: 380.0 },
    ];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0, 1],
        node_j: vec![1, 2],
        materials: vec![(0, 0, 0), (1, 0, 1)],
        friction: vec![],
        orient: vec![],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![2],
        dof: vec![0],
        value: vec![1.0],
        stage: vec![0],
    };
    let mut arc = arc_spec(0.05);
    arc.stop = stop;
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "push".into(),
            steps: 250,
            integrator: IntegratorSpec::ArcLength(arc),
            algorithm: AlgorithmSpec::NewtonRaphson,
            convergence,
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 2, dof: 0 },
            RecorderSpec::Reaction { node: 0, dof: 0 },
        ],
    };
    input
}

fn x_stop() -> Option<ArcStopSpec> {
    Some(ArcStopSpec {
        displacement: Some(DisplacementTargetSpec {
            node: 1,
            dof: 0,
            value: 0.028,
            exact: true,
        }),
        ..Default::default()
    })
}

/// Runs a session to completion with a fixed step budget per `advance`,
/// concatenating each recorder's batches and checking that batch sample
/// indices are contiguous and monotonic.
fn run_with_budget(input: CarapaceInputV1, budget: u32) -> (Vec<Vec<(f64, f64)>>, Vec<u32>) {
    let mut session = decode(input).expect("arc-length input decodes");
    let mut samples = vec![Vec::new(); 3];
    let mut steps_per_call = Vec::new();
    loop {
        let outcome = session.advance(budget);
        assert!(outcome.error.is_none(), "{outcome:?}");
        steps_per_call.push(outcome.steps_taken);
        for batch in &outcome.recorder_batches {
            let history = &mut samples[batch.recorder_index];
            assert_eq!(
                batch.first_sample as usize,
                history.len(),
                "contiguous sample indices"
            );
            history.extend_from_slice(&batch.samples);
        }
        if outcome.steps_taken > 0 {
            let continuation = outcome.continuation.expect("arc steps report diagnostics");
            assert!(continuation.arc_error <= 1e-9);
        }
        if outcome.done {
            break;
        }
    }
    (samples, steps_per_call)
}

#[test]
fn snap_back_session_orders_samples_by_index_and_stops_early_on_target() {
    let (samples, steps) = run_with_budget(snap_back_input(x_stop(), None), 1000);
    let total_steps: u32 = steps.iter().sum();
    assert!(
        total_steps < 250,
        "the stop criterion ends the stage before its step cap"
    );

    let (x, y, reaction) = (&samples[0], &samples[1], &samples[2]);
    assert_eq!(x.len(), total_steps as usize);
    // Load factors (the samples' first component) rise, then fall: they are
    // not monotonic, and storage order is by sample index regardless.
    let lambdas: Vec<f64> = x.iter().map(|s| s.0).collect();
    let peak = lambdas
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap()
        .0;
    assert!(peak > 0 && peak < lambdas.len() - 10);
    for window in lambdas[peak..].windows(2) {
        assert!(window[1] < window[0]);
    }
    for i in 0..x.len() {
        let (lambda, xi) = x[i];
        assert!((lambda - backbone(xi)).abs() < 1e-6);
        assert!((y[i].1 - (xi + backbone(xi) / 380.0)).abs() < 1e-9);
        assert!((reaction[i].1 + lambda).abs() < 1e-6);
    }
    // Landed exactly on the target.
    assert!((x.last().unwrap().1 - 0.028).abs() < 1e-9);
}

#[test]
fn advance_budget_does_not_change_accepted_samples() {
    let (one_at_a_time, steps) = run_with_budget(snap_back_input(x_stop(), None), 1);
    assert!(steps.iter().all(|&s| s <= 1));
    let (all_at_once, _) = run_with_budget(snap_back_input(x_stop(), None), 1000);
    assert_eq!(one_at_a_time, all_at_once);
}

#[test]
fn stop_detail_and_default_combined_convergence_are_reported() {
    let mut session = decode(snap_back_input(x_stop(), None)).unwrap();
    let outcome = session.advance(1000);
    assert!(outcome.done && outcome.stage_complete && outcome.error.is_none());
    let stop = outcome.continuation.unwrap().stop.expect("stop reported");
    assert_eq!(stop.reason, StopReasonDetail::DisplacementTarget);
    assert!(stop.landed_exactly);

    // An explicit combined test decodes too.
    let combined = Some(ConvergenceSpec::Combined {
        force_tol: 1e-9,
        moment_tol: None,
        relative_tol: Some(1e-9),
        displacement_tol: None,
        max_iter: 30,
    });
    let outcome = decode(snap_back_input(x_stop(), combined))
        .unwrap()
        .advance(1000);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
}

#[test]
fn spatial_profile_decodes_arc_settings_and_lands_on_a_load_target() {
    let k = 80.0;
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0; 6],
        fixed: vec![0b111111, 0b111011], // node 1 free only in uz
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e: k }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 2, 0)],
        friction: vec![],
        orient: vec![],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![2],
        value: vec![1.0],
        stage: vec![0],
    };
    let mut arc = arc_spec(0.3);
    arc.scales = ArcScalesSpec::Auto { load: 1.0 };
    arc.seed = Some(ArcSeedSpec {
        components: vec![ArcSeedComponentSpec {
            node: 1,
            dof: 2,
            value: 1.0,
        }],
        load: 0.0,
    });
    arc.stop = Some(ArcStopSpec {
        load_factor: Some(LoadFactorTargetSpec {
            value: 2.0,
            exact: true,
        }),
        ..Default::default()
    });
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "push".into(),
            steps: 100,
            integrator: IntegratorSpec::ArcLength(arc),
            algorithm: AlgorithmSpec::NewtonRaphson,
            convergence: None,
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 2 }],
    };
    let outcome = decode(input).unwrap().advance(1000);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    assert!((outcome.load_factor - 2.0).abs() < 1e-9);
    let (lambda, uz) = *outcome.recorder_batches[0].samples.last().unwrap();
    assert!((k * uz - lambda).abs() < 1e-9);
    let continuation = outcome.continuation.unwrap();
    // Auto scale: u* = load * |K^-1 p| = 1/80.
    assert!((continuation.displacement_scale - 1.0 / k).abs() < 1e-15);
    assert_eq!(
        continuation.stop.unwrap().reason,
        StopReasonDetail::LoadFactorTarget
    );
}

#[test]
fn invalid_arc_length_options_are_rejected_at_decode() {
    let expect_option = |input: CarapaceInputV1, field: &str| match decode(input) {
        Err(DecodeError::InvalidAnalysisOption { field: got, .. }) => assert_eq!(got, field),
        other => panic!("expected {field}, got {other:?}"),
    };
    let with_arc = |edit: &dyn Fn(&mut ArcLengthSpec)| {
        let mut input = snap_back_input(None, None);
        let StageSpec::Static { integrator, .. } = &mut input.sequence.stages[0] else {
            unreachable!()
        };
        let IntegratorSpec::ArcLength(arc) = integrator else {
            unreachable!()
        };
        edit(arc);
        input
    };

    expect_option(
        with_arc(&|arc| arc.min_radius = Some(0.1)),
        "integrator.initialRadius",
    );
    expect_option(
        with_arc(&|arc| {
            arc.scales = ArcScalesSpec::Explicit {
                displacement: -1.0,
                rotation: None,
                load: 1.0,
            }
        }),
        "integrator.displacementScale",
    );
    expect_option(
        with_arc(&|arc| arc.arc_tolerance = Some(f64::NAN)),
        "integrator.arcTolerance",
    );
    match decode(with_arc(&|arc| {
        arc.seed = Some(ArcSeedSpec {
            components: vec![ArcSeedComponentSpec {
                node: 1,
                dof: 7,
                value: 1.0,
            }],
            load: 0.0,
        })
    })) {
        Err(DecodeError::InvalidDof { dof: 7, .. }) => {}
        other => panic!("expected InvalidDof, got {other:?}"),
    }

    let with_stage = |algorithm: AlgorithmSpec, convergence: Option<ConvergenceSpec>| {
        let mut input = snap_back_input(None, convergence);
        let StageSpec::Static { algorithm: a, .. } = &mut input.sequence.stages[0] else {
            unreachable!()
        };
        *a = algorithm;
        input
    };
    expect_option(with_stage(AlgorithmSpec::Linear, None), "algorithm");
    expect_option(
        with_stage(
            AlgorithmSpec::Config(AlgorithmConfigSpec::Newton {
                tangent: TangentStrategySpec::Initial,
                line_search: None,
            }),
            None,
        ),
        "algorithm",
    );
    expect_option(
        with_stage(
            AlgorithmSpec::NewtonRaphson,
            Some(ConvergenceSpec::NormDispIncr {
                tol: 1e-9,
                max_iter: 10,
            }),
        ),
        "convergence",
    );
    expect_option(
        with_stage(
            AlgorithmSpec::NewtonRaphson,
            Some(ConvergenceSpec::Combined {
                force_tol: 0.0,
                moment_tol: None,
                relative_tol: None,
                displacement_tol: None,
                max_iter: 10,
            }),
        ),
        "convergence.forceTol",
    );
}

#[test]
fn runtime_arc_length_errors_reach_the_session_as_structured_detail() {
    // No reference load at all: decode can't know, the first step reports it.
    let mut input = snap_back_input(x_stop(), None);
    input.nodal_loads = NodalLoadTable::default();
    let outcome = decode(input).unwrap().advance(10);
    assert!(outcome.done);
    assert_eq!(
        outcome.error,
        Some(AnalysisErrorDetail::ZeroLoadSensitivity)
    );
}
