//! Analysis stages through the wire format: configurable algorithms, transient and pushover stages,
//! `record_initial`, and `reset_stage`.

use crate::common::{empty_input, last_sample};
use carapace_wasm::input_v1::decode;
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
};
use carapace_wasm::input_v1::tables::*;

/// Nonlinear equilibrium must be the same across each configurable solver.
#[test]
fn configurable_algorithms_solve_static_and_transient_yielding() {
    use carapace_wasm::input_v1::sequence::{
        AlgorithmConfigSpec, DampingSpec, LineSearchSpec, TangentStrategySpec,
    };
    let mut algorithms = vec![AlgorithmSpec::NewtonRaphson];
    for tangent in [
        TangentStrategySpec::Current,
        TangentStrategySpec::ReuseAtStepStart,
        TangentStrategySpec::Initial,
    ] {
        for line_search in [
            None,
            Some(LineSearchSpec::Bisection {
                tol: 1e-10,
                max_iter: 30,
                max_eta: 16.0,
            }),
            Some(LineSearchSpec::RegulaFalsi {
                tol: 1e-10,
                max_iter: 30,
                max_eta: 16.0,
            }),
        ] {
            algorithms.push(AlgorithmSpec::Config(AlgorithmConfigSpec::Newton {
                tangent,
                line_search,
            }));
        }
        algorithms.push(AlgorithmSpec::Config(AlgorithmConfigSpec::KrylovNewton {
            tangent,
            max_dimension: 3,
        }));
    }
    for algorithm in algorithms {
        for transient in [false, true] {
            let mut input = empty_input(2);
            input.nodes = NodeTable {
                coords: vec![0.0, 0.0, 0.0, 0.0],
                fixed: vec![0b111, 0b110],
                mass_node_index: vec![1],
                mass: vec![1.0, 0.0, 0.0],
            };
            input.materials = vec![
                MaterialSpec::Elastic { e: 50.0 },
                MaterialSpec::ElasticPp {
                    e: 100.0,
                    eyp: 0.01,
                },
            ];
            input.zero_lengths = ZeroLengthTable {
                node_i: vec![0, 0],
                node_j: vec![1, 1],
                materials: vec![(0, 0, 0), (1, 0, 1)],
                friction: vec![],
                orient: vec![],
            };
            input.load_patterns = LoadPatternTable {
                series: vec![if transient {
                    TimeSeriesSpec::Linear { slope: 10.0 }
                } else {
                    TimeSeriesSpec::Constant
                }],
                scale_factor: vec![1.0],
            };
            input.nodal_loads = NodalLoadTable {
                pattern: vec![0],
                node: vec![1],
                dof: vec![0],
                value: vec![if transient { 1000.0 } else { 5.0 }],
                stage: vec![0],
            };
            let convergence = Some(ConvergenceSpec::NormUnbalance {
                tol: if transient { 1e-6 } else { 1e-9 },
                max_iter: 120,
            });
            input.sequence = SequenceSpec {
                stages: vec![if transient {
                    StageSpec::Transient {
                        id: "yield".into(),
                        steps: 1,
                        dt: 0.1,
                        damping: DampingSpec {
                            alpha_m: 0.0,
                            beta_k: 0.0,
                        },
                        ground_motions: vec![],
                        algorithm,
                        convergence,
                    }
                } else {
                    StageSpec::Static {
                        id: "yield".into(),
                        steps: 1,
                        integrator: IntegratorSpec::LoadControl { increment: 1.0 },
                        algorithm,
                        convergence,
                        hold_patterns_after: vec![],
                    }
                }],
                recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
            };
            let outcome = decode(input).unwrap().advance(1);
            assert!(
                outcome.done && outcome.error.is_none(),
                "{algorithm:?}, transient={transient}: {outcome:?}"
            );
            // A ramp starts at zero force; first Newmark effective inertia = 4*m/dt² = 400.
            let expected = if transient {
                (1000.0 - 1.0) / (400.0 + 50.0)
            } else {
                (5.0 - 1.0) / 50.0
            };
            assert!(
                (last_sample(&outcome, 0).unwrap().1 - expected).abs() < 1e-8,
                "{algorithm:?}, transient={transient}: {outcome:?}, expected={expected}"
            );
        }
    }
}

/// Rayleigh-damped free vibration (core/tests/dynamics/free_vibration.rs's SDOF case),
/// decoded from a `Transient` stage: matches the closed-form damped
/// free-vibration displacement after a known number of Newmark steps.
#[test]
fn decodes_a_transient_stage_and_matches_damped_free_vibration_closed_form() {
    let (mass, k, dt, steps) = (1.0_f64, 100.0, 0.001, 500);
    let alpha_m = 0.2; // Rayleigh mass-proportional damping coefficient
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![1],
        mass: vec![mass, 0.0, 0.0],
    };
    input.materials = vec![MaterialSpec::Elastic { e: k }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)],
        friction: vec![],
        orient: vec![],
    };
    // The wire format has no initial-displacement field yet, so this
    // exercises the Newmark integration a different (still closed-form)
    // way: a `Constant` reference load, already active when the transient
    // stage starts, is a suddenly-applied constant force from rest — the
    // classic damped step response, converging to static equilibrium
    // `u_static = force/k` with decaying oscillation, rather than free
    // vibration decaying to zero from a nonzero initial offset.
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Constant],
        scale_factor: vec![1.0],
    };
    let applied_force = 10.0;
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![applied_force],
        stage: vec![0],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Transient {
            id: "dynamic".to_string(),
            steps,
            dt,
            damping: carapace_wasm::input_v1::sequence::DampingSpec {
                alpha_m,
                beta_k: 0.0,
            },
            ground_motions: vec![],
            algorithm: AlgorithmSpec::Linear,
            convergence: None,
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed transient input should decode");
    let outcome = session.advance(steps);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    // Closed-form damped step response from rest (Chopra, "Dynamics of
    // Structures"): u(t) = u_static * (1 - e^{-zeta*omega*t} * (cos(omega_d
    // t) + (zeta*omega/omega_d) * sin(omega_d t))), damping ratio
    // `zeta = alpha_m/(2*omega)` (mass-proportional-only Rayleigh damping).
    let omega = (k / mass).sqrt();
    let zeta = alpha_m / (2.0 * omega);
    let omega_d = omega * (1.0 - zeta * zeta).sqrt();
    let u_static = applied_force / k;
    let t = dt * steps as f64;
    let expected = u_static
        * (1.0
            - (-zeta * omega * t).exp()
                * ((omega_d * t).cos() + (zeta * omega / omega_d) * (omega_d * t).sin()));

    // Newmark discretization error, not exact — same tolerance
    // core/tests/dynamics/free_vibration.rs's own closed-form comparison uses.
    let (_, got) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!(
        (got - expected).abs() < 1e-4,
        "expected {expected}, got {got}"
    );
}

/// The first required case: a 2D fiber-section
/// cantilever runs gravity (`LoadControl`), freezes it, then a
/// displacement-controlled lateral pushover — built the same way
/// `core/tests/elements/force_beam.rs` builds and checks a
/// `ForceBeamColumn` cantilever, but assembled entirely from decoded wire
/// tables. Elastic-range and post-yield behavior must agree with that
/// native acceptance test, and gravity's frozen axial displacement must
/// survive the pushover stage untouched (`core/tests/analysis/load_patterns.rs`'s
/// same check, here across a decoded multi-stage `AnalysisSequence`).
#[test]
fn gravity_then_pushover_matches_native_force_beam_column_behavior() {
    let (e, length) = (30000.0_f64, 100.0);
    let area = 0.25;
    let ys = [20.0_f64, 15.0, 10.0, 5.0];
    let iz = area * ys.iter().map(|y| y * y * 2.0).sum::<f64>();
    let eyp = 0.001;
    let fy = e * eyp;
    let my = fy * iz / ys[0];
    let p_yield = my / length;

    // Fiber section: 8 fibers (±20, ±15, ±10, ±5), ElasticPP.
    let mut fiber_y = Vec::new();
    let mut fiber_area = Vec::new();
    let mut fiber_material = Vec::new();
    for y in ys.iter().flat_map(|&y| [y, -y]) {
        fiber_y.push(y);
        fiber_area.push(area);
        fiber_material.push(0u32); // single shared ElasticPP material entry
    }

    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b111, 0b000],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::ElasticPp { e, eyp }];
    input.fibers = FiberTable {
        section_offsets: vec![0, fiber_y.len() as u32],
        z: vec![],
        y: fiber_y,
        area: fiber_area,
        material: fiber_material,
    };
    input.force_beam_columns_2d = FiberBeamColumn2dTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        integration: vec![IntegrationSpec::Lobatto { points: 3 }],
        corotational: vec![false],
        density: vec![0.0],
    };
    // Pattern 0: gravity (axial, dof 0), registered on stage 0. Pattern 1:
    // lateral reference load (transverse, dof 1), registered only on stage
    // 1 — *not* upfront, since stage 0's `LoadControl` shares one
    // pseudo-time across every unfrozen pattern, so an already-registered
    // pattern 1 would ramp during gravity too (the same reason
    // core/tests/elements/force_beam.rs's native two-phase test adds its
    // lateral load only between phases). `DisplacementControl` needs this
    // reference load's nonzero sensitivity on the controlled dof, same as
    // core/tests/analysis/newton.rs.
    input.load_patterns = LoadPatternTable {
        series: vec![
            TimeSeriesSpec::Linear { slope: 1.0 },
            TimeSeriesSpec::Linear { slope: 1.0 },
        ],
        scale_factor: vec![1.0, 1.0],
    };
    // Small relative to this section's axial yield capacity (area_total *
    // fy = 2.0 * 30.0 = 60.0) and, more importantly, small enough that the
    // uniform axial strain it adds doesn't push the elastic-range pushover
    // check's outer fiber (already at half its yield curvature by design,
    // mirroring core/tests/elements/force_beam.rs) past `eyp`.
    let gravity_axial = -10.0;
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0, 1],
        node: vec![1, 1],
        dof: vec![0, 1],
        value: vec![gravity_axial, -1.0],
        stage: vec![0, 1],
    };
    input.sequence = SequenceSpec {
        stages: vec![
            StageSpec::Static {
                id: "gravity".to_string(),
                steps: 2,
                integrator: IntegratorSpec::LoadControl { increment: 0.5 },
                algorithm: AlgorithmSpec::NewtonRaphson,
                convergence: Some(ConvergenceSpec::NormUnbalance {
                    tol: 1e-10,
                    max_iter: 30,
                }),
                hold_patterns_after: vec![0],
            },
            StageSpec::Static {
                id: "pushover".to_string(),
                steps: 3,
                integrator: IntegratorSpec::DisplacementControl {
                    node: 1,
                    dof: 1,
                    increment: -0.5 * p_yield * length.powi(3) / (3.0 * e * iz),
                },
                algorithm: AlgorithmSpec::NewtonRaphson,
                convergence: Some(ConvergenceSpec::NormUnbalance {
                    tol: 1e-10,
                    max_iter: 30,
                }),
                hold_patterns_after: vec![],
            },
        ],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 1, dof: 1 },
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");

    // Drain the gravity stage first so we can read its frozen axial value.
    let mut outcome = session.advance(2);
    assert!(
        outcome.stage_complete && outcome.error.is_none(),
        "gravity stage should complete cleanly: {outcome:?}"
    );
    let axial_after_gravity = last_sample(&outcome, 0).unwrap().1;
    assert!(
        axial_after_gravity.abs() > 1e-9,
        "gravity should produce nonzero axial displacement"
    );
    // Gravity took both budgeted steps in this one advance() call, so recorder 0's batch should
    // cover samples [0, 2) of stage 0 — the "exact sample/time/value ordering across multiple
    // advance() calls and stage transitions" this bounded-batch slice needs to prove.
    let gravity_batch = outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == 0)
        .expect("recorder 0 batch");
    assert_eq!(gravity_batch.first_sample, 0);
    assert_eq!(gravity_batch.stage_index, 0);
    assert_eq!(gravity_batch.samples.len(), 2);

    // First pushover step stays within the elastic range: must match the
    // closed-form cantilever deflection exactly, same as m8's elastic
    // check. The section's fiber layout is symmetric and still uniformly
    // elastic at this point, so — same reasoning as
    // core/tests/analysis/load_patterns.rs's frozen-pattern check, here for a
    // fiber section rather than `ElasticBeamColumn` — gravity's frozen
    // axial displacement must also still be completely untouched.
    outcome = session.advance(1);
    assert!(
        outcome.error.is_none(),
        "elastic-range pushover step should converge: {outcome:?}"
    );
    let (_, transverse_1) = last_sample(&outcome, 1).unwrap();
    let (_, axial_1) = last_sample(&outcome, 0).unwrap();
    let expected_elastic = -0.5 * p_yield * length.powi(3) / (3.0 * e * iz);
    assert!(
        (transverse_1 - expected_elastic).abs() < 1e-9,
        "expected {expected_elastic}, got {transverse_1}"
    );
    assert!(
        (axial_1 - axial_after_gravity).abs() < 1e-9,
        "frozen gravity pattern's axial contribution must be unchanged while the section is still elastic: {axial_after_gravity} vs {axial_1}"
    );
    // This call crossed into the pushover stage: its one sample continues the running count
    // from gravity's two (first_sample == 2), tagged with the new stage index.
    let pushover_batch = outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == 1)
        .expect("recorder 1 batch");
    assert_eq!(pushover_batch.first_sample, 2);
    assert_eq!(pushover_batch.stage_index, 1);
    assert_eq!(pushover_batch.samples.len(), 1);

    // Remaining two steps push well past yield (cumulative target 1.5x the
    // elastic-range increment, past the outer fibers' `eyp`). Checking the
    // controlled dof's own displacement here would be circular —
    // `DisplacementControl` drives it to exactly the requested cumulative
    // increment regardless of material state, converged or not. Instead
    // check the *force* (`load_factor`, since the reference lateral load is
    // a unit load) needed to reach that displacement: a softened section
    // needs measurably less force than pure elastic stiffness would predict
    // for the same displacement, the load-controlled equivalent of m8's
    // "softer than elastic" check.
    //
    // Axial/bending decoupling itself does *not* survive into the
    // post-yield region here (unlike `analysis/load_patterns.rs`'s `ElasticBeamColumn`
    // case): once the fiber layout's yielding is asymmetric, the section's
    // *tangent* couples axial and bending even though its geometry is
    // symmetric — a real property of fiber sections, not a decode bug — so
    // this test doesn't assert axial invariance past this point.
    outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "pushover stage should finish cleanly: {outcome:?}"
    );

    let (_, transverse_final) = last_sample(&outcome, 1).unwrap();
    let elastic_force_for_final_displacement =
        transverse_final.abs() * 3.0 * e * iz / length.powi(3);
    assert!(
        outcome.load_factor.abs() < elastic_force_for_final_displacement * 0.99,
        "softened section should need measurably less force than elastic stiffness predicts: needed {}, elastic would need {elastic_force_for_final_displacement}",
        outcome.load_factor.abs()
    );
}

/// A bar loaded in two static stages that both ramp pattern 0 (load stage `0`), optionally with a
/// `Reset` between them. `hold` freezes the pattern after the first stage.
fn two_stage_bar(
    record_initial: bool,
    between: Option<StageSpec>,
    hold: bool,
) -> carapace_wasm::input_v1::Session {
    let static_stage = |id: &str, hold: bool| StageSpec::Static {
        id: id.to_string(),
        steps: 2,
        integrator: IntegratorSpec::LoadControl { increment: 0.5 },
        algorithm: AlgorithmSpec::Linear,
        convergence: Some(ConvergenceSpec::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        }),
        hold_patterns_after: if hold { vec![0] } else { vec![] },
    };
    let mut input = empty_input(2);
    input.header.record_initial = record_initial;
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 100.0, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 30000.0 }];
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![2.0],
        material: vec![0],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![60.0],
        stage: vec![0],
    };
    let mut stages = vec![static_stage("first", hold)];
    stages.extend(between);
    stages.push(static_stage("second", false));
    input.sequence = SequenceSpec {
        stages,
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
    };
    decode(input).expect("well-formed input should decode")
}

/// Every sample of recorder 0, as `(stage_index, pseudo_time, value)`, running a session to the end.
fn all_samples(session: &mut carapace_wasm::input_v1::Session) -> Vec<(usize, f64, f64)> {
    let mut out = Vec::new();
    loop {
        let outcome = session.advance(8);
        assert!(outcome.error.is_none(), "{outcome:?}");
        for batch in outcome
            .recorder_batches
            .iter()
            .filter(|b| b.recorder_index == 0)
        {
            out.extend(
                batch
                    .samples
                    .iter()
                    .map(|&(t, v)| (batch.stage_index, t, v)),
            );
        }
        if outcome.done {
            return out;
        }
    }
}

/// `header.recordInitial` adds each static stage's initial conditions as its first sample: load
/// factor 0 and the state the stage inherits (zero for the first, the first stage's end for the
/// second). Without it a stage's samples are only its steps.
#[test]
fn record_initial_adds_a_step_zero_sample_per_stage() {
    let without = all_samples(&mut two_stage_bar(false, None, true));
    assert_eq!(without.len(), 4);

    let samples = all_samples(&mut two_stage_bar(true, None, true));
    assert_eq!(samples.len(), 6);
    let stage = |i: usize| samples.iter().filter(|s| s.0 == i).collect::<Vec<_>>();
    let (first, second) = (stage(0), stage(1));
    assert_eq!((first.len(), second.len()), (3, 3));
    assert_eq!((first[0].1, first[0].2), (0.0, 0.0));
    let end_of_first = first[2].2;
    assert!((end_of_first - 60.0 * 100.0 / (2.0 * 30000.0)).abs() < 1e-9);
    assert_eq!(second[0].1, 0.0);
    assert!(
        (second[0].2 - end_of_first).abs() < 1e-12,
        "stage 2 starts where stage 1 ended"
    );
}

/// A `Reset` stage reverts displacements to the as-built state (a later stage starts from zero, not
/// from where the previous one ended). Like OpenSees' `reset`, a pattern frozen by `holdPatternsAfter`
/// stays frozen at its final factor: the next stage applies it in full from its first step.
#[test]
fn reset_stage_reverts_the_state_but_keeps_held_patterns_applied() {
    let reset = Some(StageSpec::Reset {
        id: "reset".to_string(),
    });
    let samples = all_samples(&mut two_stage_bar(true, reset, true));
    // Reset records nothing: stage indices are first = 0, second = 2.
    let second: Vec<_> = samples.iter().filter(|s| s.0 == 2).collect();
    assert_eq!(second.len(), 3);
    assert_eq!(
        second[0].2, 0.0,
        "second stage starts from the as-built state"
    );
    let end_of_first = samples.iter().filter(|s| s.0 == 0).last().unwrap().2;
    assert!(
        (second[1].2 - end_of_first).abs() < 1e-12,
        "the held load is applied in full at once"
    );
    assert!((second[2].2 - end_of_first).abs() < 1e-12);
}

/// Without a hold, the same pattern ramps again from zero after a `Reset`.
#[test]
fn reset_stage_lets_an_unheld_pattern_ramp_again() {
    let reset = Some(StageSpec::Reset {
        id: "reset".to_string(),
    });
    let samples = all_samples(&mut two_stage_bar(true, reset, false));
    let second: Vec<_> = samples.iter().filter(|s| s.0 == 2).collect();
    let end_of_first = samples.iter().filter(|s| s.0 == 0).last().unwrap().2;
    assert_eq!(second[0].2, 0.0);
    assert!(
        (second[1].2 - end_of_first / 2.0).abs() < 1e-12,
        "half the load at the first of two steps"
    );
    assert!((second[2].2 - end_of_first).abs() < 1e-12);
}
