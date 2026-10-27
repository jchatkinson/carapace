//! `loadControlCutbacks` over the wire (plan `docs/csfm-plan.md` §2.5 point 1, §4 I): a nominal
//! stage whose increments cross a regime boundary is cut back, still reaches its target load
//! factor, records every accepted substep, and `advance(budget)` counts accepted substeps.
//!
//! The model: a truss (k = 50) in parallel with two elastic-perfectly-plastic springs (k = 60,
//! yield strain 0.01; k = 40, yield strain 0.012) sharing one DOF, loaded by a reference force of
//! 5. The springs yield at load factors 0.3 and 0.336. With `max_iter: 2`, a step that crosses one
//! kink converges (one iteration to find it, one to re-form the tangent), but a step that crosses
//! both needs a third and fails, so it is halved until it crosses at most one.

use crate::common::empty_input;
use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{Domain, Element, Material, Node, Truss, ZeroLength};
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
    TangentStrategySpec,
};
use carapace_wasm::input_v1::{decode, AnalysisErrorDetail, CarapaceInputV1, DecodeError};
use serde_json::json;

const INCREMENT: f64 = 0.07;
const STEPS: u32 = 20;
const MIN_INCREMENT: f64 = 0.001;
const MAX_ITER: u32 = 2;

pub fn input(increment: f64, min_increment: f64) -> CarapaceInputV1 {
    let mut input = serde_json::to_value(empty_input(2)).unwrap();
    input["nodes"] = json!({
        "coords": [0.0, 0.0, 100.0, 0.0], "fixed": [7, 6], "massNodeIndex": [], "mass": [],
    });
    input["materials"] = json!([
        { "kind": "elastic", "e": 5000.0 },
        { "kind": "elasticPp", "e": 60.0, "eyp": 0.01 },
        { "kind": "elasticPp", "e": 40.0, "eyp": 0.012 },
        { "kind": "parallel", "children": [[1, 1.0], [2, 1.0]] },
    ]);
    input["trusses"] = json!({
        "nodeI": [0], "nodeJ": [1], "area": [1.0], "material": [0], "density": [0.0],
    });
    input["zeroLengths"] = json!({ "nodeI": [0], "nodeJ": [1], "materials": [[0, 0, 3]] });
    input["loadPatterns"] = json!({
        "series": [{ "kind": "linear", "slope": 1.0 }], "scaleFactor": [1.0],
    });
    input["nodalLoads"] = json!({
        "pattern": [0], "node": [1], "dof": [0], "value": [5.0], "stage": [0],
    });
    let mut input: CarapaceInputV1 = serde_json::from_value(input).unwrap();
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "ramp".to_string(),
            steps: STEPS,
            integrator: IntegratorSpec::LoadControlCutbacks {
                increment,
                min_increment,
            },
            algorithm: AlgorithmSpec::Newton {
                tangent: TangentStrategySpec::default(),
                line_search: None,
            },
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-10,
                max_iter: MAX_ITER,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
    };
    input
}

/// The same ramp built natively: accepted `(load_factor, displacement)` per substep until the
/// nominal `STEPS * INCREMENT` is covered.
fn native_substeps() -> Vec<(f64, f64)> {
    let mut domain: Domain = Domain::new();
    let i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let j = domain.add_node(Node::new([100.0, 0.0]).fix(1).fix(2));
    domain.load_node(j, 0, 5.0);
    domain.add_element(Element::Truss(Truss::new(
        i,
        j,
        1.0,
        Material::Elastic { e: 5000.0 },
    )));
    domain.add_element(Element::ZeroLength(ZeroLength::new(i, j).with_material(
        0,
        Material::Parallel(vec![
            (Material::elastic_pp(60.0, 0.01), 1.0),
            (Material::elastic_pp(40.0, 0.012), 1.0),
        ]),
    )));
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControlCutbacks {
            increment: INCREMENT,
            min_increment: MIN_INCREMENT,
        })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: MAX_ITER as usize,
        })
        .build(domain);
    let mut out = Vec::new();
    let mut nominal_steps_done = 0;
    while nominal_steps_done < STEPS {
        let result = analysis.step().expect("native ramp converges");
        let info = result.cutbacks.unwrap();
        if info.accepted_increment >= info.nominal_increment {
            nominal_steps_done += 1;
        }
        out.push((
            result.load_factor,
            analysis.domain().node(j).displacement[0],
        ));
    }
    out
}

#[test]
fn a_cutback_stage_reaches_its_target_and_records_every_accepted_substep() {
    let expected = native_substeps();
    assert!(
        expected.len() > STEPS as usize,
        "the ramp must need cutbacks to exercise the contract"
    );
    let mut session = decode(input(INCREMENT, MIN_INCREMENT)).unwrap();
    let outcome = session.advance(1000);
    assert!(outcome.done && outcome.stage_complete && outcome.error.is_none());
    assert_eq!(outcome.steps_taken as usize, expected.len());
    assert!((outcome.load_factor - STEPS as f64 * INCREMENT).abs() < 1e-12);
    assert_eq!(outcome.load_factor, expected.last().unwrap().0);

    let samples = &outcome.recorder_batches[0].samples;
    assert_eq!(samples.len(), expected.len());
    for (k, (&(_, got), &(_, want))) in samples.iter().zip(&expected).enumerate() {
        assert_eq!(got, want, "substep {k}");
    }
}

/// `stepDiagnostics` has one entry per accepted substep, ending where the recorders do, and its
/// `cutback` block carries the native `CutbackInfo`: some substeps are cut back (`retries > 0`),
/// every one asked for at most the nominal increment, and the accepted increments cover
/// `STEPS * INCREMENT`.
#[test]
fn step_diagnostics_report_each_substeps_cutbacks() {
    let expected = native_substeps();
    let mut session = decode(input(INCREMENT, MIN_INCREMENT)).unwrap();
    let outcome = session.advance(1000);
    let diagnostics = &outcome.step_diagnostics;
    assert_eq!(diagnostics.len(), expected.len());
    for (diagnostic, &(load_factor, _)) in diagnostics.iter().zip(&expected) {
        assert_eq!(diagnostic.progress, load_factor);
        assert!(diagnostic.iterations >= 1 && diagnostic.factorizations >= 1);
    }
    let cutbacks: Vec<_> = diagnostics.iter().map(|d| d.cutback.unwrap()).collect();
    assert!(cutbacks.iter().any(|c| c.retries > 0), "{cutbacks:?}");
    assert!(cutbacks
        .iter()
        .all(|c| c.accepted_increment <= c.nominal_increment + 1e-15
            && c.nominal_increment <= INCREMENT + 1e-15));
    let covered: f64 = cutbacks.iter().map(|c| c.accepted_increment).sum();
    assert!(
        (covered - STEPS as f64 * INCREMENT).abs() < 1e-12,
        "{covered}"
    );
}

#[test]
fn advance_budget_counts_accepted_substeps() {
    let expected = native_substeps();
    let mut session = decode(input(INCREMENT, MIN_INCREMENT)).unwrap();
    let mut calls = 0;
    let mut recorded = 0;
    loop {
        let outcome = session.advance(1);
        calls += 1;
        assert_eq!(
            outcome.steps_taken, 1,
            "one accepted substep per budget unit"
        );
        recorded += outcome
            .recorder_batches
            .iter()
            .map(|b| b.samples.len())
            .sum::<usize>();
        if outcome.done {
            assert!(outcome.stage_complete);
            break;
        }
        assert!(calls < 1000);
    }
    assert_eq!(calls, expected.len());
    assert_eq!(recorded, expected.len());
}

#[test]
fn invalid_cutback_options_are_decode_errors() {
    for (increment, min_increment, field) in [
        (0.0, 0.001, "integrator.increment"),
        (f64::NAN, 0.001, "integrator.increment"),
        (INCREMENT, 0.0, "integrator.minIncrement"),
        (INCREMENT, INCREMENT * 2.0, "integrator.minIncrement"),
    ] {
        assert_eq!(
            decode(input(increment, min_increment)).unwrap_err(),
            DecodeError::InvalidAnalysisOption {
                stage: "ramp".to_string(),
                field,
            },
            "{increment} {min_increment}"
        );
    }
}

#[test]
fn exhausted_cutbacks_report_the_last_accepted_load_factor() {
    // Nominal steps this large never converge in `MAX_ITER` iterations; a floor equal to the
    // nominal size disables retries, so the first failure exhausts them.
    let mut session = decode(input(INCREMENT, INCREMENT)).unwrap();
    let outcome = session.advance(1000);
    assert!(outcome.done && !outcome.stage_complete);
    let Some(AnalysisErrorDetail::LoadControlExhausted {
        accepted_load_factor,
        ..
    }) = outcome.error
    else {
        panic!("{:?}", outcome.error);
    };
    assert!(accepted_load_factor < STEPS as f64 * INCREMENT);
}

/// The cutback model as built by this module's tests, for other modules that need a model with no CSFM limit state.
pub fn default_input() -> CarapaceInputV1 {
    input(INCREMENT, MIN_INCREMENT)
}
