//! `embeddedBars[i].bond` through the wire (bond plan, WP2): decode and validation, `modelOutput().bond`, the
//! `zeroLengths` bulk recorder, and the perfect-bond limit against a bar without `bond`.

use crate::common::empty_input;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, BulkQuantity, BulkTable, ConvergenceSpec, IntegratorSpec, RecorderSpec,
    SequenceSpec, StageSpec,
};
use carapace_wasm::input_v1::{decode, CarapaceInputV1, DecodeError};
use serde_json::json;

/// Two 100 x 100 mm quads side by side, left edge fixed, 5 kN on each right node, one bar along y = 50.
fn tie_model(bond: serde_json::Value) -> CarapaceInputV1 {
    let mut input = serde_json::to_value(empty_input(2)).unwrap();
    input["nodes"] = json!({
        "coords": [0.0, 0.0, 100.0, 0.0, 200.0, 0.0, 0.0, 100.0, 100.0, 100.0, 200.0, 100.0],
        "fixed": [3, 0, 0, 3, 0, 0],
        "massNodeIndex": [], "mass": [],
    });
    input["planeMaterials"] =
        json!([{ "kind": "isotropic", "e": 30000.0, "nu": 0.2, "state": "planeStress" }]);
    input["materials"] = json!([{ "kind": "elastic", "e": 200000.0 }]);
    input["quads"] = json!({
        "nodeIds": [0, 1, 4, 3, 1, 2, 5, 4],
        "thickness": [100.0, 100.0], "material": [0, 0], "density": [0.0, 0.0],
        "formulation": ["full", "full"],
    });
    let have_springs = !bond.is_null();
    let mut bar = json!({ "polyline": [0.0, 50.0, 200.0, 50.0], "area": 201.0, "material": 0 });
    if !bond.is_null() {
        bar["bond"] = bond;
    }
    input["embeddedBars"] = json!([bar]);
    input["loadPatterns"] =
        json!({ "series": [{ "kind": "linear", "slope": 1.0 }], "scaleFactor": [1.0] });
    input["nodalLoads"] = json!({ "pattern": [0, 0], "node": [2, 5], "dof": [0, 0], "value": [5000.0, 5000.0], "stage": [0, 0] });
    let mut input: CarapaceInputV1 = serde_json::from_value(input).unwrap();
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "load".to_string(),
            steps: 2,
            integrator: IntegratorSpec::LoadControl { increment: 0.5 },
            algorithm: AlgorithmSpec::Newton {
                tangent: Default::default(),
                line_search: None,
            },
            convergence: Some(ConvergenceSpec::Combined {
                force_tol: 1e-6,
                moment_tol: None,
                relative_tol: Some(1e-8),
                displacement_tol: None,
                max_iter: 30,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: (0..6)
            .map(|node| RecorderSpec::NodeDisp { node, dof: 0 })
            .chain(
                have_springs
                    .then_some([
                        RecorderSpec::Bulk {
                            table: BulkTable::ZeroLengths,
                            quantity: BulkQuantity::Deformation,
                        },
                        RecorderSpec::Bulk {
                            table: BulkTable::ZeroLengths,
                            quantity: BulkQuantity::Force,
                        },
                    ])
                    .into_iter()
                    .flatten(),
            )
            .collect(),
    };
    input
}

fn bond(gb: f64, fbd: f64) -> serde_json::Value {
    json!({ "perimeter": 50.0, "fbd": fbd, "gb": gb, "slipLimit": 10.0 })
}

fn run(input: CarapaceInputV1) -> carapace_wasm::input_v1::StepOutcome {
    let mut session = decode(input).expect("decodes");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    outcome
}

#[test]
fn model_output_reports_the_bond_points_and_springs() {
    let session = decode(tie_model(bond(400.0, 3.0))).unwrap();
    let group = &session.model_output().embedded_bars[0];
    assert_eq!(group.segments.len(), 2, "split at the shared edge");
    let points = &group.bond.as_ref().expect("bond group").points;
    assert_eq!(points.len(), 3);
    assert_eq!(points[0].position, [0.0, 50.0]);
    assert_eq!(
        points.iter().map(|p| p.tributary).collect::<Vec<_>>(),
        vec![50.0, 100.0, 50.0]
    );
    assert_eq!(
        points.iter().map(|p| p.spring).collect::<Vec<_>>(),
        vec![Some(0), Some(1), Some(2)]
    );
    assert!(points.iter().all(|p| !p.locked && p.end_spring.is_none()));
}

#[test]
fn a_stiff_bond_reproduces_the_bar_without_bond() {
    let plain = run(tie_model(serde_json::Value::Null));
    let bonded = run(tie_model(bond(1e7, 1e9)));
    for node in 0..6 {
        let (_, a) = crate::common::last_sample(&plain, node).unwrap();
        let (_, b) = crate::common::last_sample(&bonded, node).unwrap();
        assert!(
            (a - b).abs() <= 1e-3 * a.abs().max(1e-9),
            "node {node}: {a} vs {b}"
        );
    }
}

#[test]
fn the_zero_lengths_recorder_reports_slip_and_bond_force() {
    let (gb, fbd) = (400.0, 3.0);
    let outcome = run(tie_model(bond(gb, fbd)));
    let batch = |index: usize| {
        outcome
            .bulk_batches
            .iter()
            .rev()
            .find(|b| b.recorder_index == index)
            .expect("batch")
    };
    let (deformation, force) = (batch(6), batch(7));
    assert_eq!(
        (
            deformation.elements,
            deformation.points,
            deformation.components
        ),
        (3, 1, 1)
    );
    let frame = |values: &[f64]| values[values.len() - 3..].to_vec();
    let (s, f) = (frame(&deformation.values), frame(&force.values));
    // Tributary lengths 50, 100, 50 mm; elastic until the spring yields at fbd P l.
    for (k, l) in [50.0, 100.0, 50.0].into_iter().enumerate() {
        let (yield_force, stiffness) = (fbd * 50.0 * l, gb * 50.0 * l);
        let expected = (stiffness * s[k]).clamp(-yield_force, yield_force);
        assert!(
            (f[k] - expected).abs() <= 2e-2 * yield_force.max(1.0),
            "spring {k}: {} vs {expected}",
            f[k]
        );
    }
    assert!(s.iter().any(|v| v.abs() > 0.0), "a loaded bar slips: {s:?}");
}

#[test]
fn invalid_bond_data_is_rejected() {
    for bad in [
        json!({ "perimeter": 50.0, "fbd": -1.0, "gb": 400.0, "slipLimit": 10.0 }),
        json!({ "perimeter": 50.0, "fbd": 3.0, "gb": 0.0, "slipLimit": 10.0 }),
        json!({ "perimeter": 50.0, "fbd": 3.0, "gb": 400.0, "slipLimit": 10.0, "rf": 2.0 }),
        json!({ "perimeter": 50.0, "fbd": 3.0, "gb": 400.0, "slipLimit": 10.0, "start": { "fau": 0.0, "ka": 1.0 } }),
    ] {
        match decode(tie_model(bad.clone())) {
            Err(DecodeError::InvalidRow { table, .. }) => assert_eq!(table, "embedded_bars"),
            other => panic!("{bad}: expected InvalidRow, got {:?}", other.err()),
        }
    }
}
