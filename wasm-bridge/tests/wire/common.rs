//! Shared builders for the wire-format tests.

use carapace_wasm::input_v1::{CarapaceInputV1, StepOutcome};

/// A decodable input with only a header; tests fill in the tables they need.
pub fn empty_input(ndm: u8) -> CarapaceInputV1 {
    serde_json::from_value(serde_json::json!({
        "header": { "schemaVersion": 1, "ndm": ndm, "engineVersion": "test" }
    }))
    .unwrap()
}

pub fn last_sample(outcome: &StepOutcome, recorder_index: usize) -> Option<(f64, f64)> {
    outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == recorder_index)
        .and_then(|batch| batch.samples.last())
        .copied()
}
