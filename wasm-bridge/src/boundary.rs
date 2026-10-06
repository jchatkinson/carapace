//! The `wasm_bindgen` boundary for M10's `input_v1` (`CarapaceInputV1`
//! decode/`Session`) API. First pass: a real JS object/array, deserialized
//! via `serde-wasm-bindgen`, not yet the transferable-typed-array wire
//! format docs/pysees-handoff.md specifies for `postMessage` — see
//! `input_v1`'s module doc comment for why that's deferred. This is
//! enough to drive `decode`/`Session::advance` from real JS (or `pysees`)
//! code end to end; the zero-copy transfer optimization can replace
//! `serde-wasm-bindgen`'s JSON-ish decoding later without changing this
//! module's shape (`decodeInput` in, `WasmSession` methods out).
//! `Ts<T>` connects signatures to generated TypeScript definitions while
//! keeping data conversion fallible inside the function body.

use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;

use crate::input_v1::{self, CarapaceInputV1};
use crate::material_probe::{MaterialProbe, MaterialProbeConfig, MaterialProbeError, MaterialProbeResponse};

/// Opaque handle to a decoded, steppable analysis session — the `Session`
/// enum itself can't cross the boundary directly, since `wasm_bindgen`
/// requires an exported type to be a plain struct.
#[wasm_bindgen]
pub struct WasmSession(input_v1::Session);

/// Decodes a `CarapaceInputV1`-shaped JS value (see `input_v1`'s table
/// doc comments for the exact field names — `#[serde(rename_all =
/// "camelCase")]` throughout) into a steppable [`WasmSession`]. Rejects
/// with a `{ kind: "...", ... }`-shaped JS error object on any
/// [`input_v1::DecodeError`], the same tagged shape a Rust caller would
/// match on.
#[wasm_bindgen(js_name = decodeInput)]
pub fn decode_input(value: Ts<CarapaceInputV1>) -> Result<WasmSession, JsValue> {
    let input: CarapaceInputV1 = value
        .to_rust()
        .map_err(|error| JsValue::from_str(&format!("malformed CarapaceInputV1: {error}")))?;
    let session = input_v1::decode(input).map_err(to_js_error)?;
    Ok(WasmSession(session))
}

#[wasm_bindgen]
impl WasmSession {
    /// See `input_v1::Session::advance`. Returns a `{ done, stageComplete, stepsTaken,
    /// loadFactor, error, recorderBatches }` object — `recorderBatches` holds only the samples
    /// this call produced (results-storage-indexeddb.md's `recorderBatch`, one entry per
    /// recorder that recorded this call), not the whole run's history. There is no separate
    /// "samples so far" accessor: a caller that needs the full history accumulates these
    /// batches itself, same as the planned results-storage worker will.
    pub fn advance(&mut self, step_budget: u32) -> Result<Ts<input_v1::StepOutcome>, JsValue> {
        let outcome = self.0.advance(step_budget);
        outcome
            .into_ts()
            .map_err(|error| JsValue::from_str(&error.to_string()))
    }

    /// Frequencies, mode shapes and participation of every `Modal` stage finished so far
    /// (`{ stages: ModalStageResult[] }`); empty until a modal stage completes. Read it after
    /// `advance` reports `done`, or whenever a stage has completed.
    #[wasm_bindgen(js_name = modalResults)]
    pub fn modal_results(&self) -> Result<Ts<input_v1::ModalResultsReport>, JsValue> {
        self.0
            .modal_results()
            .into_ts()
            .map_err(|error| JsValue::from_str(&error.to_string()))
    }

    /// The `AnalysisSequence` stage id `advance` is currently in (or, once
    /// done, whichever stage stopped it) — `undefined` once every stage has
    /// completed.
    #[wasm_bindgen(js_name = currentStageId)]
    pub fn current_stage_id(&self) -> Option<String> {
        self.0.current_stage_id().map(str::to_owned)
    }
}

fn to_js_error(error: input_v1::DecodeError) -> JsValue {
    serde_wasm_bindgen::to_value(&error).unwrap_or_else(|_| {
        JsValue::from_str("decode error (failed to serialize DecodeError itself)")
    })
}

/// Opaque handle to a persistent uniaxial material probe (`material_probe`).
/// Material state persists across `applyStrain` calls; cancel between points
/// by simply not calling it again.
#[wasm_bindgen]
pub struct WasmMaterialProbe(MaterialProbe);

#[wasm_bindgen(js_name = createMaterialProbe)]
pub fn create_material_probe(
    config: Ts<MaterialProbeConfig>,
) -> Result<WasmMaterialProbe, JsValue> {
    let config: MaterialProbeConfig = config
        .to_rust()
        .map_err(|error| JsValue::from_str(&format!("malformed MaterialProbeConfig: {error}")))?;
    MaterialProbe::new(&config)
        .map(WasmMaterialProbe)
        .map_err(probe_error)
}

#[wasm_bindgen]
impl WasmMaterialProbe {
    /// Imposes the strain and returns `{ strain, stress }`; rejects with a
    /// `{ kind, ... }` `MaterialProbeError`.
    #[wasm_bindgen(js_name = applyStrain)]
    pub fn apply_strain(&mut self, target: f64) -> Result<Ts<MaterialProbeResponse>, JsValue> {
        let response = self.0.apply_strain(target).map_err(probe_error)?;
        response
            .into_ts()
            .map_err(|error| JsValue::from_str(&error.to_string()))
    }

    pub fn reset(&mut self) -> Result<(), JsValue> {
        self.0.reset().map_err(probe_error)
    }
}

fn probe_error(error: MaterialProbeError) -> JsValue {
    serde_wasm_bindgen::to_value(&error)
        .unwrap_or_else(|_| JsValue::from_str("material probe error"))
}
