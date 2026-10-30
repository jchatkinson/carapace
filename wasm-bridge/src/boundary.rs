//! The `wasm_bindgen` boundary for `input_v1` (`CarapaceInputV1`
//! decode/`Session`) API. First pass: a real JS object/array, deserialized
//! via `serde-wasm-bindgen`, not yet the transferable-typed-array wire
//! format intended for `postMessage` — see
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
use crate::material_probe::{
    MaterialProbe, MaterialProbeConfig, MaterialProbeError, MaterialProbeResponse,
};
use crate::section_probe::{
    AxialExtremes, IsolineRequest, IsolineResult, MomentCurvatureRequest, MomentCurvatureResult,
    SectionProbe, SectionProbeConfig, SectionProbeError, SectionResponse,
};

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
    /// loadFactor, error, recorderBatches, bulkBatches, ... }` object. `recorderBatches` holds only the samples
    /// this call produced (results-storage-indexeddb.md's `recorderBatch`, one entry per
    /// recorder that recorded this call), not the whole run's history. There is no separate
    /// "samples so far" accessor: a caller that needs the full history accumulates these
    /// batches itself, same as the planned results-storage worker will. `bulkBatches` carry their
    /// `progress` and `values` as `Float64Array`s, each on its own buffer, so a worker can
    /// transfer them to the main thread without a copy.
    pub fn advance(&mut self, step_budget: u32) -> Result<Ts<input_v1::StepOutcome>, JsValue> {
        let mut outcome = self.0.advance(step_budget);
        // Bulk frames are the one large payload: leave them out of the serde pass and attach them
        // as `Float64Array`s, each a fresh buffer the caller can transfer.
        let bulk_batches = std::mem::take(&mut outcome.bulk_batches);
        let js: JsValue = Tsify::into_js(&outcome)
            .map_err(|error| JsValue::from_str(&error.to_string()))?
            .into();
        let batches = js_sys::Array::new();
        for batch in &bulk_batches {
            let object = js_sys::Object::new();
            let numbers = [
                ("recorderIndex", batch.recorder_index as f64),
                ("stageIndex", batch.stage_index as f64),
                ("firstSample", f64::from(batch.first_sample)),
                ("elements", f64::from(batch.elements)),
                ("points", f64::from(batch.points)),
                ("components", f64::from(batch.components)),
            ];
            for (key, value) in numbers {
                js_sys::Reflect::set(&object, &key.into(), &value.into())?;
            }
            let arrays = [("progress", &batch.progress), ("values", &batch.values)];
            for (key, data) in arrays {
                let array = js_sys::Float64Array::from(data.as_slice());
                js_sys::Reflect::set(&object, &key.into(), &array)?;
            }
            batches.push(&object);
        }
        js_sys::Reflect::set(&js, &"bulkBatches".into(), &batches)?;
        Ok(Ts::new_unchecked(js))
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

    /// What decoding resolved: the embedded-bar segment table and each bar group's TCM/POM
    /// routing and crack spacing (`{ embeddedBars: EmbeddedBarGroup[] }`). Available right after
    /// `decodeInput`; empty when the model has no `embeddedBars`.
    #[wasm_bindgen(js_name = modelOutput)]
    pub fn model_output(&self) -> Result<Ts<input_v1::ModelOutput>, JsValue> {
        self.0
            .model_output()
            .clone()
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

/// Opaque handle to a fiber-section probe (`section_probe`). Evaluation is stateless: nothing here commits
/// material history, so the same strain state always returns the same resultants.
#[wasm_bindgen]
pub struct WasmSectionProbe(SectionProbe);

#[wasm_bindgen(js_name = createSectionProbe)]
pub fn create_section_probe(config: Ts<SectionProbeConfig>) -> Result<WasmSectionProbe, JsValue> {
    let config: SectionProbeConfig = config
        .to_rust()
        .map_err(|error| JsValue::from_str(&format!("malformed SectionProbeConfig: {error}")))?;
    SectionProbe::new(&config)
        .map(WasmSectionProbe)
        .map_err(section_error)
}

fn section_error(error: SectionProbeError) -> JsValue {
    serde_wasm_bindgen::to_value(&error)
        .unwrap_or_else(|_| JsValue::from_str("section probe error"))
}

fn ts_error(error: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&error.to_string())
}

#[wasm_bindgen]
impl WasmSectionProbe {
    /// `{ n, mz, my, m2, m3 }` at `(eps0, kz, ky)`; strain is `eps0 - y kz + z ky`.
    pub fn evaluate(&self, eps0: f64, kz: f64, ky: f64) -> Result<Ts<SectionResponse>, JsValue> {
        self.0.evaluate(eps0, kz, ky).into_ts().map_err(ts_error)
    }

    /// `(eps0, kz, ky)` triples in, `(n, mz, my)` triples out; with `withGroups` each state also returns the
    /// per-group `(n, mz, my)` after the total.
    #[wasm_bindgen(js_name = evaluateBatch)]
    pub fn evaluate_batch(&self, states: &[f64], with_groups: bool) -> Vec<f64> {
        self.0.evaluate_batch(states, with_groups)
    }

    /// `[strain, stress]` for every fiber, flattened.
    #[wasm_bindgen(js_name = fiberResponses)]
    pub fn fiber_responses(&self, eps0: f64, kz: f64, ky: f64) -> Vec<f64> {
        self.0.fiber_responses(eps0, kz, ky)
    }

    #[wasm_bindgen(js_name = fiberCount)]
    pub fn fiber_count(&self) -> usize {
        self.0.fiber_count()
    }

    #[wasm_bindgen(js_name = groupCount)]
    pub fn group_count(&self) -> usize {
        self.0.group_count()
    }

    /// `(eps0, kz, ky)` for every (angle, profile) pair, angles outermost; `profiles` is `(et, eb)` pairs.
    #[wasm_bindgen(js_name = profileStates)]
    pub fn profile_states(&self, angles: &[f64], profiles: &[f64]) -> Result<Vec<f64>, JsValue> {
        self.0
            .profile_states(angles, profiles)
            .map_err(section_error)
    }

    /// Maximum compressive and tensile axial force over uniform strain in `[-maxStrain, maxStrain]`.
    #[wasm_bindgen(js_name = axialExtremes)]
    pub fn axial_extremes(&self, max_strain: f64) -> Result<Ts<AxialExtremes>, JsValue> {
        self.0
            .axial_extremes(max_strain)
            .map_err(section_error)?
            .into_ts()
            .map_err(ts_error)
    }

    /// Samples the PMM boundary for `solveIsolines`; returns the number of sampled states.
    #[wasm_bindgen(js_name = buildGrid)]
    pub fn build_grid(&mut self, angles: &[f64], profiles: &[f64]) -> Result<usize, JsValue> {
        self.0.build_grid(angles, profiles).map_err(section_error)
    }

    /// `[min, max]` axial force over the sampled grid, or empty before `buildGrid`.
    #[wasm_bindgen(js_name = gridAxialRange)]
    pub fn grid_axial_range(&self) -> Vec<f64> {
        self.0
            .grid_axial_range()
            .map(|(lo, hi)| vec![lo, hi])
            .unwrap_or_default()
    }

    #[wasm_bindgen(js_name = solveIsolines)]
    pub fn solve_isolines(
        &self,
        request: Ts<IsolineRequest>,
    ) -> Result<Ts<IsolineResult>, JsValue> {
        let request: IsolineRequest = request
            .to_rust()
            .map_err(|error| JsValue::from_str(&format!("malformed IsolineRequest: {error}")))?;
        IsolineResult {
            points: self.0.solve_isolines(&request).map_err(section_error)?,
        }
        .into_ts()
        .map_err(ts_error)
    }

    #[wasm_bindgen(js_name = momentCurvature)]
    pub fn moment_curvature(
        &self,
        request: Ts<MomentCurvatureRequest>,
    ) -> Result<Ts<MomentCurvatureResult>, JsValue> {
        let request: MomentCurvatureRequest = request.to_rust().map_err(|error| {
            JsValue::from_str(&format!("malformed MomentCurvatureRequest: {error}"))
        })?;
        self.0
            .moment_curvature(&request)
            .map_err(section_error)?
            .into_ts()
            .map_err(ts_error)
    }
}
