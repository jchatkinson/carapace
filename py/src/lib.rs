//! Python extension for Carapace: JSON model text in, JSON report text out. The `carapace` Python package
//! wraps this with a model builder and numpy results; keeping the boundary a string keeps it one tiny,
//! version-proof surface that mirrors the CLI exactly.

use carapace_core::analysis::{set_parallelism, Parallelism};
use carapace_input::input_v1::{
    decode, run_many as run_many_inputs, CarapaceInputV1, RunOptions, Session,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// `Err(json)` is a `{valid: false, errors: [...]}` document, `Ok` a decoded session.
fn load(model_json: &str) -> Result<Session, String> {
    let mut de = serde_json::Deserializer::from_str(model_json);
    let input: CarapaceInputV1 = serde_path_to_error::deserialize(&mut de).map_err(|e| {
        serde_json::json!({ "valid": false, "errors": [{
            "stage": "parse", "path": e.path().to_string(), "detail": e.inner().to_string(),
        }] })
        .to_string()
    })?;
    decode(input).map_err(|e| {
        serde_json::json!({ "valid": false, "errors": [{
            "stage": "decode", "detail": serde_json::to_value(&e).unwrap_or_default(),
        }] })
        .to_string()
    })
}

fn parse(model_json: &str) -> Result<CarapaceInputV1, String> {
    let mut de = serde_json::Deserializer::from_str(model_json);
    serde_path_to_error::deserialize(&mut de).map_err(|e| {
        serde_json::json!({ "valid": false, "errors": [{
            "stage": "parse", "path": e.path().to_string(), "detail": e.inner().to_string(),
        }] })
        .to_string()
    })
}

/// Decodes a model; returns `{"valid": true, "ndm": n}` or `{"valid": false, "errors": [...]}`.
#[pyfunction]
fn validate(model_json: &str) -> String {
    match load(model_json) {
        Ok(session) => {
            let ndm = if matches!(session, Session::D2(_)) {
                2
            } else {
                3
            };
            serde_json::json!({ "valid": true, "ndm": ndm }).to_string()
        }
        Err(report) => report,
    }
}

/// Decodes and runs a model with the GIL released; returns the CLI's `run` document.
#[pyfunction]
#[pyo3(signature = (model_json, step_budget = 64, threads = 1))]
fn run(py: Python<'_>, model_json: &str, step_budget: u32, threads: usize) -> PyResult<String> {
    if threads != 1 {
        set_parallelism(Parallelism::Threads(threads));
    }
    py.detach(|| match load(model_json) {
        Ok(mut session) => {
            let report = session.run(
                &RunOptions {
                    step_budget,
                    cancel: None,
                },
                |_| {},
            );
            serde_json::to_string(&serde_json::json!({ "valid": true, "report": report }))
                .map_err(|e| PyValueError::new_err(e.to_string()))
        }
        Err(report) => Ok(report),
    })
}

/// Runs independent models concurrently (`threads = 0` for all cores); one result document per model, in order.
#[pyfunction]
#[pyo3(signature = (models, threads = 0))]
fn run_many(py: Python<'_>, models: Vec<String>, threads: usize) -> PyResult<Vec<String>> {
    py.detach(|| {
        let mut out = vec![String::new(); models.len()];
        let mut inputs = Vec::new();
        let mut slots = Vec::new();
        for (at, text) in models.iter().enumerate() {
            match parse(text) {
                Ok(input) => {
                    inputs.push(input);
                    slots.push(at);
                }
                Err(report) => out[at] = report,
            }
        }
        let options = RunOptions::default();
        for (at, outcome) in slots
            .into_iter()
            .zip(run_many_inputs(inputs, threads, &options))
        {
            out[at] = match outcome {
                Ok(report) => {
                    serde_json::to_string(&serde_json::json!({ "valid": true, "report": report }))
                        .map_err(|e| PyValueError::new_err(e.to_string()))?
                }
                Err(e) => serde_json::json!({ "valid": false, "errors": [{
                    "stage": "decode", "detail": serde_json::to_value(&e).unwrap_or_default(),
                }] })
                .to_string(),
            };
        }
        Ok(out)
    })
}

/// A decoded model that is stepped from Python: the stateful half behind `analyze(n)` and the result
/// queries of the OpenSees-style `Model`.
#[pyclass(name = "Session")]
struct PySession {
    inner: Session,
}

#[pymethods]
impl PySession {
    /// Decodes a model. An invalid one raises `ValueError` whose message is the `{valid: false, errors}` JSON.
    #[new]
    fn new(model_json: &str) -> PyResult<Self> {
        load(model_json)
            .map(|inner| Self { inner })
            .map_err(PyValueError::new_err)
    }

    /// Takes up to `step_budget` steps (GIL released); returns the `StepOutcome` as JSON.
    fn advance(&mut self, py: Python<'_>, step_budget: u32) -> PyResult<String> {
        let session = &mut self.inner;
        py.detach(|| serde_json::to_string(&session.advance(step_budget)))
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

#[pymodule]
fn _carapace(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(validate, m)?)?;
    m.add_function(wrap_pyfunction!(run, m)?)?;
    m.add_function(wrap_pyfunction!(run_many, m)?)?;
    m.add_class::<PySession>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
