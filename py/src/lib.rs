//! Python extension for Carapace: JSON model text in, JSON report text out. The `carapace` Python package
//! wraps this with a model builder and numpy results; keeping the boundary a string keeps it one tiny,
//! version-proof surface that mirrors the CLI exactly.

use carapace_input::input_v1::{decode, CarapaceInputV1, RunOptions, Session};
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

/// Decodes a model; returns `{"valid": true, "ndm": n}` or `{"valid": false, "errors": [...]}`.
#[pyfunction]
fn validate(model_json: &str) -> String {
    match load(model_json) {
        Ok(session) => {
            let ndm = if matches!(session, Session::D2(_)) { 2 } else { 3 };
            serde_json::json!({ "valid": true, "ndm": ndm }).to_string()
        }
        Err(report) => report,
    }
}

/// Decodes and runs a model with the GIL released; returns the CLI's `run` document.
#[pyfunction]
#[pyo3(signature = (model_json, step_budget = 64))]
fn run(py: Python<'_>, model_json: &str, step_budget: u32) -> PyResult<String> {
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

#[pymodule]
fn _carapace(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(validate, m)?)?;
    m.add_function(wrap_pyfunction!(run, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
