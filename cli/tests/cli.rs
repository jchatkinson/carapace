//! End-to-end checks of the `carapace` binary: exit codes and the JSON it prints.

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_carapace");
const EXAMPLE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../examples/models/portal_frame_2d.json"
);

fn run(args: &[&str]) -> (Output, serde_json::Value) {
    let output = Command::new(BIN).args(args).output().unwrap();
    let json = serde_json::from_slice(&output.stdout).unwrap_or(serde_json::Value::Null);
    (output, json)
}

/// The portal frame example validates and runs; the lateral drift matches a hand check of sign and size.
#[test]
fn runs_the_example_model() {
    let (output, json) = run(&["validate", EXAMPLE]);
    assert!(output.status.success());
    assert_eq!(json["valid"], true);

    let (output, json) = run(&["run", EXAMPLE]);
    assert!(output.status.success(), "{json}");
    assert_eq!(json["report"]["completed"], true);
    let drift = json["report"]["recorders"][0]["samples"][0][1]
        .as_f64()
        .unwrap();
    assert!(drift > 0.0 && drift < 100.0, "{drift}");
}

/// A schema mismatch reports the JSON path of the bad value and exits 1.
#[test]
fn parse_errors_carry_a_path() {
    let dir = std::env::temp_dir().join("carapace_cli_test");
    std::fs::create_dir_all(&dir).unwrap();
    let text = std::fs::read_to_string(EXAMPLE)
        .unwrap()
        .replace("\"ndm\": 2", "\"ndm\": \"two\"");
    let path = dir.join("bad.json");
    std::fs::write(&path, text).unwrap();
    let (output, json) = run(&["validate", path.to_str().unwrap()]);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(json["valid"], false);
    assert_eq!(json["errors"][0]["stage"], "parse");
    assert_eq!(json["errors"][0]["path"], "header.ndm");
}

/// A missing file is a usage/I/O error (exit 2), and the schema command prints the model's definitions.
#[test]
fn io_errors_and_schema() {
    let (output, _) = run(&["run", "/nonexistent/model.json"]);
    assert_eq!(output.status.code(), Some(2));
    let (output, json) = run(&["schema"]);
    assert!(output.status.success());
    assert!(json["properties"]["header"].is_object());
}
