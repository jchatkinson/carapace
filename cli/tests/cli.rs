//! End-to-end checks of the `carapace` binary: exit codes and the JSON it prints.

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_carapace");
const EXAMPLE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../examples/models/portal_frame_2d.json"
);
const STUD_WALL: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../examples/models/stud_wall_2d.json"
);
const BUILDING: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../examples/models/two_storey_box_building.json"
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

/// The wire wall and two-storey builder examples solve through the CLI and match hand shear/displacement values.
#[test]
fn runs_macro_element_examples() {
    for (path, expected, tol) in [
        (STUD_WALL, 0.00028333333333333335, 1e-12),
        (BUILDING, 1.25, 1e-8),
    ] {
        let (output, json) = run(&["run", path]);
        assert!(output.status.success(), "{json}");
        assert_eq!(json["report"]["completed"], true);
        let value = json["report"]["recorders"][0]["samples"][0][1]
            .as_f64()
            .unwrap();
        assert!(
            (value - expected).abs() <= tol,
            "{path}: {value} vs {expected}"
        );
        if path == BUILDING {
            let rotation = json["report"]["recorders"][1]["samples"][0][1]
                .as_f64()
                .unwrap();
            assert!(
                rotation.abs() < 1e-12,
                "symmetric rigid building rotation {rotation}"
            );
        }
    }
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
    for table in [
        "studWalls",
        "floorBands",
        "diaphragms",
        "floorLoads",
        "building",
    ] {
        assert!(
            json["properties"][table].is_object(),
            "schema missing {table}"
        );
    }
}

/// `batch` returns one result per model in argument order, keeps going past an invalid model, and exits 1.
#[test]
fn batch_runs_models_concurrently_in_order() {
    let missing_field = std::env::temp_dir().join("carapace_cli_test_batch.json");
    std::fs::write(&missing_field, "{}").unwrap();
    let (output, json) = run(&[
        "batch",
        EXAMPLE,
        missing_field.to_str().unwrap(),
        EXAMPLE,
        "--threads",
        "2",
    ]);
    assert_eq!(output.status.code(), Some(1));
    let results = json["results"].as_array().unwrap();
    assert_eq!(results.len(), 3);
    assert_eq!(results[0]["report"]["completed"], true);
    assert_eq!(results[1]["valid"], false);
    assert_eq!(results[2]["path"], EXAMPLE);
    assert_eq!(
        results[0]["report"]["recorders"],
        results[2]["report"]["recorders"]
    );
}

/// Solver threads do not change the printed results.
#[test]
fn run_threads_flag_is_accepted() {
    let (a, one) = run(&["run", EXAMPLE]);
    let (b, many) = run(&["run", EXAMPLE, "--threads", "4"]);
    assert!(a.status.success() && b.status.success());
    assert_eq!(one, many);
}
