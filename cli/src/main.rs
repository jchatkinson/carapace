//! `carapace`: run and validate `CarapaceInputV1` models from the command line (JSON in, JSON out).
//!
//! Exit codes: `0` success, `1` the model is invalid or the analysis failed (the JSON report on stdout
//! says why), `2` usage or I/O problem. stdout is always one JSON document; progress lines (`--progress`)
//! go to stderr as JSON lines, so a script can pipe stdout straight into a parser.

use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::AtomicBool;

use carapace_core::analysis::{set_parallelism, Parallelism};
use carapace_input::building::{expand, BuildingSpec, ExpandOptions, UpAxis};
use carapace_input::input_v1::{decode, run_many, CarapaceInputV1, RunOptions, Session};
use clap::{Parser, Subcommand};
use serde::Serialize;

#[derive(Parser)]
#[command(
    version,
    about = "Carapace finite-element engine: run and validate JSON models"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Decode a model and report every problem found, without running it.
    Validate {
        /// Model file, or `-` for stdin.
        model: PathBuf,
    },
    /// Decode and run a model to completion; print the results as JSON.
    Run {
        /// Model file, or `-` for stdin.
        model: PathBuf,
        /// Write the report here instead of stdout.
        #[arg(short, long)]
        out: Option<PathBuf>,
        /// Steps per slice between progress reports and cancellation checks.
        #[arg(long, default_value_t = 64)]
        step_budget: u32,
        /// Print one JSON line per slice to stderr.
        #[arg(long)]
        progress: bool,
        /// Pretty-print the JSON.
        #[arg(long)]
        pretty: bool,
        /// Threads for the solver and element assembly (`0` = all cores). Results do not depend on it.
        #[arg(long, default_value_t = 1)]
        threads: usize,
    },
    /// Run several independent models concurrently; one result per model, in argument order.
    Batch {
        /// Model files.
        #[arg(required = true)]
        models: Vec<PathBuf>,
        /// Models run at once (`0` = all cores).
        #[arg(long, default_value_t = 0)]
        threads: usize,
        /// Pretty-print the JSON.
        #[arg(long)]
        pretty: bool,
    },
    /// Print the JSON Schema of the model format (`--building`: of the building spec `expand` reads).
    Schema {
        #[arg(long)]
        building: bool,
    },
    /// Expand a stud-wall building spec into an ordinary model and an index map; print `{ "input", "map" }`.
    Expand {
        /// `BuildingSpec` file, or `-` for stdin.
        spec: PathBuf,
        /// Which model axis is vertical (plan x-z maps to model x-y for `z`).
        #[arg(long, value_enum, default_value_t = Up::Y)]
        up: Up,
        /// An analysis sequence (`SequenceSpec` JSON) to put in the model. Its recorders address the map's rows.
        #[arg(long)]
        sequence: Option<PathBuf>,
        /// Print only the model (`input`), ready for `carapace run`.
        #[arg(long)]
        model_only: bool,
        /// Write the JSON here instead of stdout.
        #[arg(short, long)]
        out: Option<PathBuf>,
        /// Pretty-print the JSON.
        #[arg(long)]
        pretty: bool,
    },
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum Up {
    Y,
    Z,
}

/// A problem with the input, in the shape agents can act on: where it is and what it is.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct Diagnostic {
    /// `parse` (not valid JSON or does not match the schema) or `decode` (a model error).
    stage: &'static str,
    /// JSON path of the offending value (`parse` only).
    #[serde(skip_serializing_if = "Option::is_none")]
    path: Option<String>,
    /// The message (`parse`) or the tagged `DecodeError` object (`decode`).
    detail: serde_json::Value,
}

fn read_source(path: &PathBuf) -> std::io::Result<String> {
    if path.as_os_str() == "-" {
        let mut text = String::new();
        std::io::stdin().read_to_string(&mut text)?;
        Ok(text)
    } else {
        std::fs::read_to_string(path)
    }
}

fn load(path: &PathBuf) -> Result<Result<Session, Diagnostic>, std::io::Error> {
    let text = read_source(path)?;
    let mut de = serde_json::Deserializer::from_str(&text);
    let input: CarapaceInputV1 = match serde_path_to_error::deserialize(&mut de) {
        Ok(input) => input,
        Err(e) => {
            return Ok(Err(Diagnostic {
                stage: "parse",
                path: Some(e.path().to_string()),
                detail: e.inner().to_string().into(),
            }))
        }
    };
    Ok(decode(input).map_err(|e| Diagnostic {
        stage: "decode",
        path: None,
        detail: serde_json::to_value(&e).unwrap_or_default(),
    }))
}

fn parse_json(text: &str) -> Result<serde_json::Value, Diagnostic> {
    serde_json::from_str(text).map_err(|e| Diagnostic {
        stage: "parse",
        path: None,
        detail: e.to_string().into(),
    })
}

fn parse_diagnostic(e: &serde_path_to_error::Error<serde_json::Error>) -> Diagnostic {
    Diagnostic {
        stage: "parse",
        path: Some(e.path().to_string()),
        detail: e.inner().to_string().into(),
    }
}

fn emit(value: &impl Serialize, pretty: bool, out: Option<&PathBuf>) -> std::io::Result<()> {
    let text = if pretty {
        serde_json::to_string_pretty(value)
    } else {
        serde_json::to_string(value)
    }
    .map_err(std::io::Error::other)?;
    match out {
        Some(path) => std::fs::write(path, text + "\n"),
        None => writeln!(std::io::stdout(), "{text}"),
    }
}

fn execute(cli: Cli) -> std::io::Result<bool> {
    match cli.command {
        Command::Schema { building } => {
            if building {
                emit(&schemars::schema_for!(BuildingSpec), true, None)?;
            } else {
                emit(&schemars::schema_for!(CarapaceInputV1), true, None)?;
            }
            Ok(true)
        }
        Command::Expand {
            spec,
            up,
            sequence,
            model_only,
            out,
            pretty,
        } => {
            let fail = |diagnostic: Diagnostic| -> std::io::Result<bool> {
                emit(
                    &serde_json::json!({ "valid": false, "errors": [diagnostic] }),
                    pretty,
                    None,
                )?;
                Ok(false)
            };
            let parse = |path: &PathBuf| -> std::io::Result<Result<serde_json::Value, Diagnostic>> {
                Ok(parse_json(&read_source(path)?))
            };
            let spec: BuildingSpec = match parse(&spec)? {
                Ok(value) => match serde_path_to_error::deserialize(value) {
                    Ok(spec) => spec,
                    Err(e) => return fail(parse_diagnostic(&e)),
                },
                Err(diagnostic) => return fail(diagnostic),
            };
            let options = ExpandOptions {
                up: match up {
                    Up::Y => UpAxis::Y,
                    Up::Z => UpAxis::Z,
                },
            };
            let mut expansion = match expand(&spec, &options) {
                Ok(expansion) => expansion,
                Err(e) => {
                    return fail(Diagnostic {
                        stage: "expand",
                        path: None,
                        detail: e.to_string().into(),
                    })
                }
            };
            if let Some(path) = sequence {
                match parse(&path)? {
                    Ok(value) => match serde_path_to_error::deserialize(value) {
                        Ok(sequence) => expansion.input.sequence = sequence,
                        Err(e) => return fail(parse_diagnostic(&e)),
                    },
                    Err(diagnostic) => return fail(diagnostic),
                }
            }
            if model_only {
                emit(&expansion.input, pretty, out.as_ref())?;
            } else {
                emit(&expansion, pretty, out.as_ref())?;
            }
            Ok(true)
        }
        Command::Validate { model } => match load(&model)? {
            Ok(session) => {
                let ndm = if matches!(session, Session::D2(_)) {
                    2
                } else {
                    3
                };
                emit(
                    &serde_json::json!({ "valid": true, "ndm": ndm }),
                    false,
                    None,
                )?;
                Ok(true)
            }
            Err(diagnostic) => {
                emit(
                    &serde_json::json!({ "valid": false, "errors": [diagnostic] }),
                    false,
                    None,
                )?;
                Ok(false)
            }
        },
        Command::Batch {
            models,
            threads,
            pretty,
        } => {
            let mut results = vec![serde_json::Value::Null; models.len()];
            let mut inputs = Vec::new();
            let mut slots = Vec::new();
            for (at, path) in models.iter().enumerate() {
                let text = read_source(path)?;
                let mut de = serde_json::Deserializer::from_str(&text);
                match serde_path_to_error::deserialize::<_, CarapaceInputV1>(&mut de) {
                    Ok(input) => {
                        inputs.push(input);
                        slots.push(at);
                    }
                    Err(e) => {
                        results[at] = serde_json::json!({ "valid": false, "errors": [Diagnostic {
                            stage: "parse", path: Some(e.path().to_string()), detail: e.inner().to_string().into(),
                        }] })
                    }
                }
            }
            let options = RunOptions::default();
            for (at, outcome) in slots.into_iter().zip(run_many(inputs, threads, &options)) {
                results[at] = match outcome {
                    Ok(report) => serde_json::json!({ "valid": true, "report": report }),
                    Err(e) => serde_json::json!({ "valid": false, "errors": [Diagnostic {
                        stage: "decode", path: None, detail: serde_json::to_value(&e).unwrap_or_default(),
                    }] }),
                };
            }
            for (result, path) in results.iter_mut().zip(&models) {
                result["path"] = path.display().to_string().into();
            }
            let ok = results
                .iter()
                .all(|r| r["valid"] == true && r["report"]["completed"] == true);
            emit(&serde_json::json!({ "results": results }), pretty, None)?;
            Ok(ok)
        }
        Command::Run {
            model,
            out,
            step_budget,
            progress,
            pretty,
            threads,
        } => {
            if threads != 1 {
                set_parallelism(Parallelism::Threads(threads));
            }
            let mut session = match load(&model)? {
                Ok(session) => session,
                Err(diagnostic) => {
                    emit(
                        &serde_json::json!({ "valid": false, "errors": [diagnostic] }),
                        pretty,
                        None,
                    )?;
                    return Ok(false);
                }
            };
            let cancel = AtomicBool::new(false);
            let options = RunOptions {
                step_budget,
                cancel: Some(&cancel),
            };
            let report = session.run(&options, |outcome| {
                if progress {
                    eprintln!(
                        "{}",
                        serde_json::json!({
                            "stepsTaken": outcome.steps_taken,
                            "loadFactor": outcome.load_factor,
                            "done": outcome.done,
                        })
                    );
                }
            });
            let ok = report.completed;
            emit(
                &serde_json::json!({ "valid": true, "report": report }),
                pretty,
                out.as_ref(),
            )?;
            Ok(ok)
        }
    }
}

fn main() -> ExitCode {
    match execute(Cli::parse()) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("carapace: {e}");
            ExitCode::from(2)
        }
    }
}
