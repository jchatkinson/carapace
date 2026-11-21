# Native API: CLI and run API

How to build, validate and run Carapace models without a browser. The model format is the wire format
(`CarapaceInputV1`, see [input-format.md](input-format.md)); the CLI and the Python package are thin frontends
over `carapace-input`. Plan and rationale: [native-api-plan.md](native-api-plan.md).

## CLI

```bash
cargo run -p carapace-cli -- validate examples/models/portal_frame_2d.json
cargo run -p carapace-cli -- run examples/models/portal_frame_2d.json [--out report.json] [--progress] [--pretty]
cargo run -p carapace-cli -- schema > carapace-input.schema.json
```

- Input is a file or `-` for stdin. stdout is always one JSON document.
- Exit codes: `0` success; `1` invalid model or failed analysis (the JSON says why); `2` usage or I/O error.
- `validate` and `run` report problems as `{"valid": false, "errors": [{stage, path?, detail}]}`. `stage: "parse"`
  means the JSON does not match the schema and `path` is the JSON path of the offending value
  (`header.ndm`); `stage: "decode"` means a model error and `detail` is the tagged `DecodeError`.
  Only the first problem is reported today.
- `run` prints `{"valid": true, "report": RunReport}`: `completed`, `cancelled`, `error`, `stageId`,
  `totalSteps`, `recorders` (per recorder and stage: `(pseudoTime, value)` samples), `bulk`, `limitHistory`,
  `stepDiagnostics`, `modal`. `--progress` writes one JSON line per slice to stderr.
- `schema` prints the JSON Schema (draft 2020-12) generated from the Rust types (`schema` feature of
  `carapace-input`).

Examples live in `examples/models/` and are run by the `bundled_example_models_run` test, so they stay valid.

## Threads

Native only, behind the `parallel` feature (the CLI and Python package enable it; `carapace-wasm` never does and
`carapace-core` refuses `parallel` on `wasm32` at compile time). Results are bit-identical for any thread count.

- `carapace run model.json --threads N` (and `run(model, threads=N)` in Python): threads the element loop of
  the tangent assembly (chunks evaluated in parallel, scattered in element order) and, for large factorizations
  only, faer's sparse Cholesky. Serial is the default (threads = 1). Measured on an 8-thread machine:
  - A 100x30 fiber frame (`disp_beam_column`, 9,300 DOFs, `benchmark_scaling`): 185 ms serial, 135 ms on 4 threads
    (the assembly is the part that scales).
  - A 150x150 elastic `Shell4` plate (137k DOFs, `benchmark_shell`): no gain. About 45% of the 2.3 s is the
    factorization, which is serial below 60M entries in `L` because faer only threads inside a supernode and
    threading a 38M-entry 2D factor was 3x *slower*; the element tangents (0.3 s) cost about as much to record and
    replay as to compute. A 160M-entry 3D pattern factored 2.2x faster on 8 threads, so that path stays threaded.
  Small models (< 256 elements) stay serial. The setting is process-wide (`carapace_core::analysis::set_parallelism`).
- `carapace batch a.json b.json ... [--threads N]` and `carapace.run_many([...], threads=N)`: independent models
  concurrently, one result per model in order. This is the largest win for sweeps and load cases, and it
  switches solver threading off while it runs to avoid oversubscription.

Not threaded yet: the first assembly (triplet path), element state update, mass and load assembly, and
transient/modal solves' own loops.

## Python

`py/` is the `carapace-fe` package (PyO3, abi3, maturin; its own Cargo workspace so `cargo test` needs no Python):

```bash
cd py && python -m venv .venv && .venv/bin/pip install maturin numpy pytest
.venv/bin/maturin develop --release && .venv/bin/python -m pytest tests
```

The commands follow OpenSeesPy (same names, argument order, 1-based DOFs); see [python-api.md](python-api.md) for
the supported commands and every deviation.

```python
import carapace.opensees as ops
ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.node(1, 0.0, 0.0); ops.node(2, 0.0, 3000.0); ops.fix(1, 1, 1, 1)
ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 1, 1, 2, 8000.0, 200000.0, 1e8, 1)
ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1); ops.load(2, 1e4, 0.0, 0.0)
ops.integrator('LoadControl', 1.0); ops.analysis('Static')
ops.analyze(1); print(ops.nodeDisp(2, 1))
```

Unsupported commands raise `NotImplementedError`. Invalid models raise `ModelError` (with `.errors` and the JSON
path); an analysis failure makes `analyze` return `-1` (`model.last_error`) and `carapace.run` report
`results.completed` / `results.error`. The GIL is released while the engine runs.

## Rust: `Session::run`

```rust
let mut session = carapace_input::input_v1::decode(input)?;
let report = session.run(&RunOptions::default(), |slice| eprintln!("{}", slice.load_factor));
```

`RunOptions { step_budget, cancel }` slices the run; `cancel` is an `AtomicBool` checked between slices.
The wasm boundary keeps calling `advance` itself.

## Known sharp edges (for agents writing models)

- "Every table is optional" holds for whole tables, but a table that is present needs all its columns: a
  `nodes` table without `massNodeIndex`/`mass`, or a static stage without `holdPatternsAfter`, fails with
  `missing field`. Defaulting the empty columns is a candidate wire change (it would also change the
  generated `.d.ts`, so PySees must be updated with it).
- Decoding stops at the first error; collecting all of them needs a decoder change.
