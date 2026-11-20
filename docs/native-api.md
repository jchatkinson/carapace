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
  the tangent assembly (chunks evaluated in parallel, scattered in element order) and faer's sparse
  factorization. On a 100x30 fiber-free elastic frame (9,300 DOFs) 4 threads took the solve from 566 ms to
  330 ms; the factorization alone did not help frames, the assembly did. Small models (< 256 elements) stay serial.
  The setting is process-wide (`carapace_core::analysis::set_parallelism`).
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

```python
import carapace
m = carapace.Model(ndm=2)
a = m.add_node(0, 0, fixed="all"); b = m.add_node(0, 3000)
m.add_elastic_beam(a, b, e=200000, a=8000, iz=1e8)
m.add_static_stage(); m.add_nodal_load(b, 0, 10000)
r = m.record_node_disp(b, 0)
print(carapace.validate(m))                 # {"valid": True, "ndm": 2} or structured errors with a path
print(carapace.run(m).recorder(r))          # (n, 2) numpy array of (load factor, displacement)
```

`Model` covers nodes, trusses, elastic beams (2D/3D), nodal loads, static and modal stages and node recorders;
anything else goes in `model.raw` (a dict merged over the output) or pass a plain `CarapaceInputV1` dict or JSON
text to `validate`/`run`. Invalid models raise `ModelError` (with `.errors`); an analysis failure does not raise:
check `results.completed` / `results.error`. The GIL is released while the engine runs.

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
