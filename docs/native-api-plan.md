# Native API and workspace reorganization plan

Status: Phase A (steps 1-2) done: `carapace-input` extracted, `tsify` optional, generated `.d.ts` byte-identical. Phase B done except accumulating validation errors and result tables (see `docs/native-api.md`): `Session::run`, `carapace` CLI, JSON Schema feature, example models. Phase C done: `py/` Python package (builder, numpy results, tests). The rest is plan only.

Goals:

1. AI agents can drive native carapace from scripts to automate work, for example building models from
   structural drawings, running them, and reading results back.
2. Native-only features (multithreading first) can be added without touching or bloating the wasm build.
3. The wasm build, the wire format and PySees keep working unchanged.

## 1. Where things stand

| Fact | Consequence |
|---|---|
| `carapace-core` is plain Rust with no wasm deps. Its public API is a typed, const-generic builder (`Domain<NDIM,NDOF,NId,E>`, `Domain3`, `AnalysisBuilder` typestate, closed element enums). | Excellent for the engine; poor for scripts. An agent must write and compile Rust against 2D/3D aliases and ids from slotmaps. Compile loops are slow and the error surface is large. |
| The "model from data" layer (`input_v1`: tables, decode, `Session`, stages, recorders, `ModelOutput`; about 6.5k lines) plus `section_probe` and `material_probe` (about 2.5k lines) lives in `wasm-bridge`. | The only complete, validated, string-and-index-based way to build and run a model is in the wasm crate. Native tests already use it (`wasm-bridge/tests/wire`), so it works natively, but it is named, packaged and dependency-weighted as wasm. |
| `input_v1` types derive `serde` and `tsify::Tsify` directly. | Wire types are tied to wasm-bindgen's codegen crate. |
| `Session::advance(step_budget)` is a cooperative stepping API built for a Web Worker. | Fine for native too (progress, cancellation), but there is no run-to-completion convenience. |
| The wire format has a strong error model (`DecodeError`, `ModelErrorDetail`) that PySees turns into diagnostics. | This is exactly what an agent needs: structured, machine-readable errors it can fix and retry. |
| No threading anywhere in `core`. faer is built `default-features = false` with `std`, `sparse-linalg`. Assembly is sequential. | Threading is greenfield, so it can be designed behind a feature flag. |

Conclusion: the wire format is already the native scripting API in all but name. It is declarative,
validated, serializable, versioned, and mirrored by PySees. The plan promotes it rather than inventing a
second API.

## 2. Direction for the native API

Three layers, one source of truth.

```
Layer 3  agent-facing frontends     carapace CLI (JSON in/out)      carapace Python package (PyO3)
Layer 2  model/session API          carapace-input: wire types, decode, Session, probes, run_to_completion
Layer 1  engine                     carapace-core: Domain, elements, materials, analysis
```

### 2.1 Make the wire format the scripting API (recommended)

Agents are strongest at emitting data and iterating on error messages. Generating JSON or Python dicts from a
drawing-extraction step, validating it, running it, and reading structured results is far more reliable than
generating Rust. So:

- **CLI**: `carapace run model.json [--out results.json|results.bin] [--threads N]`, plus
  `carapace validate model.json` (decode only, emit all diagnostics as JSON, nonzero exit), and
  `carapace schema` (emit JSON Schema for `CarapaceInputV1`, via `schemars` behind a feature). Streaming
  progress as JSON lines on stderr. Fast to build, language-agnostic, trivially sandboxed.
- **Python package** (`carapace-py`, PyO3 + maturin): `carapace.run(model_dict)`, `carapace.validate(...)`,
  `Session` with `advance()` and `results()` returning numpy arrays for bulk recorders, section and
  material probes. Python is the natural home for drawing-to-model automation and matches the existing
  `comparison/` harness (OpenSeesPy). It is the layer worth investing in for agents.
- **Thin builder helpers on top of the wire types** (in the Python package, optionally mirrored in Rust):
  `Model().node(...)`, `.beam(...)`, `.fix(...)` that append to the columnar tables and return stable
  string or integer handles. The columnar format is awkward to hand-write (parallel arrays, row indices),
  and the builder removes that failure mode. Generate it from the wire types where possible so it cannot
  drift.

Not recommended: a hand-written ergonomic Rust scripting API (a second model description to maintain next
to the wire format) or a Rhai/Lua embedded scripting layer. Rust scripts can still use `carapace-core`
directly for engine work; that API stays typed and unchanged.

### 2.2 Agent-friendliness requirements

These are what make the native API good for agents, regardless of frontend:

1. **Machine-readable errors, all at once.** `validate` should collect every decode/model error rather than
   stopping at the first, each with a stable code, a path into the document (`trusses.nodeJ[12]`), and a
   suggested fix when cheap. Today's `DecodeError` is a good base. Check which paths return on first error.
2. **Self-describing**: JSON Schema output plus an `examples/` set of small valid models (cantilever, portal
   frame, shear wall, fiber section) that double as docs and as agent few-shot context.
3. **Run to completion**: `Session::run(budget, on_progress)` (and cancellation token) so native callers do
   not hand-roll the `advance` loop.
4. **Deterministic**: identical results independent of thread count (see 3.3), so agents can diff runs.
5. **Units and conventions stated in the schema** (the schema carries no units today; document the
   convention and add an optional `header.units` label that is echoed in results, not converted).
6. **Results as tables**: named columns with node/element labels, not only positional arrays. Keep the
   typed-array bulk path for speed.
7. **Model-from-drawing support** lives outside the engine: extraction (PDF/DXF to geometry) is a Python
   concern. The engine's contribution is fast validate-and-run plus good diagnostics, and `model_output`
   (geometry for rendering) so the agent can plot the generated model and check it against the drawing.

## 3. Workspace reorganization

### 3.1 Target layout

```
carapace/
  core/            carapace-core      engine; no wasm, no serde in non-dev deps
  input/           carapace-input     NEW: wire types, decode, Session/stages/recorders, probes
  wasm-bridge/     carapace-wasm      only wasm-bindgen glue (boundary.rs, Ts<T>, JsValue conversions)
  cli/             carapace-cli       NEW: `carapace` binary
  py/              carapace-py        NEW: PyO3 extension module (excluded from default workspace build)
```

`carapace-input` is the extraction of everything in `wasm-bridge/src/{input_v1,section_probe,material_probe.rs}`
that does not touch `wasm_bindgen`/`JsValue`. Grep shows only `boundary.rs` is wasm-specific; the rest uses
serde and the `tsify` derive.

### 3.2 Concrete steps

Each step leaves `cargo test` and `comparison/check_refactor.sh` green and the wasm output byte-compatible.

1. **Extract `carapace-input`.** Move `input_v1/`, `section_probe/`, `material_probe.rs` and the
   `tests/wire` binary with them (they are already native tests). `carapace-wasm` depends on it and keeps
   only `boundary.rs`. No behavior change. Update the docs paths in `AGENTS.md`, `input-format.md`, and the
   test-report generator (`docs-site/gen_tests.py` reads `cargo test` output, so check its crate filters).
2. **Make `tsify` optional.** In `carapace-input`, `#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]`
   (and `#[cfg_attr(feature = "tsify", tsify(optional))]` on fields). `carapace-wasm` enables the feature;
   native crates do not pull in `wasm-bindgen`/`js-sys`. Verify the generated `.d.ts` is identical.
   Add a `schemars` feature the same way for JSON Schema.
3. **Session ergonomics.** Add `Session::run_to_completion`, a progress callback, a cancellation flag
   (`AtomicBool`), and an error-accumulating `validate`. These are native-useful and wasm-neutral
   (wasm keeps `advance`). Bump `docs/input-format.md` for any wire changes and mirror in PySees per
   `AGENTS.md`.
4. **CLI crate.** `clap` + `serde_json`; subcommands `run`, `validate`, `schema`, `probe-section`,
   `probe-material`. Output JSON, plus a compact binary/NDJSON option for bulk. Integration tests run the
   fixtures in `wasm-bridge/tests/fixtures` through the binary.
5. **Python package.** `py/` with PyO3 (abi3) and maturin; dict/JSON in (reuse serde via `pythonize`, or take
   a JSON string for the first cut), numpy out for bulk. Ship the builder helpers and examples.
   Add to the `comparison/` harness as a consumer so Python and wasm are cross-checked against OpenSeesPy by
   the same fixtures. Exclude from `default-members` so plain `cargo test` does not need a Python toolchain.
6. **Docs.** `docs/native-api.md` (CLI and Python usage, schema, examples) and an `AGENTS.md` section
   pointing agents at `carapace validate` first. Consider shipping a short `CLAUDE.md`/skill-style
   guide inside the Python package for agent users.

Order matters: 1-2 are pure moves with a mechanical check, 3 is small, and 4-5 are additive.

## 4. Native-only features: multithreading

### 4.1 Principle

Wasm and native share the same core code path; threading is a compile-time capability, never a fork of the
algorithm.

- Cargo feature `parallel` on `carapace-core` (default off; off for wasm). Gate `rayon` and faer's
  parallel paths behind it. `cfg(target_arch = "wasm32")` plus a `compile_error!` if both are enabled
  prevents an accidental wasm build with threads (wasm threads need SharedArrayBuffer and special build
  flags, and `docs/wasm-performance.md` records that the current build is deliberately simple).
- A runtime `Parallelism` setting (`Serial` | `Threads(n)`) on `AnalysisBuilder` or analysis options, and
  `--threads` on the CLI. Default serial so results and performance are unchanged unless requested.
- Core stays `Send` where it is cheap, and a test asserts `Domain: Send` so it can move across threads.

### 4.2 Where parallelism pays, in order of value

1. **Across independent analyses** (biggest, easiest, no core change): parameter sweeps, load-case
   combinations, section/material probe batches, and multiple models from a drawing-to-model run. Each
   `Session` is independent, so the CLI and Python layers run many in a rayon pool or via a `run_many`
   API. Deterministic by construction.
2. **Sparse factorization**: faer supports parallel `Par::Rayon(n)` in its sparse LU/LLT. Plumb `Par`
   through `SparseSolver` (`core/src/analysis/solver.rs`). Gains depend on model size; measure on
   `benchmark_scaling`.
3. **Element state update and assembly**: the per-element loop is the hot path for fiber and continuum
   models. Parallelize with a per-thread triplet buffer (chunked by element) merged in element order, so
   the summed matrix is bit-identical to serial. The closed-enum, no-`dyn`, no-heap-in-loop design makes
   this safe. Needs the `TangentSink` to be per-chunk, which the existing sink trait already allows.
4. **Within a section/element**: not worth it; granularity too small.

### 4.3 Determinism

Assemble per-chunk in fixed element order and reduce in a fixed order, never with atomics or
floating-point reductions that depend on scheduling. Add a test that runs the same model at 1, 2 and 8
threads and asserts bit-identical results. This is also what lets agents trust diffs between runs.

### 4.4 Other native-only options (later, not now)

Optional native-only solvers (MKL/Pardiso, GPU), memory-mapped result files for huge recorders, and file
I/O for results (Parquet/HDF5) belong in `carapace-cli`/`carapace-py` features, not core.

## 5. Risks and open questions

- **Python vs Rust for agent scripts.** Decided: Python and the CLI (JSON) are the agent surface. Rust
  scripting is not a goal; if it becomes one, add a small typed model builder in `carapace-input` (a Rust
  twin of the Python builder) rather than exposing `Domain` directly.
- **Crate rename churn.** `AGENTS.md`, docs, `docs-site` generators, CI and PySees' build script all
  name `carapace-wasm` paths. Step 1 needs a pass over those; the comparison harness and
  `scripts/build-wasm.sh` should only need path checks.
- **Wire format growth.** More frontends make the wire format a stronger contract. Keep PySees' type and
  compiler in sync (existing rule), and add a schema snapshot test so unintended wire changes fail CI.
- **faer parallelism and wasm.** Ensure the `parallel` feature never enters the wasm dependency graph
  (`cargo tree -p carapace-wasm --target wasm32-unknown-unknown` check in CI).
- **Python packaging.** abi3 wheels via maturin and CI matrix (linux/mac/win); decide whether to publish or
  keep it as an in-repo install.
- **Results size.** Large recorders as JSON are slow; prioritize the typed-array path (numpy) and a binary
  CLI output early.

## 6. Suggested sequencing

| Phase | Scope | Outcome |
|---|---|---|
| A | Steps 1-2 (extract `carapace-input`, optional tsify) | Clean native/wasm split, zero behavior change |
| B | Step 3 + CLI (step 4) | Agents can validate and run models today |
| C | Python package + examples + docs (steps 5-6) | The main agent surface |
| D | `parallel` feature: sweeps first, then faer `Par`, then chunked assembly | Native speedups with deterministic results |

Phases B and D are independent after A, so threading can start in parallel with the frontends.
