# Embedded results storage via Turso

This document evaluates replacing pysees-handoff.md's "worker owns a
separately-instantiated sqlite-wasm, fed by `postMessage`d recorder
batches" results-storage design with a single embedded engine: compile
[`turso_core`](https://github.com/tursodatabase/turso) — a from-scratch,
pure-Rust, SQLite-compatible database engine — directly into
`carapace-wasm`, so the analysis loop writes results with an ordinary Rust
function call, not a JS round trip.

This is a real alternative, not a refinement of the existing plan: it moves
results-database ownership from the worker to `carapace-wasm`, which is a
responsibility-boundary change from pysees-handoff.md's "analysis worker:
... SQLite/OPFS ownership" line. That tradeoff is made explicit below
rather than silently absorbed.

## Why this over the alternatives already considered

Three options were compared in conversation before this document:

1. **Status quo**: worker holds a separately-instantiated official
   `sqlite-wasm` (Emscripten-built), fed by `recorderBatch`s crossing the
   `wasm_bindgen` boundary via `serde-wasm-bindgen`. Every batch pays a
   JS-object-graph marshaling cost, and two independent wasm modules can
   never call each other directly — JS always brokers.
2. **Bundle real SQLite (`libsqlite3-sys`/`rusqlite`) into `carapace-wasm`**.
   Removes the two-module problem, but SQLite's C code assumes a POSIX-ish
   filesystem; getting it to run against OPFS means hand-writing a
   `sqlite3_vfs` C-ABI struct (`xOpen`/`xRead`/`xWrite`/`xSync`/...) in
   `unsafe` Rust — real, error-prone systems work, reimplementing what the
   official sqlite-wasm project already solved once in C/Emscripten.
3. **Turso**, this document: a pluggable `IO`/`File` **Rust trait**
   already exists in `turso_core` for exactly this purpose (see
   "Verified feasibility spike" below) — no C toolchain, no hand-rolled
   `sqlite3_vfs`, no second wasm module. Turso's own JS bindings already
   ship a browser/OPFS backend behind that trait
   (`bindings/javascript/src/browser.rs` in the Turso repo), proving the
   approach works, even though that specific implementation targets
   napi-rs/emnapi rather than `wasm-bindgen` and isn't directly reusable
   as-is (see "OPFS `IO` backend" below).

Turso is pre-1.0 (`turso_core` v0.8.0-pre.12 as spiked below); its own docs
recommend independent backups until 1.0 for anything that's a system of
record. Carapace's results are regenerable from the model + sequence
snapshot pysees-handoff.md already defines runs around, which meaningfully
lowers that risk relative to a primary database — but it's a real
early-adopter bet on both the engine and its wasm story, and should be
named as one, not glossed over.

## Verified feasibility spike

Before proposing this, `turso_core` was cloned and built directly against
`wasm32-unknown-unknown` (not just read about) to check the two riskiest
claims. Both held:

- **It compiles for `wasm32-unknown-unknown`**, with:
  ```toml
  turso_core = { version = "0.8.0-pre.12", default-features = false, features = [
    "fs", "uuid", "time", "json", "percentile", "autovacuum",
  ] }
  ```
  `default-features` on its own does **not** build — some VDBE opcode
  paths (JSON temp-file spill, `DatabaseFile::new`) unconditionally call
  APIs gated `#[cfg(feature = "fs")]`, so `fs` must stay on even though it
  also gates the *native-platform* `UnixIO`/`WindowsIO` backends we don't
  want. That's fine: `core/io/mod.rs`'s own `cfg_block!` already selects
  `GenericIO` as `PlatformIO` for any target that's neither
  `target_family = "unix"` nor `target_os = "windows"` — wasm32 falls
  through to that arm automatically, so enabling `fs` does not pull in
  POSIX/Windows file code on this target. Dropped: `io_uring` (Linux-only,
  `dep:io-uring`/`rustix::io_uring`, would not build for wasm32 regardless),
  `series`/`encryption` (pulled in unrelated errors — `dbsp.rs` incremental
  views, `aegis` crypto — not needed for an append-mostly results table),
  `simd` (portable-SIMD, unnecessary risk for no benefit here).
- **No dependency conflict** with the existing `wasm-bindgen = "0.2"` /
  `serde-wasm-bindgen = "0.6"` stack: `turso_core`'s own transitive deps
  already pull `wasm-bindgen 0.2.108`, `js-sys`, and `web-sys` (for its own
  `getrandom` wasm backend), which is a semver-compatible version Cargo
  will unify with `carapace-wasm`'s existing pin.

Not yet checked (explicitly out of scope for this spike, listed as real
work below): actual behavior against a real OPFS-backed `IO` impl (only
the plain build was verified — no database was opened or written to),
final linked-and-optimized binary size once merged with `carapace-wasm`
(the unstripped, non-dead-code-eliminated `.rlib` alone is ~95 MB, which
says nothing useful about the final `wasm-opt`'d artifact — that number is
reported here only so it isn't quietly skipped, not as an estimate),
and whether the `fs`-feature's `turso_ext/vfs` pulls in anything else
wasm-unfriendly under real use rather than a bare build.

## Architecture

### Where it lives

A new module inside `carapace-wasm` — e.g. `wasm-bridge/src/results_db/`
— parallel to `wasm-bridge/src/input_v1/`. **Not** `carapace-core`:
`carapace-core`'s explicit responsibility ("numerical model and analysis,
no JS, worker, or database concepts", pysees-handoff.md) stays intact.
`carapace-core` still has zero knowledge that a database exists; it just
returns per-step state the way it already does. Everything Turso-specific
— the `IO`/`File` impl, schema, insert/query calls — lives at the wasm
boundary, next to `input_v1`, not inside the solver.

### OPFS `IO` backend

`turso_core::io::{IO, File}` (`core/io/mod.rs` in the Turso repo) are
**completion-based**, not `async`/`await`: `File::pwrite(&self, pos, buf,
c: Completion) -> Result<Completion>` issues a write and returns
immediately; the `Completion` is marked finished later, out of band.
`IO::step()` is the pump that drives outstanding completions forward, and
`drain_completions`/`wait_for_completion` busy-loop `step()` until every
named `Completion` is finished.

Turso's own browser backend (`bindings/javascript/src/browser.rs`)
resolves a `Completion` via an explicit `complete_opfs(completion_no,
result)` function **called from JS**, once the real (async) OPFS operation
settles — the Rust side registers the completion and returns; JS calls
back in later. That's the shape to port to `wasm-bindgen`/`web_sys` here
(their version is napi-rs/emnapi-specific and not directly reusable, but
the pattern — register-then-callback, not block-and-poll — is exactly
right and should be kept, not reinvented).

**This is the load-bearing detail for "non-blocking".** If Rust code ever
calls `wait_for_completion`/`drain_completions` synchronously, right after
issuing a write, inside the same JS→wasm call that issued it, it will spin
`step()` forever: `complete_opfs` can only fire once JS's event loop gets a
turn, and it can't get a turn while wasm is still busy-looping inside a
single synchronous exported call. The only correct pattern is **fire and
forget**: issue the write, do not wait on its `Completion` before
returning control to JS, and check whether earlier writes finished (and
whether any errored) on a later call — by which point JS has had a chance
to call `complete_opfs` on its own schedule.

### Ownership and schema

`ResultsDb` (new type, `results_db/mod.rs`) owns one `turso_core::Database`
+ `Connection`, opened once via `Database::open(io: Arc<dyn IO>, path:
&str, options)` — this entry point is **not** gated by the `fs` feature,
unlike the convenience wrappers (`open_file`, `open_file_with_flags`)
turso's own native CLI uses; it's the right one for an embedder supplying
its own `IO`. `path` here is whatever our `IO::open_file` interprets it
as — it never touches a real filesystem, only our OPFS `File` impl.

Schema: reuse pysees-handoff.md's `runs`/`stages`/`recorders`/
`response_blocks` shape (it was already designed to be
producer-agnostic), or — given a single `WasmSession` now corresponds to
exactly one run, decided at construction, not looked up by `run_id` — drop
the `runs` table's role down to a single header row and let the caller
(worker) key its own metadata (`run_id`, `model_hash`, `sequence_hash`)
outside this database, in whatever store already tracks the run list
across runs. That's a real design fork worth resolving with the pysees
side before implementation, not something to decide unilaterally here.

### Batching: what, where, and why

**Batch boundary = one `advance(step_budget)` call.** This is not a new
concept to introduce — `PlanarSession::advance`
(`wasm-bridge/src/input_v1/session.rs:189`) already loops `analysis.step()`
up to `step_budget` times per call, and already calls
`self.record_sample()` (`session.rs:306`) after every successful step,
appending into an in-memory `Vec<Vec<(f64, f64)>>`
(`recorder_history`, `session.rs:139`). The database write replaces (or
supplements — see "Existing code that needs to change" below) that
`Vec::push`'s destination: instead of accumulating forever in memory and
only materializing when JS calls `recorderSamples`, flush the batch
accumulated during *this* `advance()` call to `ResultsDb` once, right
before returning `StepOutcome`.

This piggybacks on cooperative cancellation for free, not as a
coincidence: `advance()` already has to return control to JS after
`step_budget` steps so a worker's cancel-check loop can run between calls
(pysees-handoff.md's "each call only ever executes up to `stepBudget`
steps before returning"). That's the exact same point at which queued OPFS
writes get their chance to complete, since JS callbacks (`complete_opfs`)
can only run between wasm calls, never during one. One design serves both
purposes — no separate flush-timer or idle-callback machinery is needed.

**Do not flush per step.** A `step_budget` of, say, 50–200 gives 50–200
solver steps per SQL transaction/OPFS write instead of one write per step
— the same amortization argument from the earlier discussion of
`recorderBatch` sizing, just now happening inside `turso_core`'s own
transaction instead of the worker's SQL insert loop. `step_budget` remains
the single tuning knob for write-amortization vs. cancellation/progress
granularity; this design doesn't add a second one.

**One transaction per flush, not one `INSERT` per sample.** `ResultsDb`
should batch every recorder's accumulated samples for the just-completed
range of steps into a single `BEGIN`/multiple bound `INSERT`s (or a bulk
`INSERT ... VALUES (...), (...), ...`)/`COMMIT`, issued as one call into
`turso_core`. `turso_core`'s own WAL/page-cache buffering then decides when
the actual `pwrite`/`fsync`-equivalent `Completion`s get issued to the
`IO` backend — that's Turso's problem to schedule well, not something this
layer needs to second-guess.

**Non-blocking write dispatch, restated concretely**: `ResultsDb::flush`
executes the batch's SQL synchronously against `turso_core` (CPU-bound
query planning/execution, no I/O wait — this part *is* fine to block on,
it's plain computation) but must not call anything that waits on a
`Completion` tied to the OPFS backend's actual disk write. Any error from
an OPFS write that hasn't resolved yet by the time of a later `flush`/
`advance` call is surfaded then, not immediately — `StepOutcome` gains an
optional field (or `AnalysisErrorDetail` gains a variant, see below) for
"a previously-issued write failed", checked and cleared at the top of the
next `advance()`.

## Implementation plan

1. **OPFS `IO`/`File` backend** (`wasm-bridge/src/results_db/opfs_io.rs`).
   Port `browser.rs`'s register/callback pattern from napi-rs conventions
   to `wasm-bindgen`/`web_sys`: a `wasm_bindgen(js_name = completeOpfsWrite)`
   exported function JS calls back into, a `HashMap<u32, Completion>` of
   outstanding ops, `File::pwrite`/`sync`/`truncate` each registering a
   completion and calling an imported JS function that performs the actual
   `FileSystemSyncAccessHandle` (or async `FileSystemWritableFileStream`,
   depending on what the worker's OPFS setup exposes) operation. This is
   the largest, riskiest, and most novel piece of work in this plan — treat
   it as its own spike with its own pass/fail checkpoint (open a handle,
   write bytes, read them back, from a real Worker) before building
   anything on top of it.
2. **`ResultsDb`** (`wasm-bridge/src/results_db/mod.rs`): schema creation
   (`CREATE TABLE IF NOT EXISTS` on first open), a `flush(&mut self,
   batch: &[RecorderSample])` that runs one transaction, and a `query`
   surface for reading back a recorder's history (range by run/stage, or
   all — matches whatever the worker's `query-results` protocol message
   needs). Depends on (1) only for construction (`Database::open` needs a
   live `Arc<dyn IO>`); the SQL layer itself has no wasm-specific code.
3. **Wire into `PlanarSession`** (`wasm-bridge/src/input_v1/session.rs`):
   `PlanarSession` gains a `results: ResultsDb` field (or an `Option`, if
   a DB-less mode stays supported for tests — see below);
   `record_sample` writes into a small per-call buffer instead of (or in
   addition to) `recorder_history`; `advance` calls
   `self.results.flush(&buffer)` once, after the step loop, before
   building `StepOutcome`. `recorder_samples` (`session.rs:174`), used
   today by `boundary.rs`'s `recorderSamples` export, either stays as an
   in-memory read of the most recent unflushed buffer (for the "what just
   happened" progress display) or becomes a `ResultsDb` query — decide
   based on whether pysees wants "poll the live buffer" or "query the
   database" for its progress UI; they don't have to be the same path.
4. **Async construction handshake**: opening the actual OPFS root
   (`navigator.storage.getDirectory()`) and any sync access handle is
   async JS, but `decodeInput` (`boundary.rs:27`) is a synchronous
   `wasm_bindgen` export today. `ResultsDb`'s `IO` needs a handle the
   *worker* obtained beforehand (async, once, at worker startup) and hands
   into `decodeInput` as an extra `JsValue` argument — mirroring the
   existing "not yet the transferable typed array wire format, but enough
   to drive this from real JS" posture `boundary.rs`'s own doc comment
   already states for the rest of the input.
5. **Error surfacing**: extend `AnalysisErrorDetail`
   (`session.rs:23`) with a `ResultsWriteFailed { detail: String }`
   variant (or a sibling field on `StepOutcome`, since a storage failure
   isn't a solver failure — conflating them would make `advance()`'s
   `done`/`error` semantics say "the analysis failed" for what's actually
   "the analysis succeeded but a write didn't land," which the worker
   needs to distinguish to decide whether to retry the write vs. abort the
   run).
6. **`boundary.rs`** gains whatever `WasmSession` methods the worker needs
   beyond what advancing already triggers — most likely none for writes
   (they're implicit in `advance()`), but a `queryResults`-shaped method
   if reads go through `ResultsDb` instead of the in-memory buffer (see
   step 3).
7. **Real JS/Node round-trip test**: this whole design is unverified
   until it's driven from actual generated `wasm-bindgen` glue in a real
   (or at least `wasmtime`/Node-simulated) OPFS environment — the existing
   `wasm-bridge/tests/m10_carapace_input_v1.rs` pattern (native Rust calls
   into `carapace_wasm`, no real JS) cannot exercise the OPFS `IO` backend
   at all, since that backend's entire reason for existing is the
   JS-callback completion path. A dedicated test harness (real browser via
   `wasm-bindgen-test` with a headless browser runner, since OPFS sync
   access handles are worker-only and not available under plain Node)
   needs to exist before this can be called done, not just built.

## Existing code that needs to change

- `wasm-bridge/Cargo.toml`: add `turso_core` with the feature set
  verified above.
- `wasm-bridge/src/input_v1/session.rs`: `PlanarSession` gains a
  `ResultsDb` (or `Option<ResultsDb>`) field; `advance` (`:189`) flushes
  once per call instead of `record_sample` (`:306`) pushing forever into
  `recorder_history`; `PlanarSession::new` (`:156`) takes the opened
  `ResultsDb` (or the `IO` handle to build one) as a new constructor
  argument; `recorder_samples` (`:174`) either stays as-is (reading the
  live buffer) or becomes a `ResultsDb` query, per the decision in plan
  step 3.
- `wasm-bridge/src/input_v1/mod.rs` / `decode.rs`: `decode` needs to
  accept (or construct) the `ResultsDb`/`IO` alongside the existing
  `CarapaceInputV1` — a new parameter, not a field folded into the wire
  format itself (opening a database handle isn't part of "decode the
  model", it's a side effect the caller should control explicitly).
- `wasm-bridge/src/boundary.rs`: `decode_input` (`:27`->`:28`) gains a
  parameter for the pre-opened OPFS handle (plan step 4); `WasmSession`
  gains no new *required* method for writing (it's implicit in
  `advance`), but may gain a query method depending on the step-3
  decision.
- `StepOutcome` (`session.rs:88`): new optional field or
  `AnalysisErrorDetail` variant for a write failure (plan step 5) — this
  is the one place the `advance()`/`StepOutcome` contract actually changes
  shape, worth flagging to the pysees side early since it's wire-visible.
- New: `wasm-bridge/src/results_db/{mod.rs, opfs_io.rs}` (plan steps 1–2).
- `docs/pysees-handoff.md`: the "Results database" and "Worker protocol and
  ownership" sections describe the worker as owning SQLite/OPFS; if this
  plan is adopted, those sections need rewriting to reflect
  `carapace-wasm` owning the database and the worker instead owning "hand
  `carapace-wasm` an OPFS handle at startup, then read query results back
  through it" — a real contract change, not a wording tweak.

## Open risks, stated plainly

- The OPFS `IO` backend (plan step 1) is unproven — the spike above only
  confirms `turso_core` *compiles* for wasm32, not that a real OPFS-backed
  `File` impl works correctly under load (concurrent writes, WAL
  checkpointing behavior against a sync access handle, recovery after a
  worker crash mid-write). This is the actual risk in this whole plan; the
  rest is comparatively mechanical wiring.
- Pre-1.0 engine: schema-breaking upgrades or correctness bugs in
  `turso_core` itself are a real possibility this plan inherits.
- Final wasm binary size is unmeasured. `turso_core` is a full SQL engine
  (query planner, VDBE, JSON, incremental views) for what this project
  needs as an append-mostly time-series table — that mismatch may or may
  not matter once `wasm-opt`'d and gzipped, but it hasn't been measured
  and shouldn't be assumed away.
- This plan changes a documented responsibility boundary
  (pysees-handoff.md's worker-owns-the-database line). That needs explicit
  agreement from the pysees side, not just a Carapace-side implementation
  decision, since it changes what the worker's `query-results` handler
  actually talks to.
