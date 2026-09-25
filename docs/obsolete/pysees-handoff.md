# pysees → Carapace Handoff

This document defines the boundary between `pysees`, the browser authoring
application, and Carapace, the WebAssembly analysis engine. It is the
implementation contract for M10/M11. Read it before adding a worker protocol,
a model serializer, or a new analysis UI workflow.

## The boundary

`pysees` owns authored intent. Carapace owns one immutable execution snapshot.
The analysis worker owns the results database.

```text
pysees Model + AnalysisSequence
        │ compile + validate once on Run
        ▼
CarapaceInputV1 (metadata + transferable typed arrays)
        │ one postMessage transfer
        ▼
analysis worker ──► carapace-wasm ──► carapace-core
        │
        └──► SQLite-over-OPFS run database
```

This is deliberately not a live, per-entity RPC API and not an interpreter
for arbitrary OpenSees commands.

Since this document was first drafted, `core` moved well past the slice this
handoff originally scoped against: spatial statics, dynamics, ground motion,
and true rigid-diaphragm constraints are all done (M15-M20; see
[`spatial-architecture.md`](spatial-architecture.md)), not just the planar
static path. `CarapaceInputV1` and the worker protocol below are designed to
be profile- and stage-generic from the start — the `space` and stage-kind
discriminants exist in the wire format even though the first decoder only
implements the planar-static arm — so extending to spatial/modal/transient
handoff later means adding table/stage variants, not a breaking format
rewrite. See "Staged decode support, not a staged wire format" below.

Runs are **snapshots, not reactive bindings**. Pressing Run captures the
current model and sequence, assigns a `run_id`, and executes that snapshot.
Edits made while it runs do not modify or cancel it. A later explicit run
creates a new snapshot and a new run. Each persisted run records its model
hash, sequence hash, engine version, and timestamps, so results always show
which input produced them.

## pysees data model required for execution

Keep editor/UI state (selection, dialogs, viewport options, drafts and undo
bookkeeping) separate from the solver-facing data. The solver-facing project
has exactly two authored inputs:

```ts
interface ProjectAnalysisInput {
  model: Model
  sequence: AnalysisSequence
}
```

`Model` should use typed discriminated unions for the Carapace-supported
subset, stable numeric tags, and explicit tag references. It includes a
required `space` discriminant — exactly planar `(NDM=2, NDF=3)` or spatial
`(NDM=3, NDF=6)` — plus nodes, fixities, masses, uniaxial materials, fiber
sections, integrations, elements, constraints, time series, load patterns,
and loads. Do not make the execution model a collection of generic
`Record<string, unknown>` OpenSees call arguments. Export-only OpenSees
features may have a separate representation. The spatial architecture and its
staged availability are defined in [`spatial-architecture.md`](spatial-architecture.md).

`AnalysisSequence` is declarative and ordered. It describes phases, not a
stream of procedural `constraints`/`system`/`analyze` commands:

```ts
type AnalysisSequence = {
  version: 1
  stages: AnalysisStage[]
  recorders: RecorderSpec[]
}

type AnalysisStage =
  | {
      kind: 'static'
      id: string
      steps: number
      integrator:
        | { kind: 'load-control'; increment: number }
        | { kind: 'displacement-control'; nodeTag: number; dof: 0 | 1 | 2; increment: number }
      algorithm: 'linear' | 'newton-raphson'
      convergence?: ConvergenceSpec
      holdPatternsAfter?: number[]
    }
  | { kind: 'modal'; id: string; modes: number }
  | { kind: 'transient'; id: string; config: TransientSpec } // specified when transient UI lands
```

The first required workflow is planar and uses two static stages: gravity with load control,
freezing its selected pattern(s), followed by lateral displacement-control
pushover. This maps directly onto Carapace's existing `Domain` handoff and
`hold_pattern_constant` APIs.

## Compiler and transport

On Run, `pysees` compiles `Model + AnalysisSequence` into
`CarapaceInputV1`. The compiler is a full validation boundary; it must report
tagged diagnostics before anything enters wasm. At minimum it validates:

- one supported execution profile: planar `(2, 3)` today; spatial `(3, 6)`
  is represented by the contract but rejected until M15;
- duplicate/missing tags and every reference between entities;
- supported materials, elements, transformations, loads, constraints and
  integration rules;
- one-based OpenSees-style DOFs converted to Carapace's zero-based DOFs;
- static-stage compatibility, including a valid reference-load sensitivity
  for displacement control;
- fiber-section conventions.

`pysees` fiber children use `yloc`/`zloc`; Carapace's planar fiber section
uses one signed coordinate. For the initial x-y frame convention, compile
using `yloc` and preserve every fiber's area. Rectangular and circular patches
must be expanded to area-weighted literal fibers before handoff. Unsupported
section child kinds must be rejected explicitly, never silently dropped.

The transport contains small structured-clone metadata and newly allocated,
transferable numeric buffers: node tags/coordinates/fixities/masses, flattened
fiber tables with offsets, typed element/load/pattern tables, and the compiled
sequence. Never transfer a buffer owned by the editor store: transferring an
`ArrayBuffer` detaches it from its sender. The exact table layout is defined
below.

## CarapaceInputV1 wire format

Structured as a small header plus one flat table per (profile, element kind,
entity kind), mirroring `core`'s closed-enum catalog instead of a generic
tagged-record array. This keeps `carapace-wasm`'s decoder a fixed set of
per-table loops, each calling the same constructors native code calls
(`Element::Truss(Truss::new(..))`), rather than growing a runtime
broker/dispatch layer that `core` deliberately doesn't have
(implementation-plan.md §2.1) — the wire format should be exactly as closed
and enumerable as the Rust catalog it feeds.

**Header** (structured-clone, not typed arrays — these are all O(1) counts):

- `schemaVersion` — the wire format's own version, independent of
  `engineVersion` below. Bump it only when a table's shape changes, not when
  a new table is added (a new table is additive and absent from older
  producers' headers).
- `space: 2 | 3` — the execution profile discriminant from
  spatial-architecture.md's "Decision: two execution profiles, selected
  once." Selects planar vs. spatial *once*, at decode time, not per entity.
- `engineVersion` — `carapace-core`'s version, recorded for run provenance
  alongside the model/sequence hashes (see "Results database").
- `tables` — a directory of which per-kind tables are present. A table for a
  variant the current decoder doesn't yet support is simply absent, not
  present-but-ignored; see "Staged decode support" below.

**Bulk tables** (transferable typed arrays; never a buffer the editor store
still owns — transferring an `ArrayBuffer` detaches it from its sender):

- **Node table**: `coords: Float64Array` (stride `NDIM`), `fixed: Uint8Array`
  (bitmask, `NDOF` bits per node), optional `mass: Float64Array` (stride
  `NDOF`, addressed via a parallel node-index array rather than a dense
  zero-filled table, since most nodes carry no mass).
- **Per-(profile, element-kind) tables**: one table per concrete Rust type
  (`Truss`/`Truss3`, `ElasticBeamColumn`/`ElasticBeamColumn3`,
  `DispBeamColumn`/`DispBeamColumn3`, `ForceBeamColumn`/`ForceBeamColumn3`,
  `ZeroLength`/`ZeroLength3`), each a flat struct-of-arrays: node-id pairs,
  section/material-arena indices, and `geomTransf`/`vec_xz` params. There is
  no generic "element record" — the table shapes are the enum variants.
- **Fiber tables**: flattened per fiber section, offset-indexed, per this
  document's existing fiber-flattening rule (rectangular/circular patches
  expanded to area-weighted literal fibers before handoff; planar sections
  carry one signed `yloc`, spatial sections carry `y`/`z`).
- **Materials/sections arena**: a small indexed arena, not flattened
  SoA — materials recurse (`Parallel`/`Series`/`MinMax`) and some leaves are
  boxed in `core` (see `core/src/model/materials/mod.rs`'s porting recipe),
  so element and fiber tables reference a material by arena index, and a
  composite entry references other arena indices in turn. This mirrors
  `core`'s own arena-of-indices philosophy (implementation-plan.md §2.2)
  instead of duplicating material definitions per fiber.
- **Load patterns / time series**: small counts as structured-clone tagged
  unions (`Constant`/`Linear`/`Path`), except a `Path` series's sample array,
  which is numeric-heavy and transferable.
- **Compiled `AnalysisSequence`**: structured-clone; stage count is always
  small enough that this never needs a typed-array table.

### Staged decode support, not a staged wire format

`space` and per-stage-kind fields exist in the header from the start
(`static`/`modal`/`transient`; `2`/`3`), even though M10's first decoder only
implements the `static` + planar arms. Rejecting spatial input or a transient
stage returns a structured `unsupported_space`/`unsupported_stage` diagnostic
from *decode*, not a schema-version bump later. This is the one deliberate
departure from the original M10 scoping in implementation-plan.md, which
scoped `CarapaceInputV1` itself to the static-planar slice: given `core`
already has full spatial statics/dynamics/diaphragms (M15-M20) and
transient/modal for both profiles, a wire format narrower than that would
guarantee a breaking rework the moment spatial or transient handoff lands —
work spatial-architecture.md's own "pysees and wasm handoff" section already
flags as outstanding.

**Implementation status:** the Rust-side half of this document — the
`CarapaceInputV1` header/table types, the hand-written decoder, and the
stepped `Session`/`advance` API described below — is implemented in
`wasm-bridge/src/input_v1/` and exercised end to end (including the
gravity → frozen-gravity → pushover acceptance case, cross-checked against
`ForceBeamColumn`'s native tests) in
`wasm-bridge/tests/m10_carapace_input_v1.rs`. These are still plain Rust
value types, not yet the transferable typed arrays a real `postMessage`
would carry — the `wasm_bindgen`/`serde-wasm-bindgen` boundary is deferred
until `pysees`'s compiler exists to actually produce `CarapaceInputV1`
bytes. `pysees`'s compiler, the worker, and the SQLite-over-OPFS results
database (everything from "Compiler and transport" through "Results
database" below, on the `pysees`/worker side) have not been started.

## Decode and session model in carapace-wasm

Decoding is hand-written per table, not generic deserialization into `core`
types: small structured-clone fields decode via `serde-wasm-bindgen` into
DTOs that live in `carapace-wasm` only, then an explicit translation step
calls the same `Domain`/`Element`/`Material` constructors native code calls.
`core` gains no serde derives, no reflection, and no knowledge that wasm
exists — this is the existing responsibility split below, just made
concrete: "`carapace-core`: numerical model and analysis, no JS, worker, or
database concepts."

Decode returns `Result<Session, DecodeError>` and must not panic on
malformed input, even though `pysees`'s compiler is a full validation
boundary and is expected to catch entity-level problems first. Decode is
defense in depth against what the compiler can't see — engine/schema version
skew, a stale compiled snapshot run against a newer or older wasm build — and
a wasm panic tears down the whole worker with no unwind, which the
cooperative-cancellation design below can't recover from. `DecodeError`
follows `AnalysisError`'s style (implementation-plan.md §2.8): tagged
variants carrying context (`UnknownMaterialIndex { table, row }`,
`UnsupportedSpace { got }`, ...), not sentinel codes.

`Session` is chosen once, at decode time, from the header's `space`
discriminant, and never branches on it again:

```rust
enum Session {
    Planar(ProfileSession<Domain, Analysis, TransientAnalysis>),
    Spatial(ProfileSession<Domain3, Analysis3, TransientAnalysis3>),
}

enum StageRunner {
    Static(Analysis),          // or Analysis3, inside the Spatial arm
    Modal,                     // one-shot solve; no stepping
    Transient(TransientAnalysis), // or TransientAnalysis3
}
```

This is the wasm-boundary expression of spatial-architecture.md's "two
execution profiles, selected once": one opaque handle exposed to JS, one
`match` on `space` inside decode, nothing downstream ever branches on
profile again.

Stepping is driven by the worker, not run to completion inside one wasm
call:

```ts
session.advance(stepBudget: number): StepOutcome
// { done, stageComplete, stepsTaken, progressSnapshot, recorderBatch? }
```

The worker loop calls `advance(N)` repeatedly, posts the throttled
`progress` message from `progressSnapshot`, appends `recorderBatch` (already
shaped as the `response_blocks` blob layout — see "Results database") to
OPFS, and checks for a pending `cancel` between calls. This is the concrete
mechanism behind "static execution advances a bounded number of steps per
worker turn, then yields" — a step budget passed into wasm, not a separate
polling or interrupt channel. Stage transitions reuse `core`'s existing
multi-phase composition: `Analysis::into_domain`/`TransientAnalysis::
into_domain` hand a finished stage's `Domain` to the next stage's builder,
and `Analysis::set_integrator` swaps the integrator in place for a
same-phase leg change (a cyclic protocol's reversals) — exactly what
`core`'s own M14 multi-phase tests already exercise natively. A stage's
`AnalysisError` stops the sequence: the run is marked `failed` with the
error's structured detail in `error_json`, and later stages are not
attempted.

## Worker protocol and ownership

The main-thread API is coarse and run-oriented:

```ts
{ type: 'run'; runId: string; input: CarapaceInputV1 }
{ type: 'cancel'; runId: string }
{ type: 'query-results'; requestId: string; runId: string; query: ResultQuery }
```

The worker returns `accepted`, throttled `progress`, `complete`, `cancelled`,
and structured `error` messages. Every response carries its run/request ID;
the UI ignores responses for a run it is no longer displaying. Cancellation is
explicit and cooperative, implemented directly by the `session.advance(
stepBudget)` loop described above: each call only ever executes up to
`stepBudget` steps before returning control to the worker's event loop, where
a pending `cancel` message is checked before the next call. This applies to
every stage kind that steps (`static`, `transient`) uniformly; a `modal`
stage's single `advance` call is checked for cancellation only at its start,
since the eigensolve itself doesn't have a natural yield point.

Responsibilities are intentionally narrow:

- `carapace-core`: numerical model and analysis, no JS, worker, or database
  concepts.
- `carapace-wasm`: decode/validate the compiled input, retain tag mappings,
  own an in-progress analysis session, and expose compact observation batches.
- analysis worker: request routing, wasm scheduling, SQLite/OPFS ownership,
  recorder writes, and result queries.
- main thread: authoring UI, run list, progress display, and query rendering.

## Results database

Use actual SQLite compiled to wasm with a worker-local OPFS VFS. The main
thread never opens the database; it obtains all data through the analysis
worker. Do not retain complete response histories in Zustand or transfer them
on every progress event.

Store recorder data in typed-array-sized blocks, rather than one SQL row per
DOF per step:

```text
runs(run_id, model_hash, sequence_hash, engine_version, status,
     created_at, completed_at, error_json)
stages(run_id, stage_id, stage_index, kind, status)
recorders(recorder_id, run_id, target_kind, target_tags_json,
          response_kind, component_layout_json, sampling_spec_json)
response_blocks(recorder_id, block_index, first_sample, sample_count,
                pseudo_time_blob, values_blob)
```

M10's first recorder is selected node displacement/load-factor history for
the static pushover. The schema is intentionally general enough for later
node velocity/acceleration, reactions, element response, modal, and transient
channels. M11 hardens this layer with migrations, indexing/paging, retention,
export, recovery, and large-run tests.

`values_blob`'s row layout (`[pseudo_time, component_0..N]` per sample) is
decided first, and `carapace-wasm`'s `recorderBatch` (see "Decode and session
model" above) is produced in exactly that shape. The worker's write is then a
raw append of the batch's bytes into the blob — no JS-side reshaping loop
between wasm's output and SQLite's input, the same "internal and external
representations should be the same shape from the start" principle
implementation-plan.md §2.7 already applies to in-process result storage,
extended across the wasm boundary.

## Required first acceptance case

A user can author a 2D fiber-section cantilever, run gravity, freeze gravity,
run a displacement-controlled lateral pushover, then query and plot the roof
displacement/load-factor history from the worker database. Its elastic and
post-yield behavior must agree with Carapace's native ForceBeamColumn
acceptance test. Unsupported model input must return an entity-level diagnostic
before a solver is created.
