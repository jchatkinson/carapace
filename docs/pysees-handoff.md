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
`ArrayBuffer` detaches it from its sender.

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
explicit and cooperative: static execution advances a bounded number of steps
per worker turn, then yields so a cancel message can be handled.

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

## Required first acceptance case

A user can author a 2D fiber-section cantilever, run gravity, freeze gravity,
run a displacement-controlled lateral pushover, then query and plot the roof
displacement/load-factor history from the worker database. Its elastic and
post-yield behavior must agree with Carapace's native ForceBeamColumn
acceptance test. Unsupported model input must return an entity-level diagnostic
before a solver is created.
