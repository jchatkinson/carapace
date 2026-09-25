# Batched results storage with IndexedDB

**Status: implemented**, on both sides, as of 2026-09-24 — this document
now describes the built system, not a proposal. See "Implementation status"
below for exactly what's built vs. still open, and README.md's "Results /
persistence" checklist for the current one-line-per-item summary.

This design keeps result persistence in the pysees analysis-worker layer and
uses IndexedDB for local browser storage. Carapace emits bounded numeric
batches from `carapace-wasm`; a separate storage worker writes those batches
to IndexedDB and serves result queries. The solver does not own a database
and does not retain a complete recorder history in memory.

This extends the ownership and worker protocol in
[`pysees-handoff.md`](obsolete/pysees-handoff.md) (moved to `docs/obsolete/`
as a superseded planning doc, but still the source for the run
lifecycle/worker-protocol concepts this design builds on). It does not
require SQLite, OPFS, or changes to `carapace-core`.

**Scope carried over from that superseded doc's first slice, still true
today:** only node-displacement recorders and static analysis stages flow
through this pipeline. There is no `responseKind`/target-kind dimension
anywhere in this schema — every recorder is implicitly "node displacement."
Velocity/acceleration, reactions, element/section response, and modal/
transient results are not recorded, not because this design excludes them,
but because nothing upstream of it (the wire format, the decoder) produces
them yet. See README.md's "Results / persistence" checklist for the exact
list of what's missing and why each is a different size of gap.

## Goals and invariants

- Keep storage calls and IndexedDB transaction waits off the solver worker.
- Move results in typed numeric buffers, not arrays of per-sample JS objects.
- Bound memory in Carapace, the analysis worker, and the storage worker.
- Persist run provenance and partial results for completed, failed, and
  cancelled runs.
- Report a run as saved only after every accepted batch has been acknowledged
  by a completed IndexedDB transaction.
- Keep result data separate from pysees model/history state and from Carapace's
  numerical core.
- Allow result queries by run, recorder, and sample range without loading every
  run into Zustand.
- Never leave a run permanently stuck in a non-terminal status after a page
  reload or tab close interrupts it mid-run.

## Ownership and flow

```text
pysees main thread
  ├─ compiles model + sequence; allocates runId (provenance hashes: placeholder only, see below)
  ├─ starts analysis worker and results-storage worker
  └─ displays progress and requests paged result queries

analysis worker
  ├─ calls carapace-wasm advance(stepBudget)
  ├─ transfers bounded recorder batches over MessagePort
  └─ observes acknowledgements / applies backpressure

carapace-wasm → carapace-core
  └─ emits only samples produced during the current advance call

results-storage worker
  ├─ owns the IndexedDB connection and schema upgrades
  ├─ batches writes into short readwrite transactions
  └─ answers query / delete requests (no export request type exists yet)
```

Use a `MessageChannel` between the two workers. The main thread creates the
storage worker and passes one port to the analysis worker when starting a run.
The storage worker remains reusable across runs in the page session. All
messages carry `runId`; write and query requests also carry a unique request
or batch ID. The main thread may receive progress and state changes, but it
does not receive complete sample histories as progress messages.

## IndexedDB layout

Database name: `pysees-results`; schema version: `2`.

Every recorded node always carries its full displacement vector — there is no
per-recorder DOF selection at the results-storage layer (a RecorderSpec's own
`dofs` subset still governs the exported openseespy script; it just doesn't
shape this cache). That gives every recorded node the same fixed width
(`dofsPerNode`, the model's `ndf`), which is what makes a single dense,
run-wide row layout possible: `responseBlocks` holds one row per chunk
covering *every* recorded node together, not one row per node or per
recorder. This favors the dominant interactive read (deformed-shape
scrubbing: one row fetch gives every node's position at a time step) over
the single-node time-history read (which scans the run's blocks and slices
its own columns out of each row) — deliberately, per "Delivery stages and
acceptance gates" below: add a transposed, per-node cache only if that scan
turns out to be too slow in practice, not up front.

| Object store | Key | Contents |
|---|---|---|
| `runs` | `runId` | Model and sequence hashes, schema/engine versions, start/end time, status, final stage, error detail, `dofsPerNode`, `nodeCount`, shared `sampleCount`, and storage failure detail if any |
| `stages` | `[runId, stageIndex]` | Stage ID/kind, order, and stage status |
| `recorders` | `[runId, recorderId]` | One row per recorded node (`recorderId` is that node's tag as a string): its `nodeIndex` (position in the dense row) and component labels |
| `responseBlocks` | `[runId, blockIndex]` | `firstSample`, `sampleCount`, optional stage index, and one packed `Float64Array` covering every recorded node |

Each sample is stored row-major as `[pseudoTime, node0.dof0..dof(D-1),
node1.dof0.., ...]` — stride `1 + nodeCount * dofsPerNode`. A given node's
columns start at `1 + nodeIndex * dofsPerNode`. Keep the numeric values in a
binary typed array. Do not serialize samples as JSON or create one object
store entry per sample/node/DOF.

Block indices and sample ranges are monotonic within a run. Make writes
idempotent using `(runId, blockIndex)` plus a batch ID, so a retry cannot
append the same block twice. Update the run's `sampleCount` in the same
transaction as the response blocks it describes. Keep each transaction
short; do not hold a transaction open while waiting for another worker
message. Use IndexedDB's `durability: 'relaxed'` transaction option for these
frequent chunk writes — losing the last unflushed chunk on a crash is an
acceptable trade for not fsync-ing every one — but keep `beginRun`/`finishRun`
at the default durability. IndexedDB supports structured-cloneable binary
data in workers and asynchronous transactions. The browser clones values
passed to `put()`, so coarse batches and transfer of the `ArrayBuffer`
between workers are important to avoid unnecessary copies. Chunk by target
byte size (roughly 1–4 MiB, buffering across several `advance()` calls, not
one block per call), flushing early at stage boundaries and run
completion/cancellation/error so no already-produced data is lost. See
[IndexedDB API](https://developer.mozilla.org/en-US/docs/Web/API/IndexedDB_API)
and [`IDBObjectStore.put()`](https://developer.mozilla.org/en-US/docs/Web/API/IDBObjectStore/put).

## Cross-worker protocol

Initial message shapes; keep these in a shared TypeScript protocol module in
pysees and document the matching Rust/WASM batch fields in Carapace:

```ts
type StorageRequest =
  | { type: 'beginRun'; requestId: string; run: RunMetadata; stages: StageMetadata[]; recorders: RecorderMetadata[] }
  | { type: 'writeBlocks'; requestId: string; runId: string; batchId: string; blocks: ResultBlock[] }
  | { type: 'query'; requestId: string; runId: string; recorderId: string; firstSample?: number; limit?: number }
  | { type: 'finishRun'; requestId: string; runId: string; status: 'complete' | 'failed' | 'cancelled' | 'storage-failed'; detail?: string }
  | { type: 'deleteRun'; requestId: string; runId: string }

// 'interrupted' is never sent by finishRun — it is written only by the
// startup reconciliation pass described in "Interrupted runs" below, for
// runs a prior page session left in a non-terminal status.
type RunStatus = 'running' | 'complete' | 'failed' | 'cancelled' | 'storage-failed' | 'interrupted'

// One row per chunk, covering every recorded node — not one row per node or per recorder (see
// "IndexedDB layout" above for why). recorderId/blockIndex addressing lives in RecorderMetadata's
// nodeIndex and RunMetadata's dofsPerNode/nodeCount instead of on the block itself.
type ResultBlock = {
  blockIndex: number
  firstSample: number
  sampleCount: number
  stageIndex: number
  data: ArrayBuffer // row-major Float64, stride 1 + nodeCount*dofsPerNode; transferred, not cloned
}
```

Storage replies with the same `requestId`, and for writes includes the
acknowledged batch IDs / block keys and committed sample counts. A write
acknowledgement means that the IndexedDB transaction completed; it does not
claim that data is immune to browser storage eviction or device failure.

### Batching and backpressure

Carapace returns samples produced during each bounded `advance()` call as
one entry per recorder that sampled this call (`StepOutcome.
recorderBatches`) — today still `serde-wasm-bindgen`-serialized `(time,
value)` pairs, not packed binary buffers (see "Batch wire shape" above; the
packed, node-major/DOF-minor buffer described here is the target shape, not
yet what crosses the wasm boundary). The analysis worker (`carapaceWorker.
ts`'s `StorageStream`) does the actual packing: it transposes those
per-recorder entries into the dense run-wide row layout and buffers them
across multiple
`advance()` calls into a chunk target (roughly 1–4 MiB, not one block per
call — small per-call blocks defeat the point of batching `put()`s), but it
must flush early on a stage-index change and at run completion/cancellation/
error so no already-produced data is lost or two stages get mixed into one
block. Pair this with a small bounded in-flight-bytes budget (for example,
8 MiB) so a slow storage worker eventually throttles the solver instead of
letting an unbounded queue build up. These are starting values for
measurement, not API guarantees.

The analysis worker can continue solving while writes are in flight. It must
stop advancing when the number or total bytes of unacknowledged blocks reaches
the configured bound, then resume on acknowledgements. This keeps memory
bounded and ensures a slower device eventually applies backpressure instead
of silently dropping results or allowing an unbounded queue. Always flush the
final partial block before sending `finishRun`.

Storage failure is distinct from solver failure. The storage worker returns a
structured failure for the affected batch; the analysis worker stops advancing
at the next yield point, asks Carapace to stop if needed, and marks the run
`storage-failed` with the last committed sample range. Never mark such a run
complete. Cancellation and solver errors still flush already-produced data
and persist their final status after outstanding writes settle.

## Implementation status by repository

### pysees

1. **Storage schema and service — done.**
   `src/app/lib/resultsStorage/db.ts` opens/upgrades IndexedDB, and
   implements `beginRun`, idempotent block writes (`put()` on the
   `(runId, blockIndex)` key plus a `Math.max`-based `sampleCount`
   reconciliation, so a retried batch can't double-count or duplicate),
   paged `query`, `finishRun`, and `deleteRun`. Reconciles any non-terminal
   run to `interrupted` on database open (see "Interrupted runs").
2. **Storage worker and protocol — done**, minus one piece.
   `src/app/workers/resultsStorageWorker.ts` is a dedicated worker with the
   typed request/reply shapes below, reachable both directly and via a
   handed-off `MessagePort`. **Not done:** no timeouts or structured errors
   for worker-startup failure — a storage worker that fails to start has no
   distinct failure path from an ordinary request timing out.
3. **Run orchestration — done.**
   `src/app/lib/carapace/carapaceWorkerClient.ts`'s `runCarapaceOnWorker`
   creates the storage worker, opens the `MessageChannel`, builds run/stage/
   recorder metadata from the compiler's output, and hands the analysis
   worker its port. Model compilation and run-ID creation stay on the
   pysees side, as designed. One gap: `runMetadataFor`'s `modelHash`/
   `sequenceHash` are the placeholder string `'unhashed'` — provenance
   hashing was never implemented.
4. **UI state — done.**
   `CarapaceRunResult` (`src/app/types/carapaceRun.ts`) is metadata/counts/
   status only (`sampleCount`, `recordedNodeTags`, `error`) — the old
   complete-in-memory-history shape this item names is gone. `AnalysisPanel.
   tsx` queries `resultsStorage` on demand for the values it shows.
5. **Imported recorder files** — unaffected by this work either way; not
   verified as part of this audit.
6. **User controls — not done.** `deleteRun`/`queryResults` exist in
   `resultsStorageClient.ts`, but nothing in the UI calls them for run
   management: there is no saved-runs list, no delete-run control, and no
   export path. `AnalysisPanel.tsx`'s results view is a single run's debug
   panel (status, diagnostics, last sample per recorded node), not the
   browsing/export surface this item describes.

Relevant files: `src/app/workers/carapaceWorker.ts`,
`src/app/workers/resultsStorageWorker.ts`,
`src/app/lib/carapace/carapaceWorkerClient.ts`,
`src/app/lib/resultsStorage/db.ts`,
`src/app/lib/resultsStorage/resultsStorageClient.ts`,
`src/app/types/carapaceRun.ts`, `src/app/types/resultsStorage.ts`.

**No automated tests exist for any of the above** — pysees has no test
runner configured at all (`package.json` has no test script). This is the
biggest rigor gap relative to Carapace's Rust side, which has native tests
for the equivalent batching behavior (see below). Treat the "duplicate
retry and transaction failure" and "worker integration" verification named
in "Delivery stages and acceptance gates" below as unverified, not just
undocumented.

### Carapace

1. **Batch wire shape — done, partially.** `StepOutcome::recorder_batches`
   (`wasm-bridge/src/input_v1/session.rs`) returns only the samples one
   `advance()` call produced, grouped by recorder, with stable
   `recorder_index`, `stage_index`, and a deterministic `first_sample`.
   **Not done:** this still crosses the `wasm_bindgen` boundary as a
   structured-clone JS object via `serde-wasm-bindgen`
   (`wasm-bridge/src/boundary.rs`'s own doc comment calls this out), not a
   transferable typed array — `carapaceWorker.ts`'s `StorageStream` does
   the packing into a real `Float64Array`/`ArrayBuffer` itself, on the JS
   side, only for the next hop (to the storage worker).
2. **Bound memory — done.** `PlanarSession` no longer holds
   `recorder_history: Vec<Vec<(f64, f64)>>`; `current_batch` is cleared at
   the top of every `advance()` call and drained into that call's
   `recorder_batches` at the end, with only a running
   `recorder_sample_counts: Vec<u32>` persisting across calls.
   `recorder_samples()` was removed entirely (not narrowed) — there is no
   accessor that can materialize a whole run's history anymore.
3. **No persistence dependency — upheld.** No IndexedDB/OPFS/SQLite/browser
   storage API exists in `carapace-core` or `carapace-wasm`.
4. **Session metadata — not needed, turned out.** No new `Session`/
   `WasmSession` accessor for recorder/stage metadata was added, but none
   was needed: pysees's compiler (`compileInputV1.ts`) already holds
   `recordedNodeTags`/`dofsPerNode` from compiling the model, before decode
   ever runs, so `carapaceWorkerClient.ts` builds `RunMetadata`/
   `RecorderMetadata` from the compiler's output directly. Provenance
   hashes remain unimplemented on the pysees side (see above), not blocked
   on anything here.
5. **Failure and cancellation contract — done.** Every `RecorderBatch`
   carries a deterministic `first_sample`/`stage_index`, verified by
   `wasm-bridge/tests/m10_carapace_input_v1.rs`'s multi-stage,
   multi-`advance()`-call assertions.

Relevant files: `wasm-bridge/src/input_v1/session.rs`,
`wasm-bridge/src/boundary.rs`, `wasm-bridge/tests/m10_carapace_input_v1.rs`.

## Delivery stages and acceptance gates

1. **Agree on the protocol — done.** Block row layout, IDs, range semantics,
   stage-boundary behavior, error/status vocabulary, and retry/idempotency
   rules match across both repositories (verified by reading both, not just
   by the protocol type files agreeing).
2. **Carapace bounded-batch slice — done and tested.** Native session tests
   (`wasm-bridge/tests/m10_carapace_input_v1.rs`) confirm exact sample/time/
   value ordering across multiple `advance()` calls and a stage transition.
3. **pysees storage worker — implemented, not verified.** IndexedDB schema
   and lifecycle operations exist (`db.ts`), but there is no test exercising
   begin/write/reopen/query/finalize/delete with synthetic blocks, duplicate
   retry, or transaction failure — pysees has no test runner at all (see
   above). The idempotency logic reads correctly by inspection; it has not
   been exercised by a test that actually retries a write.
4. **Worker integration — implemented, not verified.** The analysis and
   storage workers are connected (`carapaceWorkerClient.ts`), with transfer,
   ordering, bounded-queue backpressure, and cancellation all present in the
   code (`StorageStream` in `carapaceWorker.ts`). Same caveat as above: no
   automated test drives this end to end.
5. **Results UI — not done.** `AnalysisPanel.tsx` has a minimal debug panel
   (status, diagnostics, last sample per node), not paged queries feeding
   real plots or deformed-shape scrubbing. The query adapter
   (`resultsStorageClient.ts`'s `queryResults`) works and is exercised by
   that debug panel, but nothing consumes it for actual visualization yet.
6. **Performance check — not done.** No measurement has been taken in any
   browser. `FLUSH_BYTE_TARGET`/`IN_FLIGHT_BYTE_BUDGET` in `carapaceWorker.
   ts` are still the doc's original starting guesses (1 MiB / 8 MiB), not
   tuned from data, and their own code comment says so.

## Storage limits and recovery

**Not implemented.** None of this section is built: there is no
`navigator.storage.estimate()` call anywhere in pysees, no explicit
quota-error handling (a quota error would surface only as a generic
`storageError`/`storage-failed`, indistinguishable from any other IndexedDB
failure), no persistent-storage request, and no run export. The design
intent stands as written below.

IndexedDB and OPFS are browser-managed origin storage, subject to quotas and
possible eviction. Query `navigator.storage.estimate()` for a useful
approximation, handle quota errors as storage failures, and offer explicit
run export. Do not describe browser-local persistence as a backup or as
guaranteed archival storage. The app may request persistent storage through
the browser where appropriate, but must continue to handle refusal and
storage loss.

### Interrupted runs

**Done** — `db.ts`'s `reconcileInterruptedRuns`, run once per database open
before any request is served, matches this section exactly.

A page reload or tab close during a run leaves the in-memory analysis worker
and its `carapace-wasm` session gone, but the run's `runs` row in IndexedDB
stays wherever it last committed — there is no `finishRun` to mark it
terminal. Without a check, these accumulate as permanently "running" entries
in the run list.

On storage-worker startup (i.e. each time the app opens the `pysees-results`
database), scan the `runs` store for any row not already in a terminal
status (`complete`, `failed`, `cancelled`, `storage-failed`) and reconcile it
to `interrupted`, preserving whatever stages/blocks/sample counts were
already committed. `interrupted` is a distinct terminal status from
`storage-failed`: the storage layer itself did not fail, the run simply never
finished. Its partial data is queryable today (the reconciliation itself is
done); surfacing it in a run list alongside the other terminal-but-incomplete
states, and exporting it, both wait on §User controls above (not done).

## Open decisions

All four of these are still genuinely open — none have been decided or
resolved by the implementation work above.

- Confirm whether every run is persisted by default or whether small/temporary
  runs may remain memory-only. (Today, every run is unconditionally
  persisted — `runCarapaceOnWorker` always calls `beginRun`.)
- Set initial batch and in-flight byte targets after measuring expected result
  sizes and current `stepBudget` behavior. (Today's `FLUSH_BYTE_TARGET`/
  `IN_FLIGHT_BYTE_BUDGET` are still this doc's original unmeasured guesses.)
- Decide the export package format and whether export includes only selected
  recorder channels or the complete run. (No export exists yet at all.)
- Set a numeric throughput regression threshold for the performance gate.
  (No performance measurement has been taken yet.)

And one this audit surfaced, not in the original plan: **decide the
response-kind/target-kind schema extension** — a `responseKind` dimension
on `RecorderMetadata` and a non-node-indexed row layout for element
response — before any of velocity/acceleration/reaction/element-response
recording can be added. See README.md's "Results / persistence" checklist.
