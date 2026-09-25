# Batched results storage with IndexedDB

This plan keeps result persistence in the pysees analysis-worker layer and
uses IndexedDB for local browser storage. Carapace emits bounded numeric
batches from `carapace-wasm`; a separate storage worker writes those batches
to IndexedDB and serves result queries. The solver does not own a database
and does not retain a complete recorder history in memory.

This extends the ownership and worker protocol in
[`pysees-handoff.md`](pysees-handoff.md). It does not require SQLite, OPFS,
or changes to `carapace-core`.

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
  ├─ compiles model + sequence; allocates runId and provenance hashes
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
  └─ answers query / export / delete requests
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

Carapace should return samples produced during each bounded `advance()` call
as packed per-recorder buffers (one per wire recorder, node-major/DOF-minor
order — see "IndexedDB layout" above). The analysis worker transposes those
into the dense run-wide row layout and buffers them across multiple
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

## Work by repository

### pysees

1. **Storage schema and service** — add a results-storage module for opening
   and upgrading IndexedDB, `beginRun`, idempotent block writes, paged queries,
   run finalization, and deletion. Keep transaction creation and upgrade logic
   inside the storage worker. On open, reconcile any non-terminal run left
   over from a prior page session to `interrupted` (see "Interrupted runs").
2. **Storage worker and protocol** — add a dedicated worker plus typed request
   and reply definitions. Transfer result buffers through the channel. Add
   timeouts/structured errors for worker startup and transaction failures.
3. **Run orchestration** — update `runCarapace` and
   `carapaceWorkerClient.ts` to create/reuse the storage worker, establish the
   channel, create the immutable run metadata, and pass the port to the
   analysis worker. Keep model compilation and run-ID creation on the pysees
   side.
4. **UI state** — replace `CarapaceRunResult.recorderSamples` as a complete
   in-memory history with run metadata, counts, and status. Keep only current
   progress and the latest sample in transient UI state. Query data on demand
   for plots and deformed-shape views, requesting only the visible sample range
   or downsampled series.
5. **Imported recorder files** — keep the existing external-file import path
   available. It may share the Results view/query adapter later, but this
   implementation should not require converting imported files into the
   IndexedDB run format.
6. **User controls** — expose saved runs, delete-run behavior, and an export
   path for a run's metadata and numeric blocks. Show distinct complete,
   cancelled, solver-failed, storage-failed, and interrupted states.

Relevant current files include `src/app/workers/carapaceWorker.ts`,
`src/app/lib/carapace/carapaceWorkerClient.ts`,
`src/app/store/useAppStore.ts`, and `src/app/types/carapaceRun.ts`.

### Carapace

1. **Batch wire shape** — finalize the recorder-block contract with pysees.
   Extend the current `StepOutcome`/`advance()` surface to return only the
   samples generated in that call, grouped by recorder, with stable recorder
   identity, stage index, sample range, component count, and packed numeric
   data. Use transferable buffers at the JS worker boundary.
2. **Bound memory** — replace `recorder_history: Vec<Vec<(f64, f64)>>` in
   `PlanarSession` with a per-advance buffer (or recorder-local bounded
   buffers drained on each `advance`). Keep only latest sample/count data
   needed for progress. Remove or narrow `recorderSamples()` so it cannot
   materialize the entire run.
3. **No persistence dependency** — do not add IndexedDB, OPFS, SQLite, or
   browser storage APIs to `carapace-core` or `carapace-wasm`. Carapace's
   numerical results remain independent of storage backend and UI lifecycle.
4. **Session metadata** — expose the compact recorder/stage metadata needed
   to label batches and finalize partial runs. The run ID and model/sequence
   hashes remain pysees-owned metadata and need not be folded into the solver
   input unless a concrete diagnostic requires them.
5. **Failure and cancellation contract** — make each `advance()` batch
   deterministic and identify its first sample/count so the caller can safely
   retry persistence. Preserve solver status separately from storage status;
   storage failures are handled by the pysees worker orchestration.

Relevant current files include `wasm-bridge/src/input_v1/session.rs`,
`wasm-bridge/src/boundary.rs`, and `wasm-bridge/src/input_v1/` in Carapace.

## Delivery stages and acceptance gates

1. **Agree on the protocol.** Freeze block row layout, IDs, range semantics,
   stage-boundary behavior, error/status vocabulary, and retry/idempotency
   rules across both repositories.
2. **Carapace bounded-batch slice.** Add the per-advance packed output and
   confirm native session tests show exact sample/time/value ordering across
   multiple `advance()` calls and stage transitions.
3. **pysees storage worker.** Implement IndexedDB schema and lifecycle
   operations, then verify begin/write/reopen/query/finalize/delete behavior
   with synthetic blocks, including duplicate retry and transaction failure.
4. **Worker integration.** Connect the analysis and storage workers. Verify
   transfers, ordering, bounded queue behavior, backpressure, cancellation,
   partial-result retention, and that `complete` waits for all write acks.
5. **Results UI.** Replace full-history state with paged queries and preserve
   the existing visualizations using the query adapter.
6. **Performance check.** In target browsers, compare end-to-end analysis
   time with result persistence enabled and disabled for representative
   recorder counts and run lengths. Record solver steps/second, write queue
   high-water mark, bytes/second, transaction duration, and peak memory. Tune
   block size and queue budget from these measurements. Acceptance is no
   sustained queue growth and no material loss in solver throughput on the
   representative workloads; set a numeric threshold before implementation
   based on the expected run sizes.

## Storage limits and recovery

IndexedDB and OPFS are browser-managed origin storage, subject to quotas and
possible eviction. Query `navigator.storage.estimate()` for a useful
approximation, handle quota errors as storage failures, and offer explicit
run export. Do not describe browser-local persistence as a backup or as
guaranteed archival storage. The app may request persistent storage through
the browser where appropriate, but must continue to handle refusal and
storage loss.

### Interrupted runs

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
finished. Surface `interrupted` runs in the UI the same way as the other
terminal-but-incomplete states (§User controls), with their partial data
still queryable and exportable.

## Open decisions

- Confirm whether every run is persisted by default or whether small/temporary
  runs may remain memory-only.
- Set initial batch and in-flight byte targets after measuring expected result
  sizes and current `stepBudget` behavior.
- Decide the export package format and whether export includes only selected
  recorder channels or the complete run.
- Set a numeric throughput regression threshold for the performance gate.
