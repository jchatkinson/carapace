# Performance audit: large models

Audit date: 2026-10-13. Trigger: a ~30-story 3D elastic frame stalled after 64 transient steps, and a
pushover took ~30 s after Run completion to display the deformed shape.

Two separate problems: the engine re-factors a constant matrix every step and the recorder payload is very
wide (transient freeze), and the results reader scans the whole stored run before drawing (slow deformed
shape).

## Measured: engine (native release, 3D elastic frame, Newmark transient)

Throwaway benchmark (deleted): stories x bays grid, 3D `ElasticBeamColumn3`, nodal masses, Rayleigh
damping, uniform ground motion, `dt = 0.01`.

| Frame | Free DOFs | Time per step |
|---|---|---|
| 30 stories x 3 bays | 2.9k | 20 ms |
| 30 stories x 5 bays | 6.5k | 115 ms |
| 30 stories x 8 bays | 14.6k | 550 ms |

Doubling the DOFs makes a step ~5x slower (roughly quadratic). At 14.6k DOFs one step breaks down as:

| Phase | Time |
|---|---|
| Numeric LU factorization | 0.47 s (~90%) |
| `Domain` clone (rollback snapshot) | 0.057 s |
| Sparse-matrix build from triplets | 0.033 s |
| Element assembly | 0.018 s, run twice per step |
| Triangular solve | 0.007 s |

Wasm is typically several times slower than native. Not profiled in wasm or the browser.

## Engine findings (by payoff)

1. **Linear transient re-factors an unchanged matrix every step.** `stiffness_coeff*K + mass_coeff*M` is
   constant for elastic material and fixed `dt`. The `Algorithm::Linear` arm in
   `core/src/analysis/transient.rs` reassembles it and calls `solver.solve()`, which factors from scratch.
   Factor once and cache: step drops from ~0.55 s to ~0.1 s.
2. **`multiply_stiffness` reassembles all element tangents just to compute K*x**
   (`core/src/model/domain.rs`), called every Linear step. Cache K for linear elements or accumulate K*x in
   the same assembly pass.
3. **Pushover likely has the same problem (not benchmarked).** `TangentStrategy::Current` in
   `core/src/analysis/state.rs` factors every Newton iteration, and a step assembles the tangent extra times
   for the predictor. Skip refactoring when the tangent is unchanged, or default to a cached tangent when no
   material is nonlinear.
4. **General LU on an SPD matrix.** `core/src/analysis/solver.rs` uses `faer` `SymbolicLu` (COLAMD +
   partial pivoting). Sparse Cholesky/LDL^T with AMD should give less fill and perhaps ~2x faster factoring
   (estimate; benchmark first). Keep LU as fallback for indefinite/nonlinear stages.
5. **Every `step()` clones the whole `Domain`** for rollback. Snapshot only mutable state (node vectors,
   element history), or restore lazily on failure.
6. **`Domain::reaction` and `lumped_mass_at` scan every element per call.** ~60 ms per recorded step for 81
   base nodes x 6 DOFs on the largest frame. Build a node -> elements index once.
7. **Release profile untuned.** No `[profile.release]` (`lto`, `codegen-units = 1`) and the build script does
   not run `wasm-opt -O3`.

## Recorder / boundary findings (estimated from code, not measured)

8. **Default recorders record everything.** `pysees/src/app/lib/analysisBlocks.ts` `toRecorders` records
   disp, reaction, and (dynamic) vel + accel for all nodes, plus force for all elements (12 columns each in
   3D). The 30x8 frame is ~127k scalar recorders, ~1 MB per step.
9. **Every recorder scalar is its own `(time, value)` pair array.** `record_sample` pushes per-scalar samples,
   `advance()` clones them into `RecorderBatch.samples`, and `serde-wasm-bindgen` turns each pair into a JS
   array: ~8M tiny arrays per 64-step advance plus a ~130 MB Rust clone, inside wasm32 memory. Likely cause
   of the freeze at step 64 (progress is only reported every `STEP_BUDGET = 64` steps, so a pause may be
   memory pressure rather than solve time). Fix: one typed array per `advance()`.
10. **`element_end_force` is recomputed once per component** (12x per element, 11 discarded). Only ~13 ms per
    step natively, but compute each element once.

## Results reader findings (estimated from code)

11. **Deformed view blocks on a full-run scan.** `ResultsDriver.tsx` awaits `source.extents(...)` in
    `Promise.all` before drawing. `queryRunExtents` (`pysees/src/app/lib/resultsStorage/db.ts`) opens a
    cursor over every block, structured-cloning each, and walks every recorder component of every sample.
12. **The scan runs twice on a fresh run.** First with no case range (whole run), then again once `cases()`
    resolves and the range changes. `cases()` -> `timeline()` -> `queryColumns(runId, [])` is a third full
    read just to get pseudo-time and stage columns.
13. **Deformed shape only needs displacements,** but extents are computed for forces, vel, accel,
    reactions and shell stresses too.
14. **Dense rows make single-column reads expensive.** Store per-step `pseudoTime`/`stage` in run metadata to
    avoid the full read in `cases()`; compute max |disp| incrementally during the run and store it.

## Suggested order

1. Factor once for linear transient and cache K (items 1, 2). Biggest, simplest win; existing tests stay valid.
2. Narrow default recorders and pack samples into typed arrays (items 8, 9; touches the wire format, so
   update `docs/input-format.md` and PySees `CarapaceInputV1`).
3. Make the deformed view not wait on extents (items 11-13).
4. Adaptive `advance` budget so progress reports every few hundred ms regardless of model size.
5. Items 3-7, benchmarking each.

Open: build a wasm benchmark to confirm items 8 and 9 before changing the wire format.
