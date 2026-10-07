# Domain refactor: baseline and comparison runbook

Temporary document (deleted in Phase 5 of [domain-refactor-plan.md](domain-refactor-plan.md)).
Records what "before" looks like and how to compare against it.

## Setup

| Item | Value |
|---|---|
| Baseline tag | `pre-domain-refactor` (commit `8c6b823`) |
| Working branch | `domain-refactor` |
| Baseline worktree | `../carapace-baseline` (detached at the tag, own `target/`; has an untracked copy of `core/examples/baseline_dump.rs` because the example did not exist at the tag) |
| Recorded results | `comparison/baseline/dump_pre-domain-refactor.json`, `comparison/baseline/csv/*_carapace.csv` |

## Correctness checks

Run on the branch at any time (they take about a minute):

```sh
comparison/check_refactor.sh            # result dump diff + cyclic material CSVs
cargo test --workspace --release        # baseline: 236 passed, 0 failed, 1 ignored
```

- `core/examples/baseline_dump.rs` runs 8 model groups through the public API
  (2D two-story P-Delta frame with diaphragms and two-phase load patterns, modal; displacement-
  and force-based fiber cantilevers with displacement control and a reversal; oriented and
  friction zero-length springs; linear and nonlinear transient with ground motion; arc length
  through a degrading peak with a series spring; 3D frame with a rigid diaphragm, static and
  modal, plus a 3D spring; truss with `equal_dof` and a prescribed settlement). It prints 1318
  values (displacements, reactions, element forces, mode frequencies, iteration and
  factorization counts) as JSON. `comparison/compare_dumps.py` diffs two dumps
  (default `rtol 1e-9`, `atol 1e-12`; counters must match exactly in effect).
- To use the dump on the baseline build itself (to regenerate or sanity check):
  `cd ../carapace-baseline && cargo run --release -p carapace-core --example baseline_dump > dump.json`.
- Every dump model fixes its unused DOFs explicitly, so it is valid before and after DOF
  activation. Models that depend on behavior the refactor deliberately changes do not belong in it.
- Expected differences from the baseline: none for these models (Phase 1 is
  behavior-preserving). Any mismatch is investigated, not loosened.
- Not yet in the dump: continuum elements (new in Phase 3, no baseline exists) and the wire
  format (changes in Phase 2; convert golden inputs and compare results instead).

## Performance baseline (native, `--release`, this machine, medians of 3)

`benchmark_frame --model M --stories S --bays B --steps 50`:

| Model | Free DOFs | solve_ms |
|---|---|---|
| elastic 10x3 | 120 | 17.8 |
| fiber 10x3 | 120 | 24.5 |
| fiber 15x4 | 225 | 67.9 (single run) |
| fiber 20x5 | 360 | 142.3 (single run) |
| elastic 30x6 | 630 | 141.6 |

The fiber 30x6 case does not converge at step 14 on the baseline (benchmark limitation, not a
regression); use 20x5.

## Performance baseline (wasm, Node v24, `--warm-runs 2`, 3 runs)

Build: `wasm-pack build wasm-bridge --target nodejs --out-dir ../pkg`; run
`node comparison/wasm_frame.mjs --model M --stories S --bays B --steps 50 --warm-runs 2`.

| Model | solve_ms (3 runs) |
|---|---|
| elastic 30x6 | 164.5, 163.9, 171.2 |
| fiber 10x3 | 30.3, 33.2, 30.3 |
| fiber 20x5 | 169.4, 169.4, 195.4 |

`carapace_wasm_bg.wasm` size (that build): 1,495,510 bytes. Clean release build of
`carapace-core` plus the two examples: about 48 s; `wasm-pack` build: about 70 s.

Wasm elastic 30x6 is within about 15% of native, so native timing is a usable proxy, but the plan
still requires re-measuring in wasm because allocation behavior differs.

## Gate for each phase

1. `comparison/check_refactor.sh` passes.
2. `cargo test --workspace --release` passes (test counts may grow; none may be removed without
   a replacement).
3. Re-run the benchmarks above on the branch (native and wasm) and compare to the tables.
   Phase 1 target: within a few percent of baseline.
4. Rebuild the wasm package and compare its size.

## Results by step

| Step | Dump vs baseline | Tests (`cargo test --workspace --release`) | Native perf vs baseline | Wasm perf / size vs baseline |
|---|---|---|---|---|
| 1.1 `DofTable` | identical (max relative diff 0) | 240 passed | elastic 30x6 about +2% (141.5 vs 138.5 ms) | not measured |
| 1.2 element interface | identical (max relative diff 0) | 243 passed (3 new multi-node tests) | elastic 30x6 about +5% (147.6 vs 141.5 ms); elastic 10x3, fiber 10x3/15x4/20x5 within noise | elastic 30x6 166 ms (+1%), fiber 10x3 30.0 ms, fiber 20x5 172 ms (+2%, noisy); wasm 1,482,444 bytes (-0.9%) |

Notes on 1.2:
- The elastic 30x6 native gap (about 4-6 ms over 50 steps) is the one measurable cost of the sink
  indirection and table lookup. Two attempts to remove it (caching the per-node table lookup within
  an element; `#[inline]` on the hot helpers) had no effect and were reverted. Wasm shows no
  significant gap, and Phase 4 (pattern-cached assembly, benchmark-gated) targets assembly cost
  directly, so this is carried forward rather than chased now.
- `cargo fmt` is not clean on the baseline for files unrelated to this work (for example
  `benchmark_frame.rs`, `decode.rs`); commits only reformat files they otherwise change.
