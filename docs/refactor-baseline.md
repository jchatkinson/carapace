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

| 1.3 DOF activation | identical (max relative diff 0) | 256 passed, debug and release (13 new: 10 activation core tests, 2 mask-conformance, 1 wasm) | elastic 30x6 about +10% (156 vs 141 ms), fiber 20x5 about +10% (160 vs 146 ms) in the default release profile; **parity** with `lto = "fat"` + `codegen-units = 1` (142-144 vs 139-141; 143-145 vs 146-147) | not measured (gate at 1.5) |

| 1.4 general constraints | identical (max relative diff 0) | 278 passed, debug and release (19 new constraint/session tests) | see the Phase 1 gate | see the Phase 1 gate |

Notes on 1.2:
- The elastic 30x6 native gap (about 4-6 ms over 50 steps) is the one measurable cost of the sink
  indirection and table lookup. Two attempts to remove it (caching the per-node table lookup within
  an element; `#[inline]` on the hot helpers) had no effect and were reverted. Wasm shows no
  significant gap, and Phase 4 (pattern-cached assembly, benchmark-gated) targets assembly cost
  directly, so this is carried forward rather than chased now.
- `cargo fmt` is not clean on the baseline for files unrelated to this work (for example
  `benchmark_frame.rs`, `decode.rs`); commits only reformat files they otherwise change.

Notes on 1.3:
- Activation adds no work per step (identical equation counts for these models); the native gap
  appears only in the default release profile (16 codegen units, no LTO) and disappears with
  `lto = "fat"`, `codegen-units = 1`, on the branch, with the baseline essentially unchanged by the
  same settings. So it is an inlining artifact of the generic sink path, not extra computation.
  Timings drift by a few percent over a session on this machine, so only same-session A/B runs are
  meaningful. A `[profile.release]` change would also affect the wasm build and is left as an explicit
  decision for the Phase 1 gate (1.5), not made silently.
- `cargo test` (debug) is where the assembly guard against under-declared `dof_mask`s is active
  (`#[cfg(debug_assertions)]`); the gate runs both profiles.

## Phase 1 gate (step 1.5, after 1.4)

Correctness: `comparison/check_refactor.sh` reports all 1318 dump values and all 11 cyclic material
CSVs identical to the baseline (largest relative difference 0); `cargo test --workspace` passes in
both profiles (278 passed, 0 failed, 1 ignored; baseline 236, so 42 new tests and none removed).

Performance, same session, medians of 6 runs (50-step frame benchmark, `solve_ms`):

| Build | elastic 30x6 (630 DOF) | fiber 20x5 (360 DOF) |
|---|---|---|
| native, default release profile: baseline | 138.6 | 147.9 |
| native, default release profile: branch | 150.4 (+8.5%) | 152.6 (+3%) |
| native, `lto = "fat"` + `codegen-units = 1`: baseline | 139.0 | 148.0 |
| native, `lto = "fat"` + `codegen-units = 1`: branch | 142.0 (+2%) | 144.5 (-2%) |
| wasm (Node v24, 2 warm runs): baseline | 167.5 | 169.7 |
| wasm: branch | 161.1 | 168.3 |

Wasm fiber 10x3: 29.4 ms branch vs 30.0 ms baseline. `carapace_wasm_bg.wasm`: 1,511,474 bytes vs
1,495,510 baseline (+1.1%).

Reading: the browser target (wasm) is at parity or slightly faster. The native gap in the default
profile is an inlining artifact (it disappears under LTO with one codegen unit), not extra work per
step. Whether to adopt `[profile.release] lto = "fat", codegen-units = 1` is an open decision (it
changes native and wasm builds; wasm-pack uses the same profile).

## Phase 2 gate (wire-format unification, steps 2.1-2.4)

Correctness: `comparison/check_refactor.sh` reports all 1318 dump values and the cyclic material
CSVs identical to the baseline (largest relative difference 0). `cargo test --workspace` passes in
both profiles: 289 passed, 0 failed, 1 ignored (Phase 1: 278; 11 new wire-format tests in
`wasm-bridge/tests/unified_wire.rs`, none removed). `node wasm-bridge/tests/boundary-smoke.ts` and its
`tsc --strict` check pass against the regenerated `pkg/`.

| Step | Dump vs baseline | Tests | Wasm perf (medians of 3 runs, 2 warm runs) | Wasm size |
|---|---|---|---|---|
| 2.1-2.4 unified wire format | identical (max relative diff 0) | 289 passed | elastic 30x6 159-167 ms (Phase 1: 161), fiber 10x3 28.6-29.4 ms (29.4), fiber 20x5 164 ms (168) | 1,442,744 bytes (Phase 1: 1,511,474; -4.5%) |

The smaller binary comes from the single shared decoder (one generic implementation for both
profiles instead of two) and the dropped `*3` serde types. Core is unchanged apart from two
read-only accessors, so native timing was not re-measured.

### `[profile.release] lto = "fat", codegen-units = 1` for the wasm build (not adopted)

Measured with environment overrides (`CARGO_PROFILE_RELEASE_LTO=fat
CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 wasm-pack build ...`); `Cargo.toml` is unchanged. Two
interleaved rounds, medians of 3 runs each:

| Build | elastic 30x6 | fiber 10x3 | fiber 20x5 | wasm size | `wasm-pack` build |
|---|---|---|---|---|---|
| default release | 166.6 / 159.2 ms | 29.4 / 28.6 ms | 163.7 / 163.7 ms | 1,442,744 bytes | about 32 s |
| fat LTO, 1 codegen unit | 158.3 / 158.4 ms | 28.2 / 27.9 ms | 162.7 / 165.0 ms | 1,381,355 bytes (-4.3%) | about 78 s |

Wasm timing is within noise (at most about 3%, and not consistently in one direction); the size
saving is real but small, and the build is about 2.4 times slower. Native is where LTO closed the
8% elastic gap (Phase 1 gate). The decision is still open and belongs to the project owner.
