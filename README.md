# Carapace

A thin, curated, WebAssembly-native finite element analysis engine, written
in Rust, covering the subset of structural element formulations, material
models, and analysis types actually used in practice by
[pysees](../pysees), Carapace's companion frontend. It's designed from day
one to run efficiently inside a browser Web Worker as a `.wasm` module.

Carapace is heavily inspired by [OpenSees](https://opensees.berkeley.edu/) /
[Xara](https://xara.so) but is not a port of either — it's a from-scratch
reimplementation. See [`docs/xara-feasibility.md`](docs/xara-feasibility.md)
for why a direct C++-to-wasm port wasn't the path taken.

## Features

Both a 2D (`ux, uy, rz`) and 3D (`ux, uy, uz, rx, ry, rz`) execution
profile exist side by side as one generic implementation (const-generic
over `NDIM`/`NDOF`/`ELEMENT_DOF`), not a hand-duplicated second engine —
see each element/section entry below for which profile(s) it covers.

**Elements**
- [x] Truss (2D + 3D)
- [x] ZeroLength — independent per-DOF materials (2D + 3D; 3D springs are
  axis-aligned only, no arbitrary orientation)
- [x] ElasticBeamColumn — axial, bending, torsion, `Linear`/`PDelta`
  transforms (2D + 3D)
- [x] DispBeamColumn — fiber-discretized, displacement-based (2D uniaxial,
  3D biaxial)
- [x] ForceBeamColumn — fiber-discretized, force-based/flexibility method
  (2D uniaxial, 3D biaxial)
- [x] 2D co-rotational transform (large-displacement geometry)
- [ ] 3D co-rotational transform
- [ ] Nonlinear torsion (3D frame torsion is a decoupled elastic `G*J` term
  today, not fiber-derived)

**Materials**
- [x] Elastic
- [x] ElasticPP (elastic-perfectly-plastic, with permanent set)
- [x] Gap (one-sided gap) / Ent (no-tension)
- [x] Steel01 (bilinear kinematic + isotropic hardening)
- [x] Steel02 (Menegotto-Pinto smooth transition curve)
- [x] Concrete01 (Kent-Scott-Park envelope, no tension)
- [x] Concrete02 (Concrete01 + linear tension softening)
- [x] Hysteretic (multi-linear pinching + damage backbone)
- [x] Pinching4 (four-point pinching + independent cyclic-damage rules)
- [x] Parallel / Series / MinMax (material combinators — composite
  stress/strain coupling and strain-bound failure wrapping over any of the
  above)

**Analysis**
- [x] Static analysis: load control, displacement control; `Linear` and
  `NewtonRaphson` algorithms
- [x] Multi-phase/staged analysis with frozen (held-constant) load patterns
- [x] Modal analysis (Lanczos eigensolver)
- [x] Transient analysis (Newmark-β, Rayleigh damping, ground motion)
- [x] Sparse direct solver
- [ ] Additional algorithms: line search, Krylov acceleration, event-to-event
  stepping (design in [`docs/algorithms.md`](docs/algorithms.md))

**Constraints**
- [x] `equal_dof` (2D + 3D)
- [x] Rigid diaphragm — 2D (identity alias) and 3D (true affine constraint
  with lever-arm coupling, not identity aliasing)
- [ ] Chained/nested rigid diaphragms

**wasm / browser boundary**
- [x] `CarapaceInputV1` wire format: panic-free decode plus a stepped
  `Session`/`advance` API, budget-driven for cooperative cancellation
- [x] `wasm_bindgen` boundary exposing `decodeInput`/`advance` to JS (first
  pass: structured-clone transfer via `serde-wasm-bindgen`)
- [ ] Zero-copy transferable-typed-array wire format for `postMessage` (the
  wasm→JS leg is still structured-clone; the analysis-worker→storage-worker
  leg, entirely on the `pysees` side, is already a transferred `ArrayBuffer`)

**Results / persistence** — see
[`docs/results-storage-indexeddb.md`](docs/results-storage-indexeddb.md) for
the full design (IndexedDB-based, implemented on the `pysees` side, not
SQLite/OPFS)
- [x] Bounded, per-`advance()`-call recorder batches — the solver never
  retains a whole-run history in memory
- [x] IndexedDB storage worker: dense run-wide row layout, idempotent
  writes, interrupted-run reconciliation on reload
- [x] Analysis-worker ↔ storage-worker `MessageChannel` wiring with
  backpressure (bounded in-flight-byte budget)
- [x] Node displacement recording
- [ ] Node velocity/acceleration recording (the underlying data already
  exists on `Node` in `core`; just not wired through the wire format or
  storage schema)
- [ ] Reaction-force recording (no reaction computation exists in `core` at
  all yet — this needs new engine capability, not just plumbing)
- [ ] Element/section response recording — basic forces, fiber stress/strain
  (no public `core` API exposes an element's force/section state yet)
- [ ] Modal results (mode shapes/frequencies) or transient time-histories
  through this pipeline — `decode()` currently rejects `Modal`/`Transient`
  stages outright, even though `core` fully implements both
- [ ] Results UI beyond a minimal debug panel — paged queries work, but
  nothing yet consumes them for real plots or deformed-shape scrubbing
- [ ] Saved-run browsing, delete, and export UI (the backend calls exist —
  `deleteRun`, `queryResults` — but no UI surface calls them for run
  management)

Superseded planning docs (the original milestone-by-milestone implementation
plan, the 3D-profile architecture rationale, the pysees handoff contract,
and a rejected Turso-based storage alternative) live in
[`docs/obsolete/`](docs/obsolete/) — useful for the reasoning behind past
decisions, no longer maintained as living specs.

## Workspace layout

```
carapace/
├── core/           # carapace-core: all FE logic. No wasm-bindgen dependency —
│                    compiles and tests identically on native and wasm32 targets.
├── wasm-bridge/     # carapace-wasm: thin wasm-bindgen layer exposing core's API
│                    to JS. Should contain no numerical logic of its own.
└── docs/
    ├── algorithms.md               # design for not-yet-built algorithm extensions
    ├── results-storage-indexeddb.md # current results-persistence design (pysees side)
    ├── xara-feasibility.md          # why this is a from-scratch reimplementation, not a port
    └── obsolete/                    # superseded planning docs, kept for history
```

## Building

```bash
# native (fast iteration, real debuggers, no browser/Node round-trip)
cargo test -p carapace-core

# wasm32
cargo build -p carapace-wasm --target wasm32-unknown-unknown --release

# JS glue (wasm-bindgen-cli version MUST match the wasm-bindgen crate version
# in wasm-bridge/Cargo.toml — mismatches are the most common footgun here)
wasm-bindgen target/wasm32-unknown-unknown/release/carapace_wasm.wasm \
  --out-dir pkg --target nodejs

node -e "console.log(require('./pkg/carapace_wasm.js').axial_displacement(50, 100, 2, 30000))"
```

Always verify a change on the native target first, wasm32 second — it
isolates logic bugs from platform/toolchain bugs faster than debugging
inside a wasm runtime directly.

## Read next

Start with the **Features** checklist above for what's implemented today.

For the frontend/worker boundary and results persistence, read
**[`docs/results-storage-indexeddb.md`](docs/results-storage-indexeddb.md)**
— the current run lifecycle, transport, and IndexedDB results-storage design
(implemented on the `pysees` side).

For planned-but-not-built algorithm work (line search, Krylov acceleration,
event-to-event stepping), read
**[`docs/algorithms.md`](docs/algorithms.md)**.

For the reasoning behind building a from-scratch engine instead of porting
Xara/OpenSees to WebAssembly, read
**[`docs/xara-feasibility.md`](docs/xara-feasibility.md)**.

For historical design rationale — the original milestone plan, the
3D-profile architecture decisions, the pysees handoff contract, and a
rejected Turso-based storage alternative — see
**[`docs/obsolete/`](docs/obsolete/)**. These are no longer maintained as
living specs; treat them as an explanation of *why*, not a current *what*.
