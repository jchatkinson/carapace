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
- [x] ZeroLength — per-DOF materials and Coulomb friction coupling (2D + 3D), with optional
  OpenSees-style `-orient` local axes (default: global axes)
- [x] ElasticBeamColumn — axial, bending, torsion, `Linear`/`PDelta`
  transforms (2D + 3D)
- [x] DispBeamColumn — fiber-discretized, displacement-based (2D uniaxial,
  3D biaxial)
- [x] ZeroLengthSection — fiber section with optional independent springs (2D + 3D)
- [x] ForceBeamColumn — fiber-discretized, force-based/flexibility method
  (2D uniaxial, 3D biaxial)
- [x] 2D co-rotational transform (large-displacement geometry)
- [ ] 3D co-rotational transform — scoped, deliberately shelved: unlike planar
  rotations (a single commuting scalar angle), 3D rotations don't commute, so
  an objective large-rotation formulation needs persistent per-element
  quaternion-tracked nodal state (updated at `commit`, not just a `Corotational2d`
  -style pure function of current position) plus a dense, easy-to-mis-sign
  analytic geometric-stiffness Hessian (Xara/OpenSees's `CorotCrdTransf3d` — its
  `T` matrix and five `ksigma1`–`ksigma5` blocks). More complexity than currently
  justified; revisit only if a real model needs large 3D frame rotations.
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
  `Newton` (full/modified/initial-tangent, optional line search) algorithms
- [x] Arc-length continuation (`Integrator::ArcLength`): scaled spherical
  constraint, coupled displacement/load-factor Newton through load limits,
  softening and snap-back, adaptive radius with transactional cutbacks, stop
  criteria with exact landing, bifurcation diagnostics (design and
  verification in [`docs/arclength.md`](docs/arclength.md))
- [x] `ConvergenceTest::Combined`: per-equation force/moment equilibrium,
  with every Newton convergence check made at the accepted (post-update)
  state
- [x] Multi-phase/staged analysis with frozen (held-constant) load patterns
- [x] Modal analysis (Lanczos eigensolver)
- [x] Transient analysis (Newmark-β, Rayleigh damping, ground motion),
  including a Newton corrector (`TransientAnalysis::with_algorithm`) for
  nonlinear response driven dynamically — `Algorithm::Linear` (one
  effective-system solve per step, exact only for linear-material
  response) remains the default
- [x] `TangentStrategy` (current/reuse-at-step-start/initial-tangent),
  optional `LineSearch` (bisection/regula falsi), and `KrylovNewton`
  (stale-tangent + Krylov-subspace acceleration) — one shared iteration
  driver reused by both static and transient analysis (design in
  [`docs/algorithms.md`](docs/algorithms.md) §0–§6.1)
- [x] Sparse direct solver, with tangent-factorization caching
  (`SparseSolver::factor`/`SparseFactorization`) reused across iterations
  by `TangentStrategy::ReuseAtStepStart`/`Initial`
- [ ] Event-to-event stepping (design in
  [the archived algorithm design](docs/obsolete/algorithm-design.md) §9–§17,
  not yet implemented)

**Constraints**
- [x] `equal_dof` (2D + 3D)
- [x] Rigid diaphragm — 2D (identity alias) and 3D (true affine constraint
  with lever-arm coupling, not identity aliasing)
- [ ] Chained/nested rigid diaphragms

**wasm / browser boundary**
- [x] `CarapaceInputV1` wire format: input decoding plus a stepped
  `Session`/`advance` API, budget-driven for cooperative cancellation
- [x] `wasm_bindgen` boundary exposing `decodeInput`/`advance` as JS objects
  and arrays via `serde-wasm-bindgen`, with generated TypeScript definitions
- [ ] Zero-copy transferable-typed-array wire format for `postMessage` (the
  wasm→JS leg still converts to JS objects; the analysis-worker→storage-worker
  leg, entirely on the `pysees` side, is already a transferred `ArrayBuffer`)
- [x] Configurable `TangentStrategy`/`LineSearch`/`KrylovNewton` through
  `CarapaceInputV1`, for static and transient stages. Existing `"linear"` and
  `"newtonRaphson"` strings remain valid; object forms carry richer settings.
  Transient stages default to linear analysis when `algorithm` is omitted.
  See [the JS solver configuration guide](docs/algorithms.md).
- [x] Arc-length stages (`integrator: { kind: "arcLength", ... }`) and
  `combined` convergence through `CarapaceInputV1`, with per-step
  continuation diagnostics in `StepOutcome.continuation`

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
- [x] Engine recorders: node displacement, velocity/acceleration, reactions,
  element local forces, and fiber stress/strain (2D + 3D)
- [x] Engine modal frequencies/shape components and transient time histories
  through per-call recorder batches
- [ ] pysees compiler/storage integration for recorders beyond displacement,
  modal/transient stages, and spatial models. The engine supports these;
  the frontend compiler currently emits planar static displacement runs only.
- [ ] Results UI beyond a minimal debug panel — paged queries work, but
  nothing yet consumes them for real plots or deformed-shape scrubbing
- [ ] Saved-run browsing, delete, and export UI (`deleteRun` and `queryResults`
  backend calls exist; export and a run-management UI remain unimplemented)

Superseded planning docs (the original milestone-by-milestone implementation
plan, the 3D-profile architecture rationale, the pysees handoff contract,
algorithm and friction derivations, and a rejected Turso storage alternative) live in
[`docs/obsolete/`](docs/obsolete/) — useful for the reasoning behind past
decisions, no longer maintained as living specs.

## Workspace layout

```
carapace/
├── core/           # carapace-core: all FE logic. No wasm-bindgen dependency —
│                    compiles and tests identically on native and wasm32 targets.
├── wasm-bridge/     # carapace-wasm: input decoding, session API, and JS adapter
│   └── src/verification.rs # milestone examples for wasm/Node verification
└── docs/
    ├── algorithms.md               # current solver API and JS configuration
    ├── arclength.md                # arc-length continuation design record
    ├── results-storage-indexeddb.md # current results-persistence design (pysees side)
    ├── xara-feasibility.md          # why this is a from-scratch reimplementation, not a port
    └── obsolete/                    # superseded planning docs, kept for history
```

## Building

```bash
# native (fast iteration, real debuggers, no browser/Node round-trip)
cargo test --workspace

# one-time tool setup
cargo install wasm-pack --locked

# build wasm, generate JS glue and TypeScript types, and package for Node
wasm-pack build wasm-bridge --target nodejs --out-dir ../pkg

# check the generated package through its application API (Node 22.18+)
node wasm-bridge/tests/boundary-smoke.ts

# check the generated declarations and their use by the smoke test
npm exec --package=typescript -- tsc --strict --noEmit --target es2022 --lib esnext,dom \
  --module nodenext wasm-bridge/tests/boundary-smoke.ts
```

For a browser or browser worker, build with `--target web --out-dir ../pkg-web`
instead. This produces ES modules; initialize the wasm module with its default
`init()` export before calling `decodeInput`.

Always verify a change on the native target first, wasm32 second — it
isolates logic bugs from platform/toolchain bugs faster than debugging
inside a wasm runtime directly.

`wasm-pack` manages the matching `wasm-bindgen` tool and generates everything
in `pkg/`, including `carapace_wasm.d.ts`. Do not edit generated files.
The application interface is `decodeInput(input)` → `WasmSession`, with
`advance(stepBudget)` → `StepOutcome` and `currentStageId()` methods.
The milestone functions remain exported for compatibility; their implementation
lives in `wasm-bridge/src/verification.rs`.

Wire-format types in `wasm-bridge/src/input_v1/` derive `tsify::Tsify` alongside
Serde, so fields, camelCase names, and tagged enums generate TypeScript
definitions from the Rust source. The adapter uses `Ts<T>` to connect those
definitions to function signatures, with explicit fallible conversion inside
each call. Malformed inputs still throw an error, and decode errors retain
their structured `{ kind, ... }` shape (also exported as `DecodeError`).
The current payload uses ordinary JS objects and arrays; these declarations
do not describe the planned transferable typed-array format.

The full model/input format (tables, materials, stages, recorders, results,
and errors) is documented in [docs/input-format.md](docs/input-format.md).

When extending the engine, update the Rust input types and decoder as needed,
then rebuild to regenerate the JS types. A new binding is only needed for a
new operation exposed to JavaScript; adding a material or solver behind the
existing session interface does not require one.

## Read next

Start with the **Features** checklist above for what's implemented today.

For the frontend/worker boundary and results persistence, read
**[`docs/results-storage-indexeddb.md`](docs/results-storage-indexeddb.md)**
— the current run lifecycle, transport, and IndexedDB results-storage design
(implemented on the `pysees` side).

For implemented solver settings, JS examples, and remaining event-to-event
stepping design, read
**[`docs/algorithms.md`](docs/algorithms.md)**.

For the reasoning behind building a from-scratch engine instead of porting
Xara/OpenSees to WebAssembly, read
**[`docs/xara-feasibility.md`](docs/xara-feasibility.md)**.

For historical design rationale — the original milestone plan, the
3D-profile architecture decisions, the pysees handoff contract, friction and algorithm design rationale, and a
rejected Turso-based storage alternative — see
**[`docs/obsolete/`](docs/obsolete/)**. These are no longer maintained as
living specs; treat them as an explanation of *why*, not a current *what*.
