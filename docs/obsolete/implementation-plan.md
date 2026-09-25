# Carapace Implementation Plan

This document is the source of truth for Carapace's design and roadmap. It exists
so a coding session with no prior context can pick up the project and continue
correctly. If something here conflicts with the code, the code wins for *what
exists*, but this document should still be updated to match — don't let it rot.

Detailed rationale for individual design decisions and ported algorithms lives
in doc comments on the relevant source/test files, not here — this document
stays at the architecture/scope/roadmap level and points into the code for
specifics. See [`docs/xara-feasibility.md`](xara-feasibility.md) for the
investigation into compiling Xara/OpenSees's C++ directly to wasm that this
project's architecture reacts against.

---

## 1. Goals and non-goals

**Goal:** a small, fast, WebAssembly-native FE analysis engine covering the
subset of OpenSees/Xara functionality that `pysees` (the companion frontend)
actually uses, running inside a browser Web Worker.

**Explicit non-goals** (each was a real decision made during design, not an
oversight — don't "fix" these without revisiting the reasoning):

- **Not full OpenSees/Xara parity.** No PFEM, no fluids, no thermal, no
  reliability analysis, no arbitrary user-supplied element plugins. See §3 for
  exactly what *is* in scope.
- **Not dynamically extensible at runtime.** The element/material catalog is
  fixed at compile time (closed `enum`s, not a broker/factory). Adding a new
  element type is a recompile, not a runtime registration. This is a deliberate
  simplification enabled by "thin, curated" — see §2.1.
- **Not multi-tab / multi-writer.** A single worker owns the model and any
  persisted results (SQLite-over-OPFS). Multiple tabs querying the same
  analysis concurrently is explicitly out of scope (decided when designing the
  `pysees` results-storage architecture).
- **Not a live/interactive re-solve loop.** Workflow is "build model once,
  then analyze," not "drag a node, re-solve every frame." This shapes the
  worker/main-thread messaging design in §4.
- **Not parallel (MPI) or database-backed.** No `sendSelf`/`recvSelf`,
  no `Channel`, no `ObjectBroker`. Single-threaded, in-memory, one process.

---

## 2. Architecture principles

### 2.1 Dispatch: closed `enum`s for catalogs, not broker/factory patterns

Elements and materials are a **fixed, compile-time-known set** (§3). Represent
each as a Rust `enum`:

```rust
enum Element { Truss(Truss), ZeroLength(ZeroLength), ElasticBeamColumn(..), .. }
enum Material { Elastic(Elastic), Steel01(Steel01), Parallel(Vec<Box<Material>>), .. }
```

`match`-based dispatch, no `Box<dyn Trait>`, no runtime registration — see
`docs/xara-feasibility.md` for the broker-eager-construction problem this
sidesteps by construction (there's no broker to eagerly construct).

**Exception — analysis strategy types are different and get enums too, but for
a different reason.** `ConvergenceTest`, `Integrator`, `Algorithm`,
`ConstraintHandler` are also closed enums, but not because a broker would be
"toxic" for them (their interfaces are narrow and homogeneous, unlike the
heavy, heterogeneous `Element`/`Material` catalog) — it's simply that the set
is small, known, and user-selected per analysis, so an enum is the simplest
correct representation. `SparseSolver` (§2.6) is not an enum or trait at all —
just a concrete struct, since Carapace commits to exactly one solver
implementation.

### 2.2 Ownership via the type system

No raw pointers with ownership documented only in a comment. Model state lives
in an arena/`SlotMap`-style store keyed by generational indices, not `Box`/`Rc`
graphs mirroring a C++ pointer structure. Where a type must own another
(e.g. an `Analysis` owning its `Integrator`), plain Rust ownership (`Integrator`
by value, not `&mut Integrator` with an ownership convention) makes the
question unrepresentable rather than a runtime footgun.

### 2.3 No serialization / broker / interpreter machinery

No `sendSelf`/`recvSelf`, no `Channel`, no `ObjectBroker`, no embedded
interpreter of any kind. These exist in OpenSees/Xara for MPI domain
decomposition and database persistence — irrelevant to a single-threaded wasm
worker, and a major source of dependency-closure bloat in that architecture
(see `docs/xara-feasibility.md`).

### 2.4 Fixed-size, stack-allocated element-local linear algebra

A 2-node truss's stiffness matrix is always 4×4 (2D) or 6×6 (3D) — known at
compile time. Use const-generic fixed-size types for element-local
stiffness/state (`nalgebra::SMatrix<f64, N, N>`) rather than a generic
runtime-sized heap-allocated matrix. Zero heap allocation in the per-element
evaluation hot loop.

### 2.5 Typestate-enforced analysis composition

An `AnalysisBuilder<Unwired>` can only become `AnalysisBuilder<Ready>` once
every required piece is set, and only `Ready` exposes `.build()` →
`Analysis`, which exposes `.step()`. Illegal sequencing is a compile error,
not a runtime ordering bug. See §4.4 for the sketch, and
`core/src/analysis/builder.rs` for the real implementation and its doc
comment.

### 2.6 One fixed sparse solver, no solver abstraction layer

Carapace commits to **one** sparse solver implementation — no trait, no enum,
just a concrete `SparseSolver` struct (`core/src/analysis/solver.rs`). Modern
sparse LU libraries typically do their own fill-reducing ordering internally,
which also removes the need for a separate DOF-numbering abstraction — see
§5 decision #1.

### 2.7 Data-oriented (SoA) result storage

Element response history (strain/stress/tangent per step) is stored in flat,
contiguous, per-response-type arrays, not one heap object per element per
timestep. This shape matches what the results pipeline needs anyway (typed
arrays over `postMessage`, rows into SQLite — §4.3), so the internal and
external representations should be the same shape from the start rather than
converted at the boundary.

### 2.8 Errors are `Result<T, E>`, not sentinel integers

A real error enum (`AnalysisError::FailedToConverge{ step }`,
`::SingularSystem`, `::InvalidConstraint`, ...) carrying context — see
`core/src/analysis/error.rs`.

---

## 3. Feature scope

Scope is exactly what `pysees` uses today. Anything not listed here is out of
scope until a concrete need arises — don't speculatively add
elements/materials "while you're in there."

### 3.1 Elements — rough milestone order, see §6 for status

The implemented catalog is planar (`NDM = 2`, `NDF = 3`) today. Spatial
(`NDM = 3`, `NDF = 6`) is planned as a separate concrete execution profile,
not a mode that pads every planar model with unused DOFs; see
[`spatial-architecture.md`](spatial-architecture.md).

| Element | Complexity | Notes |
|---|---|---|
| Truss | Trivial | M1. |
| ZeroLength | Simple | No geometry/integration; direct material evaluation per DOF. |
| ElasticBeamColumn | Simple | Closed-form stiffness, no iteration. |
| DispBeamColumn | Moderate | Needs `BeamIntegration` (Gauss-Legendre/Lobatto) + fiber sections. One Newton loop total (domain-level). |
| ForceBeamColumn | **High — own milestone (M8)** | Nested element-level equilibrium iteration, architecturally distinct from every other element — see `docs/xara-feasibility.md` #5. |

Fiber sections: mechanically simple (loop over fibers, accumulate
stress-resultant + tangent from each fiber's uniaxial material response,
weighted by area/position) but evaluated many times per element per
iteration — keep allocation-free (§2.4/§2.7).

### 3.2 Materials (uniaxial)

**Leaf materials** (flat enum variants, given strain → stress + tangent from
internal state): Elastic, ElasticPP (EPP), Gap, ENT ("no tension"),
Hysteretic, Pinching4, Concrete01, Concrete02, Steel01, Steel02.

**Composite/wrapper materials** (recursive — hold other materials):
Parallel, Series, MinMax — see `core/src/model/materials/mod.rs`'s doc
comment for the trial/commit design these (and every leaf) share, and its
"Porting a new leaf material" section for the recipe to follow when adding
another.

### 3.3 Analyses

- **Static / quasi-static** — same code path (`LoadControl` or
  `DisplacementControl` integrator with small increments); not a distinct
  architecture. Multi-phase composition (gravity → pushover, a cyclic
  quasi-static protocol) and load-controlled load-path protocols are
  caller-side composition of this same code path, not separate analysis
  types — see M14 (§6).
- **Modal** — eigenvalue analysis; see §5 decision #2.
- **Time history (transient)** — `Newmark` integration, element/nodal mass
  matrices, Rayleigh damping (linear combination of mass/stiffness), and
  (M14) time-varying load patterns plus ground-motion excitation.

### 3.4 Supporting infrastructure

- **Mass:** lumped (diagonal) mass matrices per element, or nodal mass.
- **Loads:** nodal loads and element loads (distributed loads on
  beam-columns).
- **`geomTransf`:** `Linear` and `PDelta` (small-displacement paths), plus
  the planar co-rotational M9 path for large rotations. Spatial corotation
  remains separate (M18) because it needs independent rotation/objectivity
  treatment.
- **Rayleigh damping:** needed for time-history analysis.

The baseline support above is implemented through M6; planar corotational
geometry is M9 (see §6).

- **Multi-point constraints:** `Domain::equal_dof`/`rigid_diaphragm`
  (Xara/OpenSees's `equalDOF`/`rigidDiaphragm`), resolved by
  `ConstraintHandler::Transformation` — see its doc comment
  (`core/src/analysis/constraint.rs`) for the DOF-equation-aliasing scheme
  and why it's sufficient here without a general coefficient/transformation
  matrix. Deliberately narrower than Xara's `rigidDiaphragm`: this is a 2D
  (`NDF = 3`: x, y, rz) model, so there's no out-of-plane axis for a
  "diaphragm perpendicular to it" to mean anything — `rigid_diaphragm` here
  is the standard 2D-frame simplification of tying constrained nodes'
  x-translation DOF to a retained node's, not true rigid-body (lever-arm)
  kinematics. Done — see `core/tests/m_constraints.rs`.

---

## 4. System architecture (worker / JS integration)

### 4.1 Model description: batch, not per-call RPC

The frontend (`pysees`) builds the full model description as a plain
TypeScript object/typed-array bundle on the main thread (cheap, synchronous,
no wasm involved), then transfers it **once** into the worker (structured
clone for metadata, transferable `ArrayBuffer`s for numeric-heavy data —
node coordinates, element connectivity). The worker walks it once inside
wasm. This avoids mirroring a Python-style `model.node(...)`/
`model.element(...)` object-method API 1:1 across the worker boundary, which
would mean one `postMessage` round-trip per model entity — fine in-process
for Python/pybind11, expensive across a worker boundary for models with
thousands of entities.

### 4.2 Command surface: coarse, not per-node

After the initial batch transfer: `runAnalysis(recipe)`, `getResults(query)`,
not per-entity calls. Matches §1's "no live re-solve loop" non-goal.

### 4.3 Results: typed arrays + SQLite-over-OPFS, decoupled from live progress

Two separate, deliberately decoupled mechanisms:

- **Live progress (for UI, e.g. a progress bar / running plot):** the worker
  `postMessage`s a small throttled snapshot every N steps or N milliseconds
  (not every step) — step index, time, a transferred `Float64Array` of the
  current response vector. No database involved; this is the whole
  mechanism.
- **Full/persisted results:** SQLite compiled to wasm, backed by OPFS,
  **owned entirely by the analysis worker** (main thread never opens its own
  OPFS handle — it queries *through* the worker via messages). Because
  there's deliberately no multi-tab requirement (§1), the `opfs-sahpool` VFS
  (no special headers required) is sufficient — the COOP/COEP-requiring
  `opfs` VFS with `OPFSWriteAheadVFS` concurrent-read-during-write would only
  matter for a multi-reader scenario this project explicitly doesn't have.

### 4.4 Sketch of the typestate analysis API

```rust
struct AnalysisBuilder<S> { /* .. */ }
impl AnalysisBuilder<Unwired> {
    fn constraint_handler(self, h: ConstraintHandler) -> AnalysisBuilder<HasHandler> { .. }
}
// ...chained until fully wired...
impl AnalysisBuilder<Ready> {
    fn build(self, model: &Model) -> Analysis {
        // numbers DOFs, builds sparsity pattern, factors once if linear —
        // done ONCE here, not re-checked every step (§1 non-goals: no live
        // re-solve means no need for mid-loop re-number/re-factor machinery).
    }
}
impl Analysis {
    fn step(&mut self) -> Result<StepResult, AnalysisError> {
        self.integrator.new_step(&mut self.model);
        loop {
            let (k, r) = self.model.form_tangent_and_residual();
            let du = self.solver.factor_and_solve(&k, &r); // concrete, no dispatch
            self.model.update_state(&du);
            match self.test.check(&r, &du) {
                ConvergenceStatus::Converged => break,
                ConvergenceStatus::Continuing if self.iter < self.test.max_iter() => self.iter += 1,
                _ => return Err(AnalysisError::FailedToConverge { step: self.step_count }),
            }
        }
        self.integrator.commit(&mut self.model);
        Ok(StepResult { /* .. */ })
    }
}
```

(This is a sketch for orientation — `core/src/analysis/{builder,state}.rs`
are the real, current implementation.)

---

## 5. Open decisions

1. **Linear algebra / sparse solver crate: resolved.** `SparseSolver`
   (`core/src/analysis/solver.rs`) uses `faer`'s sparse LU (COLAMD/AMD
   fill-reducing ordering, pure Rust, verified on wasm32 + Node) over a dense
   `nalgebra` solve or `nalgebra-sparse` — see that file's doc comment for
   the full reasoning and `core/tests/sparse_solver.rs` for the scale
   verification (150-free-DOF chained-beam system). `Domain` assembles
   stiffness directly into sparse triplets rather than a dense buffer that's
   sparsified after. **Future optimization, not yet done:** refactoring the
   symbolic factorization across Newton iterations (same sparsity pattern,
   only values change) — `SparseSolver::solve` currently re-factors from
   scratch every call.

2. **ARPACK strategy for modal analysis: resolved.** `core/src/analysis/
   modal.rs`'s `modal_analysis` runs a hand-rolled shift-invert Lanczos
   (shift = 0, operator `K^-1 * M`), not a dense full-spectrum eigensolve or
   an ARPACK FFI binding — see that file's doc comment for why (a dense
   solve is wrong beyond just slow: O(n²) memory regardless of model
   sparsity, and it always computes every eigenpair when only the lowest
   `num_modes` are ever needed). The two numerically hard pieces (sparse LU,
   small dense eigendecomposition of the projected tridiagonal matrix) are
   both library code; the outer three-term-recurrence loop with full
   M-orthogonal re-orthogonalization is the only hand-rolled part. Only
   `shift = 0` is implemented (lowest frequencies); a nonzero shift
   (targeting a frequency band, or handling near-singular `K` for buckling)
   is a natural extension, not built until something needs it. The ARPACK
   f2c feasibility spike (`reference/xara-spike/`) remains available
   reference material for an eventual real ARPACK comparison.

3. **`cmx.h`-equivalent small-matrix inversion.** Not actually a Carapace
   decision — `nalgebra`/`faer` both handle small fixed-size matrix inversion
   natively; Xara's missing-implementation bug in this area (see
   `docs/xara-feasibility.md`) has no equivalent risk here.

---

## 6. Milestones

Status summary. Design rationale for each milestone's implementation lives in
the corresponding source file's doc comment; acceptance-test rationale lives
in the corresponding test file's doc comment (`core/tests/m*.rs`,
`wasm-bridge/src/lib.rs`'s exported functions). Verification order for every
milestone: native `cargo test` first, then `wasm32-unknown-unknown` build +
`wasm-bindgen` + Node execution (§7).

- **M0 — Workspace scaffold.** ✅ Done. `carapace-core` (no wasm deps) +
  `carapace-wasm` (thin `wasm-bindgen` layer) workspace; native and wasm32
  build/test verified against a known closed-form answer
  (`axial_displacement`, `core/src/lib.rs`).

- **M1 — Real `Domain`/`Node`/`Truss`/`Elastic` + one real solve.** ✅ Done.
  `carapace-core::model` (`Domain`/`Node`/`Element::Truss`/`Material::Elastic`)
  and `carapace-core::analysis` (typestate `AnalysisBuilder`,
  `Integrator::LoadControl`, `Algorithm::Linear`,
  `ConvergenceTest::NormUnbalance`, `ConstraintHandler::Plain`,
  `SparseSolver`) replace the M0 closed-form placeholder. See
  `core/tests/m1_truss.rs`.

- **M2 — ZeroLength + EPP/Gap/ENT materials.** ✅ Done. `Element::ZeroLength`
  (`core/src/model/elements/zero_length.rs`) evaluates a material
  independently per DOF direction. `Material` gained `ElasticPP`, `Gap`,
  `Ent` variants. See `core/src/model/materials/mod.rs`,
  `core/tests/m2_zero_length.rs`.

- **M3 — ElasticBeamColumn + geomTransf (Linear, PDelta) + element loads.**
  ✅ Done. `Element::ElasticBeamColumn`
  (`core/src/model/elements/elastic_beam_column.rs`): a prismatic 2D
  Euler-Bernoulli beam-column, `GeomTransf` (`Linear`/`PDelta`), and a
  uniform transverse element load. Bumped `NDF` from 2 to 3 (in-plane
  rotation). See `core/tests/m3_beam.rs`.

- **M4 — Analysis composition generalization.** ✅ Done.
  `Integrator::DisplacementControl` (`core/src/analysis/integrator.rs`) and
  `Algorithm::NewtonRaphson`, factored into a predictor/corrector split in
  `Analysis::step`. Added `AnalysisError::InvalidConstraint`. See
  `core/tests/m4_analysis.rs`.

- **M5 — Modal analysis.** ✅ Done. `Node` gained a lumped `mass` field;
  `analysis::modal_analysis` (`core/src/analysis/modal.rs`) solves
  `K*phi = omega^2*M*phi` for the lowest `num_modes`. See §5 decision #2,
  `core/tests/m5_modal.rs`.

- **Post-M5 — `SparseSolver`/`modal_analysis` switched to real sparse
  implementations.** ✅ Done. Not a numbered milestone — see §5 decisions
  #1 and #2. Every pre-existing M1-M5 test passed unmodified against both
  new solvers.

- **M6 — Mass, Rayleigh damping, Newmark, time-history analysis.** ✅ Done.
  Element-consistent lumped mass (`Truss`/`ElasticBeamColumn::form_mass`),
  `RayleighDamping` (`core/src/analysis/damping.rs`), and
  `TransientAnalysis` (`core/src/analysis/transient.rs`, Newmark-beta
  average-acceleration integration) — a type distinct from `Analysis`, see
  that file's doc comment for why. See `core/tests/m6_dynamics.rs`.

- **Pre-M7 — material trial/commit state + `Domain` snapshot/restore.** ✅
  Done. Landed ahead of M7 because fiber sections and path-dependent
  materials need it. See `core/src/model/materials/mod.rs`'s "Trial vs.
  commit" doc comment for the design and `core/tests/material_state.rs` for
  verification. `Domain` gained `Clone` to back `Analysis`/
  `TransientAnalysis`'s snapshot-and-restore-on-failed-step.

- **M7 — DispBeamColumn + fiber sections + full standard material catalog.**
  ✅ Done, landed in four stages given its unusual size relative to every
  other milestone (agreed with the project owner):

  - **Stage 1** — `BeamIntegration::{Legendre,Lobatto}`
    (`core/src/model/integration.rs`), `FiberSection`
    (`core/src/model/fiber_section.rs`), `DispBeamColumn`
    (`core/src/model/elements/disp_beam_column.rs`). See `core/tests/
    m7_disp_beam_column.rs`.
  - **Stage 2** — `Steel01`, `Concrete01`, `Parallel`/`Series`/`MinMax`
    composites (`core/src/model/materials/`).
  - **Stage 3** — `Steel02`, `Concrete02`.
  - **Stage 4** — `Hysteretic`, `Pinching4` (boxed payloads — see
    `core/src/model/materials/mod.rs`'s porting recipe for why).

  Each material's port-specific design notes (state-machine shape,
  deviations from the OpenSees source, damage-tracking subtleties) live in
  that material's own module doc comment
  (`core/src/model/materials/{steel01,steel02,concrete01,concrete02,
  hysteretic,pinching4}.rs`); acceptance-test reasoning lives alongside each
  material's own unit tests in that same file (and `composite.rs` for
  `Parallel`/`Series`/`MinMax`), plus
  `core/src/model/fiber_section.rs`'s composite-section test.

- **M8 — ForceBeamColumn.** ✅ Done.
  `Element::ForceBeamColumn` (`core/src/model/elements/force_beam_column.rs`):
  a flexibility-method 2D fiber beam-column, its `state_determination`
  Newton iteration ported directly from
  `xara/SRC/element/Frame/Other/Force/ForceBeamColumn2d.cpp::update()`
  (the reference C++ checkout, not re-derived from the cited papers by
  hand — see that file's doc comment for the port's scope and a
  bisection-based subdivision fallback for large load increments,
  simplified from the source's own `numSubdivide`/algorithm-ladder
  mechanism). Verified against `ElasticBeamColumn`'s closed-form elastic
  stiffness (symmetric section) and a hand-derived closed form for an
  asymmetric (axial-bending-coupled) section, plus an `ElasticPP`
  past-yield/permanent-set case — see `core/tests/m8_force_beam_column.rs`.

- **M9 — Corotational geomTransf (planar large-displacement).** Implemented;
  dedicated acceptance tests remain to be added. `GeomTransf::Corotational`
  now uses the current chord to form axial extension and both end rotations
  relative to the chord, then maps basic resistance/tangent to global DOFs
  with the exact kinematic Jacobian and its geometric Hessian. The option is
  available on `ElasticBeamColumn`, `DispBeamColumn`, and `ForceBeamColumn`;
  the nonlinear fiber elements opt in with `.with_corotational()`. The
  elastic element's uniform transverse load follows the current local frame.
  Rigid joint offsets and spatial corotation are outside this milestone; the
  latter remains M18.

- **M10 — pysees handoff and worker execution.** Partially done: the
  Rust-side `CarapaceInputV1` header/table types, decode, and stepped
  `Session` exist and are tested (`wasm-bridge/src/input_v1/`,
  `wasm-bridge/tests/m10_carapace_input_v1.rs`); the `wasm_bindgen`/
  `serde-wasm-bindgen` boundary, the `pysees` compiler, the worker, and the
  results database are not started. `pysees` hands off an immutable,
  validated `Model + AnalysisSequence` snapshot on an explicit Run action —
  not a live model binding and not a replay of arbitrary OpenSees commands.
  Compile it once to a versioned `CarapaceInputV1`: a small structured-clone
  header plus one flat table per (profile, element-kind) pair, mirroring
  `core`'s closed-enum catalog rather than a generic tagged-record array
  (§2.1) — today these are plain Rust `Vec`-backed structs
  (`wasm-bridge/src/input_v1/tables.rs`), not yet the transferable
  `Float64Array`/`Uint8Array` values that will actually cross `postMessage`;
  that binding is deferred until `pysees`'s compiler exists to produce the
  bytes. `carapace-wasm` decodes it once, panic-free (`Result<Session,
  DecodeError>`, §2.8, since the compiler's own validation can't see
  engine/schema version skew), into a `Session` chosen once from the `space`
  discriminant, and steps it via an `advance(stepBudget)` call — the
  mechanism behind cooperative cancellation, exercised today by direct Rust
  calls rather than a worker loop. Both `space` (2/3) and stage kind
  (`static`/`modal`/`transient`) are real discriminants in the wire format
  from the start, even though the first decoder only implements the
  static-planar arm (`space == 3` and `modal`/`transient` stages fail decode
  with a named `DecodeError`, not a panic): `core` already has full spatial
  statics/dynamics/diaphragms (M15-M20) and transient/modal for both
  profiles, so scoping the format itself to the first vertical slice would
  guarantee a breaking rework once those land. The first vertical slice is a
  2D fiber-column gravity → frozen-gravity → displacement-controlled
  pushover, verified against `ForceBeamColumn`'s native elastic/post-yield
  behavior; loads are grouped by which `AnalysisSequence` stage registers
  them on the `Domain` (`tables::NodalLoadTable::stage`), not all applied
  upfront, since a stage's shared pseudo-time would otherwise ramp a later
  stage's reference load too early (the same reason
  `core/tests/m8_force_beam_column.rs`'s native two-phase test defers its
  lateral pattern's load). Selected node response history is currently kept
  in an in-memory `Vec` per recorder (`PlanarSession::recorder_samples`);
  writing it into the worker-owned SQLite-over-OPFS run database in the same
  row layout `carapace-wasm` produces it in (§2.7, extended across the wasm
  boundary — no JS-side reshaping between wasm output and the SQLite blob)
  is still to come. Progress messages, the worker protocol, and run
  provenance persistence are all still to come too. See
  [`pysees-handoff.md`](pysees-handoff.md) for the exact contract, wire
  format, session/decode model, worker protocol, and first acceptance case.

- **M11 — Results database hardening.** Not started. Build on M10's initial
  SQLite-over-OPFS recorder/write path: schema migrations, efficient
  range-query paging and indexes, retention, export, interruption/recovery,
  and large-run storage tests. Extend result channels beyond M10's static
  nodal pushover history (node velocity/acceleration and reactions, element
  and fiber responses, modal and transient channels) only when their core
  response APIs exist and a UI query needs them. The analysis worker remains
  the database's sole owner; the main thread queries it by run ID.

- **M12 — Algorithm richness (line search, initial/secant tangent, Krylov
  acceleration).** Not started. `Algorithm` is currently just
  `Linear`/`NewtonRaphson` (full Newton, current tangent every iteration,
  no globalization) — see [`docs/algorithms.md`](algorithms.md) for the
  full design: a `TangentStrategy` axis (current/initial/reuse-first,
  needing `SparseSolver` to support factor-once-solve-many first), a
  `LineSearch` modifier (Bisection/RegulaFalsi to start), and a
  Krylov-subspace-accelerated modified-Newton variant, each ported from a
  specific `xara/SRC/analysis/algorithm/equiSolnAlgo/` source file rather
  than re-derived. Scoped to static `Analysis` only —
  `TransientAnalysis`'s own already-flagged Newton-corrector gap (§6 M6
  note) is separate follow-on work that can reuse whatever comes out of
  this milestone.

- **M13 — Event-to-event (EtE) stepping.** Not started. A genuinely
  different `Algorithm` variant from M12's Newton family — for a model
  built from piecewise-linear materials, the response between "events"
  (a fiber yielding, a gap closing, a hysteretic control point crossed)
  is exactly linear, so the load-factor distance to the next event can be
  solved for directly (one linear solve, zero iteration, zero
  convergence tolerance) rather than approximated by Newton — the
  numerical basis of PERFORM-3D, which OpenSees/Xara has no equivalent
  of (confirmed: nothing under `xara/SRC/analysis/`), so unlike every
  other milestone here there is no source tree to port from — designed
  instead from the originating paper (Karamchandani & Cornell 1992) and
  CSI's own public terminology for the same technique in SAP2000/ETABS.
  See [`docs/algorithms.md`](algorithms.md)'s Part II for the full
  design: a new `Material::distance_to_event` query (closed-enum match
  arm per variant, §12), its aggregation up through `FiberSection`/
  element/`Domain` (§13, including the harder `ForceBeamColumn` and
  `Series`/`MinMax`-composite sub-problems flagged as open rather than
  hand-waved), event lumping's precise exactness guarantee (§14), and an
  explicit open decision on smooth (non-piecewise-linear) materials
  (`Steel02`/`Concrete01`/`Concrete02`, §13.4) that needs resolving
  before implementation, not during it.

- **M14 — Load patterns, multi-phase analysis composition, ground
  motion.** ✅ Done. Generalizes load application from one implicit
  pattern (`Node::load` + a baked-in element-load field, both scaled by a
  single global `load_factor`) to `Domain`-owned, independently-scaled
  `LoadPattern`s (`core/src/model/load_pattern.rs`), each driven by a
  `LoadSeries` (`Constant`/`Linear`/`Path`, ported from Xara's
  `TimeSeries` family — `ConstantSeries`/`LinearSeries`/`PathSeries`) and
  freezable in place (`Domain::hold_pattern_constant`, Xara's
  `LoadPattern::setLoadConstant`). Multi-phase analyses (gravity →
  pushover, a cyclic quasi-static protocol, gravity → earthquake) are
  caller-side composition of existing types, not a new orchestration
  layer — matches how Xara itself works (Tcl/Python scripts call
  `analyze()` repeatedly and swap integrators; `Domain` carries state
  across the swap): `Analysis::into_domain`/`TransientAnalysis::into_domain`
  hand a finished phase's `Domain` to the next phase's builder, and
  `Analysis::set_integrator` swaps the integrator in place (no rebuild) for
  cheap same-phase-type leg changes, e.g. a cyclic protocol's many
  reversals. A load-controlled load-path protocol (Xara's `Integrator::
  LoadPath`) falls out of `LoadSeries::Path` + the existing `LoadControl`
  integrator — no new `Integrator` variant needed.
  `TransientAnalysis` gained real time-varying load
  (`assemble_reference_load` now takes the current pseudo-time) and
  `GroundMotion` (`core/src/analysis/ground_motion.rs`, Xara's
  `UniformExcitation`): the effective inertial force `-M·ι·ag(t)` via a
  precomputed direction-influence vector (`Domain::direction_incidence`).
  Breaking migration: `Node::load` and `ElasticBeamColumn::w_transverse`
  are gone — loads are added via `Domain::add_nodal_load`/
  `add_element_load` (or the `default_pattern()`-scoped convenience
  `Domain::load_node`) after `add_node`/`add_element`, not as node/element
  builder fields. See `core/tests/m_load_patterns.rs`,
  `m_cyclic_protocol.rs`, `m_ground_motion.rs`, `m_load_path.rs`, and the
  doc comments on `LoadSeries`/`LoadPattern`
  (`core/src/model/load_pattern.rs`), `Domain::hold_pattern_constant`, and
  `GroundMotion` for the design rationale.

- **M15 — Spatial foundation (3D).** Done, except transient/modal.
  `Node<NDIM, NDOF>` (generic, `core/src/model/node.rs`) is specialized as
  `Node3 = Node<3, 6>`. `Domain`, `Analysis`, `AnalysisBuilder`,
  `Integrator`, and `LoadPattern` are each a single generic implementation
  (an `ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>` trait — `core/src/model/
  elements/mod.rs` — abstracts the per-profile element catalog out of
  `Domain`'s bookkeeping), not hand-duplicated planar/spatial pairs;
  `Domain3`/`Analysis3` are type aliases for the spatial instantiation.
  `Truss3` and `ZeroLength3` (`core/src/model/elements/`) are the spatial
  elements: fixed-size `[ux, uy, uz, rx, ry, rz]` layout, direction-cosine
  truss stiffness, independent per-global-direction zero-length materials,
  three-translation lumped mass. A full static solve (`Domain3` →
  `AnalysisBuilder` → `Analysis3::step`) is verified end to end in
  `core/tests/m15_spatial_truss.rs` against a closed-form space truss.
  Element formulas themselves (`Truss` vs `Truss3`) stay separate concrete
  types — see `ElementOps`'s doc comment for why that part isn't unified
  the same way.
  Remaining: spatial modal/transient plumbing (`TransientAnalysis`/
  `modal_analysis` still only accept the planar `Domain`), a real spatial
  `equal_dof`/`rigid_diaphragm` beyond planar identity ties, and the M10
  wasm handoff (still rejects spatial input entirely — unrelated to this
  core-side progress). See [`spatial-architecture.md`](spatial-architecture.md).

- **M16 — Spatial elastic frame and constraints (3D).** Partially done.
  `GeomTransf3` (`core/src/model/transform.rs`) is a first-class, reusable
  3D geometric-transform type — local axes built from a `vec_xz`
  orientation vector, plus the block-diagonal local↔global rotation — the
  spatial counterpart of planar `GeomTransf`, so future spatial frame
  elements share it rather than each inlining their own copy. It has two
  variants: `Linear3`, and `PDelta3` (a linearized geometric-stiffness
  correction from axial force in both bending planes — the `Iy`-plane
  block is a hand-derived sign flip of the `Iz`-plane one, verified by
  exercising both planes in `core/tests/m16_pdelta3.rs`, not just asserted
  by analogy). `ElasticBeamColumn3`
  (`core/src/model/elements/elastic_beam_column.rs`), which holds one, is a
  12×12 Euler-Bernoulli frame (axial, torsion, biaxial bending via `E, G,
  A, J, Iy, Iz`). Verified against Xara's `ElasticBeam3d`/
  `LinearCrdTransf3d` and closed-form cantilever deflections, at the
  transform level (`transform.rs`'s own tests), element level
  (`spatial_tests` in `elastic_beam_column.rs`), and through the full
  `Domain3`/`Analysis3` stack (`core/tests/m16_elastic_beam3.rs`,
  `core/tests/m16_pdelta3.rs`).
  Remaining: local-y/local-z beam loads and a real constraint
  transformation for spatial rigid diaphragms (identity DOF aliasing cannot
  express lever-arm kinematics). See [`spatial-architecture.md`](spatial-architecture.md).

- **M17 — Spatial fiber/nonlinear frame elements (3D).** Done.
  `FiberSection3` (`core/src/model/fiber_section.rs`) has local y/z fiber
  locations and coupled `[N, Mz, My]` response with a 3×3 tangent
  (product-of-inertia coupling included); torsion is explicit `G`/`J` on
  the element (decoupled elastic `G*J/L`), not a fiber-derived quantity.
  `DispBeamColumn3` and `ForceBeamColumn3`
  (`core/src/model/elements/disp_beam_column.rs`,
  `.../force_beam_column.rs`) are independently verified formulations, not
  3D constructor flags on their planar counterparts —
  `ForceBeamColumn3`'s basic system generalizes `ForceBeamColumn`'s to five
  dof (axial plus chord-relative end rotations in both bending planes),
  with the y-bending rows' rotation-coefficient sign flipped relative to
  z-bending's (`ry = -dw/dx` vs `rz = +dv/dx`, the same substitution used
  throughout the spatial bending code). Verified against
  `ElasticBeamColumn3`'s exact closed-form biaxial elastic stiffness
  (`core/tests/m17_disp_beam_column3.rs`,
  `core/tests/m17_force_beam_column3.rs`) and, for `ForceBeamColumn3`, a
  doubly-symmetric `ElasticPP` section yielding independently in both
  bending planes with permanent set on unload, plus (closing this
  milestone's last open item) an asymmetric, non-doubly-symmetric section
  with nonzero product of inertia against a hand-derived closed form
  (`asymmetric_biaxial_section_matches_hand_derived_closed_form`), the
  spatial analogue of `m8_force_beam_column.rs`'s
  `asymmetric_elastic_section_matches_hand_derived_closed_form`. See
  [`spatial-architecture.md`](spatial-architecture.md).

- **M18 — Spatial corotational geometry (3D), if confirmed in scope.** Not
  started. Revisit after M15–M17; spatial and planar corotational transforms
  need different rotation handling and independent objectivity tests. See
  [`spatial-architecture.md`](spatial-architecture.md).

- **M19 — Spatial dynamics: modal/transient analysis (3D).** Done.
  `modal_analysis` and `TransientAnalysis`
  (`core/src/analysis/modal.rs`, `core/src/analysis/transient.rs`)
  generalized to the same `<NDIM, NDOF, ELEMENT_DOF, NId, E>` profile as
  `Domain`/`Analysis`, with `TransientAnalysis3` as the spatial type alias
  (`modal_analysis` infers its profile from the `domain` argument, needing
  no separate spatial entry point). No dynamics formula changed — Lanczos,
  Newmark-beta, the initial-acceleration solve, and `GroundMotion`'s
  effective inertial force already went through `Domain`'s generic free-DOF
  interface. Verified through the full `Domain3`/`Analysis3` stack
  (`core/tests/m19_spatial_dynamics.rs`): a spatial `Truss3` SDOF natural
  frequency, undamped/Rayleigh-damped spatial `ZeroLength3` SDOF free
  vibration on two different translational DOFs, and two simultaneous,
  independent `GroundMotion`s (`ux`/`uz`) driving an uncoupled spatial mass
  with no cross-talk — the 3D ground-motion case. See
  [`spatial-architecture.md`](spatial-architecture.md).

- **M20 — Spatial beam loads and true rigid-diaphragm constraints.** Done.
  `ElementLoad3::UniformTransverse { wy, wz }`
  (`core/src/model/load_pattern.rs`) is `ElasticBeamColumn3`'s biaxial
  local transverse load, `ElementOps::Load` for the spatial profile in
  place of the uninhabited `Infallible` (scope matches the planar profile
  exactly — only `ElasticBeamColumn3`, not `DispBeamColumn3`/
  `ForceBeamColumn3`, same as their planar counterparts). Verified against
  Xara/OpenSees's `ElasticBeam3d::addLoad` fixed-end-force derivation and
  through the full `Domain3`/`Analysis3` stack
  (`core/tests/m20_spatial_beam_loads.rs`).
  `Domain3::rigid_diaphragm`/`rigid_diaphragm_about` replace identity
  aliasing with a genuine `AffineConstraint` transformation
  (`core/src/model/domain.rs`): each constrained node's two in-plane
  translations (the axes other than an explicit `normal: Axis3`, default
  `Y`) tie to the retained node's same-axis translations plus a lever-arm
  term from its rotation about `normal` — every rotational dof and the
  out-of-plane translation stay completely free, narrower in scope than
  Xara/OpenSees's `rigidDiaphragm` by design. Resolved once at
  `Domain::number_dofs` time into a `dof_transform` map that every
  assembly hot loop (stiffness, mass, load) and state-scatter path reads
  through one shared `dof_terms` helper, rather than three separate code
  paths for free/identity-tied/affine-tied dofs. `assemble_mass_diagonal`
  panics rather than silently dropping physics if nonzero mass ever lands
  on an affine-tied dof (the lumped-mass diagonal can't represent the
  resulting rotational inertia at the retained node). Verified through the
  full `Domain3`/`Analysis3` stack (`core/tests/m20_rigid_diaphragm3.rs`):
  the classic asymmetric-diaphragm torsion problem against an
  independently-derived closed-form stiffness matrix, rigid-rotation
  invariance of every untied dof, and an explicit non-default `normal`
  cross-checked against Xara/OpenSees's own documented `RigidDiaphragm`
  equations. See [`spatial-architecture.md`](spatial-architecture.md).

---

## 7. Verification methodology

Every milestone with numerical output gets checked against a known-correct
value — either closed-form (where one exists) or against **Xara's own native
build as an oracle** for cases without a simple closed form (see
`docs/xara-feasibility.md`). When in doubt about expected behavior for a
formulation, run and diff against the Xara spike rather than relying on
documentation alone.

**Verification order for every milestone:** native `cargo test` first, then
`wasm32-unknown-unknown` build + `wasm-bindgen` + Node execution. Don't skip
the native step to save time — it isolates logic bugs from
platform/toolchain bugs (see `docs/xara-feasibility.md`'s `s_copy` note for
why this discipline exists).

---

## 8. Reference material location

The Xara C++ spike and the standalone ARPACK f2c feasibility spike that
grounded this plan's early decisions are preserved (source-level artifacts
and dependency-closure lists, not build outputs) at `reference/xara-spike/`
— see that directory's README for what's there, and
`docs/xara-feasibility.md` for the design conclusions drawn from it.

---

## 9. Known follow-ups / performance backlog

Not scheduled against any milestone — things noticed along the way that are
deliberately *not* being acted on without a profile or a concrete need,
recorded here so they don't have to be rediscovered from scratch. Don't
implement anything in this section speculatively; it exists so a future
session (or this one, later) knows the question was already thought through
once.

### 9.1 `trial_stress_tangent` unnecessarily allocates for boxed materials

**The situation:** `Material`'s large leaf variants (`Hysteretic`,
`Pinching4`) are boxed specifically so one huge variant doesn't force every
`Material` value in the system — including a trivial `Elastic { e: f64 }` —
up to its size (closed-enum dispatch means `size_of::<Material>()` is always
the size of the largest variant). That boxing has a real, separate cost:
`Material::evaluate` (the shared core behind both `trial_stress_tangent` and
`commit`) always constructs the full next-state `Material`, box included,
even when the caller is `trial_stress_tangent` and immediately discards it.
Since `trial_stress_tangent` runs once per Newton iteration (far more often
than `commit`), every iteration touching a boxed material allocates and
immediately frees a small heap block for no reason.

**Why this isn't being fixed reflexively:**

- *Scope*: the easy part — giving a boxed leaf material's own `evaluate_xxx`
  a variant that returns the raw, unboxed fields struct instead of a boxed
  `Material`, letting `trial_stress_tangent`/`commit` each decide separately
  whether to box the result — is clean and wouldn't change the trial/commit
  purity design at all. The hard part: `Parallel`/`Series`/`MinMax`
  recursively call `child.evaluate(strain)` through a single generic
  dispatcher that has no way to know whether the *outer* caller wanted a
  trial or a commit, so a boxed material inside a composite would still box
  on every trial unless composites also got a trial/commit-aware recursive
  path — a materially bigger change than the leaf-material case alone.
- *Impact*: unmeasured. Plausibly noise next to the cost already paid per
  Newton iteration (global tangent assembly, sparse LU solve) unless a model
  has enough boxed-material fibers, in enough elements, over enough
  iterations, for it to add up. Nobody has profiled this against a real
  model yet.

**If this becomes worth doing:** benchmark first (allocation count and/or
wall time per step, before/after, on a model with many boxed-material
fibers under `Algorithm::NewtonRaphson`) to confirm it's actually
load-bearing before spending the implementation effort, and scope the first
pass to the leaf-material case only.
