# Carapace Implementation Plan

This document is the source of truth for Carapace's design and roadmap. It exists
so a coding session with no prior context can pick up the project and continue
correctly. If something here conflicts with the code, the code wins for *what
exists*, but this document should still be updated to match — don't let it rot.

Cross-references to `xara/...` refer to the OpenSees/Xara C++ repository at
`/home/jchat/projects/xara`, which Carapace uses as a design reference and a
numerical verification oracle. Carapace does not depend on Xara at build time.

---

## 1. Goals and non-goals

**Goal:** a small, fast, WebAssembly-native FE analysis engine covering the
subset of OpenSees/Xara functionality that `pysees` (the companion frontend)
actually uses, running inside a browser Web Worker.

**Explicit non-goals** (each was a real decision made during design, not an
oversight — don't "fix" these without revisiting the reasoning):

- **Not full OpenSees/Xara parity.** No PFEM, no fluids, no thermal, no
  reliability analysis, no arbitrary user-supplied element plugins. See §4 for
  exactly what *is* in scope.
- **Not dynamically extensible at runtime.** The element/material catalog is
  fixed at compile time (closed `enum`s, not a broker/factory). Adding a new
  element type is a recompile, not a runtime registration. This is a deliberate
  simplification enabled by "thin, curated" — see §3.1.
- **Not multi-tab / multi-writer.** A single worker owns the model and any
  persisted results (SQLite-over-OPFS). Multiple tabs querying the same
  analysis concurrently is explicitly out of scope (decided when designing the
  `pysees` results-storage architecture).
- **Not a live/interactive re-solve loop.** Workflow is "build model once,
  then analyze," not "drag a node, re-solve every frame." This shapes the
  worker/main-thread messaging design in §5.
- **Not parallel (MPI) or database-backed.** No `sendSelf`/`recvSelf`,
  no `Channel`, no `ObjectBroker`. Single-threaded, in-memory, one process.

---

## 2. Why this exists / what we learned porting Xara directly

A parallel investigation (see `xara/` conversation history, and the spike
artifacts referenced in §7) attempted the more obvious path: compile Xara's
actual C++ to WebAssembly via `emcc`, translating its Fortran BLAS/LAPACK/ARPACK
dependencies with `f2c`. **This worked** — a real, running proof of concept
(`Domain`/`Node`/`Truss`/`ElasticMaterial`/`BandGenLinSOE`, static analysis,
correct displacement output) is proof the numerics and toolchain are sound.

But it surfaced concrete architectural problems worth designing around rather
than inheriting:

1. **`XaraClassBroker` eager construction.** `xara/SRC/runtime/runtime/
   ModelRegistry.h` holds a `ProcessContext` *by value*, whose constructor
   unconditionally builds an `XaraClassBroker` (`xara/SRC/runtime/runtime/
   XaraClassBroker.cpp`, 230 `#include`s — Shell elements, MVLEM, Fortran
   `feap` materials). Since `BasicAnalysisBuilder` depends on `ModelRegistry`,
   *any* use of Xara's own analysis-composition class drags in nearly the
   whole element/material catalog at compile time, regardless of what the
   model actually uses. → Carapace uses closed `enum` dispatch instead of a
   factory/broker (§3.1).

2. **Ownership by convention, not enforcement.** `DOF_Numberer` assumes
   ownership of its `GraphNumberer&` and deletes it in its destructor — this
   is a comment-level convention, not compiler-checked. Passing a
   stack-allocated `RCM` where heap ownership was expected caused a real
   double-free in the spike, caught only by AddressSanitizer. → Carapace
   relies on Rust's ownership model to make this class of bug a compile error
   (§3.2).

3. **Tcl coupling leaks through headers that don't need it.**
   `xara/SRC/logging/logging.cpp` unconditionally `#include <tcl.h>` despite
   using zero `Tcl_*` symbols — purely incidental coupling. → Carapace has no
   interpreter dependency anywhere; there is nothing to leak.

4. **A pre-existing, unrelated bug was found along the way:**
   `xara/SRC/matrix/routines/cmx.h` declares `cmx_inv2..cmx_inv6` with no
   implementation anywhere in the checkout (git history shows the header
   deleted in one commit and reintroduced later without its source). Not a
   wasm-specific issue — worth a real upstream fix in Xara, unrelated to
   Carapace's own work.

5. **Force-based elements do nested equilibrium iteration.**
   `xara/SRC/element/Frame/Force/ForceFrame3d.h` carries its own `tol`/
   `max_iter` for internal "local iterations" — a force-based beam-column
   element solves its own internal equilibrium (section forces consistent
   with element-end forces) *inside* the domain-level Newton loop. This is a
   genuinely different shape from every other element and is its own
   milestone (M8), not a variant of the displacement-based element work.

6. **Some materials are recursive compositions, not leaves.**
   `xara/SRC/material/wrapper/{ParallelMaterial,SeriesMaterial,MinMaxMaterial}.h`
   hold pointers to *other* `UniaxialMaterial` instances. A flat material enum
   doesn't work for these — see §3.1 and §4.2.

---

## 3. Architecture principles

### 3.1 Dispatch: closed `enum`s for catalogs, not broker/factory patterns

Elements and materials are a **fixed, compile-time-known set** (§4). Represent
each as a Rust `enum`:

```rust
enum Element { Truss(Truss), ZeroLength(ZeroLength), ElasticBeamColumn(..), .. }
enum Material { Elastic(Elastic), Steel01(Steel01), Parallel(Vec<Box<Material>>), .. }
```

`match`-based dispatch, no `Box<dyn Trait>`, no runtime registration. This is
the direct fix for problem #1 above: there is no broker to eagerly construct,
because there is no broker.

**Exception — analysis strategy types are different and get enums too, but for
a different reason.** `ConvergenceTest`, `Integrator`, `Algorithm`,
`ConstraintHandler` are also closed enums, but not because a broker would be
"toxic" for them (their interfaces are narrow and homogeneous, unlike the
heavy, heterogeneous `Element`/`Material` catalog) — it's simply that the set
is small, known, and user-selected per analysis, so an enum is the simplest
correct representation. `SparseSolver` (§3.6) is not an enum or trait at all —
just a concrete struct, since Carapace commits to exactly one solver
implementation.

### 3.2 Ownership via the type system

No raw pointers with ownership documented only in a comment. Model state lives
in an arena/`SlotMap`-style store keyed by generational indices, not `Box`/`Rc`
graphs mirroring the C++ pointer structure. Where a type must own another
(e.g. an `Analysis` owning its `Integrator`), plain Rust ownership (`Integrator`
by value, not `&mut Integrator` with an ownership convention) makes the
question unrepresentable rather than a runtime footgun — directly addressing
problem #2 above.

### 3.3 No serialization / broker / interpreter machinery

No `sendSelf`/`recvSelf`, no `Channel`, no `ObjectBroker`, no embedded
interpreter of any kind. These exist in Xara for MPI domain decomposition and
database persistence — irrelevant to a single-threaded wasm worker, and a
major source of the "everything is transitively reachable from everything"
coupling that made the Xara C++ dependency closure balloon (348 files for a
one-element static analysis, before hand-wiring around it).

### 3.4 Fixed-size, stack-allocated element-local linear algebra

A 2-node truss's stiffness matrix is always 4×4 (2D) or 6×6 (3D) — known at
compile time. Use const-generic fixed-size types for element-local
stiffness/state (e.g. `nalgebra::SMatrix<f64, 4, 4>`, or hand-rolled — see §6
open decision) rather than a generic runtime-sized heap-allocated `Matrix`
(Xara's `SRC/matrix/Matrix.cpp`, which also carries a manual `fromFree`
ownership flag as a workaround for not having Rust-style ownership). Zero heap
allocation in the per-element evaluation hot loop.

### 3.5 Typestate-enforced analysis composition

Xara's `BasicAnalysisBuilder` requires calling `setLinks` on five-plus objects
in a specific order before use — undocumented in the type system, verified
only by reading `xara/SRC/runtime/runtime/BasicAnalysisBuilder.cpp` line by
line (which is literally what happened when hand-wiring the spike's driver).
Carapace uses the typestate pattern: an `AnalysisBuilder<Unwired>` can only
become `AnalysisBuilder<Ready>` once every required piece is set, and only
`Ready` exposes `.build()` → `Analysis`, which exposes `.step()`. Illegal
sequencing is a compile error. See §5 for the concrete sketch.

### 3.6 One fixed sparse solver, no solver abstraction layer

Xara has `LinearSOE`/`LinearSOESolver`/`DOF_Numberer` as three separate
interchangeable-backend abstractions (to support banded/sparse/parallel/PETSc
backends). Carapace commits to **one** sparse solver implementation — no
trait, no enum, just a concrete `SparseSolver` struct. This also likely
removes the need for a separate DOF-numbering abstraction: modern sparse LU
libraries typically do their own fill-reducing ordering internally (verify
this against whichever solver crate is chosen — see §6).

### 3.7 Data-oriented (SoA) result storage

Element response history (strain/stress/tangent per step) is stored in flat,
contiguous, per-response-type arrays, not one heap object per element per
timestep. This shape matches what the results pipeline needs anyway (typed
arrays over `postMessage`, rows into SQLite — §5.3), so the internal and
external representations should be the same shape from the start rather than
converted at the boundary.

### 3.8 Errors are `Result<T, E>`, not sentinel integers

Xara's analysis methods return `-1`..`-5` with meaning defined only in
scattered comments (had to be replicated by hand in the spike driver).
Carapace uses a real error enum (`AnalysisError::FailedToConverge{ step }`,
`::SingularSystem`, `::InvalidConstraint`, ...) carrying context.

---

## 4. Feature scope

Scope is exactly what `pysees` uses today, mapped from a direct conversation
with the project owner. Anything not listed here is out of scope until a
concrete need arises — don't speculatively add elements/materials "while
you're in there."

### 4.1 Elements — rough milestone order, see §7 for details

| Element | Complexity | Notes |
|---|---|---|
| Truss | Trivial | Done in the Xara C++ spike; M1 target for Carapace. |
| ZeroLength | Simple | No geometry/integration; direct material evaluation per DOF. |
| ElasticBeamColumn | Simple | Closed-form stiffness, no iteration. |
| DispBeamColumn | Moderate | Needs `BeamIntegration` (Gauss-Legendre/Lobatto) + fiber sections. One Newton loop total (domain-level). |
| ForceBeamColumn | **High — own milestone (M8)** | Nested element-level equilibrium iteration (`xara/SRC/element/Frame/Force/ForceFrame3d.h` has its own `tol`/`max_iter`), distinct architecture from every other element. |

Fiber sections: mechanically simple (loop over fibers, accumulate
stress-resultant + tangent from each fiber's uniaxial material response,
weighted by area/position) but evaluated many times per element per
iteration — keep allocation-free (§3.4/§3.7).

### 4.2 Materials (uniaxial)

**Leaf materials** (flat enum variants, given strain → stress + tangent from
internal state): Elastic, ElasticPP (EPP), Gap, ENT ("no tension"),
Hysteretic, Pinching4, Concrete01, Concrete02, Steel01, Steel02.

**Composite/wrapper materials** (recursive — hold other materials, see
`xara/SRC/material/wrapper/{ParallelMaterial,SeriesMaterial,MinMaxMaterial}.h`):
Parallel, Series, MinMax. These require `Material::Parallel(Vec<Box<Material>>)`
/ `Material::Series(Vec<Box<Material>>)` / `Material::MinMax(Box<Material>, ..)`
— design the `Material` enum and its evaluation trait/interface to support
recursive calls from the start (M7), don't bolt it on after leaf materials are
done.

### 4.3 Analyses

- **Static / quasi-static** — same code path (`LoadControl` or
  `DisplacementControl` integrator with small increments); not a distinct
  architecture.
- **Modal** — eigenvalue analysis. See §6 open decision (FFI-bind the
  f2c'd ARPACK from the Xara spike, vs. pure-Rust Lanczos/Arnoldi).
- **Time history (transient)** — `Newmark` integration, element/nodal mass
  matrices, Rayleigh damping (linear combination of mass/stiffness).

### 4.4 Supporting infrastructure

- **Mass:** lumped (diagonal) mass matrices per element, or nodal mass —
  needed for modal and time-history analyses.
- **Loads:** nodal loads (done in M1) and element loads (distributed loads on
  beam-columns — needed once `ElasticBeamColumn`/`DispBeamColumn` land).
- **`geomTransf`:** `Linear` and `PDelta` are cheap (small correction terms).
  **Corotational is a separate, high-complexity milestone (M9)** if/when
  large-displacement analysis is actually needed — it requires element-local
  frame tracking and large-rotation updates, comparable in complexity to
  force-based elements. Confirm this is actually in scope before starting it;
  it was flagged as "decide if you need it" during scoping, not committed.
- **Rayleigh damping:** needed for time-history analysis (M6).

---

## 5. System architecture (worker / JS integration)

Decided in earlier design conversation, restated here for continuity:

### 5.1 Model description: batch, not per-call RPC

The frontend (`pysees`) builds the full model description as a plain
TypeScript object/typed-array bundle on the main thread (cheap, synchronous,
no wasm involved), then transfers it **once** into the worker (structured
clone for metadata, transferable `ArrayBuffer`s for numeric-heavy data —
node coordinates, element connectivity). The worker walks it once inside
wasm. This avoids the alternative — mirroring a Python-style
`model.node(...)`/`model.element(...)` object-method API 1:1 across the
worker boundary — which would mean one `postMessage` round-trip per model
entity; fine in-process for Python/pybind11, expensive across a worker
boundary for models with thousands of entities.

### 5.2 Command surface: coarse, not per-node

After the initial batch transfer: `runAnalysis(recipe)`, `getResults(query)`,
not per-entity calls. Matches §1's "no live re-solve loop" non-goal.

### 5.3 Results: typed arrays + SQLite-over-OPFS, decoupled from live progress

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

### 5.4 Sketch of the typestate analysis API (from design discussion)

```rust
struct AnalysisBuilder<S> { /* .. */ }
impl AnalysisBuilder<Unwired> {
    fn constraint_handler(self, h: ConstraintHandler) -> AnalysisBuilder<HasHandler> { .. }
}
// ...chained until fully wired...
impl AnalysisBuilder<Ready> {
    fn build(self, model: &Model) -> Analysis {
        // numbers DOFs, builds sparsity pattern, factors once if linear —
        // done ONCE here, not re-checked every step (see §1 non-goals:
        // no live re-solve means no need for Xara's mid-loop
        // hasDomainChanged()/re-number/re-factor machinery).
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

---

## 6. Open decisions — resolve early, before M2

These were identified but deliberately not settled during design, because
they're better answered empirically than by discussion:

1. **Linear algebra / sparse solver crate: `nalgebra` vs `faer`. Resolved
   for real, post-M5.** `SparseSolver` (`core/src/analysis/solver.rs`) uses
   `faer`'s sparse LU (`SparseColMat::sp_lu()`), not a dense `nalgebra`
   solve — the M1-M5 dense placeholder was deliberately deferred (see the
   superseded note this replaced) until it actually mattered, and "real
   problems are not small" was the trigger to stop deferring rather than
   wait for a specific milestone to force the issue. `faer` over
   `nalgebra-sparse`: a built-in sparse LU with COLAMD/AMD fill-reducing
   ordering (answers the original "does the solver do its own fill-reducing
   ordering" question — yes, so no separate `DOF_Numberer` abstraction is
   needed, confirming §3.6's original bet), pure Rust (no C dependency to
   fight through `wasm32-unknown-unknown`), and confirmed by hand to build
   and solve correctly under `wasm32-unknown-unknown` + Node with a minimal
   feature set (`default-features = false, features = ["std",
   "sparse-linalg"]` — the defaults pull in `rand`/`rayon`/`npy`, and
   `getrandom` fails the wasm32 build without the `wasm_js` backend those
   defaults don't request). `Domain` assembles stiffness directly into
   sparse triplets (`core/src/model/domain.rs`'s
   `assemble_stiffness_triplets`) rather than a dense buffer that's
   sparsified after — the whole point was avoiding the O(n²) memory/compute
   a dense assembly-then-convert would still cost on a real model.
   `core/tests/sparse_solver.rs` exercises a 150-free-DOF chained-beam
   system (exact against the closed-form cantilever tip deflection) as a
   scale sanity check beyond the single- and few-DOF systems every other
   milestone's tests use. Refactoring the symbolic factorization across
   Newton iterations (same sparsity pattern, only values change) is a real
   future optimization this pass didn't need — `SparseSolver::solve`
   currently re-factors from scratch every call.

2. **ARPACK strategy for modal analysis (M5). Resolved for real, post-M5:
   a hand-rolled shift-invert Lanczos — not because "no mature crate
   exists" (the original framing), but because once `SparseSolver` went
   sparse (decision #1), a dense full-spectrum `SymmetricEigen` on the
   whole stiffness matrix turned out to be wrong beyond just slow: it
   still cost O(n²) memory regardless of how sparse the model actually
   was (could exhaust a wasm worker's memory well before compute time
   mattered), and it always computed *every* eigenpair when the standard
   use — including the modal-superposition damping M6+ will want — only
   ever needs the lowest `num_modes`.** `core/src/analysis/modal.rs`'s
   `modal_analysis` now takes a `num_modes` argument and runs a
   shift-invert (shift = 0, i.e. operator `K^-1 * M`) Lanczos: each
   iteration's shift-invert solve reuses `SparseSolver`'s sparse LU
   directly, and the only dense step is `nalgebra::SymmetricEigen` on the
   small `m×m` *projected* tridiagonal matrix (`m` = Lanczos subspace
   size, always small — `min(2*num_modes+8, n)` — never `n` itself). This
   is not "rolling our own ARPACK": the two numerically hard pieces
   (sparse LU, small dense eigendecomposition) are both library code; what's
   ours is the well-understood outer three-term-recurrence loop with full
   M-orthogonal re-orthogonalization (cheap, since the subspace is small).
   Only `shift = 0` is implemented (lowest frequencies — the standard
   structural-dynamics case); a nonzero shift (targeting a frequency band,
   or handling near-singular `K` for buckling) is a natural extension, not
   built until something needs it. `Domain::assemble_mass_diagonal` also
   replaced the old `assemble_mass` (which built a full dense N×N matrix
   for what's always diagonal lumped mass — the same class of waste as a
   dense stiffness matrix). **Verified:** `core/tests/m5_modal.rs`'s
   `requesting_fewer_modes_than_free_dofs_returns_the_lowest_ones` checks a
   partial-spectrum request (2 of 5 modes) against a full-spectrum solve of
   the same system; a 500-free-DOF chain solved for its lowest 5 modes in
   ~2ms in manual testing (not a committed benchmark, just a scale sanity
   check). The ARPACK spike (§9) remains available reference material for
   an eventual real ARPACK comparison, not wasted work.

3. **`cmx.h`-equivalent small-matrix inversion.** Not actually a Carapace
   decision — `nalgebra`/`faer` both handle small fixed-size matrix inversion
   natively; this is only listed to note that Xara's missing-implementation
   bug (§2, problem #4) has no equivalent risk here.

---

## 7. Milestones

Each milestone: scope, acceptance criteria, verification order (native
target first, always — this caught a real wasm-only bug, an f2c/wasm-ld
function-signature mismatch, during the Xara spike; the Rust equivalent might
be a numeric or codegen difference between targets, not a memory-safety bug,
but the discipline of isolating "is my logic wrong" from "is something
wasm-specific wrong" still applies).

- **M0 — Workspace scaffold.** ✅ Done. `carapace-core` (no wasm deps) +
  `carapace-wasm` (thin `wasm-bindgen` layer) workspace; native build/test
  and wasm32 build/JS-glue/Node-execution all verified against a known
  closed-form answer (`axial_displacement`, matches the Xara C++ spike's
  oracle value: L=100, A=2, E=30000, P=50 → dx=0.0833333333).

- **M1 — Real `Domain`/`Node`/`Truss`/`Elastic` + one real solve.** ✅ Done.
  `carapace-core::model` (`Domain`/`Node`/`Element::Truss`/`Material::Elastic`,
  generational `NodeId`/`ElementId` via `slotmap`) and `carapace-core::analysis`
  (typestate `AnalysisBuilder<Unwired..Ready>` per §5.4, `Integrator::LoadControl`,
  `Algorithm::Linear`, `ConvergenceTest::NormUnbalance`, `ConstraintHandler::Plain`,
  `SparseSolver` — a dense `nalgebra` LU at the time, later replaced by
  `faer`'s sparse LU, see §6 decision #1) replace the closed-form
  placeholder. **Acceptance verified:** same
  2-node truss case as M0 (L=100, A=2, E=30000, P=50 → dx=0.0833333333),
  computed through the real architecture — `core/tests/m1_truss.rs` (native)
  and `wasm-bridge`'s `axial_displacement_via_analysis` (wasm32 + Node, via
  `wasm-bindgen`) both match to 1e-9. The M0 closed-form path is kept
  alongside as a fixed oracle constant, not removed.

- **M2 — ZeroLength + EPP/Gap/ENT materials.** ✅ Done. `Element::ZeroLength`
  (`core/src/model/element.rs`) evaluates a material independently per DOF
  direction (global axes only, no orientation vectors — out of scope until
  needed), no geometry/integration. `Material` gained `ElasticPP`, `Gap`,
  `Ent` variants (`core/src/model/material.rs`), all still stateless
  (pure functions of current strain) — implemented as **reversible**
  bilinear/gap envelopes, not true path-dependent plasticity with permanent
  set. This is an intentional simplification, not the real OpenSees
  `ElasticPPMaterial`/`ElasticPPGap` semantics (which track plastic strain
  and unload elastically from wherever they last yielded); it's exact for
  monotonic loading that stays in one regime, and diverges once loading
  reverses direction or crosses a yield/gap boundary mid-step. Real
  path-dependent state (trial vs. committed, `&mut self`) is deferred
  alongside the M7 stateful materials (Steel01, Concrete01, ...), once M4
  brings `Algorithm::NewtonRaphson` and there's an iteration loop that
  actually needs to distinguish trial from committed state. **Acceptance
  verified:** `core/src/model/material.rs` unit tests exercise each
  variant's regimes directly; `core/src/model/element.rs` unit tests drive
  `ZeroLength` across regimes via manually-set node displacements
  (bypassing `Analysis`, since resolving equilibrium across a material's
  nonlinear regimes needs the Newton iteration M4 adds); `core/tests/
  m2_zero_length.rs` + `wasm-bridge`'s `zero_length_ent_displacement` run a
  `ZeroLength`+`Ent` connector through the real `Domain`/`Analysis` pipeline
  end-to-end (native and wasm32+Node), matching hand calc for a load that
  stays within `Ent`'s single linear (engaged) regime — the same
  `Algorithm::Linear`-is-exact condition M1's truss case relied on.

- **M3 — ElasticBeamColumn + geomTransf (Linear, PDelta) + element loads.**
  ✅ Done. `Element::ElasticBeamColumn` (`core/src/model/beam.rs`): a
  prismatic 2D Euler-Bernoulli beam-column with closed-form local stiffness
  (no `Material` dispatch — unlike `Truss`/`ZeroLength`, its response is
  fully determined by `e`/`a`/`iz`, no nonlinear stress-strain law), a
  `GeomTransf` (`Linear` / `PDelta`) resolving local↔global via the
  standard 2D block-rotation transform, and a uniform transverse
  `w_transverse` element load converted to consistent (virtual-work)
  equivalent nodal loads. **This bumped `NDF` from 2 to 3** (added in-plane
  rotation) — `Truss`/`ZeroLength` only populate translational DOF entries
  of the now-6×6 element matrices, so any node connected only to those
  needs its rotation DOF fixed explicitly (M1/M2's tests and wasm exports
  were updated accordingly).

  `GeomTransf::PDelta` adds a linearized geometric-stiffness correction
  from the current axial force to the *tangent* only — the resisting-force
  recovery stays purely elastic. This is a **first-order, non-path-
  following** approximation: within a single `Algorithm::Linear` solve from
  a zero initial state the axial force it reacts to is read from the
  pre-step (usually zero) displacement, so the correction only does
  anything across multiple `LoadControl` increments (see
  `core/tests/m3_beam.rs`'s two-increment test). Full path-consistent
  P-Delta needs M4's Newton iteration; revisit then.

  **Acceptance verified:** `core/tests/m3_beam.rs`
  (native) + `wasm-bridge`'s `simply_supported_beam_end_rotation` (wasm32 +
  Node) — a single-element simply-supported beam under a uniform transverse
  load matches the closed-form end rotation `theta = w*L^3/(24*E*I)`
  *exactly* (not approximately: a single 2-node beam element with
  consistent equivalent nodal loads reproduces the continuous beam's nodal
  DOFs exactly whenever there's no interior point load — the same identity
  behind the classical fixed-end-moment method). A second native test
  confirms `PDelta` matches `Linear` exactly at zero axial force and
  softens a cantilever's tip deflection under compressive axial force, in
  the intended direction.

- **M4 — Analysis composition generalization.** ✅ Done.
  `Integrator::DisplacementControl` (`core/src/analysis/integrator.rs`) picks
  each step's load factor via the standard unit-load technique (solve the
  current tangent against the reference load pattern once, scale so the
  controlled node/DOF's response matches the target increment) rather than
  a fixed increment. `Algorithm::NewtonRaphson` re-forms the tangent and
  residual every iteration at a load factor held fixed for the step,
  reusing `ConvergenceTest` (now also `NormDispIncr`, `EnergyIncr` alongside
  `NormUnbalance`) to decide when to stop; `Algorithm::Linear` is
  unchanged, still a single unconditional solve. `Integrator` and
  `Algorithm` were factored into a clean predictor/corrector split
  (`Analysis::step`: integrator predicts the step's load factor once, then
  the algorithm iterates — or doesn't — at that fixed factor), which is why
  this needed touching both rather than just adding enum variants in
  isolation. Added `AnalysisError::InvalidConstraint` for a
  `DisplacementControl` target DOF that turns out to be fixed.

  **Acceptance verified:** `core/tests/m4_analysis.rs` (native) +
  `wasm-bridge`'s `newton_raphson_elastic_plastic_displacement` (wasm32 +
  Node) — a `Truss` (elastic) in parallel with a `ZeroLength`+`ElasticPP`
  spring, loaded past the EPP spring's yield point *within a single step*
  (something `Algorithm::Linear` structurally cannot resolve, since its
  one-shot solve is only exact within a single material regime — see M1-M3
  notes above), matches the closed-form post-yield displacement
  `u = (F - fy) / k_truss` under all three `ConvergenceTest` variants. A
  separate test confirms `DisplacementControl` on M1's truss case recovers
  the same P=50 load factor M1 verified directly, when given the equivalent
  target displacement instead.

- **M5 — Modal analysis.** ✅ Done. `Node` gained a `mass` field (lumped,
  per-DOF — element-consistent mass matrices are still M6's job) and
  `Domain::assemble_mass_diagonal` builds the free-DOF mass diagonal from
  it (a `DVector`, not a dense or sparse N×N matrix — `M` is diagonal by
  construction, so anything more is pure waste). `analysis::modal_analysis`
  (`core/src/analysis/modal.rs`) solves the generalized eigenproblem
  `K*phi = omega^2*M*phi` at the domain's current state, requesting the
  lowest `num_modes` via a shift-invert Lanczos — see §6 decision #2 for
  the full resolution (a dense full-spectrum solve was the initial version
  of this milestone, later replaced once it turned out to be wrong beyond
  just slow — real projects should read decision #2's writeup, not just
  this summary, before assuming "dense" anywhere in this codebase). A free
  DOF with zero mass, or a `num_modes` outside `[1, num_free_dofs]`, is
  reported as `AnalysisError::SingularSystem` /
  `AnalysisError::InvalidModeCount` respectively — not silently wrong.
  Frequencies are returned sorted ascending.

  **Acceptance verified:** `core/tests/m5_modal.rs` (native) +
  `wasm-bridge`'s `mass_spring_chain_frequencies` (wasm32 + Node) — the
  classic 2-DOF "1-1-1-1" mass-spring chain (unit masses, unit spring
  stiffnesses, `ZeroLength`+`Elastic` in series) matches its exact
  closed-form natural frequencies `1/phi` and `phi` (golden ratio, from the
  characteristic polynomial `lambda^2 - 3*lambda + 1 = 0`) to 1e-9. A
  second test confirms the massless-free-DOF case is reported as an error
  rather than silently dividing by zero. A third
  (`requesting_fewer_modes_than_free_dofs_returns_the_lowest_ones`) checks
  the Lanczos partial-spectrum request itself: asking a 5-free-DOF chain
  for only its lowest 2 modes matches a full 5-mode solve of the same
  system.

- **Post-M5 — `SparseSolver` switched to a real sparse solver (`faer`), and
  `modal_analysis` to a shift-invert Lanczos.** Not a numbered milestone
  (see §6 decisions #1 and #2's full writeups for the reasoning): the
  M1-M5 dense placeholders were deliberately deferred until they mattered,
  and "real problems are not small" was the trigger to stop deferring —
  for both the linear solve and, once that went sparse, the eigensolve
  that had been quietly relying on it staying dense-compatible.
  `core/tests/sparse_solver.rs` adds a 150-free-DOF chained-beam scale
  check for the linear solve; the Lanczos partial-spectrum test above
  covers the eigensolve. Every pre-existing milestone test (M1-M5) passed
  unmodified against both new solvers — a real correctness signal, not
  just "it compiles."

- **M6 — Mass, Rayleigh damping, Newmark, time-history analysis.** First
  transient analysis. Needs element/nodal mass matrices (lumped, to start)
  and Rayleigh damping. **Acceptance:** a simple SDOF or 2-DOF dynamic
  problem with a known closed-form or well-established numerical response
  matches within tolerance.

- **M7 — DispBeamColumn + fiber sections + full standard material catalog.**
  `BeamIntegration` (Gauss-Legendre/Lobatto), fiber section stress-resultant
  integration, and the remaining leaf materials (Hysteretic, Pinching4,
  Concrete01/02, Steel01/02) plus the **recursive** composite materials
  (Parallel, Series, MinMax — §4.2). This is the milestone where the
  `Material` enum's recursive-composition design (decided at M2, not
  deferred) gets exercised for real.

- **M8 — ForceBeamColumn.** Its own milestone per §2/§4.1 — nested
  element-level equilibrium iteration, architecturally distinct from every
  prior element. Don't attempt to generalize the M3/M7 element pattern to
  cover this; it needs its own `form_tangent_and_residual` shape (an
  internal iterative solve, not a direct evaluation).

- **M9 — Corotational geomTransf (large-displacement).** Only if confirmed
  still in scope (§4.4) — comparable complexity to M8. Revisit scope with
  the project owner before starting.

- **M10 — JS/TS worker integration.** Expand `carapace-wasm` to the real
  command-buffer interface (§5.1/§5.2): batch model-description ingestion,
  `runAnalysis`/`getResults` commands, typed-array boundaries via `js-sys`.
  This is where `pysees` actually starts calling Carapace instead of Xara.

- **M11 — Results persistence.** SQLite-over-OPFS integration (§5.3),
  throttled progress-snapshot `postMessage` streaming, decoupled from each
  other as designed.

---

## 8. Verification methodology

Every milestone with numerical output gets checked against a known-correct
value — either closed-form (where one exists, e.g. `δ = PL/(AE)` for M1) or
against **Xara's own native build as an oracle** for cases without a simple
closed form. The Xara C++ spike (native binary, not the wasm build — the
native one is faster to run and equally valid as ground truth) is the
reference implementation; when in doubt about expected behavior for a
formulation, that's the thing to run and diff against, not the OpenSees
documentation alone (docs can describe intent; the actual compiled behavior
is ground truth for "does my Rust port match").

**Verification order for every milestone:** native `cargo test` first, then
`wasm32-unknown-unknown` build + `wasm-bindgen` + Node execution. Don't skip
the native step to save time — it isolates logic bugs from
platform/toolchain bugs, which was the single most useful debugging
discipline throughout the Xara C++ investigation (it's how the `s_copy`
wasm-only signature-mismatch trap and the `RCM` double-free were each
identified as belonging to a specific layer, rather than being debugged
blind inside a wasm runtime).

---

## 9. Reference material location

The Xara C++ spike and the standalone ARPACK f2c feasibility spike that this
whole plan is grounded in are preserved (source-level artifacts and
dependency-closure lists, not build outputs) at `reference/xara-spike/` in
this repo — see `reference/xara-spike/README.md` for what's there and how it
maps to the sections above. The original work happened in an ephemeral
session scratchpad under `/tmp`; the durable copy in this repo is now the
canonical reference — don't go looking for the `/tmp` path, it may not exist.

Key Xara source files referenced throughout this plan, for direct inspection:

- `xara/SRC/runtime/runtime/BasicAnalysisBuilder.cpp` — the real analysis
  composition sequence (`setLinks`/`number`/`domainChanged`/`initialize`/
  `analyzeStatic`) that §5.4's typestate design is modeled on.
- `xara/SRC/runtime/runtime/{ModelRegistry,ProcessContext,XaraClassBroker}.h/.cpp`
  — the broker eager-construction problem (§2, problem #1).
- `xara/SRC/element/Frame/Force/ForceFrame3d.h` — nested iteration evidence
  for M8.
- `xara/SRC/material/wrapper/{ParallelMaterial,SeriesMaterial,MinMaxMaterial}.h`
  — recursive material composition evidence for M7.
- `xara/SRC/matrix/routines/cmx.h` — the pre-existing missing-implementation
  bug (§2, problem #4; not a Carapace concern, just documented context).
- `xara/OTHER/ARPACK/` (legacy tree) vs `xara/OTHER/arpack-ng/` — the
  f2c-compatibility distinction relevant to §6 decision #2.
