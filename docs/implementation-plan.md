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
  needed for modal and time-history analyses. ✅ Done at M5 (nodal) / M6
  (element).
- **Loads:** nodal loads (done in M1) and element loads (distributed loads on
  beam-columns — needed once `ElasticBeamColumn`/`DispBeamColumn` land). ✅
  Done at M3.
- **`geomTransf`:** `Linear` and `PDelta` are cheap (small correction terms).
  ✅ Done at M3. **Corotational is a separate, high-complexity milestone
  (M9)** if/when large-displacement analysis is actually needed — it
  requires element-local frame tracking and large-rotation updates,
  comparable in complexity to force-based elements. Confirm this is
  actually in scope before starting it; it was flagged as "decide if you
  need it" during scoping, not committed.
- **Rayleigh damping:** needed for time-history analysis (M6). ✅ Done.

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

- **M6 — Mass, Rayleigh damping, Newmark, time-history analysis.** ✅ Done.

  **Element mass.** `Truss`/`ElasticBeamColumn` gained a `density` field
  (`with_density`, default 0.0 — massless, so no existing model changes
  behavior) and `form_mass` (`core/src/model/element.rs`,
  `core/src/model/beam.rs`): half the element's total mass
  (`density*area*length`) lumped to each node's translational DOFs, zero
  rotational contribution — the simplest standard lumped-mass model, per
  §4.4's "lumped, to start." `Domain::assemble_mass_diagonal` (from M5,
  now additive over `Node::mass` *and* every element's `form_mass`) stays
  a `DVector`, not a matrix — `M` is diagonal by construction regardless
  of how many contributors feed it.

  **`RayleighDamping`** (`core/src/analysis/damping.rs`): a concrete
  `{ alpha_m, beta_k }` struct (`C = alpha_m*M + beta_k*K`), not an enum —
  it's the one form of damping in scope, per §4.4.

  **`TransientAnalysis`** (`core/src/analysis/transient.rs`): a type
  distinct from `Analysis`/`AnalysisBuilder`, not a third `Integrator`
  variant — static analysis (scalar load factor, optionally Newton-iterated
  displacement at fixed time) and dynamic analysis (displacement/velocity/
  acceleration propagated through actual time) are different enough
  physical processes that unifying them would blur both, the same reasoning
  behind `modal_analysis` being a free function rather than an `Analysis`
  mode. Fixed at Newmark's "average acceleration" parameters
  (`beta=1/4, gamma=1/2` — unconditionally stable, OpenSees' own default),
  one linear solve per step (exact for this milestone's linear-elastic
  scope, mirroring `Algorithm::Linear`'s identical caveat for the static
  case — a Newton-iterated corrector for nonlinear dynamic response is a
  natural extension, not built until needed). `TransientAnalysis::new`
  computes a consistent initial acceleration from equilibrium at t=0
  (`M*a0 = F0 - K*u0 - C*v0`) so free-vibration problems (the standard
  validation case) don't need a fabricated starting acceleration.
  `Node` gained `velocity`/`acceleration` fields and
  `with_initial_displacement`/`with_initial_velocity` builders (bypassing
  the usual equation-driven displacement path, since an initial condition
  is given, not solved for). `Domain` gained the gather/scatter plumbing
  (`gather_displacement`/`velocity`/`acceleration`, `scatter_state`) plus
  `assemble_newmark_system` (one traversal producing both the effective
  sparse stiffness `K_eff = (1+a4*beta_k)*K + (a1+a4*alpha_m)*diag(mass)`
  and the `K*damp_vector` term the effective load needs, reusing the same
  raw triplets `assemble_tangent_and_resistance` builds from, and
  `multiply_stiffness`, a single sparse mat-vec product used once for the
  initial-acceleration solve) — no new `faer` sparse-matrix-arithmetic API
  needed; scaling/adding was done directly on the triplets.

  **Acceptance verified:** `core/tests/m6_dynamics.rs` (native) +
  `wasm-bridge`'s `damped_sdof_free_vibration_displacement` (wasm32 +
  Node). Two element-mass checks (`Truss` and `ElasticBeamColumn`) against
  the exact SDOF closed form `omega = sqrt(k/m)` via `modal_analysis`
  (solver precision, 1e-9 — a static/modal check, no time discretization
  involved). Two Newmark + Rayleigh-damping checks — undamped and 5%
  mass-proportional-damped SDOF free vibration — against the classical
  closed forms `u(t) = u0*cos(omega*t)` and `u(t) = exp(-xi*omega*t)*u0*
  [cos(omega_d*t) + (xi*omega/omega_d)*sin(omega_d*t)]`. These match to
  ~7e-6 (checked against a 1e-4 tolerance) — real, small Newmark
  discretization error at this `dt`/period ratio (~1/630), not solver
  precision; the first milestone where a "closed form" acceptance check
  is honestly approximate rather than exact, and the test comments say so.

- **Pre-M7 — material trial/commit state + `Domain` snapshot/restore.** ✅
  Done. Landed ahead of M7's fiber-section/full-catalog work because that
  work needs real path-dependent materials (Steel01, Concrete01/02,
  Hysteretic, Pinching4) to be correct, and *how* material state is
  organized needed settling first — see the design discussion this session
  that led here.

  **Not** OpenSees's parallel `C*`/`T*` field-pair pattern (verified
  directly against `xara/SRC/material/uniaxial/{UniaxialMaterial,Steel01}.h`
  before deciding against it) — that duplication solves a C++ problem Rust
  doesn't have. Instead: a `Material` holds only committed history (e.g.
  `ElasticPP`'s new `ep`, real permanent plastic strain now, not the old
  reversible envelope — construct via `Material::elastic_pp`, not the
  struct literal, since `ep` is internal history, not a parameter). Two
  pure (`&self`, non-mutating) views over that state, sharing one private
  `evaluate`: [`trial_stress_tangent`](core/src/model/material.rs) — called
  every Newton iteration, always relative to the same fixed committed
  baseline, which is the real correctness property this design is built
  around (not just tidiness: if a material mutated on every trial call,
  iteration N's result would depend on iteration N-1's discarded trial,
  compounding error every time a step needs more than one iteration) — and
  [`commit`](core/src/model/material.rs), called once, only after a step
  actually converges, producing the new committed `Material`.
  `Element::commit`/`Domain::commit` are the only things that ever call it,
  and only from `Analysis::step`/`TransientAnalysis::step`'s success path.

  This also means materials need no rollback machinery of their own —
  nothing mutates until commit, so a failed step never left one in an
  inconsistent trial state to revert. Nodal displacement/velocity/
  acceleration still do (mutated eagerly every Newton iteration, correctly,
  but not something a failed step should leave half-applied) — fixed by
  giving `Domain` `#[derive(Clone)]` and having `Analysis::step`/
  `TransientAnalysis::step` snapshot it up front, restoring on any `Err`
  path via a `step`/`try_step` split. This closes a real, pre-existing gap:
  before this, a `FailedToConverge` step already left displacement
  partway through its last (discarded) Newton iterations, for every
  milestone back to M4 — it just never mattered while every material was
  a pure function of current displacement.

  A useful side effect for M7: recursive composite materials
  (Parallel/Series/MinMax) get correct commit/revert for free once they
  exist — nothing needs to manually propagate `commit()` down through a
  nested tree the way OpenSees's wrapper materials forward `commitState()`/
  `revertToLastCommit()` to each child by hand; `Domain::commit`'s single
  top-level call reaches the whole nested structure uniformly.

  **Acceptance verified:** every M1-M6 test passed either unmodified or
  with only the `Material::elastic_pp` constructor rename (a real
  regression check that the trial/commit split doesn't change single-step
  behavior, which it shouldn't — see `core/src/model/material.rs`'s
  `elastic_pp_clamps_symmetrically_on_first_loading`). New coverage in
  `core/src/model/material.rs` (`elastic_pp_remembers_permanent_set_after_
  commit` — a partial unload after yielding shows real elastic unloading
  the old reversible envelope couldn't) and `core/tests/material_state.rs`:
  a *two-step* `Analysis` (load past yield, then partially unload) matches
  the closed-form path-dependent answer and explicitly diverges from the
  path-*independent* (wrong) answer a stateless model would give; a
  separate test confirms a step with an unsatisfiable `ConvergenceTest`
  (`max_iter: 0`) leaves the domain's displacement completely unchanged
  rather than partially applied.

- **M7 — DispBeamColumn + fiber sections + full standard material catalog.**
  Being landed in stages (agreed with the project owner given this
  milestone's unusual size relative to every other one — six new
  materials, several genuinely complex hysteretic state machines, plus the
  recursive composites, all bundled together in the original scoping).

  **Stage 1 — `BeamIntegration` + fiber sections + `DispBeamColumn`.** ✅
  Done. `BeamIntegration::{Legendre, Lobatto}` (`core/src/model/
  integration.rs`) — point/weight tables copied directly from OpenSees's
  own source (`xara/SRC/quadrature/Frame/{Legendre,Lobatto}
  BeamIntegration.cpp`), not computed via a general-purpose quadrature
  crate: domain-specific numerics are copied from the reference
  implementation, the same policy already applied to element/material
  formulations — heavily-optimized linear algebra (`nalgebra`/`faer`) is
  the one place a library is preferred over hand-rolling. `FiberSection`
  (`core/src/model/fiber_section.rs`): given section deformation
  `(eps0, kappa)`, each fiber's strain is `eps0 - y*kappa`; the
  work-conjugate stress resultants `(N, M)` and section tangent are
  derived directly from virtual work (documented in the module), not
  assumed. `DispBeamColumn` (`core/src/model/disp_beam_column.rs`): the
  same cubic-Hermite/linear-axial local shape functions
  `ElasticBeamColumn` uses give the strain-displacement field directly
  from nodal displacements (a single evaluation, no Newton loop of its
  own — unlike M8's `ForceBeamColumn`), integrated over `BeamIntegration`
  points, each owning an independent `FiberSection` (prismatic-member
  assumption: same fiber layout everywhere, but independent material
  history per point, since each point's strain path differs). `Linear`
  transform only, no element loads yet — both explicitly deferred, not
  silently unsupported (see the type's doc comment).

  **Acceptance verified:** a 2-fiber elastic section
  (`y = ±sqrt(Iz/A)`, area `A/2` each) reproduces `EA`/`EI` exactly
  (`FiberSection`'s unit test), and — because the underlying integrand is
  then a low-degree polynomial any >=2-point Gauss rule integrates
  exactly — a cantilever `DispBeamColumn` built from it matches
  `ElasticBeamColumn`'s closed-form tip deflection *exactly* (1e-9,
  `core/tests/m7_disp_beam_column.rs`, native + wasm32/Node via
  `disp_beam_column_cantilever_tip_deflection`), and Legendre vs. Lobatto
  integration give the identical exact answer — not an approximation
  that happens to be close, a real equivalence.

  **Stage 2 — `Steel01`, `Concrete01`, `Parallel`/`Series`/`MinMax`.** ✅
  Done. `Material` (`core/src/model/material.rs`) gained five new
  variants, all following the trial/commit purity recipe from the type's
  doc comment — every `evaluate` is still a pure `&self` function of the
  committed state, never mutating anything. `Steel01`: bilinear kinematic
  hardening with isotropic hardening on load reversal (Filippou 1983),
  ported from `Steel01::determineTrialState` — the `Steel01Loading`
  history field is a real 3-state enum, not OpenSees' sentinel `-1/0/1`.
  `Concrete01`: Kent-Scott-Park compression envelope with the
  Karsan-Jirsa degrading unload/reload rule and zero tensile strength,
  ported from `Concrete01::setTrialStrain` (the version actually driving
  the state machine — `determineTrialState` is dead code in the OpenSees
  source, never called). `Parallel(Vec<(Material, f64)>)`: same strain to
  every child, factor-weighted stress/tangent sum. `Series`: same stress
  across children, total strain = sum of individual strains — not
  directly invertible, so `evaluate` runs its own small internal
  flexibility-based Newton loop entirely within one call (a deliberate,
  documented deviation from OpenSees' `SeriesMaterial`, which instead
  carries a warm-started trial state across many *external* Newton calls
  — see the type's doc comment for why applying the stress correction
  every inner iteration, instead of once at the end, is necessary here).
  `MinMax`: wraps one material with strain bounds; once a committed
  strain exits them, `failed` permanently latches `true`, stress reads
  `0.0`, and tangent reads `1e-8 * inner.initial_tangent()` (not exactly
  zero, matching OpenSees — a structurally zero tangent risks a singular
  global tangent). New `Material::initial_tangent()` method supports
  `MinMax`'s post-failure tangent and the composites' propagation of it.

  **Structural ripple:** `Material` dropped `#[derive(Copy)]` (kept
  `Clone`) now that `Parallel`/`Series`/`MinMax` hold a `Vec`/`Box` of
  children — `Fiber` (`core/src/model/fiber_section.rs`) lost `Copy` too,
  and `ZeroLength`'s `[Option<Material>; NDF]` array-repeat initializer
  became `std::array::from_fn(|_| None)`. No other call sites needed
  changes (every other site already used `.clone()` or read-only
  references).

  **Acceptance verified:** unit tests per variant in `material.rs`
  (bilinear hardening, isotropic-hardening shift on reversal, the
  Kent-Scott-Park envelope, Karsan-Jirsa unload degradation, `Parallel`
  of two `Elastic`s exactly equaling one `Elastic` with the summed
  stiffness, `Series` of two equal `Elastic`s halving the stiffness and
  splitting strain for equal child stress, `MinMax` latching failure and
  its `1e-8`-scaled tangent). Plus the suggested composite-section
  acceptance test (`core/src/model/fiber_section.rs`): a `FiberSection`
  built from one `Steel01` rebar layer and two `Concrete01` layers,
  pushed past first yield/cracking, deeper into the nonlinear range, then
  partially unloaded — checked at each state against a hand computation
  (each fiber's `Material` evaluated independently at its own strain and
  summed into `N`/`M`, outside `FiberSection`'s own code path). No closed
  form exists for a nonlinear composite section, so this is what
  "verified" means for this stage, unlike stage 1's exact elastic match.

  **Stage 3 — `Steel02`, `Concrete02`.** ✅ Done. `Steel02`: Menegotto-
  Pinto smooth transition curve between an elastic and a strain-hardening
  asymptote, with Filippou isotropic hardening on reversal — genuinely
  more involved than `Steel01`'s piecewise-linear rule, as flagged in the
  original stage-3 handoff notes, but the same trial/commit recipe still
  applies cleanly: `kon` (`Steel02Kon`, mirroring `Steel02::konP` minus
  its `sigini`-only pinned state — `sigini` itself isn't supported, same
  scope decision as `Steel01`) tracks which asymptote the curve is
  heading toward; `asymptote_strain`/`asymptote_stress` and
  `reversal_strain`/`reversal_stress` are the two curve anchors, recomputed
  (with the isotropic-hardening shift) on every direction reversal.
  `Concrete02`: `Concrete01`'s compression envelope plus linear tension
  softening (`ft`/`ets`) and a stiffer EERC-report reloading rule —
  smaller increment over `Concrete01` than `Steel02` is over `Steel01`,
  confirming the original handoff note's expectation. `tension_strain`
  (`Concrete02`'s `deptP`) is the one genuinely new history field: the
  peak tensile-excursion strain, with no `Concrete01` analogue since
  `Concrete01` carries no tension branch at all. One real (not
  simplified-away) divergence between the two concrete ports: the flat
  crushed/ruptured branches' tangent floor is `0.0` in `Concrete01` but
  `1e-10` in `Concrete02`, matching each one's own OpenSees source exactly
  rather than reconciling them to look alike.

  **Acceptance verified:** unit tests per variant in `material.rs` — for
  `Steel02`, the Menegotto-Pinto curve gives the *exact* initial tangent
  `E0` at zero strain (eps_ratio is exactly `0` there, independent of
  `R`), approaches the hardening slope `Esh` far past yield, starts a
  reversal near the elastic slope again, and correctly latches
  `max_strain`/`min_strain` across a reversal; for `Concrete02`, the
  compression envelope matches `Concrete01`'s formula exactly (as
  expected, since it's unchanged), tension is linear-elastic below `ft`
  and softens at `-Ets` above it (verifying real tensile behavior
  `Concrete01` has none of), and reload stiffness degrades below `Ec0`
  after a compression excursion.

  **Stage 4a — `Hysteretic`.** ✅ Done (`Pinching4`, stage 4's other
  material, is not — see below). Multi-linear (up to 3 points per side)
  backbone with pinching, deformation/energy-based damage, and ductility-
  degraded unloading stiffness, ported from `HystereticMaterial`. Genuinely
  the most-fields leaf material yet (~30), which motivated a structural
  change beyond the usual recipe: `Material::Hysteretic(Box<
  HystereticFields>)` — the payload is a boxed, `Copy` struct rather than
  inlined enum fields like every other leaf. Inlining it would have made
  *every* `Material` (including a plain `Elastic { e: f64 }`) pay the
  largest variant's stack size everywhere a `Material` is stored (`Fiber`,
  `ZeroLength`'s per-DOF array, every composite's `Vec<Material>`) —
  confirmed by `clippy::large_enum_variant` firing on `Element` (which
  embeds `Material` three deep via `ZeroLength`) once `Hysteretic` landed
  inline; boxing it made the warning disappear. `Pinching4` should follow
  the same boxed-payload pattern — it will be larger still.

  Two of OpenSees' early-return guards weren't ported (documented in
  `evaluate_hysteretic`'s own comment): one reproduces exactly via the
  general formula regardless (verified by a test); the other — skipping
  re-derivation when `dStrain==0` — has no interior-branch fallback at all
  in the original source (a real gap, since the C++ only reaches that
  branch after the guard already intercepted the zero-increment case).
  Resolved by treating `dstrain>=0.0` (not `>0.0`) as the "positive
  increment" branch, which is both well-defined at exactly zero and
  consistent with this same file's own tie-break for the very first-ever
  increment.

  **Acceptance verified:** unit tests in `material.rs` against a shared,
  fully hand-computed symmetric bilinear-then-flat backbone (`beta=0`,
  `damfc1=damfc2=0` to keep the arithmetic tractable) — exact initial
  tangent, exact envelope stress on first loading through two segments,
  and (the real test of the state machine) reloading after a full
  positive-then-negative reversal: the pinching reference point `rot_nu`
  and the resulting pinched-plateau stress/tangent both matched a from-
  scratch hand derivation to `1e-9`, including confirming the pinched
  path is genuinely softer than straight elastic reload (the whole point
  of pinching) rather than just structurally different.

  **Stage 4b — `Pinching4` — not landed.** Read in full
  (`xara/SRC/material/uniaxial/Pinching4Material.cpp`, ~1725 lines) but
  deliberately not ported this session: it's a 5-state cyclic-damage state
  machine (states 0-4, not a simple positive/negative/interior split) with
  an intricate ~250-line geometric correction cascade (`getState3`/
  `getState4`) that repairs a trilinear unload-reload path when it comes
  out geometrically invalid (reload point behind the target point, wrong-
  sign slopes, etc.) through a sequence of special-cased fallback
  constructions — exactly the kind of logic where a rushed, unverified
  translation risks a subtly wrong structural-engineering result. See
  §7.1 below for the full distilled handoff (state machine shape, the key
  insight that its "damaged envelope" arrays are cheap to treat as
  *derived*, not stored, committed state, and where the real risk is).

### 7.1 M7 handoff: remaining stages

Written for a session with no memory of the ones that did stages 1-3 —
read this before opening any OpenSees source yourself; it's the distilled
result of already having read the relevant files once.

**The general porting recipe** (applies to every material below): OpenSees
implements each as `setTrialStrain`/`determineTrialState` mutating
parallel `T*`/`C*` field pairs (see `xara/SRC/material/uniaxial/
UniaxialMaterial.h` and e.g. `Steel01.h`). Carapace's `Material` (§
"Pre-M7" entry above) does **not** carry that duplication — every method
you port becomes a pure `fn evaluate(&self, strain: f64) -> (stress,
tangent, next_self)` where `self` is *always* the committed state (there
is only one copy of each history field, not a `T`/`C` pair). The
mechanical translation:

1. Read the C++ `determineTrialState`/`setTrialStrain`. Every place it
   reads a `C*` field, read `self`'s corresponding field instead. Every
   place it reads/writes a `T*` field, use a local variable instead (never
   write to `self` — `evaluate` takes `&self`, not `&mut self`).
2. Wherever the algorithm references `Cstrain`/`Cstress` directly (not
   just `CminStrain`/`Cloading`-style history), your `Material` variant
   needs to carry the **last committed strain and stress as explicit
   fields**, not just the "obvious" history variables — `Steel01` and
   `Concrete01` both do this (see below) and it's easy to miss if you
   only skim for the fields declared as "History Variables" in the header,
   since `Cstrain`/`Cstress` are filed separately as "State Variables" in
   `Steel01.h` but are just as load-bearing for the algorithm.
3. The function's return value is `(Tstress, Ttangent, Self { ..committed
   fields updated.. })`. `trial_stress_tangent` discards the third
   element; `commit` discards the first two — exactly the same split
   `ElasticPP` already uses (`core/src/model/material.rs`).
4. Port the *math* faithfully (it's validated, published research — see
   each material's `// References` comment in its `.cpp`), but the
   *state-machine mechanics* (trial/commit, no field duplication) always
   follow the recipe above, not the C++ structure.

Stage 2 already went through the `Material`-loses-`Copy` ripple (see the
M7 entry in §7 above); stages 3-4 only add more leaf variants and
shouldn't need any further structural changes on this front. Stage 3
(`Steel02`/`Concrete02`, also done — see §7) confirmed the recipe scales
past `Steel01`/`Concrete01` cleanly, including a genuinely more involved
state machine (`Steel02`'s Menegotto-Pinto curve): the trial/commit purity
split was never the hard part, faithfully tracing each C++ source's
branches is. Stage 4a (`Hysteretic`, also done — see §7) added one more
structural pattern worth reusing for `Pinching4` below: once a variant's
field count gets large (`Hysteretic` has ~30), box the payload as a
`Copy` struct (`Material::Hysteretic(Box<HystereticFields>)`) instead of
inlining it into the enum — otherwise *every* `Material` pays that
variant's stack size everywhere a `Material` is stored, not just the
values that are actually that variant. `clippy::large_enum_variant`
(checked via `cargo clippy --workspace --all-targets`, already clean
otherwise) is the trip-wire for this — if it fires on `Element` after
`Pinching4` lands, that's the fix.

#### Stage 4b: `Pinching4` (the only piece of M7 not yet landed)

Read in full this session (`xara/SRC/material/uniaxial/
Pinching4Material.cpp`, ~1725 lines) but deliberately not ported — see the
M7 entry in §7 for why. This is the distilled result of that reading;
still budget the most time of anything in M7 for actually writing it,
the porting recipe notwithstanding.

**Shape of the state machine — genuinely different from every other
material here.** Not a positive/negative/interior split. Five states
(`Tstate`/`Cstate`, an `int` in the source — model as a real 5-variant
enum): `0` (near-origin elastic, before the first real excursion in
either direction), `1` (on the positive envelope), `2` (on the negative
envelope), `3` (a trilinear unload/reload path connecting a negative
excursion back toward — but not necessarily reaching — the positive
side), `4` (the mirror: a trilinear path connecting a positive excursion
back toward the negative side). `getstate()` (`Pinching4Material.cpp`
~874) is the transition function; state changes only when the trial
strain exits the current state's `[lowTstateStrain, hghTstateStrain]`
bracket (tracked explicitly, not re-derived) or the strain-rate sign
flips (`du*CstrainRate<=0.0`, the `cid` flag) — read this function first,
before anything else, since every other piece exists to compute what a
transition *into* a given state initializes.

**The envelope is 6 points, not 4** — `SetEnvelope()` (~810) synthesizes
point 0 (a tiny near-origin point at `1e-4` of the first real breakpoint,
giving the elastic state-0 slope `k = max(stress1p/strain1p,
stress1n/strain1n)`) and point 5 (a far extrapolated point at `1e6 *
strain4p`/`strain4n`, extending the 3-4 segment's slope if it's still
rising, or flat at `1.1x` the point-4 stress if it's already softening)
around the 4 real backbone points — `posEnvlpStress`/`negEnvlpStress`/
`*Tangent` (~1038-1121) piecewise-search this 6-point array, not a
hand-unrolled 4-branch `if` chain like `Hysteretic`'s envelope — port as
a small loop or `Vec`/array of breakpoints, not 6 named fields per side.

**The key implementation insight, easy to miss on a first read:**
`envlpPosDamgdStress`/`envlpNegDamgdStress` (the "damaged" — strength-
degraded — envelope actually used by `posEnvlpStress`/`negEnvlpStress`)
*look* like they need to be stored, mutable, per-side 6-element arrays
mirrored in the committed state. They don't. Every write site computes
them as `envlpPosStress(i)*(1.0-gammaFUsed)` (or `...Neg...`), and
`gammaFUsed`/`kElasticPosDamgd`/`kElasticNegDamgd` only change at a state
*transition* (inside `getstate()`, using `gammaFUsed = CgammaF` — the
*committed*, not trial, damage factor) — never on every `evaluate` call.
So the actual state that needs to persist is just two scalars per side
(`gamma_k_used`, `gamma_f_used` — `Hysteretic`-recipe committed fields,
point 2 of the general recipe above, distinct from the `TgammaK`/
`TgammaD`/`TgammaF` computed by `updateDmg` below), and the "damaged
envelope" is a cheap on-demand computation (`envlp_stress[i] * (1.0 -
gamma_f_used)`) wherever the C++ reads `envlpPosDamgdStress`/
`envlpNegDamgdStress` — never its own stored field. Confirm this by
checking every read site in `getstate()`, `posEnvlpStress`/
`negEnvlpStress`, and `getState3`/`getState4` actually only ever wants
the *last-transition* value, not something recomputed mid-call — it does
(re-verify this claim before relying on it; it's the single highest-
leverage simplification for this material, so worth double-checking
against the source rather than taking the handoff's word for it).

**`getState3`/`getState4`** (~1124-1395, mirror images of each other,
same structure as `posEnvlpStress`/`negEnvlpStress` and `Envlp3`/
`Envlp4Stress`/`Tangent` below) build the 4-point path
`(lowTstateStrain, lowTstateStress) → (state3/4[1]) → (state3/4[2]) →
(hghTstateStrain, hghTstateStress)` for states 3/4's trilinear
unload-reload. Points 0 and 3 are given (the state's bracket); points 1
and 2 are constructed from the pinching parameters (`rDispN/P`,
`rForceN/P`, `uForceN/P`) and then run through a cascade of geometric
sanity checks — reload point behind the target point, reload stiffness
exceeding unload stiffness, wrong-sign or decreasing segments — each
with its own fallback (usually "redraw points 1/2 as a straight
33%/67% split between 0 and 3"). This cascade is the highest-risk part
of the whole material to mistranslate: port it as a literal, checkable
line-by-line trace of the C++ (variable names and all, adapted only for
the trial/commit purity split), not a "cleaned up" reimplementation from
first principles — the fallback conditions encode real degenerate-
geometry cases discovered empirically, not something to rederive from
the docstring-level description above. `Envlp3Stress`/`Envlp4Stress`/
`*Tangent` (~1397-1496) then piecewise-interpolate this 4-point path,
structurally identical to `posEnvlpStress` but over `state3`/`state4`
instead of the 6-point backbone.

**`updateDmg`** (~1498-1541) computes `TgammaK`/`TgammaD`/`TgammaF` from
the ductility demand (`umaxAbs`/`uultAbs`) and either the dissipated-
energy ratio or an explicit cycle count (`DmgCyc` — port as a 2-variant
enum, `EnergyBased`/`CycleBased`, not a raw `int`), clamped against
`gammaKLimit`/`gammaDLimit`/`gammaFLimit` and an *envelope-derived*
stiffness limit (`gammaKLimEnv`, from how much the current
maximum-ductility-demand point has already softened relative to the
undamaged elastic stiffness) so damage can never make the tangent
increase. This runs unconditionally at the end of every `setTrialStrain`
(**not** gated behind a state transition, unlike `gammaKUsed`/
`gammaFUsed` above) — it's the source of `TgammaK`/`TgammaD`/`TgammaF`,
which then become next step's `CgammaK`/`CgammaD`/`CgammaF`, which then
become `gammaKUsed`/`gammaFUsed` *only if* that next step happens to
transition state. Three separate damage "clocks" (`used` vs. `committed`
vs. `this step's fresh computation`) — keep them straight; this is the
easiest place in the whole material to introduce a subtle lag/no-lag bug
that only shows up several cycles into a real hysteretic loop, not in a
single-step unit test.

**Suggested acceptance tests**, in increasing order of how much of the
state machine they exercise: (1) state-0 exact elastic slope at zero
strain (`envlpPosStress(0)/envlpPosStrain(0)`, i.e. the synthetic
near-origin point's slope — analogous to every other material's
"exact initial tangent" test); (2) state-1/2 envelope match on first
loading through all 4 real breakpoints plus the extrapolated point 5,
same spirit as `Hysteretic`'s bilinear-envelope test but now against the
6-point piecewise search; (3) a full reversal into states 3/4, checked
against a hand-walked `getState3`/`getState4` trace for a case that does
*not* hit any of the degenerate-geometry fallbacks (pick backbone/
pinching parameters that keep the trilinear path well-behaved, to keep
the hand computation tractable — the same tactic `Hysteretic`'s
`beta=0`/`damfc1=damfc2=0` test material used); (4) at least one test
that deliberately drives a fallback branch in `getState3`/`getState4`
(e.g. pinching parameters that push the reload point behind the target
point), since that cascade is exactly the risk called out above and
deserves its own direct check, not just incidental coverage from (3).

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

---

## 10. Known follow-ups / performance backlog

Not scheduled against any milestone — things noticed along the way that are
deliberately *not* being acted on without a profile or a concrete need,
recorded here so they don't have to be rediscovered from scratch. Don't
implement anything in this section speculatively; it exists so a future
session (or this one, later) knows the question was already thought through
once.

### 10.1 `trial_stress_tangent` unnecessarily allocates for boxed materials

**The situation:** `Material`'s large leaf variants (currently `Hysteretic`;
`Pinching4` will be the same — see §7.1's stage-4b handoff) are boxed
(`Material::Hysteretic(Box<HystereticFields>)`) specifically so one huge
variant doesn't force every `Material` value in the system — including a
trivial `Elastic { e: f64 }` — up to its size (§3.1's closed-enum dispatch
means `size_of::<Material>()` is always the size of the *largest* variant;
confirmed empirically: `Steel02`, at ~20 `f64` fields, currently sets it at
160 bytes for every `Material` value that exists, regardless of which
variant it actually is).

That boxing has a real, separate cost: `Material::evaluate` (the shared
core behind both `trial_stress_tangent` and `commit` — see `material.rs`'s
top-level doc comment for the trial/commit design) always constructs the
*full* next-state `Material`, box included, even when the caller is
`trial_stress_tangent` and immediately discards it. Since
`trial_stress_tangent` runs once per Newton iteration (far more often than
`commit`, which runs once per converged step), every Newton iteration
touching a boxed material allocates and immediately frees a small heap
block for no reason.

**Why this isn't being fixed reflexively:** two reasons, one about scope and
one about actual impact.

- *Scope*: a real fix has an easy part and a hard part. The easy part —
  giving a boxed leaf material's own `evaluate_xxx` a variant that returns
  the raw, unboxed fields struct (e.g. `(f64, f64, HystereticFields)`
  instead of `(f64, f64, Material)`), and letting `trial_stress_tangent`
  and `commit` each decide separately whether to box the result — is clean
  and contained, and wouldn't change the trial/commit purity design at all
  (see §10 intro and the "why not OpenSees' mutable `T*`/`C*` pattern"
  discussion this section is transcribing — same conclusion: keep the pure
  design, it already gets step rollback "for free" via `Domain` cloning,
  documented at `core/src/analysis/state.rs:36` and
  `core/src/model/domain.rs:9`). The hard part is that `Parallel`/
  `Series`/`MinMax` all recursively call `child.evaluate(strain)` through
  the single generic dispatcher, which has no way to know whether the
  *outer* caller wanted a trial or a commit — so a boxed material sitting
  inside a composite would still box on every trial unless composites also
  got a trial/commit-aware recursive path, which means duplicating (or
  generically parameterizing) the entire `evaluate` dispatch tree, not
  just patching one material's call sites. That's a materially bigger,
  more invasive change than the leaf-material case alone.
- *Impact*: unmeasured. A handful of small heap allocations per Newton
  iteration per boxed-material fiber is plausibly noise next to the cost
  already paid per iteration (global tangent assembly, sparse LU solve) —
  unless a model has enough `Hysteretic`/`Pinching4` fibers, in enough
  elements, over enough iterations, for it to add up. Nobody has profiled
  this against a real model yet.

**If this becomes worth doing:** benchmark first (allocation count and/or
wall time per step, before/after, on a model with many boxed-material
fibers under `Algorithm::Newton` — something like a multi-story frame with
`Hysteretic` or `Pinching4` sections would exercise it realistically) to
confirm it's actually load-bearing before spending the implementation
effort, and scope the first pass to the leaf-material case only (skip
composites unless the profile specifically implicates a boxed material
used inside a `Parallel`/`Series`/`MinMax`).
