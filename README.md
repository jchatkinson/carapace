# Carapace

A thin, curated, WebAssembly-native finite element analysis engine, written in Rust,
heavily inspired by [OpenSees](https://opensees.berkeley.edu/) / [Xara](https://xara.so).

Carapace is **not** a port of OpenSees/Xara. It's a from-scratch reimplementation of
a deliberately small subset of element formulations, material models, and analysis
types — the ones actually used in practice by [pysees](../pysees), Carapace's
companion frontend — designed from day one to run efficiently inside a browser
Web Worker as a `.wasm` module, with an architecture that fixes several structural
issues found while investigating a more literal "compile OpenSees/Xara's C++ to
wasm" approach (see below).

## Why not just compile Xara to wasm?

We tried. A working spike proved it's *possible*: a hand-wired `Domain`/`Node`/
`Truss`/`ElasticMaterial`/`BandGenLinSOE` static analysis, using Xara's actual
production C++ classes, was compiled to WebAssembly (via `f2c` for the Fortran
BLAS/LAPACK/ARPACK dependencies, and `emcc` for the C++) and produced correct
results under Node. That spike is the numerical oracle Carapace verifies against
(see `docs/implementation-plan.md`, "Verification methodology").

But getting there surfaced real architectural friction that a clean-sheet design
avoids entirely:

- **`XaraClassBroker`** (Xara's runtime element/material factory) is eagerly
  constructed by `ModelRegistry`, which `BasicAnalysisBuilder` depends on — so
  even a "just one truss element" build transitively pulls in ~230 classes
  (Shell elements, MVLEM, Fortran `feap` materials) via a single `#include`.
- **Ownership is convention, not compiler-enforced.** A stack-allocated `RCM`
  numberer passed where `DOF_Numberer` expected heap ownership caused a real
  double-free, caught only by AddressSanitizer, not the compiler.
- **Tcl leaks through headers that don't need it** — e.g. `SRC/logging/logging.cpp`
  unconditionally `#include <tcl.h>` despite using zero `Tcl_*` symbols.
- The Fortran numerics (BLAS/LAPACK/ARPACK) work fine via `f2c`, but add a whole
  translation-and-runtime-shimming step that a Rust rewrite of the *formulations*
  (not the framework) sidesteps by construction.

None of this is a knock on Xara — it's 90s/2000s C++ architected for arbitrary
extensibility, MPI parallelism, and database persistence, none of which a
single-threaded browser worker needs. Carapace's architecture (enum-dispatched
closed catalogs, typestate-enforced analysis wiring, no broker/serialization
layer) is designed around what we're actually building, not what OpenSees needed
to be in 1997.

## Status

M1–M6 done. Real `Domain`/`Node`/
`Element::{Truss,ZeroLength,ElasticBeamColumn}`/
`Material::{Elastic,ElasticPP,Gap,Ent}` types and a typestate-composed
`Analysis` replace the M0 closed-form placeholder, verified against known
values on both native and wasm32 + Node. `ElasticBeamColumn` (M3) brought
`NDF` up to 3 (2D translation + rotation), a `GeomTransf` (`Linear`/`PDelta`)
local↔global transform, and element loads (uniform transverse).
`Algorithm::NewtonRaphson` + `Integrator::DisplacementControl` (M4) generalize
analysis composition past the M1 baseline of `LoadControl` + `Linear` —
Newton iteration is what lets a single step now cross a material's nonlinear
regime boundary correctly, which `Linear`'s one-shot solve structurally
can't do. `modal_analysis` (M5) adds lumped nodal mass and a generalized
eigenproblem solve. Post-M5, both solvers in the codebase moved from dense
to real sparse implementations, once "real problems are not small" made the
M1-M5 dense placeholders indefensible: `SparseSolver` now uses `faer`'s
sparse LU (COLAMD/AMD fill-reducing ordering, pure Rust, verified on wasm32
+ Node), with `Domain` assembling stiffness directly into sparse triplets
rather than a dense buffer; `modal_analysis` now runs a shift-invert Lanczos
that reuses `SparseSolver`'s sparse LU for each iteration's solve and only
asks for the `num_modes` actually wanted (a dense full-spectrum eigensolve
can't do partial-spectrum extraction at all — it was the wrong shape of
computation independent of raw size, since real modal analysis and
mode-superposition damping only ever want the lowest handful of modes out
of a model with many more DOFs than that). Every pre-existing milestone test
passed unmodified against both new solvers. `TransientAnalysis` (M6) adds
Newmark-beta time-history integration and `RayleighDamping`
(`C = alpha_m*M + beta_k*K`) — a type distinct from `Analysis`, not a third
`Integrator` variant, since propagating displacement/velocity/acceleration
through actual time is different enough from static equilibrium at a fixed
load factor to warrant its own type (the same reasoning behind
`modal_analysis` being a free function). `Truss`/`ElasticBeamColumn` gained
an optional `density` for element-consistent lumped mass, additive with
M5's nodal `Node::mass`. Verified against classical closed-form SDOF free
vibration (undamped and 5%-damped) to ~7e-6 — real Newmark discretization
error at typical `dt`/period ratios, the first milestone where a
closed-form check is honestly approximate rather than exact to solver
precision. Ahead of M7, `Material` gained a real trial/commit state design:
not OpenSees's parallel `C*`/`T*` field duplication (a C++ workaround Rust
doesn't need), but a single committed-state value with two pure views —
`trial_stress_tangent` (called every Newton iteration, always relative to
the same fixed baseline — a correctness requirement, not just tidiness) and
`commit` (called once, only after a step converges). `ElasticPP` is now
genuinely stateful (real permanent plastic strain, not the old reversible
envelope). `Domain` gained `Clone` and `Analysis`/`TransientAnalysis` now
snapshot-and-restore on a failed step — closing a real pre-existing gap
where a `FailedToConverge` step left displacement partway through its
discarded Newton iterations, for every milestone back to M4. See
`docs/implementation-plan.md`'s "Pre-M7" entry for the full design
reasoning.

M7 (`DispBeamColumn` + fiber sections + full standard material catalog) is
large enough — six new hysteretic materials plus three recursive
composites — that it's being landed in stages rather than one pass.
**Stage 1 done:** `BeamIntegration::{Legendre,Lobatto}` (point/weight
tables copied directly from OpenSees's own source, not a quadrature
crate — domain-specific FEM numerics are copied from the reference
implementation, the same policy as element/material formulations;
`nalgebra`/`faer` remain the exception, since heavily-optimized linear
algebra is the one place a library beats hand-rolling), `FiberSection`
(stress resultants and section tangent derived directly from virtual
work), and `DispBeamColumn` itself. A 2-fiber elastic section built to
reproduce `EA`/`EI` exactly makes a cantilever `DispBeamColumn` match
`ElasticBeamColumn`'s closed-form deflection *exactly*, not
approximately — a real equivalence check, verified native + wasm32/Node.
Remaining stages: `Steel01`/`Concrete01` + `Parallel`/`Series`/`MinMax`
(which get correct commit/revert for free from the trial/commit
mechanism just landed), then `Steel02`/`Concrete02`, then
`Hysteretic`/`Pinching4`.

## Workspace layout

```
carapace/
├── core/           # carapace-core: all FE logic. No wasm-bindgen dependency —
│                    compiles and tests identically on native and wasm32 targets.
├── wasm-bridge/     # carapace-wasm: thin wasm-bindgen layer exposing core's API
│                    to JS. Should contain no numerical logic of its own.
└── docs/
    └── implementation-plan.md   # comprehensive design doc + milestones — read this first
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

Always verify a change on the native target first, wasm32 second. That ordering
caught a real wasm-specific bug during the Xara C++ spike (an f2c/wasm-ld
function-signature mismatch that only trapped under wasm, not natively) — it'll
catch the Rust equivalent (e.g. numeric behavior that differs between targets)
faster than debugging inside a wasm runtime directly.

## Read next

**`docs/implementation-plan.md`** — architecture principles and rationale, full
feature scope, milestone-by-milestone plan with acceptance criteria, and open
decisions that need resolving early (linear algebra crate, ARPACK strategy).
Any new coding session on this repo should start there.
