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

M1 done: real `Domain`/`Node`/`Element::Truss`/`Material::Elastic` types and a
typestate-composed `Analysis` (`LoadControl` + `Linear` + a real linear solve)
replace the M0 closed-form placeholder, verified against the same oracle value
on both native and wasm32 + Node. See `docs/implementation-plan.md` for the
milestone roadmap, starting at M2 (`ZeroLength` + EPP/Gap/ENT materials).

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
