# Why not just compile Xara to WebAssembly?

Before starting Carapace as a from-scratch Rust engine, we investigated the
more obvious path: compile [Xara](https://xara.so) (the actively-maintained
fork of [OpenSees](https://opensees.berkeley.edu/)) directly to WebAssembly.

## The spike

A hand-wired `Domain`/`Node`/`Truss`/`ElasticMaterial`/`BandGenLinSOE` static
analysis, using Xara's actual production C++ classes, was compiled to
WebAssembly — `f2c` for the Fortran BLAS/LAPACK/ARPACK dependencies, `emcc`
for the C++ — and produced correct results under Node. This worked: it's a
real, running proof that the numerics and toolchain are sound, and it's the
numerical oracle Carapace's own tests verify against (see
`docs/implementation-plan.md` §8, "Verification methodology"). The spike
itself is preserved at `reference/xara-spike/` (source-level artifacts and
dependency-closure lists, not build outputs) — see that directory's README
for what's there.

## What it surfaced

Getting the spike running exposed architectural friction that a clean-sheet
design avoids entirely:

1. **`XaraClassBroker` eager construction.** `SRC/runtime/runtime/
   ModelRegistry.h` holds a `ProcessContext` by value, whose constructor
   unconditionally builds an `XaraClassBroker` (`SRC/runtime/runtime/
   XaraClassBroker.cpp`, ~230 `#include`s — Shell elements, MVLEM, Fortran
   `feap` materials). Since `BasicAnalysisBuilder` depends on `ModelRegistry`,
   any use of Xara's own analysis-composition class drags in nearly the whole
   element/material catalog at compile time, regardless of what the model
   actually uses.
2. **Ownership by convention, not enforcement.** `DOF_Numberer` assumes
   ownership of its `GraphNumberer&` and deletes it in its destructor — a
   comment-level convention, not compiler-checked. Passing a stack-allocated
   `RCM` where heap ownership was expected caused a real double-free in the
   spike, caught only by AddressSanitizer.
3. **Tcl coupling leaks through headers that don't need it.**
   `SRC/logging/logging.cpp` unconditionally `#include <tcl.h>` despite using
   zero `Tcl_*` symbols — purely incidental coupling.
4. **A pre-existing, unrelated bug found along the way:** `SRC/matrix/
   routines/cmx.h` declares `cmx_inv2..cmx_inv6` with no implementation
   anywhere in the checkout (git history shows the header deleted in one
   commit and reintroduced later without its source). Not wasm-specific —
   worth an upstream fix in Xara, unrelated to Carapace.
5. **Force-based elements do nested equilibrium iteration.**
   `SRC/element/Frame/Force/ForceFrame3d.h` carries its own `tol`/`max_iter`
   for internal "local iterations" — a force-based beam-column element
   solves its own internal equilibrium (section forces consistent with
   element-end forces) *inside* the domain-level Newton loop. Genuinely
   different from every other element; Carapace treats this as its own
   milestone (M8) rather than a variant of the displacement-based element
   work.
6. **Some materials are recursive compositions, not leaves.**
   `SRC/material/wrapper/{ParallelMaterial,SeriesMaterial,MinMaxMaterial}.h`
   hold pointers to *other* `UniaxialMaterial` instances — a flat material
   enum can't represent these directly.
7. **A real wasm-specific toolchain bug.** An `f2c`/`wasm-ld` function-signature
   mismatch (in `s_copy`) only trapped under the wasm build, not natively.
   Debugging it is why Carapace's own verification order is "native `cargo
   test` first, `wasm32-unknown-unknown` second, always" (implementation-plan
   §8) — isolating logic bugs from platform/toolchain bugs is the discipline
   that found this one, and it'll catch the Rust equivalent (a numeric or
   codegen difference between targets) faster than debugging inside a wasm
   runtime directly.

None of this is a knock on Xara — it's 90s/2000s C++ architected for
arbitrary extensibility, MPI parallelism, and database persistence, none of
which a single-threaded browser worker needs. Carapace's architecture
(enum-dispatched closed catalogs, typestate-enforced analysis wiring, no
broker/serialization layer — see `docs/implementation-plan.md` §3) is
designed around what's actually being built, not what OpenSees needed to be
in 1997.

## Key Xara source files referenced during design

For direct inspection, if revisiting any of the above:

- `SRC/runtime/runtime/BasicAnalysisBuilder.cpp` — the real analysis
  composition sequence (`setLinks`/`number`/`domainChanged`/`initialize`/
  `analyzeStatic`) that Carapace's typestate `AnalysisBuilder` is modeled on.
- `SRC/runtime/runtime/{ModelRegistry,ProcessContext,XaraClassBroker}.h/.cpp`
  — the broker eager-construction problem (#1 above).
- `SRC/element/Frame/Force/ForceFrame3d.h` — nested iteration evidence for M8.
- `SRC/material/wrapper/{ParallelMaterial,SeriesMaterial,MinMaxMaterial}.h`
  — recursive material composition evidence for M7.
- `SRC/matrix/routines/cmx.h` — the pre-existing missing-implementation bug
  (#4 above; not a Carapace concern, just documented context).
- `OTHER/ARPACK/` (legacy tree) vs `OTHER/arpack-ng/` — the f2c-compatibility
  distinction relevant to the eigensolver decision (implementation-plan §6).
