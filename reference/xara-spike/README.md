# Xara wasm spike (reference material)

Source-level artifacts from the C++ investigation that led to Carapace's
design (see `docs/implementation-plan.md` §2, §6, §7, §9). Build artifacts
(`.o` files, logs) were intentionally not copied — regenerate them by
following the steps in the implementation plan / conversation history if
needed; they're not useful to keep around.

- **`arpack_spike/`** — proved `f2c`-translating Xara's *legacy* `OTHER/ARPACK`
  tree (not `arpack-ng`, which fails to translate — see implementation plan
  §6) to WebAssembly produces correct results for `dsaupd`/`dseupd`, verified
  against the known eigenvalues of a diagonal test operator.
  - `src/test_dsaupd.cpp` — the verification harness.
  - `filelist_legacy.txt` — the resolved 38-routine BLAS/LAPACK/ARPACK
    dependency closure, traced from `dsaupd`/`dseupd`'s Fortran call graph.
  - `f2c_runtime_lite/` — the ~7-function minimal libf2c runtime subset
    actually needed (not the full 160-file libf2c), plus `io_stubs.c`, hand
    -written stand-ins for the 3 formatted-I/O symbols ARPACK's debug-print
    paths reference at link time but never execute at default verbosity.

- **`poc/`** — the full Truss+ElasticMaterial static-analysis wasm POC, using
  Xara's actual production C++ classes (not stubs), hand-wired to bypass the
  `XaraClassBroker` eager-construction problem (implementation plan §2).
  - `src/main.cpp` — the driver. This is what `docs/implementation-plan.md`
    §5.4's typestate `AnalysisBuilder` sketch is informed by: it replicates
    `BasicAnalysisBuilder`'s real `setLinks`/`number`/`domainChanged`/
    `initialize`/`analyzeStatic` sequence by hand.
  - `patched/BasicAnalysisBuilder.{h,cpp}` and `ProcessContext.h` — **unmodified
    copies of the real Xara files**, kept here as a stable reference for the
    composition sequence, since `main.cpp` deliberately bypasses using this
    class directly (see implementation plan §2, problem #1).
  - `patched/cmx.cpp` — a Gauss-Jordan stand-in for Xara's missing
    `cmx_inv2..6` implementation (a real, pre-existing Xara bug, unrelated to
    wasm — see implementation plan §2, problem #4).
  - `patched/logging.cpp` — Xara's `SRC/logging/logging.cpp` with its
    unnecessary `#include <tcl.h>` removed (unused in that file).
  - `cpplist2.txt` — the resolved 97-file dependency closure for
    Domain/Node/Truss/ElasticMaterial + the hand-wired analysis chain,
    traced by header-closure BFS *after* excluding `ModelRegistry`/
    `BasicAnalysisBuilder`/`XaraClassBroker` (which alone would balloon this
    to 348 files — see implementation plan §2).
  - `lapack_extra.txt` — the additional LAPACK/BLAS closure needed for
    `Matrix::Solve()`/`Invert()` (`dgesv`/`dgetrf`/`dgetri` + deps), beyond
    what `BandGenLinLapackSolver`'s `dgbsv` closure alone requires.
  - `incdirs_raw.txt` — the `-I` search path used to compile against Xara's
    `SRC/` tree directly (classic flat-include style), excluding
    reliability/mpm/PFEM/executable/runtime-interpreter directories.

**Verified results at time of writing** (native and wasm32, both matching):
the 2-node truss case (L=100, A=2, E=30000, P=50) computed `dx =
0.0833333333`, matching `δ = PL/(AE)` exactly — this is the oracle value
Carapace's M0/M1 milestones check against.
