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

## Status

Milestones M0–M7 are done: `Domain`/`Node`, the full element catalog through
`DispBeamColumn` (`Truss`, `ZeroLength`, `ElasticBeamColumn`,
`DispBeamColumn`), the full standard uniaxial material catalog (`Elastic`,
`ElasticPP`, `Gap`, `Ent`, `Steel01`/`02`, `Concrete01`/`02`, `Hysteretic`,
`Pinching4`, plus the `Parallel`/`Series`/`MinMax` composites), a
typestate-composed `Analysis` (`LoadControl`/`DisplacementControl` ×
`Linear`/`NewtonRaphson`), modal analysis (shift-invert Lanczos), and
`TransientAnalysis` (Newmark-beta time-history with Rayleigh damping) are
all in place and verified on both native and `wasm32` + Node targets, backed
by a real sparse solver (`faer`'s sparse LU).

Remaining work: `ForceBeamColumn` (M8, nested element-level equilibrium
iteration), optional corotational geometry (M9), and the JS/TS worker
integration + results persistence that let `pysees` actually call Carapace
(M10–M11).

3D (spatial) support is in progress, separately from M8–M11: `Domain`,
`Analysis`, `AnalysisBuilder`, `Integrator`, and `LoadPattern` are each one
generic implementation covering both a planar and a spatial (six-DOF-per-
node) profile, with `Domain3`/`Analysis3` as the spatial instantiation.
`Truss3`, `ZeroLength3`, `ElasticBeamColumn3` (linear-elastic 3D
Euler-Bernoulli, axial/torsion/biaxial bending, `Linear3` or `PDelta3`
geometry via the reusable `GeomTransf3` transform), and the nonlinear
fiber-section pair `DispBeamColumn3`/`ForceBeamColumn3` (biaxial bending via
a new `FiberSection3`, decoupled elastic `G*J` torsion) are implemented and
verified end to end (`Domain3` → `Analysis3::step`) against closed-form
cases, elastic equivalence to `ElasticBeamColumn3` in both bending planes,
and, for the beam, against Xara's `ElasticBeam3d`/`FiberSection3d` reference
formulas — see `core/tests/m15_spatial_truss.rs`,
`core/tests/m16_elastic_beam3.rs`, `core/tests/m16_pdelta3.rs`,
`core/tests/m17_disp_beam_column3.rs`, and
`core/tests/m17_force_beam_column3.rs`.
A true spatial rigid diaphragm and spatial modal/transient analysis are not
yet implemented. See
[`docs/spatial-architecture.md`](docs/spatial-architecture.md) for the full
status and remaining milestones (M15–M18).

See [`docs/implementation-plan.md`](docs/implementation-plan.md) for the
full design, scope, and milestone-by-milestone plan.

## Workspace layout

```
carapace/
├── core/           # carapace-core: all FE logic. No wasm-bindgen dependency —
│                    compiles and tests identically on native and wasm32 targets.
├── wasm-bridge/     # carapace-wasm: thin wasm-bindgen layer exposing core's API
│                    to JS. Should contain no numerical logic of its own.
└── docs/
    ├── implementation-plan.md   # comprehensive design doc + milestones — read this first
    ├── pysees-handoff.md        # pysees model/sequence → worker contract for M10/M11
    └── spatial-architecture.md  # 3D execution-profile architecture and milestones M15–M18
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

**`docs/implementation-plan.md`** — architecture principles and rationale, full
feature scope, milestone-by-milestone plan with acceptance criteria, and open
decisions that need resolving early (linear algebra crate, ARPACK strategy).
Any new coding session on this repo should start there.

For the frontend/worker boundary, then read
**[`docs/pysees-handoff.md`](docs/pysees-handoff.md)** — the required pysees
model and analysis-sequence contract, run lifecycle, transport, and results
database design for M10/M11.

For the 3D (spatial) profile — in progress, see Status above — read
**[`docs/spatial-architecture.md`](docs/spatial-architecture.md)** before
changing node DOFs, element transformations, fiber sections, or constraints.
