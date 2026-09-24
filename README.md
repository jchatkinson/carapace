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

- [x] M0–M8: Planar analysis core.
- [x] M9: Planar co-rotational transform implementation.
- [ ] M9: Dedicated co-rotational acceptance tests.
- [ ] M10: pysees handoff and worker execution.
- [ ] M11: Results database hardening.
- [ ] M12–M13: Algorithm extensions and event-to-event stepping.
- [x] M14: Load patterns, multi-phase analysis, and ground motion.
- [x] M15: Spatial modal and transient analysis.
- [x] M16: Spatial beam loads and diaphragm constraints.
- [ ] M17: Asymmetric biaxial section coverage.
- [ ] M18 (optional): Spatial co-rotational geometry.

See [`docs/implementation-plan.md`](docs/implementation-plan.md) and
[`docs/spatial-architecture.md`](docs/spatial-architecture.md) for milestone
details.

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
