# Carapace

A from-scratch Rust finite-element engine compiled to WebAssembly, covering the subset of OpenSees/Xara
features that its frontend, [PySees](../pysees), needs. Runs in a browser Web Worker. Not a port of OpenSees.
Feature list: `README.md`.

## Layout

Cargo workspace of two crates:

- `core/` (`carapace-core`): the native engine. `model/` (domain, nodes, elements, materials, sections,
  constraints, load patterns) and `analysis/` (solvers, integrators, algorithms, arc-length, modal, transient).
  Plain Rust, no wasm dependencies. Tests in `core/tests/`, benchmarks in `core/examples/`.
- `wasm-bridge/` (`carapace-wasm`): the wasm-bindgen boundary. `input_v1/` decodes the JSON wire format
  (`CarapaceInputV1`) into a core `Domain` and analysis sequence; also `material_probe.rs` (unit zero-length
  spring for PySees' Material Preview) and result/recorder streaming. Tests in `wasm-bridge/tests/`.

Other directories: `docs/` (design and wire-format docs), `comparison/` (Python harness vs real OpenSeesPy,
uses `comparison/.venv`), `pkg/` and `pkg-web/` (generated wasm-bindgen output, gitignored).

## Docs to read first

- `docs/architecture.md`: Domain, element interface, DOF activation, constraints, how to add an element or a
  wire table. Follow its checklists.
- `docs/input-format.md`: the wire format. Keep it in sync with any change to `input_v1`.
- `docs/algorithms.md`, `docs/arclength.md`: solution algorithms.

## Rules

- 2D (ndm=2, 3 DOF) and 3D (ndm=3, 6 DOF) share one const-generic implementation (`Domain` / `Domain3`).
  Share assembly, numbering, constraints and the analysis layer; split leaf elements whose kinematics differ.
  See "When to share and when to split" in `architecture.md`.
- Element catalogs are closed enums: no `dyn`, no heap allocation in the assembly loop.
- `dof_mask` must declare exactly the slots an element stiffens. Never read state through `equation_of`; use
  `Domain::dof_terms` / `Domain::value_at`.
- Unsupported or invalid input must surface as an error (`ModelError`, `DecodeError`), never be silently
  dropped or degraded. PySees turns these into user-facing diagnostics.
- Any change to the wire format (`input_v1`) must be mirrored by PySees: its `CarapaceInputV1` type
  (`../pysees/src/app/types/carapaceInputV1.ts`) and compiler (`src/app/lib/carapace/compileInputV1.ts`).
- Results must agree with OpenSees. Cross-check new elements and materials against OpenSeesPy where one
  exists; PySees also runs `carapaceVsOpenSees.test.ts` against the bundled wasm.
- 2D/3D, not "planar"/"spatial", in prose and docs. (Some existing file names use `spatial`; leave them.)
- Keep code compact and match the surrounding style. Run `cargo fmt` and `cargo clippy`.

## Commands

```
cargo test                                   # all native tests (core + wasm-bridge)
cargo test -p carapace-core --test elements   # one integration test binary (elements, analysis, dynamics, constraints, continuum, materials; wasm-bridge: wire)
cargo build -p carapace-wasm --target wasm32-unknown-unknown --release
```

Integration tests are one binary per area (`tests/<area>/main.rs` plus a module per feature):
`core/tests/{elements,analysis,dynamics,constraints,continuum,materials}` and `wasm-bridge/tests/wire`
(shared builders in `wire/common.rs`). Add a test to the module for its feature and give 2D and 3D cases
to the same module. Each file's `//!` header says what it verifies and against which reference.

The wasm bundle PySees uses is built from here by `../pysees/scripts/build-carapace.sh` (cargo build, then
`wasm-bindgen --target web`, then copied into `src/app/carapace/wasm`). Rebuild and copy it after any
change PySees should see. `comparison/check_refactor.sh` checks numerical output against recorded baselines
(`comparison/baseline/`) and should show no drift after a refactor.
