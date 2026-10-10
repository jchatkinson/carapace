# 2D membranes and 3D shells

## Status

| Piece | State |
|---|---|
| Variable-node assembly (`ElementOps`, sinks, `DofMask`, `DofTable`, DOF activation) | **Done** (see [architecture.md](architecture.md)) |
| 2D membranes: `Tri3` (constant strain) and `Quad4` (bilinear, `Full` or `Enhanced` Wilson-Taylor) | **Done** |
| `PlaneMaterial` (isotropic plane stress/strain, orthotropic, packed elastic `D`) | **Done** |
| Continuum element loads (body, edge traction, edge pressure), lumped mass, Gauss-point strain/stress recorders | **Done** |
| Wire: `planeMaterials`, `triangles`, `quads`; modal and transient verified with membranes; OpenSees cross-check | **Done** |
| PySees compile of `tri3`/`quad`, export of `quad` | Not done (PySees compiler has no surface support) |
| 3D shell (MITC4), `ElasticMembranePlateSection`, `Element3` surface variants | **Not started** |
| Zero-mass DOFs in transient analysis | Not done (modal already condenses them) |

Departures from the original plan: there is no separate `(ndm=2, ndf=2)` profile and no header `ndf`. Membranes run in
the ordinary 2D `Domain` (`ux, uy, rz`); `dof_mask` leaves `rz` inactive, so membrane-only nodes need no rotation fixed
and a beam and a membrane can share a node. The wire header carries only `ndm`; `ndf` is 3 in 2D and 6 in 3D.
Membranes are 2D only: `Domain3` rejects the tables.

## What exists (membranes)

- `model/continuum/`: `tri3`/`quad4` shape functions, Gauss rules (1 point for Tri3, 2x2 for Quad4), `Jacobian`,
  `physical_derivatives`, scale-relative `validate_triangle`/`validate_quad` (counter-clockwise, positive Jacobian at
  all corners).
- `PlaneMaterial` is a trial/commit seam on `[eps_x, eps_y, gamma_xy]` (engineering shear); each Gauss point owns a copy.
- `Quad4Formulation::Enhanced` (incompatible modes, condensed) accepts linear materials only and makes no
  volumetric-locking guarantee.
- `density` gives lumped translational mass; mass on a multi-term constraint slave is a `ModelError`.
- Recorders: `gaussPoint` with `elementKind` `tri3`/`quad4`, `quantity` strain/stress, `component` 0..3, static and
  transient stages only. Nodal force recording uses `local_force_width`.
- Wire layouts: see [input-format.md](input-format.md).

## What a 3D shell needs

Carapace has no 3D surface element. The pieces, in dependency order:

### 1. Catalog and interface changes

- Add shell variants to `Element3`. `Element3::nodes()` currently returns `[Node3Id; 2]`; it must return a `NodeList`
  like `Element::nodes()`, and the two-node `pair()` path must exclude the new variants.
- `ElementLoad3` is a single-variant enum. Decide whether shells take loads in v1 (recommended: nodal only, reject with
  `IncompatibleElementLoad` via `accepts_load`).
- `GaussResponse` is `{ strain: [f64; 3], stress: [f64; 3] }`, which fits plane elements only. A shell reports
  generalized strains and resultants (membrane `n_xx n_yy n_xy`, bending `m_xx m_yy m_xy`, shear `q_xz q_yz`: 8
  components). Make the response width element-defined (or add a separate `shell_responses`) and extend
  `gauss_point_count`/recorder validation to match.
- `MAX_ELEMENT_NODES = 4` already covers a quad shell; `DofMask::all(6)` and the 24-component `SMatrix` sink path need
  no interface change.

### 2. Formulation (`ShellMITC4`, 24 DOFs)

- Reuse the Quad4 shape functions, Gauss rule and `Jacobian`, but evaluate them in a local 2D frame because the nodes
  are 3D: build an orthonormal local basis from the element geometry (OpenSees `ShellMITC4` convention: normal from the
  node ordering) and project nodal coordinates into it.
- Membrane (plane-stress `D`, thickness `h`), bending (`D h^3/12`), MITC4 assumed transverse shear with the 5/6 shear
  correction, and a drilling (`rz`) stiffness derived from membrane stiffness (OpenSees penalty), all in the local frame;
  transform with a 24x24 block-diagonal rotation (translations and rotations) to global.
- Section: `ElasticMembranePlateSection` (`E`, `nu`, `h`, `rho`, optional plate-modulus multiplier). A new wire table; it
  is not a `PlaneMaterial`, which has no bending.
- Validation (structured `InvalidElement`): four distinct finite nodes, `h > 0`, `E > 0`, `-1 < nu < 0.5`, density `>=
  0`, convex, positive Jacobian, and **flatness** (warp tolerance relative to element size). Today's `validate_quad` is
  2D only; there is no warp check yet.
- Mass: translational lumped mass from density; rotations carry none.
- Local-axis conventions and Gauss-point order must be documented because PySees will draw resultants from them.

### 3. Analysis support

- Modal: already condenses massless DOFs (`modal_analysis` accepts non-negative mass with at least `num_modes` positive
  entries). Check that returned shapes give the rotational components sensible values, since shells make most
  rotations massless.
- Transient: `TransientAnalysis::new` returns `SingularSystem` if any free DOF has `mass <= 0`, and the initial
  acceleration divides by mass. Shell models fail here unless rotations are fixed or given artificial mass. Needs
  consistent initialization of the massless equilibrium rows (keep Newmark average acceleration).

### 4. Wire and bridge

- New tables (for example `sectionsShell` and `shells`) in `tables.rs`, decoded in `decode/elements_3d.rs`;
  `reject_tables` makes the 2D profile refuse them, and `Domain3` already refuses `planeMaterials`/`triangles`/`quads`.
- Extend `elementKind` and the `gaussPoint` recorder (currently `tri3`/`quad4`, 3 components) for shell resultants;
  `local_force_width` is 24 for the force recorder.
- Add the shell to the mask-conformance test in `elements/mod.rs`, update [input-format.md](input-format.md),
  regenerate TypeScript declarations and rebuild the wasm.

### 5. PySees

- `compileInputV1.ts` has no surface compilation, for membranes or shells. It needs: node/fix/load compile for 3D shells
  (Z-up), schema argument names for `ShellMITC4` and the section (the OpenSees schema already contains both),
  diagnostics for unsupported options, and element-specific recorder labels (the current ones assume beams/trusses).
- Rendering and result display (surface meshes, contours) are separate work.
- Export: `ShellMITC4`/`quad` must round-trip through `commands/decode.ts`/`encode.ts` and `tables.ts`, and the exports
  must be accepted by real OpenSees.

### 6. Verification

- Stiffness symmetry and exactly six rigid-body modes (zero strain energy under rigid translation and rotation).
- Membrane and bending patch tests; rotated and distorted valid quads; thin-shell locking checks (MITC4 should not
  lock); cantilever and simply supported plate under **out-of-plane** pressure against closed-form solutions
  (nodal-equivalent loads first, since element loads are deferred).
- Drilling DOF: no spurious energy, and the stiffness scale does not affect results.
- Compare displacements, reactions and Gauss-point resultants with OpenSees `ShellMITC4`; for dynamics, assign matching
  lumped masses in OpenSees (it uses consistent mass otherwise).
- Mixed shell and frame models (shared nodes carry all six DOFs).

## Deferred

Pressure and body-force loads on shells, triangular shells (`ShellDKGT`), warped or non-convex quads, nonlinear or
layered sections, geometric nonlinearity, 3D membranes, and surface rendering.
