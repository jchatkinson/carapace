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
| 3D shells: `Shell4` (MITC4) and `Shell3` (DKT/Allman triangle), `ShellSection` (`elasticMembranePlate`), `Element3` surface variants | **Done** |
| Shell loads (`shellPressure`, `shellBody`), lumped mass, 8-component Gauss-point resultant recorders, mixed shell and beam models | **Done** |
| Wire: `shellSections`, `shell3s`, `shell4s`; OpenSees cross-check (`ShellMITC4`, `ShellDKGT`) | **Done** |
| Zero-mass DOFs in transient analysis | **Done** (massless rows get zero initial acceleration; no shell-specific transient test yet) |
| Shell stiffness modifiers (`f11 f22 f12 m11 m22 m12 v13 v23`) on `elasticMembranePlate` | **Done** (core and wire; PySees mirror pending) |

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

## What exists (shells)

3D only. `Domain3` accepts them; the 2D profile rejects the tables.

- **`Shell4`** (`elements/shell4.rs`): the OpenSees `ShellMITC4` formulation, 24 DOFs. Membrane, bending, MITC4
  assumed transverse shear with the 5/6 correction, and a drilling penalty, reproduced term by term.
- **`Shell3`** (`elements/shell3.rs`): the OpenSees `ShellDKGT` formulation, 18 DOFs. DKT bending plus Allman membrane,
  4-point rule `(1/3,1/3,1/3)` with weight `-9/16` and the `(1/5,3/5,1/5)` permutations with `25/48`. No drilling
  penalty and no transverse shear (`gamma_xz`, `gamma_yz`, `Qx`, `Qy` are zero).
- **`ShellSection`** (`materials/shell_section.rs`): `elasticMembranePlate` (`e`, `nu`, `h`, `rho`, optional
  `ShellModifiers`: ETABS-style factors `f11 f22 f12 m11 m22 m12 v13 v23` on the individual terms of the tangent; see
  [input-format.md](input-format.md) for the rules and the OpenSees equivalence). The trial/commit
  seam for the section, one copy per Gauss point, with OpenSees' generalized-strain ordering and sign convention.
  It is not a `PlaneMaterial`, which has no bending.
- **Local frame** follows OpenSees: `e1` from the edge midpoints, `e3 = e1 x e2` (node order sets the sign). A quad is
  treated as flat in that frame; warp is tolerated as OpenSees does, so there is no flatness check. A triangle's local
  x axis is the edge from node 0 to node 1.
- **Interface**: shells go through the sink path like `Quad4`; `Element3::nodes()` returns a `NodeList`. Every node
  carries all six DOFs, so a beam can share a shell node with no extra constraint.
  `ElementOps::shell_responses` reports eight strains and resultants per Gauss point and `gauss_component_count`
  bounds the recorder component.
- **Loads**: `shellPressure` (per unit area, along `e3`) and `shellBody` (body acceleration, force per area
  `rho h b`), both consistent nodal forces on translations only. A shell takes no `uniform` load. Nodal loads and
  self-weight work as for any element.
- **Mass**: row-sum lumped translational mass `rho h`; rotations carry none. OpenSees uses a consistent mass, so
  modal frequencies agree with OpenSees only approximately (closed-form plate frequencies are the check).
- **Recorders**: `gaussPoint` on a shell has 4 points and 8 components per quantity. Strain is
  `[eps_x, eps_y, gamma_xy, kappa_x, kappa_y, 2 kappa_xy, gamma_xz, gamma_yz]`; stress is the resultants
  `[Nx, Ny, Nxy, Mx, My, Mxy, Qx, Qy]` in the local axes, in OpenSees' order and sign. The force recorder returns the
  24 (quad) or 18 (triangle) global nodal forces and moments. `principalStrain`/`principalStress` record nothing for a
  shell.
- **Analysis**: static, modal (massless rotations are condensed) and transient (massless
  equilibrium rows start at zero acceleration).
- **Wire**: `shellSections`, `shell3s`, `shell4s`, shell element loads and shell Gauss-point recorders; see
  [input-format.md](input-format.md).

### OpenSees differences to know

- `ShellMITC4` has no pressure load and applies `-selfWeight` with the opposite sign; `ShellDKGT` applies twice the
  self-weight (an integration-weight bug). Carapace applies the correct weight, so exporters expand pressure, and a
  triangle's self-weight, into nodal loads and negate the quad's body vector.
- Lumped versus consistent mass (above).
- OpenSees' bending resultants have the opposite sign to `M = D kappa`; Carapace reports them as OpenSees does.

### Verification

`core/tests/shells/shell4.rs` and `shell3.rs` (tolerance 1e-8 relative to the largest value against OpenSeesPy):

- Distorted, tilted quad and triangle patches: displacements, clamp reactions and Gauss-point resultants.
- Self-weight against OpenSees (sign and the `ShellDKGT` double weight), and pressure resultant `p A e3`.
- Simply supported square plate: centre deflection `0.004062 p a^4 / D` and fundamental frequency
  `2 pi^2 / a^2 sqrt(D / (rho h))` (closed form, since the mass differs from OpenSees).
- Warped quad against OpenSees, a column framing into a shell corner, and the Scordelis-Lo barrel vault
  (8x8 mesh against OpenSees; about 0.30 against the published 0.3024).

`core/tests/shells/modifiers.rs` checks the modified tangent and resultants term by term, the drilling stiffness
(including extreme factors), validation, and exact element scalings (bending factors on `Shell3`, a uniform factor
on `Shell4` equal to scaling `E`, the `f11` axial stiffness, shear factors acting on `Shell4` only). There is no
OpenSees reference for non-uniform modifiers.

`wasm-bridge/tests/wire/shells.rs` covers the wire tables end to end, serde round trips and decode errors
(including `invalidShellSection`).

### PySees

`compileInputV1.ts` now has shell support (a `shell` recorder kind with 24 force labels). Treat the PySees repo as the
source of truth for what it compiles and exports.

## Deferred

- Nonlinear or layered shell sections, geometric nonlinearity for shells.
- Consistent mass for shells.
- A flatness (warp) check, and non-convex quads.
- A shell-specific transient test.
- 3D membranes.
- The PySees compiler and type for shell stiffness modifiers.
- 2D surface support in the PySees compiler (`tri3`/`quad`).
