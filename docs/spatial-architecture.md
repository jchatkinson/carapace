# Spatial (3D) Architecture Plan

Carapace is currently a planar frame engine: global coordinates `[x, y]`,
node DOFs `[ux, uy, rz]`, and two-node element matrices of order 6. `pysees`
also authors spatial models, with global **y-up** coordinates. This document
defines how Carapace grows to full spatial analysis without turning planar
analysis into a slow, ambiguous special case.

## Decision: two execution profiles, selected once

The public model declares one of two supported spaces:

```text
Planar:  NDM = 2, NDF = 3, coordinates [x, y], DOFs [ux, uy, rz]
Spatial: NDM = 3, NDF = 6, coordinates [x, y, z],
         DOFs [ux, uy, uz, rx, ry, rz]
```

`y` is the global vertical/up axis in both profiles. It is a global authoring
and rendering convention, not an element-local-axis rule.

Use distinct concrete `Domain2`/`Node2`/`Element2` and
`Domain3`/`Node3`/`Element3` families, with generic analysis/assembly plumbing
only where it genuinely removes duplication. `Analysis<S>` and
`TransientAnalysis<S>` are parameterized by a sealed space marker; the wasm
boundary selects the profile once and holds an enum of the two session types.
This preserves fixed-size hot-loop matrices:

```text
planar two-node element:  6 × 6
spatial two-node element: 12 × 12
```

Do not represent planar models as spatial models with restrained
out-of-plane DOFs, and do not move all element algebra to heap-sized matrices.
That would add artificial equations, conceal model mistakes, and regress the
browser workload optimized by the current implementation.

The sparse solver, load-series machinery, trial/commit material model,
analysis algorithms, and worker protocol are dimension-independent concepts.
Element kinematics, coordinate transforms, section response, and constraints
are not.

## Coordinate frames and orientation

All spatial input coordinates use global `[x, y, z]` with y-up. Every 3D frame
member has a right-handed local basis and must carry an explicit
`vec_xz: [f64; 3]` (or a named orientation resolving to one), matching
Xara/OpenSees's 3D `geomTransf` convention:

```text
x_local = normalize(p_j - p_i)
y_local = normalize(vec_xz × x_local)
z_local = x_local × y_local
```

Reject a zero-length member or a `vec_xz` parallel to its axis. Never silently
choose an alternate roll direction: roll is mechanically meaningful for
unsymmetric sections, local loads, and recorder output. Local x is along the
member; local y/z are section axes and are unrelated to global y-up.

For planar members, local x follows the chord in global x-y and local
`y = [-x.y, x.x]`, giving positive local z out of plane. Existing planar signs
remain unchanged.

## Required core changes

### Nodes, assembly, state, and loads

- Split global `NDF = 3` and `ELEMENT_DOF = 6` by profile. Node coordinates,
  fixities, mass, state, and equation maps become profile-sized fixed arrays.
- Keep sparse assembly and the solver shared, but make scatter/gather, mass,
  ground-motion incidence, and recorder layout profile-aware.
- Spatial nodal loads have six components. Uniform excitation applies in
  translational directions 0–2; it does not imply rotational excitation.
- Generalize element loads to named local components. The first spatial beam
  load is uniform local `wy`/`wz`, never a reinterpretation of planar
  `UniformTransverse`.
- Modal/transient algorithms remain shared, but their state gathering and
  lumped mass must handle three translations and explicit rotational inertia.

### Elements and transforms

- `Truss3`: three direction cosines, 12-slot vectors with rotational rows
  zero, and mass over three translations. It is the first spatial element.
- `ZeroLength3`: six global directions initially. Arbitrarily oriented local
  springs are a later explicit transform feature.
- `ElasticBeamColumn3`: a 12×12 Euler-Bernoulli frame with `E, A, G, J, Iy,
  Iz`; its local response includes axial force, torsion, and bending about
  both axes. Implement `Linear3` first; add spatial P-Delta only once force
  and tangent recovery are consistent.
- Split `GeomTransf` into planar/spatial transform types. Corotational
  geometry is correspondingly two implementations, not a shared switch.
- `DispBeamColumn3` and `ForceBeamColumn3` are new formulations. The spatial
  force-based basic system includes axial, torsion, and independent end
  rotations about both bending axes; it needs its own flexibility/state-
  determination implementation and reference verification.

### Fiber sections

The current `Fiber { y, area, material }` and `[N, Mz]` response are planar.
A spatial fiber is `{ y, z, area, material }`; its strain depends on centroid
strain and both curvatures. Its section response is `[N, My, Mz]` with a 3×3
tangent, including coupling from unsymmetric layouts. Torsion does not arise
from uniaxial fibers: supply explicit `GJ` (or a future torsion material) in
the spatial section type.

Keep compact `FiberSection2`; do not give planar sections a meaningless
`z = 0` coordinate or singular spatial tangent merely for API uniformity.

### Constraints

`equal_dof` generalizes to selected DOFs. The current transformation handler
is equation aliasing and is correct only for identity ties. It cannot implement
a true spatial rigid diaphragm: slave translation includes retained-node
rotation through the lever arm.

Before exposing spatial `rigidDiaphragm`, implement a sparse constraint
transformation `u_full = C * u_reduced` (or equivalent verified reduction).
Retain the identity-tie fast path for planar `equal_dof`. The spatial API names
the retained node and diaphragm normal, rejecting unsupported affine/mixed
constraints instead of silently aliasing them.

## pysees and wasm handoff

`Model` and `CarapaceInputV1` gain a `space` discriminant, valid only for
`(2, 3)` or `(3, 6)`. Spatial element records include `vec_xz`; spatial
fibers include local y/z. Result metadata records the profile and component
layout rather than assuming three response values per node.

M10's first runnable vertical slice remains planar. Its schema must reserve
the discriminant but report `unsupported_space` for Spatial until M15. The
worker and SQLite response-block schema must remain dimension-neutral.

## Delivery sequence and acceptance tests

### M15 — Spatial foundation

Add the sealed profile architecture; `Node3`, `Domain3`, profile-aware
assembly/state/load handling, `Truss3`, `ZeroLength3`, six-DOF nodal loads,
identity `equal_dof`, and spatial modal/transient plumbing. Verify a skew 3D
truss against analytical direction-cosine stiffness, a six-DOF spring system,
spatial mass/modal and ground-motion cases, and unchanged results for every
existing planar test.

### M16 — Spatial elastic frame and constraints

Add `Linear3`, `ElasticBeamColumn3`, local `wy`/`wz` loads, correct 3D lumped
mass, and a true spatial diaphragm transformation. Verify axial, torsional,
and both bending deflections of rotated cantilevers; rigid rotation invariance;
diaphragm lever-arm kinematics; and diagnostics for invalid orientation hints.

### M17 — Spatial fiber and nonlinear frames

Add `FiberSection3`, `DispBeamColumn3`, then `ForceBeamColumn3` as separately
verified sub-stages. Start with biaxial axial/bending and elastic `GJ`; add
nonlinear torsion only when a real model requires it. Verify elastic
equivalence to `ElasticBeamColumn3`, asymmetric biaxial coupling,
yielding/unloading under combined bending, and a reference force-based case.

### M18 — Spatial corotational geometry, if needed

Revisit only after M15–M17. Planar and spatial corotational transforms need
different rotation parameterizations and independent objectivity tests.

No M15+ implementation belongs inside M10. M10 must cleanly reject spatial
input, but its public handoff contract must preserve the profile discriminant.
