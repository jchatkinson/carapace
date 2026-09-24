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

**As implemented** (`Domain`, `Analysis`, `AnalysisBuilder`, `Integrator`,
`LoadPattern` in `core/src/model`/`core/src/analysis`): one generic
implementation per type, parameterized by const generics (`NDIM`, `NDOF`,
`ELEMENT_DOF`) and the node/element types, not a sealed space-marker type
parameter or a pair of hand-duplicated concrete structs. `Domain3`/
`Analysis3` are plain type aliases for the spatial instantiation
(`Node3Id`/`Element3`); bare `Domain`/`Analysis` default to the planar one,
so existing planar code is source-compatible. An `ElementOps<NDIM, NDOF,
ELEMENT_DOF, NId>` trait (`core/src/model/elements/mod.rs`) is what lets
`Domain`'s assembly/state bookkeeping stay a single implementation: it
abstracts exactly the shape that bookkeeping needs from an element catalog
(nodes, tangent/mass/load-vector formation, commit), without needing either
a shared enum across profiles or `Box<dyn Trait>`. `Element`/`Element3`
themselves — and every concrete formulation inside them (`Truss` vs
`Truss3`, and any future `DispBeamColumn3`) — stay separate concrete types
implementing that trait; their actual formulas differ (different direction-
cosine counts, biaxial bending, torsion) in ways no amount of genericity
unifies, so this is real, irreducible duplication, not an oversight.
`TransientAnalysis`/`modal_analysis` have **not** been generalized this way
yet — they still only accept the planar `Domain` (see the M15 status in
`implementation-plan.md`).

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

- **Done.** `NDF`/`ELEMENT_DOF` are split by profile via const generics
  (`NDOF`/`ELEMENT_DOF` on `Domain`/`Analysis`/etc., not a runtime field);
  node coordinates, fixities, mass, state, and equation maps are
  profile-sized fixed arrays (`Node<NDIM, NDOF>`).
- **Done** for static assembly/scatter/gather (`Domain::assemble_*`,
  `apply_displacement_increment`, `gather_*`, `scatter_state` are one
  generic implementation, profile-aware via `NDOF`/`ELEMENT_DOF`). **Not
  done** for ground-motion incidence/recorder layout — `TransientAnalysis`
  and `modal_analysis` haven't been generalized yet (see below).
- **Done.** Spatial nodal loads have six components
  (`Domain::add_nodal_load`/`load_node` take a `dof` index up to `NDOF`,
  generic over profile).
- **Not done.** Element loads are still the planar `ElementLoad` enum
  (`UniformTransverse`); the spatial profile's `ElementOps::Load` is
  `Infallible` (uninhabited — no spatial element load exists yet, so
  `add_element_load` is uncallable there until a real one lands, e.g. a
  beam's local `wy`/`wz`).
- **Not done.** `TransientAnalysis`/`modal_analysis` still only accept the
  planar `Domain`; three-translation-plus-rotational-inertia state handling
  for the spatial profile hasn't been built.

### Elements and transforms

- **Done.** `Truss3`: three direction cosines, 12-slot vectors with
  rotational rows zero, and mass over three translations, verified through
  the full `Domain3`/`Analysis3` stack (`core/tests/m15_spatial_truss.rs`),
  not just in isolation.
- **Done.** `ZeroLength3`: independent per-DOF materials along the six
  global directions `[ux, uy, uz, rx, ry, rz]`, no orientation vectors —
  the same simplification `ZeroLength` already makes relative to OpenSees'
  general `ZeroLength`. Arbitrarily oriented local springs (an explicit
  orientation transform, matching OpenSees' arbitrary-direction materials)
  remain a later feature.
- **Done.** `ElasticBeamColumn3`: a 12×12 Euler-Bernoulli frame with `E, G,
  A, J, Iy, Iz`, holding a `GeomTransf3` (`Linear3` or `PDelta3`). Verified
  against Xara's `ElasticBeam3d`/`LinearCrdTransf3d` (basic-stiffness
  terms, local DOF order, local-axis construction from `vec_xz`) and
  against textbook closed-form axial/biaxial-bending/torsion cantilever
  deflections, both at the element level (`core/src/model/elements/
  elastic_beam_column.rs`'s `spatial_tests`) and through the full
  `Domain3`/`Analysis3` stack (`core/tests/m16_elastic_beam3.rs`).
- **Done.** `GeomTransf3` (`core/src/model/transform.rs`): a first-class 3D
  geometric-transform type, the spatial counterpart of `GeomTransf` — this
  doc's `x_local`/`y_local`/`z_local` construction from `vec_xz` and the
  block-diagonal local↔global rotation live here, not inlined per element,
  so `DispBeamColumn3`/`ForceBeamColumn3` can reuse the same `Linear3`/
  `PDelta3` variants once they exist (unlike planar `GeomTransf`, which
  stays a behavior-only flag with each planar element inlining its own
  cheap 2×2 rotation — see `GeomTransf`'s doc comment for why the spatial
  case can't get away with that).
- **Done.** Spatial P-Delta (`GeomTransf3::PDelta3`): the analogue of
  planar `GeomTransf::PDelta` — a linearized geometric-stiffness correction
  from axial force, applied to *both* bending planes independently
  (`ElasticBeamColumn3::geometric_stiffness`). The `Iy`-plane block is a
  hand-derived sign flip of the `Iz`-plane one (the standard consistent
  geometric-stiffness block doesn't depend on `EI`, only on axial force and
  the shape functions, but `ry = -dw/dx` while `rz = +dv/dx` flips its
  off-diagonal terms the same way `local_elastic_stiffness`'s bending
  blocks flip) — verified, not just asserted, by exercising both planes in
  `core/tests/m16_pdelta3.rs` (zero-axial-force-matches-`Linear3`, and
  compressive-softening, in each plane independently). Corotational
  geometry remains its own, later milestone (M18).
- **Done.** `DispBeamColumn3` and `ForceBeamColumn3`: the spatial (biaxial)
  counterparts of `DispBeamColumn`/`ForceBeamColumn`, built directly on
  `FiberSection3`. Both are `GeomTransf3::Linear3`-only (a `vec_xz` field,
  not a full `GeomTransf3`) — same scope limit as their planar counterparts,
  for the same reason (a fiber section's state isn't a single scalar axial
  force, so it doesn't plug into the closed-form geometric-stiffness formula
  directly). Torsion is excluded from both the fiber loop and (for
  `ForceBeamColumn3`) the basic system entirely, folded in afterward as a
  decoupled elastic `G*J/L` term — see "Fiber sections" below.
  `ForceBeamColumn3`'s basic system (`v1` axial, `v2`/`v3` z-bending
  chord-relative rotations, `v4`/`v5` y-bending chord-relative rotations) is
  the direct biaxial generalization of `ForceBeamColumn`'s, with `v4`/`v5`'s
  rotation-coefficient sign flipped relative to `v2`/`v3`'s (the same
  `ry = -dw/dx` vs `rz = +dv/dx` substitution used throughout the spatial
  bending code) — verified, not just asserted by analogy, against
  `ElasticBeamColumn3`'s exact closed-form biaxial elastic stiffness in
  `core/tests/m17_disp_beam_column3.rs` and
  `core/tests/m17_force_beam_column3.rs`, plus a doubly-symmetric
  `ElasticPP` section yielding/unloading independently in both bending
  planes for `ForceBeamColumn3`.

### Fiber sections

**Done.** `Fiber3 { y, z, area, material }` and `FiberSection3` (biaxial):
strain is `eps0 - y*kappa_z + z*kappa_y`, section response is `[N, Mz, My]`
with a 3×3 tangent including product-of-inertia coupling from unsymmetric
layouts — verified against Xara's `FiberSection3d::addFiber`/
`getStressResultant` (`matData`/`sData` layout, `strain = e0 - y*e1 + z*e2`,
`N = sum(stress*A)`, `Mz = sum(-y*stress*A)`, `My = sum(z*stress*A)`) and,
independently, a four-corner-fiber elastic-tensor check
(`core/src/model/fiber_section.rs`'s
`four_corner_elastic_fibers_reproduce_ea_eiy_eiz_with_zero_cross_coupling`).
Torsion does not arise from uniaxial fibers, as planned: `DispBeamColumn3`/
`ForceBeamColumn3` supply explicit `G`, `J` fields and add torsion as a
decoupled elastic term, rather than the section carrying a torsion material.

`FiberSection` (planar) is unchanged — no meaningless `z = 0` coordinate was
added to it.

### Constraints

**Partially done.** `Domain::equal_dof` (identity ties on selected DOFs) is
already generic and works unchanged for `Domain3` — the reasoning below
about *why* a real spatial rigid diaphragm can't reuse this still holds, so
`rigid_diaphragm` itself stays planar-only (a separate, non-generic `impl
Domain` block; see `Domain`'s doc comment in `core/src/model/domain.rs`).
The transformation matrix described below is **not started**.

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

**Done.** The generic profile architecture (not a sealed marker type — see
"Decision" above), `Node3`, `Domain3`, profile-aware static
assembly/state/load handling, `Truss3`, `ZeroLength3`, six-DOF nodal loads,
identity `equal_dof`. Verified: a symmetric space-truss tripod against
closed-form vertical stiffness through the full `Domain3`/`Analysis3` stack
(`core/tests/m15_spatial_truss.rs`), a skew-truss direction-cosine check at
the `Domain3` assembly level (`core/src/model/domain.rs`'s
`domain3_tests`), element-level `Truss3`/`ZeroLength3` tests, and unchanged
results for every existing planar test.

Not done: spatial modal/transient plumbing (`TransientAnalysis`/
`modal_analysis` still only accept the planar `Domain`), and a
ground-motion acceptance case (no spatial transient support to verify it
against yet).

### M16 — Spatial elastic frame and constraints

**Partially done.** `GeomTransf3` (`Linear3` and `PDelta3`,
`core/src/model/transform.rs`) and `ElasticBeamColumn3` (Euler-Bernoulli,
holding one) are implemented and verified against Xara's `ElasticBeam3d`/
`LinearCrdTransf3d` and closed-form cantilever deflections (`transform.rs`'s
own tests, `core/src/model/elements/elastic_beam_column.rs`'s
`spatial_tests`, `core/tests/m16_elastic_beam3.rs`,
`core/tests/m16_pdelta3.rs`) — "axial, torsional, and both bending
deflections of cantilevers" from the acceptance list below, on an
axis-aligned member (isolating the local-stiffness formula) and a skew one
(isolating the transform), plus P-Delta's zero-axial-force/compressive-
softening checks in both bending planes.

Not done: local `wy`/`wz` element loads and the true spatial diaphragm
transformation (`Domain::equal_dof` on `Domain3` only does the identity
ties it already could; see "Constraints" above). Not yet verified: rigid
rotation invariance, diaphragm lever-arm kinematics, or diagnostics for
invalid orientation hints (`vec_xz` parallel to the member axis, or a
zero-length member, currently just panics via a plain `assert!` in
`GeomTransf3::local_axes` — not the structured `AnalysisError`-style
diagnostic a real API would want).

### M17 — Spatial fiber and nonlinear frames

**Done.** `FiberSection3`, `DispBeamColumn3`, and `ForceBeamColumn3`,
biaxial axial/bending with elastic (decoupled) `G*J` torsion — nonlinear
torsion deferred until a real model requires it, as planned. Verified:
elastic equivalence to `ElasticBeamColumn3` in both bending planes plus
axial and torsion (`core/tests/m17_disp_beam_column3.rs`,
`core/tests/m17_force_beam_column3.rs`), and — for `ForceBeamColumn3` — a
doubly-symmetric `ElasticPP` section yielding independently under a
`y`-only or `z`-only tip load and showing permanent set on unload
(`core/tests/m17_force_beam_column3.rs`'s `elastic_pp_section_softens_
past_yield_in_both_planes_and_shows_permanent_set_on_unload`), the spatial
analogue of `m8_force_beam_column.rs`'s equivalent planar test.

Not done: asymmetric (non-doubly-symmetric) biaxial section coupling isn't
exercised by a dedicated test the way `m8_force_beam_column.rs`'s
`asymmetric_elastic_section_matches_hand_derived_closed_form` covers the
planar element — the four-corner and doubly-symmetric-eight-fiber sections
used so far have zero product-of-inertia coupling by construction.

### M18 — Spatial corotational geometry, if needed

Revisit only after M15–M17. Planar and spatial corotational transforms need
different rotation parameterizations and independent objectivity tests.

No M15+ implementation belongs inside M10. M10 must cleanly reject spatial
input, but its public handoff contract must preserve the profile discriminant.
