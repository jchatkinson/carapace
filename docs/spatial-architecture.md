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
`TransientAnalysis`/`modal_analysis` are now generalized the same way (see
the M19 status below): `pub fn modal_analysis<const NDIM: usize, const NDOF:
usize, const ELEMENT_DOF: usize, NId, E>(domain: &mut Domain<NDIM, NDOF,
ELEMENT_DOF, NId, E>, ...)` and `TransientAnalysis<NDIM, NDOF, ELEMENT_DOF,
NId, E>` (with `TransientAnalysis3` the spatial type alias), inferred from
`domain`'s concrete type at each call site exactly like
`AnalysisBuilder::build`. Neither needed a single DOF-count-specific line:
Lanczos, Newmark, and ground-motion excitation only ever touch `Domain`
through its already-generic free-DOF interface.

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
  generic implementation, profile-aware via `NDOF`/`ELEMENT_DOF`). **Done**
  for ground-motion incidence too — `direction_incidence` was already
  generic; `GroundMotion::direction` in `0..NDOF` selects `ux`/`uy`/`uz` (and
  is simply out of the planar 0..3 range) with no change needed. **Not
  done**: recorder/result-layout code outside `core` (worker/pysees
  handoff — see "pysees and wasm handoff" below) still assumes a
  three-value-per-node planar layout.
- **Done.** Spatial nodal loads have six components
  (`Domain::add_nodal_load`/`load_node` take a `dof` index up to `NDOF`,
  generic over profile).
- **Done (M20).** `ElementLoad3` (`core/src/model/load_pattern.rs`) —
  `ElasticBeamColumn3`'s biaxial local `wy`/`wz` uniform transverse load,
  the spatial counterpart of planar `ElementLoad::UniformTransverse`, is
  now `ElementOps::Load` for the spatial profile (no longer the
  uninhabited `Infallible`). Scope matches the planar profile exactly:
  only `ElasticBeamColumn3` supports it (`DispBeamColumn3`/
  `ForceBeamColumn3` don't, same as their planar counterparts). Verified
  against Xara/OpenSees's `ElasticBeam3d::addLoad`'s `Beam3dUniformLoad`
  fixed-end-force derivation and through the full `Domain3`/`Analysis3`
  stack (`core/tests/m20_spatial_beam_loads.rs`): a simply-supported
  spatial beam under simultaneous biaxial UDL matches the same closed-form
  end rotation `theta = w*L^3/(24*E*I)` independently in both bending
  planes.
- **Done (M19).** `TransientAnalysis<NDIM, NDOF, ELEMENT_DOF, NId, E>` and
  `modal_analysis<const NDIM, const NDOF, const ELEMENT_DOF, NId, E>` are
  generic over the same kinematic profile as `Domain`/`Analysis`
  (`TransientAnalysis3` is the spatial type alias, mirroring `Analysis3`).
  Newmark stepping, the Lanczos eigensolver, initial-acceleration
  equilibrium, and ground-motion excitation only ever went through
  `Domain`'s already-generic free-DOF interface (mass diagonal,
  gather/scatter, tangent assembly, `direction_incidence`), so no
  DOF-count-specific logic needed writing — this was a mechanical
  generalization, not new dynamics code. Verified through the full
  `Domain3`/`Analysis3` stack: a spatial SDOF truss's natural frequency
  (`core/tests/m19_spatial_dynamics.rs`'s
  `spatial_truss_mass_matches_sdof_closed_form_frequency_via_modal_
  analysis3`), undamped and Rayleigh-damped spatial SDOF free vibration on
  two different translational DOFs, and two simultaneous, independent
  `GroundMotion`s (`ux` and `uz`) driving an uncoupled spatial mass with no
  cross-talk between directions — the spatial "3D ground motion" case.

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

**Done (M20).** `Domain::equal_dof` (identity ties on selected DOFs) stays
exactly as before — still generic, still resolved by plain equation
aliasing (the "fast path" the plan below calls for keeping). Planar
`Domain::rigid_diaphragm` is unchanged too (its own non-generic `impl
Domain` block, still an identity alias on `x` only).

A genuine affine multi-point constraint now exists alongside it:
`AffineConstraint<NId>` (`core/src/model/domain.rs`) expresses `u_c[dof] =
sum(coeff * u_r[retained_dof])` — a linear combination of *several*
retained dofs with real coefficients, not just an identity alias. Resolved
at `Domain::number_dofs` time into `dof_transform: HashMap<(NId, usize),
Vec<(usize, f64)>>` (retained-dof references become actual free-DOF
equation numbers there; a retained dof that's fixed contributes nothing
and is dropped, one that's itself another diaphragm's slave is rejected
with a panic — chained diaphragms aren't supported). Every assembly hot
loop (`assemble_stiffness_triplets`, `assemble_mass_diagonal`,
`assemble_load_with`) and state-scatter path (`apply_displacement_
increment`, `scatter_state`) reads a shared `dof_terms` helper instead of
`node.equation[dof]` directly, so ordinary free dofs, identity-tied dofs,
and affine-tied dofs are all handled by one code path — not three.
`assemble_mass_diagonal` is the one deliberate exception: it panics if
nonzero nodal/element mass ever lands on an affine-tied dof with more than
one term, since the lumped-mass *diagonal* it returns can't represent the
resulting rotational inertia at the retained node (a real limitation,
documented and tested, not silently wrong physics — see
`core/tests/m20_rigid_diaphragm3.rs`'s panic test).

`Domain3::rigid_diaphragm`/`rigid_diaphragm_about` (its own non-generic
`impl Domain3` block, spatial-only — the kinematics need real 3D
coordinates and a rotation-about-normal dof neither the trait bound nor
the planar profile has) build one `AffineConstraint` per constrained node:
its two in-plane translations (the two axes other than `normal`, an
`Axis3`, default `Axis3::Y` — this crate's "up" axis) each tie to the
retained node's same-axis translation plus a lever-arm term from the
retained node's rotation about `normal`. Per this milestone's explicit
scope (**narrower than Xara/OpenSees's `rigidDiaphragm`**): only those two
translations are ever tied — every rotational dof of the constrained node
(including rotation about `normal` itself) and its out-of-plane
translation stay completely free, never constrained by this call. A floor
diaphragm is rigid in its own plane; this doesn't impose a shared rotation
on whatever's attached to it or claim any out-of-plane stiffness.

Verified through the full `Domain3`/`Analysis3` stack
(`core/tests/m20_rigid_diaphragm3.rs`): the classic asymmetric-diaphragm
torsion problem (two offset lateral springs reproduce the standard
`[[sum(k), sum(k*z)], [sum(k*z), sum(k*z^2)]]` rigid-diaphragm torsional
stiffness matrix, derived independently in the test, not copied from the
implementation) with both the closed-form translation/rotation solution
and each wall's raw kinematic displacement checked; rigid-rotation
invariance (a constrained node's vertical translation and every rotational
dof survive a large master-node rotation completely untouched, while its
two in-plane translations respond exactly per the lever-arm formula); and
an explicit non-default `normal` (`Axis3::Z`) cross-checked against
Xara/OpenSees's own documented `RigidDiaphragm` constraint equations for
their default z-normal case.

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

Spatial modal/transient plumbing and a ground-motion acceptance case are now
covered by M19 below, not part of this milestone's original scope.

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

Local `wy`/`wz` element loads and the true spatial diaphragm transformation
are now done — see M20 below and the "Constraints" section above. Still
not done: structured diagnostics for invalid orientation hints (`vec_xz`
parallel to the member axis, or a zero-length member, currently just
panics via a plain `assert!` in `GeomTransf3::local_axes` — not the
structured `AnalysisError`-style diagnostic a real API would want).

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

### M19 — Spatial dynamics (modal/transient, 3D)

**Done.** `modal_analysis` and `TransientAnalysis` (`core/src/analysis/
modal.rs`, `core/src/analysis/transient.rs`) generalized to the same
`<const NDIM, const NDOF, const ELEMENT_DOF, NId, E>` profile as `Domain`/
`Analysis` — `TransientAnalysis3` is the new spatial type alias;
`modal_analysis` is a free function whose profile is inferred from its
`domain` argument, so no separate `modal_analysis3` is needed. This closed
out the one piece M15's "Required core changes" table had flagged as not
done for ground-motion/dynamic state handling. No dynamics formula needed
to change: the Lanczos shift-invert eigensolver, Newmark-beta stepping, the
initial-acceleration equilibrium solve, and `GroundMotion`'s effective
inertial force all already went through `Domain`'s generic free-DOF
interface (mass diagonal, gather/scatter, tangent assembly,
`direction_incidence`) rather than anything hardcoded to three planar DOFs.

Verified through the full `Domain3`/`Analysis3` stack
(`core/tests/m19_spatial_dynamics.rs`): a spatial `Truss3`'s lumped-mass
natural frequency against the closed-form SDOF `omega = sqrt(k/m)` via
`modal_analysis`; undamped and Rayleigh-damped `ZeroLength3` SDOF free
vibration against Chopra's closed forms via `TransientAnalysis3`, each
isolating a different one of the three spatial translational DOFs; and two
simultaneous, independent `GroundMotion`s (`ux` and `uz`) driving an
uncoupled two-direction spatial mass, each direction matching the same
single-direction step-response closed form `m_ground_motion.rs` already
verifies for the planar case, with no cross-talk between directions — this
is the spatial "3D ground motion" acceptance case referenced above.

Not covered by this milestone (tracked separately, see M20 below and the
top-level summary): spatial beam-load assembly, the true rigid-diaphragm
constraint transformation, and the pysees/WASM handoff for spatial
models/results.

### M20 — Spatial beam loads and true rigid-diaphragm constraints

**Done.** Both pieces M19 explicitly left for later:

- `ElementLoad3::UniformTransverse { wy, wz }` (`core/src/model/
  load_pattern.rs`) — `ElasticBeamColumn3`'s biaxial local transverse load,
  `ElementOps::Load` for the spatial profile in place of `Infallible`. See
  "Required core changes" above for the acceptance detail.
- `Domain3::rigid_diaphragm`/`rigid_diaphragm_about`, backed by a genuine
  `AffineConstraint` transformation (not identity aliasing) generalized
  through every assembly hot loop via a shared `dof_terms` helper. See
  "Constraints" above for the full design and acceptance detail.

Both are verified through the full `Domain3`/`Analysis3` stack
(`core/tests/m20_spatial_beam_loads.rs`, `core/tests/
m20_rigid_diaphragm3.rs`), independently of each other.

Not covered by this milestone: axial (`wx`) spatial beam loads,
`DispBeamColumn3`/`ForceBeamColumn3` element-load support (neither exists
for their planar counterparts either — no scope creep beyond parity),
chained rigid diaphragms (a diaphragm's retained node can't itself be
another diaphragm's slave — rejected with a panic, not silently wrong),
and dynamic (nodal-mass) analysis of a diaphragm-tied translation (the
lumped mass *diagonal* can't represent the rotational inertia a real
tributary mass there would induce at the retained node — also rejected
with a panic; see "Constraints" above).

No M15+ implementation belongs inside M10. M10 must cleanly reject spatial
input, but its public handoff contract must preserve the profile discriminant.
