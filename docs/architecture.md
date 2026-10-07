# Architecture

How the core engine is put together and how to extend it. Design history and the reasoning behind
these choices is in [domain-refactor-plan.md](domain-refactor-plan.md); the wire format is in
[input-format.md](input-format.md).

## Domain

`Domain<NDIM, NDOF, NId, E>` owns nodes, elements, constraints and load patterns. The two profiles are
`Domain` (2D: `ux, uy, rz`, element catalog `Element`) and `Domain3` (3D: six DOFs, catalog `Element3`).
The catalogs are closed enums: no `dyn`, no heap in the assembly loop, fixed-size `SMatrix` element math.
The analysis layer (`analysis/`) is generic over the profile.

## Element interface

Every catalog implements `ElementOps`. An element reports its nodes (`nodes()`, at most
`MAX_ELEMENT_NODES` = 4) and pushes its own contribution through sinks that are generic over the
element's matrix size:

- `TangentSink::add::<N>(&[DofRef; N], &SMatrix<N,N>, &SVector<N>)` receives stiffness and resistance.
- `VectorSink::add::<N>(&[DofRef; N], &SVector<N>)` receives loads and mass.
- `DofRef` is `(node, slot)`, so an element's local layout is its own business.

The same element code path serves tangent assembly, load assembly, mass, and reactions. Other trait
members: `dof_mask` (which node slots the element stiffens), `accepts_load`, `validate`, `prepare` (cache
immutable reference geometry; called from `Domain::validate`), `commit`, `local_force_width`,
`local_force`, `local_load_force`, `fiber_responses`, `gauss_point_count`, `gauss_responses`.
Concrete elements (`Truss`, `Quad4`, ...) keep their own typed math; only the catalog impl builds `DofRef`
arrays.

## DOF activation and state mapping

A node slot becomes an equation only if an element's `dof_mask` touches it or a constraint uses it as a
master; then fixed slots get no equation and constraint slaves are eliminated. Truss-only and
membrane-only nodes therefore need no rotation fixed by hand, and a beam and a membrane can share a node
(the membrane does not resist `rz`). A nodal load or mass on an inactive slot is a `ModelError`
reported by `Domain::validate`, never dropped. In debug builds the scatter path asserts that a nonzero
contribution lands only on declared slots.

State lives in a dense `DofTable`. Because a slave can have several terms, **never read state through
`equation_of`** (it returns only a DOF's own equation). Use `Domain::dof_terms` and
`Domain::value_at(q, node, slot)`.

## Constraints

One mechanism: `LinearConstraint`, `u_slave = sum coeff * u_master`, masters on any nodes
(`Domain::add_constraint`). `equal_dof`, `rigid_diaphragm`, `rigid_diaphragm_about` and `rigid_link` are
helpers that generate them. Chains are resolved when equations are numbered; cycles, duplicate slaves,
fixed slaves, prescribed masters (`can_prescribe`) and inconsistent initial slave state are rejected by
`validate`/`can_prescribe`. Inhomogeneous constraints are not implemented. Lumped mass or element mass
on a multi-term slave is a `ModelError::MassOnConstrainedDof`.

## Element loads

The effective load on an element is the sum of every pattern's load times its factor, so loads are an
additive accumulator. 2D `ElementLoad` is a struct (`beam_uniform`, `body`, `edge_traction[4]`,
`edge_pressure[4]`) with a fixed 16-component layout for the `elementLoad` recorder; `ElementLoad3` is
still a single-variant enum. `accepts_load` decides compatibility per element and `Domain::validate`
reports `IncompatibleElementLoad`; the wire decoder checks the same at decode time.

## Continuum elements

`model/continuum/` has Tri3/Quad4 shape functions, Gauss rules, the Jacobian and shape validation
(quads are validated at the four corners; Gauss-point positivity is not enough). `PlaneMaterial` is a
trial/commit seam on `[eps_x, eps_y, gamma_xy]` with an independent copy per Gauss point.
`Quad4Formulation::Enhanced` adds Wilson-Taylor incompatible modes with static condensation and accepts
linear materials only. `prepare` caches per-point `B` and `detJ * w`.

## When to share and when to split the 2D/3D code

**Share (one implementation, parameterized):** assembly, numbering, constraints, sinks, `DofMask`,
recorders and the whole analysis layer; the wire format's nodes, materials, fibers, loads, constraints,
sequence and recorders; isoparametric machinery (shape functions, Gauss rules, Jacobian, `B`-matrix
scatter loops) written over `NDIM`-sized arrays.

**Split (separate types):** formulations whose kinematics differ, e.g. 2D vs 3D beams (biaxial bending,
torsion, orientation vector), plane-stress/strain vs 3D solid constitutive reduction, and
friction/orientation in `ZeroLength`.

**Do not force genericity onto leaf elements.** Stable Rust cannot express `2 * NDIM` in a const
generic, so a dimension-generic `Truss` would need a redundant const parameter to size its matrices. That
costs more than the roughly 100 duplicated lines it saves, so `Truss`/`Truss3` stay as two types.

## Adding an element

1. Write the typed element (`model/elements/<name>.rs`); add a variant to `Element`/`Element3` and the
   catalog arms.
2. `dof_mask`: declare every node slot the element stiffens, no more.
3. `accepts_load`: which load kinds it takes (default: none).
4. `validate` (geometry and parameters, returning `ModelError::InvalidElement`) and `prepare` (cache
   immutable geometry).
5. Call `sink.add::<N>` from `assemble_tangent`/`assemble_load`/`assemble_mass`; implement `commit`.
6. `local_force_width` and `local_force` (and `gauss_*`/`fiber_responses` if it has them) so recorders are
   bounds-checked at decode time.
7. Tests: add the element to the mask-conformance test in `elements/mod.rs`; rigid-body and patch tests;
   closed-form check; mass test if it carries density; a cross-check against an independent solver where
   one exists.
8. Wire: see below.

## Adding a wire table

Add the table struct and an `ElementKind` variant in `wasm-bridge/src/input_v1/tables.rs` (a per-formulation
table, optional via `#[serde(default)]`, with a tsify optional marker); decode it in
`decode/elements_2d.rs` or `elements_3d.rs` and make the other profile reject it through `reject_tables`;
extend recorder validation if it has new responses; add a `DecodeError` variant for any new failure; add a
native-vs-wire test in `wasm-bridge/tests/` and a case in `boundary-smoke.ts`; document it in
`input-format.md`.
