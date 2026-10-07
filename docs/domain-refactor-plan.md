# Domain refactor and elastic Q4/T3 plan

Status: plan only, nothing implemented. Revision 3: integrates the first review
(section 8) and records the resolved decisions (section 5). Scope note: 3D support
means frames and 3D planar (shell) elements; 3D solid elements are not planned.

Scope: migrate the core model/assembly architecture so it is not tied to
two-node elements or a fixed per-profile DOF layout, unify the wire format across
2D and 3D, then add general elastic Q4 and T3 continuum elements on top. CSFM is
out of scope (see [csfm-integration.md](csfm-integration.md) for why this comes
first).

Ground rules from the project owner: nothing is in production, so there is no
backward-compatibility requirement and breaking changes are fine; 2D and 3D stay
two profiles; DOF activation is the chosen approach; the wire format is unified
now.

Terminology: 2D/3D throughout (matches pysees). Older docs and some code say
"planar"/"spatial"; those names go away where they are renamed in this plan.

## 1. Why: what the current design assumes

Findings from reading `core/src/model/domain.rs`, `elements/mod.rs`, the
analysis layer and `wasm-bridge`:

| # | Assumption | Where it lives |
|---|---|---|
| A1 | An element has exactly two nodes | `ElementOps::nodes() -> [NId; 2]`; `[id_i, id_j] = element.nodes()` at 9 sites in `domain.rs` (L250, 270, 294, 333, 359, 669, 706, 813, 880); every `form_*`/`commit` takes `(node_i, node_j)`; the row-to-node map `if a < NDOF { i } else { j }` appears 3 times |
| A2 | Element DOF count is a const generic threaded everywhere | `ELEMENT_DOF` on `Domain`, `Analysis`, `AnalysisBuilder::build`, `iterate_to_equilibrium`, `arclength`, `transient`, `modal`, `ModelSession` (about 130 occurrences across 20 files). Every element in a catalog must share one size. |
| A3 | Every node owns every profile DOF | `Node.equation: [Option<usize>; NDOF]`, numbering loop over `0..NDOF`. Unused DOFs (a truss's `rz`) must be fixed by hand or the system is singular. |
| A4 | Constraints are identity ties or a rigid diaphragm | `MpConstraint` (identity), `AffineConstraint` (single retained node, `MAX_AFFINE_TERMS = 2`), no chaining (panics), `HashMap` lookup in the assembly hot loop |
| A5 | Local force vectors are fixed `SVector<_, ELEMENT_DOF>` | `element_local_force`, `element_end_force`, recorders |
| A6 | Assembly rebuilds a triplet list and a sparse matrix every iteration | `assemble_stiffness_triplets` + `try_new_from_triplets` |
| A7 | The wire format is two parallel table sets | `CarapaceInputV1` carries `nodes`/`nodes3`, `trusses`/`trusses3`, ... (about 1.9k lines across `tables`/`tables3`/`decode`/`decode3`), `RecorderSpec`/`RecorderSpec3`, `SequenceSpec`/`SequenceSpec3` |
| A8 | One equation per physical DOF is assumed by state consumers | `ModelSession` modal export reads `domain.equation_of(node, dof)` (`session.rs` L890, L916, L1044), so a multi-term slave would export as zero; `Domain::direction_incidence`, `gather_*`, `rotational_equations` assume the same |
| A9 | Element loads are one enum with a destructuring `Add` | `ElementLoad::add` destructures its single variant irrefutably; a second variant does not compile and `Uniform + Body` has no defined meaning |
| A10 | Prescribed displacements are validated only against user fixity | `Domain::can_prescribe` / `Analysis::step_prescribed` check only `fixed[dof]`; after constraints are resolved a fixed master's contribution is already eliminated |

A good sign: the analysis layer (`newton_loop`, `arclength`, `transient`,
`modal`, `integrator`) reaches `Domain` through about 25 methods
(`gather_displacement`, `assemble_tangent_and_resistance`,
`apply_displacement_increment`, `commit`, `num_free_dofs`, ...), none of which
depend on element size. The refactor is therefore concentrated in `domain.rs`
and the element catalog, plus mechanical signature plumbing and an audit of the
state consumers in A8.

How OpenSees differs (for reference): variable node count per element
(`getNumExternalNodes`, `getNumDOF`), per-node DOF count at runtime, `NDMaterial`
strain vectors typed by material, and a general `MP_Constraint` (retained node +
constraint matrix). It does not have DOF activation: a `FourNodeQuad` requires
`ndf = 2`, so frame and quad elements cannot share a node. This plan does better
on that point.

## 2. Target architecture

### 2.1 Element interface: variable nodes, element-driven scatter

Keep the closed enum catalogs (`Element` for 2D, `Element3` for 3D) and
fixed-size `SMatrix` element math (no `dyn`, no heap in the hot loop). Replace
the fixed two-node, fixed-size-return trait with one where the element pushes
its own contribution through a sink that is generic over the element's matrix
size:

```rust
pub struct DofMask(u8);                       // bit s set => node slot s is touched
pub type DofRef<NId> = (NId, u8);             // (node, slot)
pub const MAX_ELEMENT_NODES: usize = 4;   // frames, Tri3/Quad4, 3D shell tri/quad; 3D solids are out of scope
pub type NodeList<NId> = ArrayVec<NId, MAX_ELEMENT_NODES>;

pub trait TangentSink<NId> {
    fn add<const N: usize>(&mut self, dofs: &[DofRef<NId>; N],
                           k: &SMatrix<f64, N, N>, r: &SVector<f64, N>);
}
pub trait VectorSink<NId> {                   // loads, mass
    fn add<const N: usize>(&mut self, dofs: &[DofRef<NId>; N], v: &SVector<f64, N>);
}

pub trait ElementOps<const NDIM: usize, const NDOF: usize, NId: Copy> {
    type Load: Copy + ElementLoadComponents + Add<Output = Self::Load> + Mul<f64, Output = Self::Load>;
    type Id: Key;
    fn nodes(&self) -> NodeList<NId>;
    fn dof_mask(&self) -> DofMask;            // which node slots this element stiffens
    fn accepts_load(&self, load: &Self::Load) -> bool;       // see 2.6
    fn validate(&self, nodes: &NodeView<..>) -> Result<(), ModelError> { Ok(()) }
    fn prepare(&mut self, nodes: &NodeView<..>) {}           // cache immutable reference geometry
    fn assemble_tangent<S: TangentSink<NId>>(&self, nodes: &NodeView<..>, load: Option<&Self::Load>, sink: &mut S);
    fn assemble_load<S: VectorSink<NId>>(&self, nodes: &NodeView<..>, load: Option<&Self::Load>, sink: &mut S);
    fn assemble_mass<S: VectorSink<NId>>(&self, nodes: &NodeView<..>, sink: &mut S);
    fn commit(&mut self, nodes: &NodeView<..>, load: Option<&Self::Load>);
    fn local_force_width(&self) -> usize;                    // recorder component bound
    fn local_force(&self, nodes: &NodeView<..>) -> ElementForce;
    fn local_load_force(&self, nodes: &NodeView<..>, load: Option<&Self::Load>) -> ElementForce;
    fn fiber_responses(&self, nodes: &NodeView<..>) -> Option<Vec<Vec<(f64, f64)>>>;
    fn gauss_responses(&self, nodes: &NodeView<..>) -> Option<Vec<GaussResponse>>;   // new, continuum
}
```

Key points:

- `ELEMENT_DOF` disappears from `Domain`/`Analysis`/builder/session generics
  (A2). The profile is `(NDIM, NDOF, NId, E)`.
- Each variant arm in the catalog impl knows its own `N` and calls
  `sink.add::<N>(..)`. A 6x6 truss and an 8x8 Q4 coexist in one catalog with no
  padding.
- **Concrete elements keep their existing math and signatures**
  (`Truss::form_tangent_and_resistance(&Node, &Node) -> (SMatrix<6,6>, SVector<6>)`
  etc.). Only the two catalog impls (`Element`, `Element3`) change; they do the
  node lookups and build the `DofRef` arrays. This keeps the diff, and the
  risk, small.
- `ElementForce` is a small newtype over `Vec<f64>` (indexable, `len()`, `iter()`)
  used for recorder paths only, which are not hot. `local_force_width` is static
  per element so the wire decoder can validate a recorder's `component` at decode
  time instead of assuming `2 * NDOF`.
- Elements stay free to use whatever local layout they like because the sink
  takes explicit `(node, slot)` references. The 2-node frame elements keep
  `[node i all slots, node j all slots]`; Q4 uses `[n1 ux, n1 uy, n2 ux, ...]`.
- Loads and mass go through `VectorSink`, so `Domain::assemble_load_with`,
  `assemble_mass_diagonal` and `reaction` all use the same element code path
  instead of three hand-rolled node-i/node-j loops.
- `prepare` is called from `Domain::validate()` (idempotent). It lets continuum
  elements cache reference geometry (per-Gauss-point `B`, `detJ * w`) once. This
  relies on node coordinates being immutable after `add_node`, which is true today
  (no `node_mut`); if coordinates ever become mutable (updated Lagrangian), `prepare`
  must be re-run.

### 2.2 DOF activation and state mapping

A node keeps profile slots (3 in 2D, 6 in 3D). A slot becomes an equation only
if something uses it:

```text
active(slot) = touched by an element's dof_mask
             OR referenced as a master by a constraint
and then:     fixed(slot)      -> no equation (user boundary condition, unchanged)
              slave(slot)      -> eliminated by constraint
              otherwise        -> free, numbered in node-insertion order
```

Consequences:

- A truss-only node no longer needs `rz` fixed by hand (the README and
  `Truss`/`Truss3` doc comments call this out today).
- A Q4/T3 node has no `rz`. A node shared by a beam and a membrane gets `rz`
  from the beam only: the beam end rotation is free, the membrane does not
  resist it. This is the standard behavior and is what makes mixed frame and
  continuum models work (OpenSees cannot do this).
- Constraint masters must count as active. The 2D diaphragm retained node usually
  has no element at all; without this rule its DOFs would be silently dropped.
  This is the trap in a naive implementation, so it is called out and tested.
- **Silent-drop guard.** A DOF that is inactive, unfixed and unconstrained but
  receives a nonzero element contribution means an element under-declared its
  `dof_mask`. A debug-build assertion in the scatter path catches this (and a
  test per element kind checks `dof_mask` covers every nonzero row of `k`/`r`).
- A nodal load or nodal mass on an inactive DOF is a model error
  (`ModelError::LoadOnInactiveDof` / `MassOnInactiveDof`), reported by
  `Domain::validate()`, never silently ignored. `prescribe_displacement`
  continues to require the DOF to be user-fixed (and see 2.3 for constraint
  masters).
- Equation order is unchanged for fully active models (node insertion order,
  then slot), so mode-shape indices and existing numerics are preserved.

DOF state moves out of `Node` into a dense `DofTable` owned by `Domain`
(`node ordinal -> slot states`, plus a terms arena for constrained DOFs). This
replaces `Node.equation` and the `HashMap<(NId, usize), Vec<..>>` in the hot
loop.

**State-mapping API (A8).** Once a constraint slave has several terms, "the
equation of a DOF" no longer exists, so every consumer that goes from a
free-DOF-indexed vector to a node/slot value must use one mapping instead of
`equation_of`:

- `Domain::dof_terms(node, slot) -> DofTerms` where `DofTerms` is
  `None | Single { eq, coeff } | Many(&[(eq, coeff)])`. `Single` carries the
  coefficient: `u_slave = 2 * u_master` is representable, and coefficient 1 is just
  the common case, not a special type.
- `Domain::value_at(q: &DVector<f64>, node, slot) -> f64`: the transformed value
  (`sum coeff * q[eq]`; `0` for fixed or inactive DOFs).
- `equation_of` stays, documented as "the DOF's own equation if it is a free
  unknown, `None` for a slave or fixed or inactive DOF". It must no longer be used
  to read state.

Consumers to audit and convert (each gets a test in step 1.4):

| Consumer | Today | Needed |
|---|---|---|
| Modal shape export (`session.rs` L916, L1044) | `equation_of` -> multi-term slave reads as 0 | `value_at` |
| Ground-motion influence (`session.rs` L890; `Domain::direction_incidence`) | `1.0` at each node's own equation for that slot | must still be a rigid translation of the constrained model; verify it equals `T * 1` for rigid links and diaphragms |
| Participation factors and total mass (`session.rs` L393-L395, L921-L935) | reduced-space `phi^T M iota` | stays valid in reduced space; verify against full-space result on a constrained model |
| `gather_displacement/velocity/acceleration`, `scatter_state`, `apply_displacement_increment` | per-node equations plus affine map | read/write via the same `DofTable` |
| `rotational_equations` | classifies by the slot index of the equation's owner | unchanged rule (an equation's classification is its own slot), but confirm for constrained models |
| `Integrator::DisplacementControl` / arc-length seed (`integrator.rs`, `arclength.rs` L920) | `equation_of(..).ok_or(..)` | slave or inactive DOF must return a clear error, not a wrong DOF |
| `Domain::reaction` | element-driven sum at one `(node, slot)` | unchanged (does not need equations) |

### 2.3 General linear multi-point constraints

One mechanism replaces `MpConstraint` and `AffineConstraint`:

```rust
pub struct LinearConstraint<NId> {
    slave: (NId, usize),
    terms: Vec<(NId, usize, f64)>,     // u_slave = sum coeff * u_master, masters on any node
}
```

- `equal_dof`, planar `rigid_diaphragm`, 3D `rigid_diaphragm_about` become
  helpers that generate `LinearConstraint`s (same public API, same results).
- New helper `rigid_link(master, slave)`: full rigid-body coupling with the lever
  arm (2D: `ux, uy, rz`; 3D: `ux..rz`). It is what joins a shell/continuum edge to
  a beam and what an eventual embedded/offset feature needs.
- Resolution at numbering time substitutes chains (a slave whose master is itself
  a slave) with cycle detection (`ModelError::ConstraintCycle`). This lifts the
  existing "no chained diaphragms" limitation.
- The terms live in an arena indexed from the `DofTable`, so the hot loop reads a
  slice with no `HashMap` and no `MAX_AFFINE_TERMS` cap.

**Normalization and conflicts (applied when constraints are added or resolved).**

- Combine repeated masters in one constraint (sum coefficients), drop terms whose
  coefficient is exactly zero or below a scale-aware tolerance relative to the
  largest coefficient, reject an empty result unless it means `u_slave = 0`
  (which is then a valid fixity-equivalent constraint and is kept).
- Reject: a slave that appears in two constraints (`DuplicateSlave`); a slave that
  is also user-fixed (`SlaveIsFixed`); a self-referencing term; a node or slot out
  of range; a cycle.
- Masters that are user-fixed are permitted (their contribution is zero) **only if
  they stay at zero**; see the prescription rule below.

**Prescribed masters (A10).** The check cannot live only at build/decode time,
because `Analysis::step_prescribed` can change a fixed DOF after numbering, when
its contribution has already been eliminated. Rule:

- `Domain::can_prescribe(node, slot, value)` additionally returns `false` if
  `(node, slot)` appears as a master term in any constraint. Chains need no extra
  handling: a chain is built from direct references, so a DOF that feeds a chain
  is a direct master of its first link.
- `Analysis::step_prescribed` already validates every target before changing
  anything (`InvalidConstraint`, state untouched); the new rule is enforced by
  that same call, so a rejected step leaves the domain exactly as it was.
- Static check too: `Domain::validate()` rejects a constraint master that carries a
  nonzero `with_initial_displacement` (`ConstraintOnPrescribedDof`).
- Inhomogeneous constraints (offset term) are not implemented; it would touch
  `apply_displacement_increment`, `scatter_state` and `reaction`. Revisit only if a
  model needs it.

**Initial state.** The solver updates slaves incrementally, so an initially
inconsistent slave state would persist. `Domain::validate()` rejects, rather than
silently projecting, any slave whose initial displacement or velocity differs from
what its masters imply (`InconsistentInitialState`, with a small scale-aware
tolerance). The default (all zero) always passes.

Lumped mass on a multi-term slave stays unsupported (it would need a full
`T^T M T`), now as a `ModelError` instead of a panic.

### 2.4 Two profiles: when to share, when to split

The 2D/3D split stays: `Element`/`Element3`, `Domain`/`Domain3`, `Node2`/`Node3`.
The refactor shares everything that is dimension-agnostic and keeps a split where
the formulation genuinely differs. Rule of thumb for future work:

**Share (one implementation, parameterized):** assembly, numbering, constraints,
sinks, `DofMask`, recorders and the whole analysis layer (already generic); the
wire format's nodes, materials, fibers, loads, constraints, sequence and
recorders; isoparametric machinery (shape functions, Gauss rules, Jacobian,
`B`-matrix scatter loops) written over `NDIM`-sized arrays.

**Split (separate types):** formulations whose kinematics differ, e.g. 2D vs 3D
beams (biaxial bending, torsion, orientation vector), plane-stress/strain vs 3D
solid constitutive reduction, and friction/orientation in `ZeroLength` (the 2D and
3D versions differ in logic, not just array length).

**Do not force genericity onto leaf elements.** Stable Rust cannot express
`2 * NDIM` in a const generic (`generic_const_exprs`), so a dimension-generic
`Truss` would need a redundant extra const parameter to size its matrices. That
costs more than the roughly 100 duplicated lines it saves, so `Truss`/`Truss3`
stay as two types. Revisit if that language feature stabilizes.

### 2.5 Unified wire format

One `CarapaceInputV1` with `header.ndm` (2 or 3; replaces `space`). Everything
dimension-agnostic is shared and written once. Tables are per **element
formulation**, not per profile, and every table is `#[serde(default)]`.

| Table | Change |
|---|---|
| `nodes` | one `NodeTable`; `coords` stride `ndm`; `fixed` bitmask (6 bits fits `u8`); `mass` stride `ndf` |
| `trusses` | one table (identical fields in 2D and 3D) |
| `elasticBeamColumns2d` / `3d` | split by formulation (different fields), 2D table must be empty in a 3D model and vice versa |
| `dispBeamColumns2d` / `3d`, `forceBeamColumns2d` / `3d` | same |
| `zeroLengths`, `zeroLengthSections` | one table each; `friction` rows carry `shearDofs: Vec<u8>` (length 1 in 2D, 2 in 3D); `orient` rows `{row, x: [f64;3], yp?: [f64;3]}` |
| `fibers` | one `FiberTable`; `z` empty for 2D sections, same length as `y` for 3D |
| `equalDofs` | one table |
| `rigidDiaphragms` | one table with optional `normal` (3D; defaults to Y), absent in 2D |
| `rigidLinks` (new) | `{ master, slave }` node pairs; core `rigid_link` |
| `linearConstraints` (new) | general `u_slave = sum coeff * u_master`: parallel arrays `slaveNode`, `slaveDof`, `termOffsets`, `termNode`, `termDof`, `termCoeff` (flattened sparse rows, like the existing sparse tables) |
| `elementLoads` | tagged `ElementLoadSpec` (see 2.6); kind-incompatible assignments are a decode error |
| `sequence` | one `SequenceSpec` / `RecorderSpec`; `ElementKind` becomes one enum covering all formulations (`truss`, `elasticBeamColumn2d`, `elasticBeamColumn3d`, ..., `quad4`, `tri3`); recorder `component` is validated against `ElementOps::local_force_width` (and the load width for `elementLoad`) |

Decode picks the session profile once from `ndm`. A table or element kind that
does not belong to the profile is a `DecodeError::TableNotInProfile { table }`
(or `ElementKindNotInProfile`), not silently ignored. `Session::{Planar, Spatial}`
become `Session::{D2, D3}`. Every new core capability gets a wire table in the same
phase that introduces it (rigid links and general constraints in Phase 2, continuum
tables and the Q4 formulation selector in Phase 3); nothing is left reachable from
core only.

### 2.6 Element-load algebra

The effective load on an element is the sum of every pattern's load, scaled by its
factor (`effective_element_load`), so loads must form an additive accumulator.
Today `ElementLoad` is a one-variant enum with a destructuring `Add`; a second
variant would not compile, and `Uniform + Body` has no meaning. Resolution:

- Make the 2D `ElementLoad` a **struct accumulator** with one field group per load
  kind, all additive:
  `beam_uniform: [wx, wy]` (local axes, per length), `body: [bx, by]` (global axes,
  per volume), `edge_traction: [[tx, ty]; 4]` (global axes, per edge length) and
  `edge_pressure: [f64; 4]` (positive = acting into the element along the inward
  normal). `Add`/`Mul<f64>` are field-wise and total; `Uniform + Body` is simply
  both fields set.
- Compatibility is decided per element by `ElementOps::accepts_load`: beam-columns
  accept only `beam_uniform`, continuum elements accept only `body`, `edge_*` (edges
  beyond the element's edge count must be zero). It is checked **before** any
  accumulation: `Domain::validate()` walks every pattern's element loads and returns
  `ModelError::IncompatibleElementLoad { element }`; the wire decoder checks at
  decode time and returns a `DecodeError`.
- `ElementLoadComponents` gets a documented fixed component layout
  (`[wx, wy, bx, by, t0x, t0y, p0, t1x, t1y, p1, ...]`) and a `component_count`, so
  the existing `elementLoad` recorder keeps working and its `component` is bounds
  checked.
- `ElementLoad3` stays an enum for now (single variant, still compiles); the rule for
  any future 3D load kind is the same accumulator-struct approach. Converting it
  now is optional.

## 3. Step-by-step plan

Each step is one reviewable commit that leaves `cargo test --workspace` green.
Steps within a phase are ordered by dependency. Behavior changes (DOF activation,
new constraint rules) are never mixed into a structural commit.

### Phase 0: safety net (before touching anything)

**0.1 Baseline.** Record, on `main`: `cargo test --workspace` pass; the OpenSees
comparison harness outputs (`comparison/run_all.py`, results in `comparison/out/`);
and `benchmark_frame` timings (elastic and fiber models, one size) as the
performance reference, **natively and under wasm** (build with wasm-pack and run the
same model through Node, as `comparison/wasm_frame.mjs` already does for the frame).
Reasoning: this refactor must be behavior-preserving for every existing element, so
the existing 36 core test files (about 7.3k lines) plus the wasm tests are the
oracle; native-only timing does not establish browser performance, which is the
actual target.
Files: none changed; store numbers in the PR description (or
`docs/refactor-baseline.md`, deleted at the end).

**0.2 Snapshot test for cross-cutting behavior.** Add
`core/tests/refactor_snapshot.rs`: a handful of representative models (2D frame with
diaphragm; 2D fiber cantilever; 3D frame with 3D diaphragm; a transient run; an
arc-length run; a modal run with a constraint) that assert displacement, reaction,
mode shapes and element-force values to 1e-12 against values captured on `main`.
Reasoning: existing tests assert physics tolerances, not bit-for-bit continuity; this
makes any numbering or assembly-order drift visible immediately. Delete or fold into
existing tests when the refactor lands.

### Phase 1: core Domain refactor

Ordered so every commit builds and passes: `DofTable` lands with every DOF active
before `Node.equation` goes away; the trait change carries its dependent signature
changes with it; activation (the behavior change) is its own commit.

**1.1 `DofTable` (structure only, every DOF still active).** *Done.*
Add the `DofTable` owned by `Domain` (`NodeView`, a thin read-only wrapper over the
node store, moved to 1.2 where it is first used, to avoid dead code in between). The table reproduces today's
numbering exactly (all slots active, node-insertion order, same identity/affine
handling), is the only source of equation numbers, and replaces `Node.equation` in
the same commit (`equation_of`, `rotational_equations`, `direction_incidence`,
`gather`, `scatter_state`, `apply_displacement_increment` read the table). Keep the
existing `MpConstraint`/`AffineConstraint` for now.
Files: `model/node.rs`, `model/domain.rs`, `model/mod.rs`.
Reasoning: introduces the numbering structure without any behavior change, and keeps
the tree green because nothing is removed before its replacement exists.
Tests: whole existing suite and `refactor_snapshot` identical; unit test that the
table reproduces the old numbering on an identity-tied and an affine-tied model.

**1.2 New `ElementOps`, sinks, and signature plumbing (one commit).** *Done; deviations below.*
Implement section 2.1: `DofMask`, `DofRef`, `NodeList`, `TangentSink`, `VectorSink`,
`ElementForce`, `ModelError`, the new trait and the two catalog impls (every element
returns the **full** mask for now, so behavior is unchanged). Rewrite in `domain.rs`:
`assemble_stiffness_triplets`, `assemble_load_with`, `assemble_mass_diagonal`,
`reaction`, `commit`, `element_local_force`, `element_end_force`,
`element_fiber_responses` to use the sinks (removes all 9 `[id_i, id_j]` sites and
the three row-to-node maps). In the **same commit** remove `ELEMENT_DOF` from
`Analysis`, `AnalysisBuilder::build`, `iterate_to_equilibrium`, the `ArcLength` impl
blocks, `TransientAnalysis`, `modal_analysis`, and `ModelSession`; remove the public
`ELEMENT_DOF`/`SPATIAL_ELEMENT_DOF` constants and `Spatial*Matrix/Vector` aliases
(or make them element-private). Aliases `Domain`, `Domain3`, `Analysis`, `Analysis3`,
`TransientAnalysis3` keep their names, so almost all tests are unaffected.
Add `arrayvec` to `core/Cargo.toml` for `NodeList` (decision in section 5).
Files: `model/elements/mod.rs` (largest change: both catalog impls, roughly 500 lines
of mostly mechanical rewrite), `elements/{truss,zero_length,elastic_beam_column,
disp_beam_column,force_beam_column,uniform_load}.rs` (private consts,
`local_force_width`), `model/{domain,mod}.rs`, new `model/error.rs`,
`analysis/{builder,state,newton_loop,arclength,transient,modal,integrator}.rs`,
`wasm-bridge/src/input_v1/session.rs` (27 occurrences), `core/tests/
arclength_continuation.rs` (4), `core/tests/m19_spatial_dynamics.rs` (1), the 7 tests
that call `element_local_force`/`element_end_force` (`element_local_force.rs`,
`disp_beam_element_load.rs`, `force_beam_element_load.rs`, `element_load_axial.rs`,
`fiber_recorder.rs`, `m21_transient_corrector.rs`, plus `wasm-bridge/tests/
m10_carapace_input_v1.rs`) for the `ElementForce` return type, `core/Cargo.toml`.
Reasoning: removes A1, A2, A5. The trait change and its plumbing cannot be split
without a non-compiling intermediate state, so they are one (large but mechanical)
commit; the numerics are identical because every element is still full-mask.
Tests: whole existing suite unchanged and passing; `refactor_snapshot` identical.
New unit tests: a mock 3-node and a mock 4-node element in a test catalog (a stand-in
`ElementOps` impl) assembled into a `Domain`, checking `K`, `r`, load and mass vectors
against hand-computed values; `ElementForce` indexing; `local_force_width` matches the
length of `local_force` for every element kind.
Docs: doc comments on the trait (replace the `ElementOps` essay in `elements/mod.rs`).
As built: `validate`, `prepare`, `accepts_load`, `ModelError` and the real `dof_mask` values moved to
1.3 / 3.x, where they are first used (every element returns the full mask in 1.2, so no behavior
changes and no dead code); `NodeView` landed here. The old per-element dispatch methods on `Element`
and `Element3` (`form_tangent_and_resistance`, ...) are now private and the `ElementOps` impls call
them, so concrete element files were not touched. `Domain::reaction` only evaluates elements that
touch the node (as before). Result: dump identical to baseline, 243 tests pass, native elastic
30x6 about 5% slower and wasm within noise (see [refactor-baseline.md](refactor-baseline.md)).

**1.3 DOF activation (behavior change).** *Done; deviations below.*
Implement section 2.2's activation rule on top of the `DofTable`: activation pass
(element masks, then constraint masters), `Domain::is_active`, `Domain::validate() ->
Result<(), ModelError>` (inactive-DOF loads/mass, element `validate`, `prepare`), the
debug assertion in the scatter sink. `AnalysisBuilder::build` calls `validate` and
panics with the error text on failure (existing convention for model-construction
errors); wasm decode calls `validate` first and maps to `DecodeError::InvalidModel`.
Each concrete element gets its real `dof_mask`:

| Element | Mask |
|---|---|
| `Truss` / `Truss3` | translations (`ux, uy` / `ux, uy, uz`) |
| `ElasticBeamColumn`, `Disp`/`ForceBeamColumn` (2D) | `ux, uy, rz` |
| 3D beam-columns | all six |
| `ZeroLength` / `ZeroLength3` | union of DOFs that carry a material or a friction coupling; with an `orient`ation, all translations (and all rotations if any rotational DOF is used), since the local frame mixes them |
| `ZeroLengthSection` / `3` | section axial and flexural DOFs plus any spring DOFs |

Files: `model/domain.rs`, `model/elements/*.rs` (masks), `analysis/builder.rs`,
`model/error.rs`; audit `analysis/integrator.rs` L77/L135 and `analysis/arclength.rs`
L920 (`equation_of(...).ok_or(...)`: an inactive DOF now also returns `None`, and the
error returned should still be correct).
Reasoning: removes A3. Models that fix unused DOFs by hand keep working (explicit
`fixed` wins). The user-visible change is that truss-only and frame-plus-membrane
models stop needing manual `fix`.
Tests (new, `core/tests/dof_activation.rs`): truss-only 2D and 3D models solve with no
rotations fixed; frame + truss joint; load or mass on an inactive DOF returns
`ModelError`; `equation_of` is `None` for inactive DOFs; constraint masters stay active
(a diaphragm retained node with no elements still carries the diaphragm DOFs);
`dof_mask` conformance test per element kind (nonzero `k`/`r` rows must be within the
mask); modal analysis of a truss with unfixed rotations; arc-length does not demand a
rotation scale when no rotation is active. Existing tests that fix unused DOFs by hand
keep their fixes.
Docs: README note about manually fixing unused rotations removed; `Truss`/`Truss3` doc
comments updated.

As built: (1) **A nonzero nodal mass also activates its DOF**, and `MassOnInactiveDof` was dropped.
The plan treated mass on an unused DOF as an error, but `core/tests/m20_rigid_diaphragm3.rs`
deliberately builds mass-only DOFs (a free mass is legal under Newmark), so treating mass as
inertia keeps that behavior exactly; only loads on unused DOFs are errors. (2) Validation reaches
callers as values, not panics: `AnalysisBuilder::try_build`, `AnalysisError::InvalidModel`,
`modal_analysis` and `TransientAnalysis::new` validate, and the wasm session maps the error to
`AnalysisErrorDetail::InvalidModel` at stage start (stage loads are registered then) while
`decode` reports geometry/mass problems as `DecodeError::InvalidModel`; a panic would abort the
wasm module. (3) `ElementOps::validate` exists with a default of "accept"; no element overrides
it yet (the first real use is Phase 3). `prepare` and `accepts_load` remain deferred. (4) The
mask conformance test lives in `elements/mod.rs` (`mask_conformance`) because it needs the
crate-private `NodeView`; it checks `mask == stiffened slots` exactly for every element kind except
oriented 3D zero-length (covered, not exact). It runs in the default test profile, which is also
where the debug-build guard is active.

**1.4 General linear constraints, state mapping, prescription rule.**
Implement sections 2.2 (state-mapping API) and 2.3: `LinearConstraint`, normalization
and conflict checks, chain substitution with cycle detection, terms arena,
`DofTerms::{None, Single{eq, coeff}, Many}`; re-express `equal_dof`,
`rigid_diaphragm`, `rigid_diaphragm_about` on top; add `rigid_link`; add
`ModelError::{ConstraintCycle, DuplicateSlave, SlaveIsFixed, ConstraintOnPrescribedDof,
InconsistentInitialState, MassOnConstrainedDof}`. Convert every consumer in the A8 table
to `dof_terms`/`value_at`. `Domain::can_prescribe` rejects constraint masters.
Files: `model/domain.rs`, `analysis/constraint.rs` (doc), `analysis/{integrator,
arclength,modal}.rs` and `wasm-bridge/src/input_v1/session.rs` (state consumers),
`model/error.rs`.
Tests: `core/tests/m_constraints.rs` and `m20_rigid_diaphragm3.rs` pass unchanged (the
regression proof that the new engine reproduces identity and diaphragm ties); new
`core/tests/linear_constraints.rs`: chained diaphragm; rigid link against the analytic
`u_s = u_m + theta x r` kinematics in 2D and 3D; **coefficient-2 tie** (`Single` with
`coeff != 1`); cycle, duplicate slave, slave-is-fixed, self-reference, zero-coefficient
removal and repeated-master merging; congruence check (`K_reduced == T^T K_full T` on a
small model); `equal_dof` vs the same tie written as a general constraint solve
identically; **constrained modal test** that exports shapes through the session and
verifies slave-node motion equals the transformed master values (and that participation
factors and total mass match a full-space computation); ground-motion influence vector
correct on a rigid-link model; `can_prescribe` rejects a constraint master (including
the first link of a chain), and `step_prescribed` with such a target returns
`InvalidConstraint` with the domain bit-for-bit unchanged afterwards; inconsistent
initial slave state rejected.
Docs: `ConstraintHandler::Transformation` doc comment (it currently describes aliasing
and says general affine ties are out of scope).

**1.5 Phase 1 gate.** `cargo test --workspace`; `refactor_snapshot` identical;
`benchmark_frame` within a few percent of the Phase 0 baseline in both native and
wasm; `comparison/` harness outputs unchanged; wasm size noted. Then: delete
`refactor_snapshot.rs` (or fold into existing tests) and the baseline note.

### Phase 2: wire-format unification

Done after Phase 1 so `session.rs` is already profile-generic without `ELEMENT_DOF`;
done before the continuum elements so their tables are added once.

**2.1 Unified tables and header.**
Rewrite `tables.rs`, delete `tables3.rs`; collapse `RecorderSpec3`/`SequenceSpec3` into
`RecorderSpec`/`SequenceSpec`; unify `ElementKind`; header `space` -> `ndm`;
`#[serde(default)]` on every table; per-formulation tables per section 2.5; add
`rigidLinks` and `linearConstraints` tables (core features from Phase 1.4); recorder
component bounds from `local_force_width`. New `DecodeError` variants:
`TableNotInProfile`, `ElementKindNotInProfile`, `InvalidModel`, `InvalidRecorderComponent`.
Files: `wasm-bridge/src/input_v1/{mod,tables,sequence,error}.rs`; delete `tables3.rs`.

**2.2 Shared decoding.**
Split `decode.rs`/`decode3.rs` into: shared helpers (nodes generic over `<NDIM, NDOF>`,
materials, fibers, equal-dof, diaphragms, rigid links, linear constraints, load
patterns, stage compilation, recorder resolution via a `(ElementKind, row) -> E::Id`
map) and two small per-profile element decoders (`elements_2d.rs`, `elements_3d.rs`).
`Session::{Planar, Spatial}` -> `Session::{D2, D3}`; `PlanarSession`/`SpatialSession` ->
`Session2`/`Session3`.
Files: `wasm-bridge/src/input_v1/{decode,decode3,session,sequence}.rs` ->
`decode/{mod,shared,elements_2d,elements_3d}.rs`; `boundary.rs` (names only).
Reasoning: the duplication that exists today is not the formulations (genuinely
different) but the bookkeeping around them.

**2.3 Tests.**
Port, do not delete, `wasm-bridge/tests/m10_carapace_input_v1.rs` (2096 lines) and
`m10_spatial_input_v1.rs` (677 lines): mostly field renames and removal of the empty
opposite-profile tables; the semantics they verify are unchanged. Update
`arclength_session.rs`, `material_probe.rs`, `boundary-smoke.ts`. New tests: a table
that belongs to the other profile is rejected with the right error; defaulted (omitted)
tables decode; serde round trip of the unified input; unknown `ndm` rejected; rigid-link
and general-constraint tables decode and produce the same results as the native core
model; a recorder `component` out of range is a decode error.

**2.4 Consumers (sibling repos and local tooling).**
- `comparison/wasm_frame.mjs` (builds input objects by hand): update.
- Regenerate `pkg/` and `pkg-web/` (wasm-pack output; the `.d.ts` is generated by
  tsify, do not hand-edit).
- pysees (sibling repo; the change is theirs to make): `src/app/types/carapaceInputV1.ts`
  (mirror of the wire shape), `src/app/lib/carapace/compileInputV1.ts` (it currently
  emits every `*3` field as an empty table; those lines are deleted, so the compiler
  gets smaller), `src/app/carapace/wasm/` (copy of the regenerated package), and a
  check of `carapaceWorker.ts`/`resultsStorage.ts` for `space` references. Land this in
  lockstep with 2.1-2.3.
Docs: rewrite the table sections of `docs/input-format.md` (Nodes, Elements,
Constraints, Loads, Recorders, "Extending the format").

### Phase 3: continuum foundations and elastic T3/Q4

**3.1 Isoparametric utilities and geometry validation.** New
`core/src/model/continuum/` module: `shape.rs` (T3 and Q4 shape functions and
derivatives, written over `NDIM`-sized arrays so a 3D shell's surface mapping can reuse the Jacobian machinery),
`quadrature.rs` (Gauss rules: 1-point triangle, 3-point triangle, 2x2 quad, 2-point
line for edges), `jacobian.rs` (`detJ`, inverse, physical derivatives).
Validation (`validate_triangle`, `validate_quad`) with a scale-aware tolerance
(relative to the element's characteristic length squared, not an absolute number):
- Reject repeated or coincident nodes (distance below tolerance).
- Triangle: area above tolerance and positive orientation (counter-clockwise).
- Quad: **`detJ` evaluated at the four reference corners, not only at Gauss points.**
  For a bilinear map `detJ` is itself linear in `(xi, eta)`, so positivity at the four
  corners is necessary and sufficient for positivity everywhere; positivity at the
  Gauss points alone is not (the quad (0,0), (1,0), (0.4,0.4), (0,1) has all four
  Gauss-point determinants positive, 0.187/0.1/0.1/0.013, but the determinant at the
  reflex corner is -0.05, so the element is invalid). Also reject a negative corner
  determinant for clockwise node ordering, with a message that names the offending
  corner.
Reasoning: this is the shareable part (section 2.4); everything dimension-specific (`B`
matrix shape, constitutive reduction) stays in the element.
Tests: partition of unity and `sum dN = 0`; derivative check against finite differences;
quadrature exactness (areas, polynomial integrals, triangle rules, edge rule); `detJ`
equals area scale for affine maps; **the reflex-corner quad above is rejected while a
convex trapezoid and a parallelogram are accepted**; repeated-node and clockwise
inputs rejected; tolerance behaves the same for a unit element and a 1e6-scaled copy.

**3.2 `PlaneMaterial` catalog (n-D material seam).** New `core/src/model/materials/plane.rs`:
a closed enum for 3-component plane strain vectors with the same trial/commit design as
the uniaxial `Material` (pure `trial_stress_tangent(&self, strain) -> (stress, tangent)`,
`commit(strain) -> Self`, plus `is_linear(&self)`).
Conventions, fixed now so nonlinear materials later do not have to change them:
- Strain vector `[eps_x, eps_y, gamma_xy]` with **engineering shear** `gamma_xy = 2 *
  eps_xy`; stress `[sigma_x, sigma_y, tau_xy]`; `D` maps one to the other.
- `ElasticMatrix { d: [d11, d12, d13, d22, d23, d33] }`: packed upper triangle, row
  by row, symmetric `D` (so the engineering-shear convention is explicit).
- `Orthotropic { ex, ey, nu_xy, g_xy, angle }` (plane stress): `angle` is the
  counter-clockwise angle in radians from global x to material axis 1; `D` is rotated
  into global axes with the stress/strain transformations that respect engineering shear.
- `PlaneStress { e, nu }`, `PlaneStrain { e, nu }` isotropic, sharing parameter helpers
  with a 3D shell/section material later rather than duplicating them. (3D solid elements and
  6-component solid materials are out of scope for this project.)
**History ownership:** an element owns one independent `PlaneMaterial` value per Gauss
point (`Tri3`: 1, `Quad4`: 4), cloned from the prototype given to its constructor and
committed individually in `commit` (exactly as the fiber elements own one material per
fiber). This is the seam where nonlinear plane materials (CSFM later, damage,
plasticity) plug in with per-point history.
Reasoning: "general elastic" means the element should not hard-wire isotropy; and
defining ownership and conventions now avoids a breaking redesign when the first
stateful material arrives.
Tests: `D` symmetric and positive definite; plane stress vs plane strain equivalence
(`E' = E/(1-nu^2)`, `nu' = nu/(1-nu)` give the same in-plane response); orthotropic
with equal moduli reduces to isotropic for any `angle`; rotation by 90 degrees swaps
`ex`/`ey` terms; `ElasticMatrix` ordering round-trips against the equivalent isotropic
`D`; invalid constants rejected (`nu` outside bounds, nonpositive moduli, non-PD
matrix) with `Result`, not panic; independent per-point commit (a stateful test
material in the unit tests).
Files: new `materials/plane.rs`; `materials/mod.rs`, `model/mod.rs` exports.

**3.3 `Tri3` (constant-strain triangle).** `core/src/model/elements/tri3.rs`: 3 nodes,
6x6, one Gauss point, thickness, `PlaneMaterial`, density, `dof_mask = {ux, uy}`,
`validate` via 3.1, `prepare` (caches `B` and `A`), `gauss_responses`, `local_force`
(nodal resisting forces in global DOF order), `local_force_width = 6`. Mass: exact
`rho t A / 3` per node, applied to **each translational DOF**. Added to the 2D `Element`
enum.
Tests: closed-form `K` for a known triangle (Cook/Hughes reference); constant-strain
patch test on arbitrary triangulations (exact to machine precision); rigid-body modes
(`K` has exactly 3 zero eigenvalues, stress-free under rigid motion); symmetry and
positive semi-definiteness; **mass: each translational direction sums to `rho t A`
separately (summing the ux and uy diagonals together gives twice the physical mass, so
the test sums per direction)**; passes the `dof_mask` conformance test; solves with no
`rz` fixed (activation); `prepare` is idempotent and cached values equal recomputation.

**3.4 `Quad4` (bilinear quadrilateral, plain isoparametric).**
`core/src/model/elements/quad4.rs`: 4 nodes, 8x8, 2x2 Gauss, same material/thickness/
density interface, `validate` via 3.1, `prepare` (caches per-Gauss-point `B` and
`detJ * w`).
**Lumped mass is the row-sum of the consistent mass:** `m_i = sum_g rho t N_i(g) detJ(g)
w(g)` with the 2x2 rule, which is exact here (bilinear `N_i` times linear `detJ` is
degree at most 2 per axis). It equals `rho t A / 4` only for parallelograms; for a
trapezoid it does not (for example (0,0), (2,0), (1.5,1), (0.5,1) with `rho t = 1` gives
nodal masses 5/12, 5/12, 1/3, 1/3, against `A/4 = 0.375`). `m_i` is applied to each
translational DOF separately, so each direction sums to `rho t A`. The same weights
`m_i / rho` drive the consistent body-force load (3.6), so mass and load integration
share one routine.
Tests: same battery as `Tri3`; patch test on distorted quad meshes (exactly reproduced);
invariance to cyclic relabeling of the nodes (identical `K` after permutation);
**individual nodal masses on a trapezoid equal 5/12, 5/12, 1/3, 1/3 and a parallelogram
gives `A/4` each; per-direction mass sum equals `rho t A`**; cantilever tip deflection
convergence against Timoshenko beam theory (documents the expected locking ratio of the
plain element and asserts mesh-refinement convergence rate); a thick-cylinder (Lame)
stress/displacement check on a coarse mesh.

**3.5 `Quad4` enhanced formulation (separate commit, own acceptance tests).**
Adds `Quad4Formulation::{Full, Enhanced}`. `Enhanced` is the Wilson-Taylor
incompatible-modes element (QM6/Q6) with this exact contract:
- *Internal modes:* for each displacement component add the two incompatible shapes
  `(1 - xi^2)` and `(1 - eta^2)`, giving four internal parameters `alpha` per element,
  discontinuous across edges (so they condense out locally).
- *Enhanced strain matrix:* `G_tilde(xi, eta) = (detJ0 / detJ(xi, eta)) * J0^{-T} *
  G_hat(xi, eta)`, with `J0` the Jacobian at the element center and `G_hat` the
  reference-coordinate derivatives of the internal shapes. The `detJ0/detJ` scaling is the
  Taylor-Beresford-Wilson correction that preserves the patch test on distorted
  elements; without it the element fails constant strain.
- *Condensation:* `K_e = K_uu - K_ua K_aa^{-1} K_au` with the 4x4 `K_aa` factored per
  element (singular `K_aa` is a `ModelError`); the element resistance is `K_e * u_e`
  (the element is linear elastic, see below).
- *Strain and stress recovery:* `alpha = -K_aa^{-1} K_au u_e`, then at each Gauss point
  `eps = B u_e + G_tilde alpha` and `sigma = D eps`. `gauss_responses` reports these
  recovered values (stress **including** the internal-mode contribution), so the
  verification `D * B * u` is only valid for `Full`; for `Enhanced` it must be
  `D * (B u + G_tilde alpha)`.
- *Scope limits:* `Enhanced` accepts only `PlaneMaterial` values with `is_linear()`
  (checked in `validate`), because for a nonlinear material `alpha` becomes element
  internal state with its own local Newton iteration; that is a later feature, not part
  of this milestone. It relieves **shear (bending) locking** only. It does **not** cure
  volumetric locking in nearly incompressible plane strain (`nu -> 0.5`); that needs a
  mixed or B-bar formulation, which is out of scope, and this is stated in the doc
  comment and a test (the locking persists and is documented, not hidden).
- *Distortion:* incompatible-mode elements can remain sensitive to trapezoidal
  distortion. Acceptance tests therefore assert what is guaranteed (exact patch test on
  distorted meshes; exactness for rectangular elements in pure bending; improvement over
  `Full` on rectangles in bending) and **record, rather than assert superiority of**,
  results on trapezoidal meshes (MacNeal-Harder-style cantilever with rectangular,
  parallelogram and trapezoid meshes), compared against the independent OpenSees
  `enhancedQuad` oracle in 3.9.
Tests: patch test on distorted meshes; pure-bending cantilever exact for rectangles with
one element through the depth; recovered-stress consistency (`gauss_responses` equals
`D (B u + G_tilde alpha)` and differs from `D B u`); condensation identity (`alpha`
recovery reproduces the condensed force); nonlinear material rejected; the
nearly-incompressible locking test.

**3.6 Element loads: body force and edge traction/pressure; recorders.**
Implement section 2.6 (accumulator `ElementLoad`, `accepts_load`, validation). New load
fields: `body: [bx, by]` and per-edge `edge_traction`/`edge_pressure`. Body and edge
loads are integrated with the element shape functions (body: the 3.4 weights; edge:
2-point Gauss on the edge with the edge's tangent length scale, edge `k` joins local node
`k` to `k+1`), producing consistent nodal loads. `Domain::element_gauss_responses(id)`
wired through `ElementOps::gauss_responses`.
Files: `model/load_pattern.rs`, `model/domain.rs`, `elements/{tri3,quad4}.rs`,
`model/error.rs`.
Tests: `Uniform + Body` accumulation yields both fields set; several body-force
patterns with different factors sum correctly; a beam rejects `body`/`edge`, a continuum
element rejects `beam_uniform`, an out-of-range edge rejected (`validate` in core,
`DecodeError` in wasm, no accumulation attempted); total applied force equals `b t A` and
`t * integral(traction) dl` for both elements; self-weight column reaction equals total
weight; **patch tests with edge loads**: a uniform traction on the boundary of a distorted
mesh reproduces the exact constant stress field (this validates load integration without
hand-computed nodal forces); edge pressure sign convention; Gauss-point responses match
`D * B * u` (`Full`) by hand evaluation.

**3.7 Mixed and dynamic verification (core).** New `core/tests/continuum_*.rs`:
- Mixed 2D model: Q4 panel with a beam along one edge sharing nodes. Axial stiffness
  equals `E A_panel / L + E A_beam / L` exactly for `nu = 0` (proves activation plus
  assembly across element kinds).
- Rigid link between a beam end and a Q4 edge; constrained modal shapes include slave
  motion (re-uses the 1.4 session path).
- Modal: cantilever strip frequencies approach beam theory under refinement; total mass
  check.
- Transient: undamped free vibration energy conservation of a Q4 block under Newmark.
- Nonlinear solver path: a continuum model through `Algorithm::Newton`, `LineSearch`, and
  `Integrator::ArcLength` converges in one iteration for a linear element (guards the whole
  analysis layer with a multi-node element).

**3.8 Wire exposure.** Add `quads`, `triangles` tables (`nodeIds` with stride 4 and 3,
`thickness`, `material` into a new `planeMaterials` arena, `density`, and for quads
`formulation: full | enhanced`); `PlaneMaterialSpec` tagged union (isotropic with `state:
planeStress | planeStrain`, `orthotropic` with `angle`, `elasticMatrix` with the packed
ordering from 3.2); `ElementKind::{Quad4, Tri3}`; `RecorderSpec::GaussPoint {
elementKind, elementIndex, point, quantity: strain | stress, component }` with `point` and
`component` bounds checked against element metadata; `ElementLoadSpec` variants for body
force and edge traction/pressure (aggregated into the accumulator at decode, wrong-kind
and wrong-profile assignments rejected); decode rejects all of these in a 3D model
(`ElementKindNotInProfile`).
Files: `wasm-bridge/src/input_v1/{tables,sequence,error,decode/elements_2d,session}.rs`.
Tests: decode round trip; bad material index, bad node index, wrong profile, wrong-kind
load errors; a `quad4` with `enhanced` and a nonlinear material rejected; end-to-end
`advance` on a 2-element panel against the same model built natively in `core`; recorder
output equals the `core` Gauss-point values; multiple body-force patterns through wasm.

**3.9 OpenSees cross-check.** Extend `comparison/` (`opensees_frame.py`, `run_all.py`,
`wasm_frame.mjs`) with a cantilever panel using OpenSees `quad` (`PlaneStress`, the plain
formulation), `tri31`, and `enhancedQuad`, comparing nodal displacements, reaction, and
Gauss-point stress to Carapace (tolerance 1e-9 relative on displacements for the plain
elements, which should agree to round-off for the same formulation). For `enhancedQuad`,
OpenSees uses its own enhanced-strain variant; agreement is expected on rectangles and is
**checked and reported, not assumed**, on distorted meshes (a difference there is
recorded with its cause, not silently loosened). This is the independent oracle; the
closed-form tests above are the analytical one.

### Phase 4: performance (benchmark-gated, native and wasm)

**4.1 Measure first.** Add a `core/examples/benchmark_continuum.rs` (an N x N Q4 panel
under load, N around 50-100, so 5k to 20k DOF) reporting time in element assembly,
sparse-matrix construction, symbolic analysis, numeric factorization and solve, and run
the same model through the wasm build under Node (extend `comparison/wasm_frame.mjs` or
add a sibling). Native numbers alone do not establish browser performance. Reasoning:
whether the triplet-rebuild path (A6) or element geometry computation dominates for 8x8
blocks on large meshes is an empirical question; decide with numbers.

**4.2 Cheap, safe first: reference geometry cache (already built in).** `prepare` (3.3,
3.4) caches per-Gauss-point `B` and `detJ * w` at validation time, so no Jacobian or
shape-derivative work happens per iteration. This costs a few hundred bytes per element
and is valid for every material, so it ships with the elements, not behind a gate.

**4.3 Gated options, chosen by the 4.1 breakdown:**
- *Pattern-cached assembly* (if matrix construction exceeds roughly a quarter of
  per-iteration time): build the CSC pattern once (element connectivity x constraint
  terms, diagonal always present) and cache it in the `Domain` (invalidated by
  `add_node`, `add_element`, `add_constraint`, renumbering); each iteration scatters
  values into the fixed pattern (per-entry lookup by binary search in the column, no
  per-element slot tables that would cost memory on large wasm models).
  `multiply_stiffness` and `assemble_newmark_system` reuse it. The existing
  `SparseSolver` pattern-equality check stays. Tests: assembled `K` equals the triplet
  path value-for-value (property test over randomized small models including
  constraints), symbolic analysis count stays 1 across iterations, pattern invalidation
  after `add_element`.
- *Cached element stiffness for linear materials* (only for `PlaneMaterial::is_linear()`;
  the stiffness is then constant, so assembly is a copy and the resistance is `K_e u`):
  gated separately because it costs 8x8 f64 (512 B) per Q4, about 50 MB at 100k
  elements, which matters in a wasm memory budget. Adopt only if the measurement shows
  element kernels dominate and the memory cost is acceptable; never for nonlinear
  materials or the enhanced element's internal state.
Tests: benchmark shows the gain; full suite unchanged.

### Phase 5: documentation and cleanup

- New `docs/architecture.md` (current, short): the Domain model, the element interface
  and sinks, DOF activation and the state-mapping API, general constraints, the load
  accumulator, the 2D/3D split rules (section 2.4 verbatim), how to add an element
  (checklist: `dof_mask`, `accepts_load`, `validate`, `prepare`, `local_force_width`, sink
  calls, tests including mask conformance and patch/rigid-body tests), how to add a wire
  table.
- `README.md`: features list (Q4/T3 elements, plane materials, general constraints, rigid
  link, edge/body loads), remove the manual-fix note, link `architecture.md`.
- `docs/input-format.md`: rewritten per Phase 2; Q4/T3 tables, loads and recorders per 3.8.
- `docs/algorithms.md`/`docs/arclength.md`: one-line note that rotational scale is only
  needed when rotational DOFs are active.
- `docs/obsolete/spatial-architecture.md` and `implementation-plan.md`: add a top-of-file
  pointer to `architecture.md` for the superseded "do not represent planar as spatial"
  statement (history is kept, not edited).
- `docs/csfm-integration.md`: replace gap 1 (separate `membranes` store) and gap 3 (MPC)
  with a pointer to this plan; the CSFM spike then depends on Phase 3 and, for per-point
  nonlinear history, on the material seam in 3.2.
- Remove `docs/refactor-baseline.md` and the temporary snapshot test.

## 4. Files touched (summary)

| Area | Files |
|---|---|
| Core model | `model/{mod,node,domain,load_pattern}.rs`, new `model/error.rs`, `model/elements/{mod,truss,zero_length,elastic_beam_column,disp_beam_column,force_beam_column,uniform_load}.rs`, new `elements/{tri3,quad4}.rs`, new `model/continuum/*`, new `materials/plane.rs` |
| Core analysis | `analysis/{builder,state,newton_loop,arclength,transient,modal,integrator,constraint,solver}.rs` (signatures, state-mapping consumers) |
| Core tests | new `dof_activation.rs`, `linear_constraints.rs`, `continuum_*.rs`, temporary `refactor_snapshot.rs`; edits to the 7 `element_local_force`-using files, `arclength_continuation.rs`, `m19_spatial_dynamics.rs`; `testkit.rs` unaffected |
| wasm-bridge | `input_v1/{mod,tables,sequence,error,session}.rs`, `decode.rs`+`decode3.rs` -> `decode/*`, delete `tables3.rs`, `boundary.rs` (names), `tests/*` (port), `tests/boundary-smoke.ts` |
| Tooling | `comparison/{wasm_frame.mjs,opensees_frame.py,run_all.py}`, regenerated `pkg/`, `pkg-web/`, new `core/examples/benchmark_continuum.rs` |
| pysees (sibling) | `types/carapaceInputV1.ts`, `lib/carapace/compileInputV1.ts`, `carapace/wasm/*`, worker check |
| Docs | `README.md`, `docs/{architecture(new),input-format,algorithms,arclength,csfm-integration}.md`, pointers in two obsolete docs |

## 5. Decisions (all resolved)

| # | Decision | Outcome |
|---|---|---|
| 1 | `NodeList` implementation and capacity | `arrayvec` crate, capacity **4** (frames, Tri3, Quad4, 3D shell triangles and quads). 3D solids (Hex8, Quad8, Tet4) are not planned; raising the constant is a one-line change if that ever changes. |
| 2 | `DofMask` granularity | One mask per element, applied to all of its nodes. |
| 3 | Enhanced Q4 | Included in this milestone as its own commit (3.5), linear materials only. |
| 4 | Element-load representation | Accumulator struct with `accepts_load` validation (2.6). |
| 5 | Initial slave state | Reject an inconsistent initial slave displacement/velocity (`InconsistentInitialState`); never project. |
| 6 | Prescribed constraint masters | Reject through `can_prescribe` and `validate`; inhomogeneous constraints are not implemented. |
| 7 | Phase 4 option selection | Decided after the 4.1 benchmark (native and wasm), using its time and memory breakdown. |

## 6. Risks

- **Hot-loop regression from the sink design.** The sink is monomorphized per element
  size, so it should be no slower than today; the Phase 0 baseline (native and wasm) and
  the Phase 1 gate exist to prove it. The old code also looked up `dof_terms` through a
  `HashMap` for affine DOFs, so the arena should be faster, not slower.
- **Silent DOF dropping from under-declared masks.** Mitigated by the debug assertion and
  the per-element mask conformance test; this is the main new way to get a wrong answer
  without an error, so the test is mandatory for every element (including future ones).
- **State consumers reading `equation_of`.** Any missed consumer returns zero for a
  multi-term slave without failing. Mitigated by converting every one in the A8 table in
  one commit (1.4), a constrained modal test through the session, and by keeping
  `equation_of` documented as "own equation only".
- **Zero-length mask with `orient`.** The orientation mixes translations (and rotations)
  in one local frame; the mask rule above must be validated against
  `core/tests/m2_zero_length.rs` and the friction tests before 1.3 merges.
- **Enhanced Q4 correctness and distortion sensitivity.** The formulation contract in 3.5
  is a published method but non-trivial; mitigations are the patch tests, the recovery
  identities, and the independent OpenSees oracle. I have not verified published
  distortion-sensitivity figures myself; the acceptance tests record results rather than
  assert a performance claim.
- **Blast radius in `wasm-bridge`.** The wire change is the biggest single breaking diff;
  it is isolated in Phase 2 so Phase 1 can merge and be validated on its own.
- **pysees lockstep.** The sibling repo must update with Phase 2; since nothing is in
  production this is a coordination cost, not a compatibility one.

## 7. Definition of done

- `cargo test --workspace` green; `refactor_snapshot` and the OpenSees `comparison/`
  outputs unchanged for all pre-existing element types after Phase 1.
- `benchmark_frame` within a few percent of baseline, native and wasm (or faster after
  Phase 4).
- No `ELEMENT_DOF` generic anywhere; no `[NId; 2]` in `ElementOps`; no `*3` wire tables; no
  manual rotation fixing needed for truss-only models; no state consumer reads state
  through `equation_of`.
- Q4 and T3 pass the patch, rigid-body, mask-conformance, closed-form, mass (per direction,
  trapezoid) and OpenSees cross-check tests; invalid quads (including the reflex-corner
  case) are rejected; a mixed frame + membrane model solves correctly with no manual DOF
  fixing; edge-load patch tests pass.
- Every core capability introduced (rigid link, general constraints, Q4 formulation
  selector, body/edge loads) is reachable and tested through the wire format.
- `docs/architecture.md` written; README, `input-format.md` and `csfm-integration.md`
  updated; pysees compiler updated and building against the regenerated package.

## 8. Review disposition

Comments from the first review, and what was done with each:

| Comment | Disposition |
|---|---|
| 1. Q4 mass lumping wrong for distorted quads | **Accepted, verified.** The trapezoid (0,0), (2,0), (1.5,1), (0.5,1) has row-sum masses 5/12, 5/12, 1/3, 1/3 against `A/4 = 0.375`. 3.4 now defines lumping as `m_i = sum rho t N_i detJ w`, shared with the body-force integration, with per-node trapezoid tests and per-direction total-mass tests (3.3 and 3.4 also now state that summing both translational diagonals doubles the mass). |
| 2. Gauss-point `detJ > 0` does not validate a Q4 | **Accepted, verified.** For the quad (0,0), (1,0), (0.4,0.4), (0,1) the Gauss-point determinants are all positive but the reflex-corner determinant is -0.05. 3.1 validates at the four reference corners (sufficient because the Q4 `detJ` is linear in `xi, eta`), rejects repeated nodes, and uses a scale-aware tolerance. |
| 3. Constraints need a full state-mapping API | **Accepted, verified.** `session.rs` reads `equation_of` at L890, L916, L1044, so a multi-term slave exports as zero. Added the `dof_terms`/`value_at` API, an audit table of consumers (2.2), a constrained modal test through the session (1.4), and `DofTerms::Single { eq, coeff }` so `u_s = 2 u_m` is representable. |
| 4. Prescribed-master rejection must hold when prescriptions change | **Accepted, verified.** `can_prescribe`/`step_prescribed` only check user fixity. The rule is now in `can_prescribe` (so `step_prescribed` rejects before touching state), with `validate` covering initial displacements, and tests for chains and for the unchanged-state-on-failure property. Chains need no extra logic because they are built from direct references. |
| 5. `Body` breaks the load algebra | **Accepted, verified.** `ElementLoad::add` destructures one variant irrefutably. Resolved with an accumulator struct plus `accepts_load` validated before accumulation (2.6, 3.6), with core and wasm tests. Alternative (fallible `add_element_load`) listed as a decision. |
| 6. Commit sequence cannot stay green | **Accepted.** Phase 1 re-ordered: `DofTable` (all DOFs active) first and `Node.equation` removed in the same commit; the trait change and its dependent signature plumbing are one commit; activation is separate. One consequence: 1.2 is large but mechanical, because splitting it would leave a non-compiling tree. |
| 7. Enhanced Q4 needs a precise contract | **Accepted.** Own commit (3.5) with internal modes, Taylor correction, condensation, strain/stress recovery including the internal-mode contribution, linear-material-only scope, explicit statement that volumetric locking is not cured, and acceptance tests that record (not assert superiority on) trapezoidal results. The reference to Abaqus documentation was not independently checked. |
| Edge traction/pressure loads | **Accepted** (3.6), including edge-load patch tests. This is added scope. |
| Constraint normalization and conflicts; initial state | **Accepted** (2.3). One modification: reject an inconsistent initial slave state rather than project it, to avoid silently overwriting user input (decision 5). |
| Future material seam (history ownership, axis angle, shear convention, packed order) | **Accepted** (3.2). |
| Wire exposure (constraints, Q4 formulation, recorder bounds from metadata) | **Accepted** (2.5, 2.1, 3.8). The plan previously added core features without wire tables; it now says every capability must be reachable through the wire in the phase that introduces it. |
| Benchmark wasm; cache element geometry | **Accepted** (0.1, 4.1, 4.2): wasm timing is part of the baseline and the geometry cache ships with the elements via `prepare`. **Modified:** caching the full element stiffness is *not* adopted by default, only gated on measurement, because 512 B per Q4 is about 50 MB at 100k elements, which matters under a wasm memory budget, and it is invalid for nonlinear materials. |

Not adopted, with reasons: none of the comments was rejected outright. The two
deliberate deviations are the initial-state rule (reject, not project) and gating
element-stiffness caching behind a measurement.
