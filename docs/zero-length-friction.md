# Handoff: Coulomb friction coupling for `ZeroLength`/`ZeroLength3`

## Goal

Give `ZeroLength`/`ZeroLength3` (core/src/model/elements/zero_length.rs) a shear
direction whose resistance is bounded by `mu * |normal force|`, where the
normal force comes from a *different* DOF's own material on the same element
(e.g. an `Elastic` axial/vertical spring). Today every DOF on these elements
is evaluated fully independently; this is the first cross-DOF coupling either
element will have.

## Why this isn't a `Material` variant

`Material::evaluate(&self, strain: f64)` (materials/mod.rs:352) is a pure
function of one scalar strain, and every call site (`Truss`, `ZeroLength`,
`Fiber`, the `Parallel`/`Series`/`MinMax` composites) depends on that
signature. Friction's yield force needs a second element's worth of state
(the normal DOF's current force) that a single `Material` has no way to see.
So the coupling has to live in the *element*, not the material catalog —
`Material` is not touched by this work at all.

This also isn't a new `Element`/`Element3` enum variant. `ZeroLength` is
already the right shape (two nodes, independently configured DOFs, no
geometry/local frame); friction just adds one coupling term between two of
its DOFs while every other DOF keeps working exactly as today. Do not add a
`FrictionZeroLength` variant to the `Element`/`Element3` enums
(elements/mod.rs:99-105, 229-235) — extend the existing `ZeroLength`/
`ZeroLength3` structs in place.

Precedent for "couple two DOFs while keeping `Material::evaluate` at one
scalar" already exists in this codebase: `FiberSection::trial`
(fiber_section.rs:50-63) computes each fiber's strain as `eps0 - y*kappa` — a
linear combination of two section DOFs — then integrates many independent
1-scalar `Material` evaluations into a 2x2 stress-resultant/tangent pair.
**Friction cannot reuse that exact trick** — there's no shared linear
kinematic field relating a `ZeroLength`'s normal and shear DOFs the way a
fiber's offset `y` relates axial strain and curvature. The coupling here
comes from a yield-surface constraint between two force resultants, not a
projection of two displacements onto one shared strain, so the resulting
tangent block is not guaranteed symmetric the way `FiberSection`'s `eq` term
is. Build it as its own small return-map construct (below), not as a second
`FiberSection`-style integrator.

## Scope

**Phase 1 (required): planar `ZeroLength`, one normal DOF + one shear DOF.**
Full uniaxial Coulomb return-map, coupled tangent, tests.

**Phase 2 (required): spatial `ZeroLength3`, one normal DOF + two shear
DOFs, treated independently** (each shear DOF gets its own copy of the
Phase-1 logic against the same normal DOF — a square, not circular,
interaction surface between the two shear directions).

**Phase 3 (explicit stretch — do not build unless asked):** true biaxial
friction cone for `Friction3` (one scalar slip magnitude, vector radial
return shared across both shear DOFs, circular interaction surface). This is
the physically correct behavior for a friction-pendulum-style bearing
sliding in two horizontal directions at once, but it's a materially harder
return-map problem than two independent copies of Phase 1 — don't attempt it
in the same pass as Phase 1/2.

**Explicit non-goals for this work:**
- Velocity-dependent `mu` (OpenSees `VelDependent`/`VelPressureDep`-style).
  `Node`/`Node3` already carry `velocity` (node.rs:69) so it's wireable
  later, but only meaningful under `TransientAnalysis` — leave it out now.
- Any change to `Algorithm` (algorithm.rs:4 — today just `Linear` and
  `NewtonRaphson`, no line search, no modified/initial-tangent option, no
  event-to-event). Convergence for this feature is handled entirely inside
  the element (see "Convergence" below), not by new solver infrastructure.
- `wasm-bridge/src/input_v1`'s `ZeroLengthTable` (tables.rs:104-110) is the
  serialized app-facing schema for building `ZeroLength` elements and does
  not currently have a slot for friction. Extending it is a separate
  follow-on task once the core element is done and tested — don't fold it
  into this pass.

## Design

Add to core/src/model/elements/zero_length.rs, next to the existing
`materials` array on both structs:

```rust
/// Couples one "shear" DOF to one "normal" DOF's own material: sticks
/// (stiffness `k0`) until |V| reaches `mu * yield-eligible normal force`,
/// then slides — `ElasticPP`'s return map (elastic_like.rs:13-36) with the
/// yield force fed in from the normal direction's current trial force each
/// call instead of a fixed material constant, plus a small post-slip
/// hardening slope `b*k0` (see "Convergence").
///
/// Sign convention: this codebase is compression-negative (see
/// `Material::Ent`/`Gap`, elastic_like.rs:44-58 — `strain <= 0.0` is the
/// compression/engaged branch). Friction is mobilized only when the normal
/// direction is in compression: yield force is `mu * (-normal_force).max(0.0)`.
/// If a given element's normal DOF uses the opposite sign convention, flip
/// that clamp — don't silently assume `.abs()`.
#[derive(Debug, Clone, Copy)]
pub struct Friction {
    pub normal_dof: usize,
    pub shear_dof: usize,
    pub mu: f64,
    pub k0: f64,
    /// Post-slip hardening ratio (fraction of `k0`); see "Convergence".
    /// Not optional-by-omission — require it at construction so callers
    /// make an active choice instead of silently inheriting a default.
    pub b: f64,
    slip: f64,
}

impl Friction {
    pub fn new(normal_dof: usize, shear_dof: usize, mu: f64, k0: f64, b: f64) -> Self {
        Friction { normal_dof, shear_dof, mu, k0, b, slip: 0.0 }
    }

    /// `(shear_force, d(force)/d(shear_disp), d(force)/d(normal_force), next_slip)`
    /// — same `(value, tangent, next_committed_state)` shape as
    /// `Material::evaluate` (materials/mod.rs:352), plus the extra coupling
    /// derivative friction needs. Pure — safe to call every Newton
    /// iteration, same contract as `Material::trial_stress_tangent`.
    fn evaluate(&self, normal_force: f64, shear_relative: f64) -> (f64, f64, f64, f64) {
        let yield_force = self.mu * (-normal_force).max(0.0);
        let trial = self.k0 * (shear_relative - self.slip);
        if trial.abs() <= yield_force {
            (trial, self.k0, 0.0, self.slip)
        } else {
            let excess = trial.abs() - yield_force;
            let force = (yield_force + self.b * self.k0 * excess).copysign(trial);
            let d_force_d_normal = -self.mu * trial.signum(); // d(yield_force)/d(normal_force) when normal_force < 0
            let next_slip = shear_relative - (force - self.b * self.k0 * excess.copysign(trial)) / self.k0; // see note below
            (force, self.b * self.k0, d_force_d_normal, next_slip)
        }
    }
}
```

> Implementer note on `next_slip`: work out the exact return-mapping algebra
> for the hardening branch yourself rather than trusting the line above
> verbatim — it's included to show the shape of the computation (elastic
> predictor minus whatever the plastic/frictional branch consumed), not
> checked against a worked example. Validate it with the finite-difference
> test specified below before trusting it.

`Friction3` is the same struct with `shear_dofs: [usize; 2]` instead of
`shear_dof: usize`, and `slip: [f64; 2]`. For Phase 2, its `evaluate` just
runs the Phase-1 algorithm independently for each entry of `shear_dofs`
against the same `normal_dof` — do not attempt to share a single scalar slip
between the two axes (that's Phase 3).

Add to both element structs:

```rust
pub struct ZeroLength { pub node_i: NodeId, pub node_j: NodeId, materials: [Option<Material>; NDF], friction: Option<Friction> }
```

(and the `Friction3` analogue on `ZeroLength3`), with a builder method
mirroring `with_material`:

```rust
pub fn with_friction(mut self, friction: Friction) -> Self { self.friction = Some(friction); self }
```

**The shear DOF's slot in `materials` must stay `None`** — `friction` owns
that DOF's resistance entirely. Don't try to make both populate the same DOF;
if you want to guard against misconfiguration, a debug assertion in
`with_friction` that `materials[friction.shear_dof].is_none()` is reasonable,
but isn't required for correctness.

## Wiring

### `form_tangent_and_resistance`

After the existing per-DOF `materials` loop (zero_length.rs:43-58), add:

```rust
if let Some(fr) = &self.friction {
    let normal_material = self.materials[fr.normal_dof].as_ref().expect("friction needs a normal-direction material");
    let normal_rel = node_j.displacement[fr.normal_dof] - node_i.displacement[fr.normal_dof];
    let (normal_force, normal_tangent) = normal_material.trial_stress_tangent(normal_rel);
    let shear_rel = node_j.displacement[fr.shear_dof] - node_i.displacement[fr.shear_dof];
    let (force, dv_dshear, dv_dnormal, _) = fr.evaluate(normal_force, shear_rel);

    let mut b_v = SVector::<f64, ELEMENT_DOF>::zeros(); b_v[fr.shear_dof] = -1.0; b_v[NDF + fr.shear_dof] = 1.0;
    let mut b_n = SVector::<f64, ELEMENT_DOF>::zeros(); b_n[fr.normal_dof] = -1.0; b_n[NDF + fr.normal_dof] = 1.0;

    k += dv_dshear * (b_v * b_v.transpose()) + (dv_dnormal * normal_tangent) * (b_v * b_n.transpose());
    resistance += force * b_v;
}
```

Note this only ever adds a `b_v * b_n.transpose()` block, not its transpose —
the resulting `k` is **not symmetric** in general. Confirmed fine for this
codebase: the solver (analysis/solver.rs:29-35) uses `faer`'s general sparse
LU, not a symmetric-only factorization, so a non-symmetric block from
friction is not a solver-compatibility problem.

For `ZeroLength3`/Phase 2, do this once per entry of `fr.shear_dofs`
(each against the same `fr.normal_dof`), accumulating into the
`SpatialElementMatrix`/`SpatialElementVector` the same way the existing
per-DOF loop does (zero_length.rs:112-124).

### `commit`

After the existing per-DOF commit loop (zero_length.rs:63-69):

```rust
if let Some(fr) = &mut self.friction {
    let normal_rel = node_j.displacement[fr.normal_dof] - node_i.displacement[fr.normal_dof];
    let normal_force = self.materials[fr.normal_dof].as_ref().unwrap().trial_stress_tangent(normal_rel).0;
    let shear_rel = node_j.displacement[fr.shear_dof] - node_i.displacement[fr.shear_dof];
    let (.., next_slip) = fr.evaluate(normal_force, shear_rel);
    fr.slip = next_slip;
}
```

Read the normal force via `trial_stress_tangent` (not `.commit()`) here —
the normal DOF's own slot in `materials` advances its own state through the
existing loop above this block; don't commit it twice. Every leaf material
in this catalog is constructed so that re-evaluating at the exact committed
strain reproduces the same force regardless of whether you read it just
before or just after that slot's own `commit()` runs, so ordering relative
to the existing loop doesn't matter — but don't call the normal material's
`.commit()` a second time from inside this block.

### `local_force`

No change needed — it's already just
`self.form_tangent_and_resistance(node_i, node_j).1` (zero_length.rs:75-77),
so it picks up the friction contribution automatically once the above is in
place. Confirm this with a test rather than assuming.

## Convergence

While sliding, `dv_dshear` is small (`b*k0`) rather than exactly zero. This
matters: a friction isolator is commonly the *only* stiffness path on its
shear DOF while fully sliding (unlike an `ElasticPP` damper, which rarely
is), so an exactly-zero post-slip tangent risks handing the solver a
structurally singular row/column — see `AnalysisError::SingularSystem`
(analysis/error.rs:7, thrown from analysis/solver.rs:32). This codebase
already accepts an analogous non-physical tangent floor for exactly this
reason: `Material::MinMax`'s failed state returns `stress=0.0` but
`tangent = 1e-8 * inner.initial_tangent()`, "not exactly zero, to avoid
handing back a structurally singular row/column" (materials/mod.rs:284-287,
composite.rs:58).

`Friction::b` (required, not defaulted — see Design) is the mechanism here,
and it's deliberately part of the actual force law rather than a tangent-only
fudge: `Steel01`-style kinematic hardening keeps the reported tangent
consistent with the real derivative of the force being returned, and rounds
the stick/slip corner slightly, which also reduces Newton bouncing right at
the transition (a pure tangent-only floor, MinMax-style, would fix the
singularity but leave the corner perfectly sharp). Suggest documenting
`1e-3`–`1e-4` as the recommended starting range in the constructor's doc
comment; a caller who wants an exact textbook Coulomb comparison can set it
to `0.0` and accept the singularity risk that comes back with it.

`k0` (stick-phase penalty stiffness) is a separate tuning knob: too soft
gives spurious pre-slip flexibility, too stiff risks ill-conditioning.
Recommend documenting "a few orders of magnitude above the softest real
stiffness elsewhere in the model" as the sizing rule of thumb, same as any
penalty-type contact stiffness.

## Testing

Add to zero_length.rs's existing `#[cfg(test)] mod tests` block, following
the style of `zero_length_dispatches_across_material_regimes` and
`zero_length3_dispatches_independently_per_global_direction` (drive the
element directly off manually-set node displacements, no `Analysis`/
`Algorithm` needed):

1. **Sticks below the yield surface**: small shear displacement, |V| well
   under `mu*|N|` → `V == k0 * shear_rel` exactly, tangent `== k0`.
2. **Slides at `mu * normal force`**: large shear displacement → `V`
   approaches `mu*|N| + b*k0*excess`, tangent `== b*k0`.
3. **Yield surface tracks a changing normal force**: same shear
   displacement, two different committed normal-direction states (e.g. two
   different `Elastic` axial displacements) → confirm the slip threshold
   scales with `N` — this is the direct test of "shear resistance depends on
   axial load level," not just a generic plasticity check.
4. **Remembers permanent slip after commit**: mirrors
   `elastic_pp_remembers_permanent_set_after_commit` (elastic_like.rs:76-94)
   — slide past the yield surface, commit, partially unload, confirm the
   response is elastic from the new committed `slip`, not from zero.
5. **Drops capacity out of compression**: normal force `>= 0` (tension,
   per this codebase's sign convention) → yield force is `0`, any nonzero
   shear displacement should slide immediately at (near-)zero force.
6. **Off-diagonal tangent matches finite difference** — the most important
   test here, and the easiest place to get a sign or chain-rule factor
   wrong: perturb the normal DOF's relative displacement by a small `h`,
   recompute the shear force, and check
   `(V(u_n + h) - V(u_n)) / h ≈ k[shear_dof, normal_dof]` (the analytical
   `dv_dnormal * normal_tangent` term) to a tight tolerance while sliding.
   Do this before trusting the hardening-branch `next_slip` algebra above.
7. **Phase 2**: `ZeroLength3` with one normal DOF and two independent
   shear DOFs — confirm each slides independently against the same normal
   force with no cross-talk between the two shear axes (that absence of
   coupling is exactly what distinguishes Phase 2 from the deferred Phase 3
   circular interaction).

## Files touched

- `core/src/model/elements/zero_length.rs` — everything above.
- `core/src/model/elements/mod.rs` — only if `Friction`/`Friction3` need to
  be constructible from outside the `elements` module; add them to the
  existing `pub use zero_length::{ZeroLength, ZeroLength3};` (mod.rs:84) if so.
- No changes to `core/src/model/materials/mod.rs` or any file under
  `core/src/model/materials/` — the `Material` catalog is untouched by this
  work, by design (see "Why this isn't a `Material` variant").
- No changes to the `Element`/`Element3` enums or their dispatch `match`
  arms in `core/src/model/elements/mod.rs` — `ZeroLength`/`ZeroLength3`
  remain the same variants, just with more capability inside them.
