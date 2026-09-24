use std::collections::HashSet;

use faer::sparse::Triplet;
use nalgebra::DVector;
use slotmap::{Key, SlotMap};

use super::load_pattern::LoadPattern;
use super::{
    Element, Element3, ElementOps, LoadPatternId, LoadSeries, Node, Node3Id, NodeId, SparseMatrix, NDF, PLANAR_NDIM,
    SPATIAL_ELEMENT_DOF, SPATIAL_NDF, SPATIAL_NDIM,
};

/// A multi-point constraint tying `dofs` of `constrained` exactly to the
/// same DOFs of `retained` (`u_c[dof] = u_r[dof]`) — identity ties only, no
/// coefficients, no cross-DOF terms. See `ConstraintHandler::Transformation`
/// for how this is resolved (DOF-equation aliasing, not a real
/// transformation matrix) and why that's sufficient for `equal_dof`/
/// `rigid_diaphragm` specifically.
#[derive(Debug, Clone)]
struct MpConstraint<NId> {
    retained: NId,
    constrained: NId,
    dofs: Vec<usize>,
}

/// Owns all nodes, elements, and load patterns. No serialization/broker
/// machinery (§2.3) — this is the whole model, in memory, for one worker.
/// `Clone` backs `Analysis`/`TransientAnalysis`'s snapshot-and-restore-on-
/// failure (a failed step shouldn't leave nodal displacement/velocity/
/// acceleration partway through a discarded Newton iteration) — see
/// `Material`'s doc comment for why materials themselves don't need this
/// protection. Also what makes multi-phase analysis composition possible:
/// `Analysis::into_domain`/`TransientAnalysis::into_domain` hand this same
/// `Domain` (nodal state, committed material history, and load-pattern
/// freeze state all intact) to the next phase's builder.
///
/// Generic over the kinematic profile: `NDIM`/`NDOF`/`ELEMENT_DOF` (const
/// generics, kept as separate parameters rather than computed from each
/// other — stable Rust can't evaluate `2 * NDOF` inside a generic item), the
/// node ID type, and `E: ElementOps<...>` (the element catalog — `Element`
/// or `Element3`; its `Domain`-element-store key comes from `E::Id`, an
/// `ElementOps` associated type rather than a further generic parameter
/// here — see `ElementOps::Id`'s doc comment for why: a truly independent
/// `EId` parameter had nothing else to constrain it and broke type
/// inference at every `Domain::new()` call site). This bookkeeping
/// (equation numbering, sparse assembly, state gather/scatter) doesn't care
/// about element physics, only about DOF counts and `ElementOps`'s uniform
/// shape, so it's genuinely one implementation rather than two — unlike
/// `Truss`/`Truss3`, whose actual formulas differ and stay separate
/// concrete types (see `ElementOps`'s doc comment for why that part can't
/// be unified the same way).
///
/// All five parameters default to the planar profile, the same trick `Node`
/// uses, so every existing bare `Domain` usage (`Domain::new()`, `fn f(d:
/// &Domain)`, ...) keeps compiling unchanged. `Domain3` is the spatial
/// instantiation.
#[derive(Debug)]
pub struct Domain<
    const NDIM: usize = PLANAR_NDIM,
    const NDOF: usize = NDF,
    const ELEMENT_DOF: usize = { super::ELEMENT_DOF },
    NId = NodeId,
    E = Element,
> where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>,
{
    nodes: SlotMap<NId, Node<NDIM, NDOF>>,
    elements: SlotMap<E::Id, E>,
    mp_constraints: Vec<MpConstraint<NId>>,
    load_patterns: SlotMap<LoadPatternId, LoadPattern<NDOF, NId, E::Id, E::Load>>,
    default_pattern: LoadPatternId,
    num_free_dofs: usize,
}

/// `Domain`'s spatial instantiation — six-DOF `Node3`/`Element3`. See
/// `Domain`'s doc comment for why this is a type alias over one generic
/// implementation rather than a hand-duplicated struct.
pub type Domain3 = Domain<SPATIAL_NDIM, SPATIAL_NDF, SPATIAL_ELEMENT_DOF, Node3Id, Element3>;

/// Manual rather than `#[derive(Clone)]`: the derive macro only adds `E:
/// Clone`, not `E::Load: Clone` (it can't see through the associated type),
/// so it under-constrains this struct and fails to compile at every call
/// site instead. `Analysis::step`'s snapshot-and-restore-on-failure needs
/// this — see `Domain`'s doc comment.
impl<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E> Clone
    for Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId> + Clone,
    E::Load: Clone,
{
    fn clone(&self) -> Self {
        Domain {
            nodes: self.nodes.clone(),
            elements: self.elements.clone(),
            mp_constraints: self.mp_constraints.clone(),
            load_patterns: self.load_patterns.clone(),
            default_pattern: self.default_pattern,
            num_free_dofs: self.num_free_dofs,
        }
    }
}

impl<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E> Default
    for Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E> Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>,
{
    /// A fresh, empty model — already carrying one `LoadPattern`
    /// (`default_pattern`, `LoadSeries::Linear { slope: 1.0 }`, unscaled)
    /// so simple single-pattern models (the common case) don't need to
    /// name a pattern at all: `domain.load_node(id, dof, value)` reaches
    /// for it implicitly. Reproduces, exactly, the single-pattern behavior
    /// every `Analysis` had before multi-pattern support existed.
    pub fn new() -> Self {
        let mut load_patterns = SlotMap::default();
        let default_pattern = load_patterns.insert(LoadPattern::new(LoadSeries::Linear { slope: 1.0 }));
        Domain {
            nodes: SlotMap::default(),
            elements: SlotMap::default(),
            mp_constraints: Vec::new(),
            load_patterns,
            default_pattern,
            num_free_dofs: 0,
        }
    }

    pub fn add_node(&mut self, node: Node<NDIM, NDOF>) -> NId {
        self.nodes.insert(node)
    }

    pub fn add_element(&mut self, element: E) -> E::Id {
        self.elements.insert(element)
    }

    pub fn node(&self, id: NId) -> &Node<NDIM, NDOF> {
        &self.nodes[id]
    }

    /// `Domain::new()`'s always-present pattern — see its doc comment.
    pub fn default_pattern(&self) -> LoadPatternId {
        self.default_pattern
    }

    /// Add a new, independently-scaled `LoadPattern` (scale factor 1.0;
    /// chain `.with_scale_factor` via `add_load_pattern_scaled` if needed).
    /// Xara/OpenSees's `pattern Plain <tag> <series> {...}`.
    pub fn add_load_pattern(&mut self, series: LoadSeries) -> LoadPatternId {
        self.load_patterns.insert(LoadPattern::new(series))
    }

    /// Same as `add_load_pattern`, with an explicit scale factor (e.g. a
    /// ground-motion PGA scale, or a load-factor-vs-displacement-target
    /// unit conversion) applied on top of the series' own factor.
    pub fn add_load_pattern_scaled(&mut self, series: LoadSeries, scale_factor: f64) -> LoadPatternId {
        self.load_patterns.insert(LoadPattern::new(series).with_scale_factor(scale_factor))
    }

    /// Add a nodal load to `pattern` — Xara/OpenSees's `load <node> <...>`
    /// inside a `pattern` block.
    pub fn add_nodal_load(&mut self, pattern: LoadPatternId, node: NId, dof: usize, value: f64) {
        self.load_patterns[pattern].add_nodal_load(node, dof, value);
    }

    /// Convenience for the common single-pattern case: adds `value` to
    /// `default_pattern()` directly, so simple models don't need to name a
    /// pattern (`domain.load_node(id, 0, force)` instead of
    /// `domain.add_nodal_load(domain.default_pattern(), id, 0, force)`).
    pub fn load_node(&mut self, node: NId, dof: usize, value: f64) {
        self.add_nodal_load(self.default_pattern, node, dof, value);
    }

    /// Add an element load (e.g. a beam-column's uniform transverse load)
    /// to `pattern` — Xara/OpenSees's `eleLoad` inside a `pattern` block.
    /// `load`'s type is this profile's `ElementOps::Load` — `ElementLoad`
    /// for the planar profile, and the uninhabited `Infallible` for the
    /// spatial one (so this is uncallable there until a real spatial
    /// element-load type exists — see `ElementOps`'s doc comment).
    pub fn add_element_load(&mut self, pattern: LoadPatternId, element: E::Id, load: E::Load) {
        self.load_patterns[pattern].add_element_load(element, load);
    }

    /// Freeze `pattern`'s contribution at its value as of `pseudo_time` —
    /// Xara/OpenSees's `loadConst -time <pseudo_time>` (there, applied to
    /// every pattern in the domain; here, one pattern at a time, so a
    /// caller composing phases can freeze only the pattern(s) that should
    /// stop ramping, e.g. gravity, while others continue). Call this with
    /// whatever pseudo-time the phase you're ending actually stopped at
    /// (`StepResult::load_factor` / `TransientStepResult::time`) — see
    /// `LoadPattern`'s doc comment for why the frozen value is captured
    /// explicitly here rather than lazily cached on every assemble call.
    pub fn hold_pattern_constant(&mut self, pattern: LoadPatternId, pseudo_time: f64) {
        self.load_patterns[pattern].hold_constant(pseudo_time);
    }

    /// Tie `dofs` of `constrained` exactly to the same DOFs of `retained`
    /// (`u_c = u_r` for each listed dof) — Xara/OpenSees's `equalDOF`.
    /// Requires `ConstraintHandler::Transformation`; see its doc comment
    /// for how this is resolved.
    pub fn equal_dof(&mut self, retained: NId, constrained: NId, dofs: &[usize]) {
        self.mp_constraints.push(MpConstraint {
            retained,
            constrained,
            dofs: dofs.to_vec(),
        });
    }

    /// Equation number of a node's DOF, or `None` if it's fixed. Used by
    /// `Integrator::DisplacementControl` to locate its controlled DOF in
    /// the free-DOF system.
    pub(crate) fn equation_of(&self, node: NId, dof: usize) -> Option<usize> {
        self.nodes[node].equation[dof]
    }

    /// Assign a sequential equation number to every free, unconstrained
    /// DOF, in node insertion order, then alias every multi-point-
    /// constrained DOF to its retained node's equation number for that DOF
    /// (`ConstraintHandler::Transformation` — see `MpConstraint`). Fixed
    /// DOFs, and constrained DOFs whose retained DOF is itself fixed, get
    /// no equation number. This is the entire numbering pass — done once,
    /// not re-checked every step (§4.4: no live re-solve means no
    /// `hasDomainChanged()`-style machinery).
    ///
    /// Constrained DOFs must be excluded from the *first* pass (not just
    /// overwritten afterwards) — otherwise the equation number allocated
    /// to them before being overwritten is never referenced by any node,
    /// leaving an all-zero row/column in the assembled system.
    ///
    /// Returns the number of free DOFs (the size of the global system).
    pub(crate) fn number_dofs(&mut self) -> usize {
        let constrained_dofs: HashSet<(NId, usize)> = self
            .mp_constraints
            .iter()
            .flat_map(|c| c.dofs.iter().map(move |&dof| (c.constrained, dof)))
            .collect();

        let mut next = 0;
        for (id, node) in self.nodes.iter_mut() {
            for dof in 0..NDOF {
                node.equation[dof] = if node.fixed[dof] || constrained_dofs.contains(&(id, dof)) {
                    None
                } else {
                    let eq = next;
                    next += 1;
                    Some(eq)
                };
            }
        }

        for constraint in &self.mp_constraints {
            let retained_eq = self.nodes[constraint.retained].equation;
            for &dof in &constraint.dofs {
                self.nodes[constraint.constrained].equation[dof] = retained_eq[dof];
            }
        }

        self.num_free_dofs = next;
        next
    }

    pub fn num_free_dofs(&self) -> usize {
        self.num_free_dofs
    }

    /// Whether any `equal_dof`/`rigid_diaphragm` constraint has been added —
    /// used by `AnalysisBuilder<Ready>::build` to reject `ConstraintHandler
    /// ::Plain` (which can't resolve them) early, at model-construction
    /// time rather than as a silently-wrong solve.
    pub(crate) fn has_mp_constraints(&self) -> bool {
        !self.mp_constraints.is_empty()
    }

    /// Number DOFs (idempotent given a fixed set of nodes) and assemble the
    /// lumped-mass matrix's diagonal over free DOFs — user-assigned
    /// `Node::mass` plus each element's own lumped mass (`Element::
    /// form_mass`, M6; zero unless a `Truss`/`ElasticBeamColumn` was given
    /// a nonzero `density`). Needed by modal analysis (M5) and time-history
    /// analysis (M6). Just the diagonal, not a full (dense or sparse) N×N
    /// matrix — `M` is diagonal by construction (lumped mass), so storing
    /// anything more is pure waste, the same class of mistake a dense
    /// stiffness matrix would be (see `SparseMatrix`'s doc comment).
    pub fn assemble_mass_diagonal(&mut self) -> DVector<f64> {
        self.number_dofs();
        let n = self.num_free_dofs;
        let mut mass = DVector::<f64>::zeros(n);
        for (_, node) in self.nodes.iter() {
            for dof in 0..NDOF {
                if let Some(eq) = node.equation[dof] {
                    mass[eq] += node.mass[dof];
                }
            }
        }

        for (_, element) in self.elements.iter() {
            let [id_i, id_j] = element.nodes();
            let node_i = &self.nodes[id_i];
            let node_j = &self.nodes[id_j];
            let mass_local = element.form_mass(node_i, node_j);

            let mut equations = [None; ELEMENT_DOF];
            equations[..NDOF].copy_from_slice(&node_i.equation);
            equations[NDOF..].copy_from_slice(&node_j.equation);

            for (a, eq_a) in equations.iter().enumerate() {
                let Some(eq_a) = eq_a else { continue };
                mass[*eq_a] += mass_local[a];
            }
        }

        mass
    }

    /// Per-element tangent-stiffness triplets and the assembled internal
    /// resisting force, over free DOFs, at the current nodal displacement
    /// state. Duplicate `(row, col)` triplets (every DOF shared by more
    /// than one element) are summed by whoever consumes them — `faer`'s
    /// triplet constructor does this automatically.
    fn assemble_stiffness_triplets(&self) -> (Vec<Triplet<usize, usize, f64>>, DVector<f64>) {
        let n = self.num_free_dofs;
        let mut triplets = Vec::new();
        let mut resistance = DVector::<f64>::zeros(n);

        for (_, element) in self.elements.iter() {
            let [id_i, id_j] = element.nodes();
            let node_i = &self.nodes[id_i];
            let node_j = &self.nodes[id_j];
            let (k_local, r_local) = element.form_tangent_and_resistance(node_i, node_j);

            let mut equations = [None; ELEMENT_DOF];
            equations[..NDOF].copy_from_slice(&node_i.equation);
            equations[NDOF..].copy_from_slice(&node_j.equation);

            for (a, eq_a) in equations.iter().enumerate() {
                let Some(eq_a) = eq_a else { continue };
                resistance[*eq_a] += r_local[a];
                for (b, eq_b) in equations.iter().enumerate() {
                    let Some(eq_b) = eq_b else { continue };
                    triplets.push(Triplet::new(*eq_a, *eq_b, k_local[(a, b)]));
                }
            }
        }

        (triplets, resistance)
    }

    /// Assemble the global tangent stiffness (sparse) and internal
    /// resisting force over free DOFs only, at the current nodal
    /// displacement state. No load-pattern contribution — see
    /// `assemble_reference_load`.
    pub(crate) fn assemble_tangent_and_resistance(&self) -> (SparseMatrix, DVector<f64>) {
        let n = self.num_free_dofs;
        let (triplets, resistance) = self.assemble_stiffness_triplets();
        let k = SparseMatrix::try_new_from_triplets(n, n, &triplets)
            .expect("equation numbers are always in [0, num_free_dofs)");
        (k, resistance)
    }

    /// Assemble the total applied load (every `LoadPattern`'s nodal +
    /// element equivalent loads, each scaled by its own factor at
    /// `pseudo_time` — frozen patterns use their frozen value regardless of
    /// `pseudo_time`, see `LoadPattern::factor`) over free DOFs. Already
    /// fully scaled — unlike the single-pattern design this replaced,
    /// callers no longer multiply the result by an external load factor
    /// (see `form_tangent_and_residual`).
    pub(crate) fn assemble_reference_load(&self, pseudo_time: f64) -> DVector<f64> {
        self.assemble_load_with(|p| p.factor(pseudo_time))
    }

    /// `d(assemble_reference_load)/d(pseudo_time)` at `pseudo_time` — the
    /// unit-load probe `Integrator::DisplacementControl` needs (how much
    /// the *total* applied load changes per unit pseudo-time increment;
    /// frozen/`Constant`-series patterns contribute zero, see `LoadPattern::
    /// sensitivity`/`LoadSeries::slope`'s doc comments). With the single
    /// default `Linear { slope: 1.0 }` pattern (every model that hasn't
    /// added a second pattern), this is numerically identical to
    /// `assemble_reference_load` — no behavior change for single-pattern
    /// models.
    pub(crate) fn assemble_reference_load_sensitivity(&self, pseudo_time: f64) -> DVector<f64> {
        self.assemble_load_with(|p| p.sensitivity(pseudo_time))
    }

    /// Shared core of `assemble_reference_load`/`_sensitivity`: sum every
    /// pattern's nodal + element reference load, each scaled by whatever
    /// `pattern_factor` computes for that pattern (its current factor, or
    /// its sensitivity — the only difference between the two callers).
    fn assemble_load_with(&self, pattern_factor: impl Fn(&LoadPattern<NDOF, NId, E::Id, E::Load>) -> f64) -> DVector<f64> {
        let n = self.num_free_dofs;
        let mut load = DVector::<f64>::zeros(n);

        for (_, pattern) in self.load_patterns.iter() {
            let factor = pattern_factor(pattern);
            if factor == 0.0 {
                continue;
            }

            for (node_id, node) in self.nodes.iter() {
                let Some(nodal_load) = pattern.nodal_load(node_id) else { continue };
                for dof in 0..NDOF {
                    if let Some(eq) = node.equation[dof] {
                        load[eq] += factor * nodal_load[dof];
                    }
                }
            }

            for (element_id, element) in self.elements.iter() {
                let Some(element_load) = pattern.element_load(element_id) else { continue };
                let [id_i, id_j] = element.nodes();
                let node_i = &self.nodes[id_i];
                let node_j = &self.nodes[id_j];
                let load_local = element.form_load_vector(node_i, node_j, Some(element_load));

                let mut equations = [None; ELEMENT_DOF];
                equations[..NDOF].copy_from_slice(&node_i.equation);
                equations[NDOF..].copy_from_slice(&node_j.equation);

                for (a, eq_a) in equations.iter().enumerate() {
                    let Some(eq_a) = eq_a else { continue };
                    load[*eq_a] += factor * load_local[a];
                }
            }
        }

        load
    }

    /// Assemble the global tangent stiffness and unbalanced-force residual
    /// (`reference_load(pseudo_time) - internal_resistance`) over free DOFs
    /// only. Called once per Newton iteration by `Analysis::step`. Unlike
    /// before multi-pattern support, `pseudo_time` is *not* an external
    /// multiplier on a single reference-load vector — each `LoadPattern`
    /// scales its own contribution internally (`assemble_reference_load`),
    /// so this is just their difference from internal resistance.
    pub(crate) fn form_tangent_and_residual(&self, pseudo_time: f64) -> (SparseMatrix, DVector<f64>) {
        let (k, resistance) = self.assemble_tangent_and_resistance();
        let residual = self.assemble_reference_load(pseudo_time) - resistance;
        (k, residual)
    }

    /// A unit "influence vector" over free DOFs: `1.0` at every DOF whose
    /// local index is `dof_direction`, `0.0` elsewhere — Chopra's ι, what
    /// `TransientAnalysis`'s ground-motion excitation needs to turn a
    /// scalar ground acceleration into an effective nodal force vector
    /// (`-M·ι·ag(t)`, mass-proportional, not reference-load-proportional —
    /// see `GroundMotion`). Computed once at `TransientAnalysis::new` time,
    /// not every step (like `mass` itself).
    pub(crate) fn direction_incidence(&self, dof_direction: usize) -> DVector<f64> {
        let mut incidence = DVector::<f64>::zeros(self.num_free_dofs);
        for (_, node) in self.nodes.iter() {
            if let Some(eq) = node.equation[dof_direction] {
                incidence[eq] = 1.0;
            }
        }
        incidence
    }

    /// Commit every element's material(s) at the current (final, converged)
    /// nodal state — called once per successful `Analysis`/
    /// `TransientAnalysis` step, never during Newton iteration itself. See
    /// `Material`'s doc comment for the trial/commit design this closes
    /// the loop on.
    pub(crate) fn commit(&mut self) {
        for (_, element) in self.elements.iter_mut() {
            let [id_i, id_j] = element.nodes();
            // `self.nodes` and `self.elements` are disjoint fields, so
            // borrowing one immutably while iterating the other mutably
            // is fine.
            element.commit(&self.nodes[id_i], &self.nodes[id_j]);
        }
    }

    /// Scatter a global free-DOF displacement increment back onto nodes.
    pub(crate) fn apply_displacement_increment(&mut self, du: &DVector<f64>) {
        for (_, node) in self.nodes.iter_mut() {
            for dof in 0..NDOF {
                if let Some(eq) = node.equation[dof] {
                    node.displacement[dof] += du[eq];
                }
            }
        }
    }

    /// Gather free-DOF displacement/velocity/acceleration into dense
    /// vectors — `TransientAnalysis` (M6) state, unlike static `Analysis`
    /// which only ever mutates displacement incrementally.
    pub(crate) fn gather_displacement(&self) -> DVector<f64> {
        self.gather(|node, dof| node.displacement[dof])
    }
    pub(crate) fn gather_velocity(&self) -> DVector<f64> {
        self.gather(|node, dof| node.velocity[dof])
    }
    pub(crate) fn gather_acceleration(&self) -> DVector<f64> {
        self.gather(|node, dof| node.acceleration[dof])
    }
    fn gather(&self, get: impl Fn(&Node<NDIM, NDOF>, usize) -> f64) -> DVector<f64> {
        let mut v = DVector::<f64>::zeros(self.num_free_dofs);
        for (_, node) in self.nodes.iter() {
            for dof in 0..NDOF {
                if let Some(eq) = node.equation[dof] {
                    v[eq] = get(node, dof);
                }
            }
        }
        v
    }

    /// Replace (not increment — every `TransientAnalysis` step computes a
    /// full new state, not a correction) each free DOF's displacement,
    /// velocity, and acceleration.
    pub(crate) fn scatter_state(&mut self, u: &DVector<f64>, v: &DVector<f64>, a: &DVector<f64>) {
        for (_, node) in self.nodes.iter_mut() {
            for dof in 0..NDOF {
                if let Some(eq) = node.equation[dof] {
                    node.displacement[dof] = u[eq];
                    node.velocity[dof] = v[eq];
                    node.acceleration[dof] = a[eq];
                }
            }
        }
    }

    /// `K * v` for an arbitrary free-DOF vector `v` — used once, to compute
    /// `TransientAnalysis`'s initial acceleration from equilibrium (`M*a0 =
    /// F0 - K*u0 - C*v0`, and `C` has a `K*v0` term). Reuses the same
    /// triplets `assemble_tangent_and_resistance` builds its sparse matrix
    /// from, computed directly (no `faer` matrix construction needed for a
    /// single matrix-vector product).
    pub(crate) fn multiply_stiffness(&self, v: &DVector<f64>) -> DVector<f64> {
        let (triplets, _resistance) = self.assemble_stiffness_triplets();
        let mut result = DVector::<f64>::zeros(self.num_free_dofs);
        for t in &triplets {
            result[t.row] += t.val * v[t.col];
        }
        result
    }

    /// Newmark's effective dynamic system for one `TransientAnalysis` step:
    /// `K_eff = stiffness_coeff*K + mass_coeff*diag(mass)`, plus `K*damp_vector`
    /// (the stiffness-proportional-damping contribution to the effective
    /// load, `beta_k*K*(a4*u_n + a5*v_n + a6*a_n)` — see `TransientAnalysis::
    /// step` for where `stiffness_coeff`/`mass_coeff`/`damp_vector` come
    /// from). One traversal shared between both outputs, rather than the
    /// two separate `assemble_tangent_and_resistance` / `multiply_stiffness`
    /// calls that would otherwise be needed every step.
    pub(crate) fn assemble_newmark_system(
        &self,
        mass: &DVector<f64>,
        mass_coeff: f64,
        stiffness_coeff: f64,
        damp_vector: &DVector<f64>,
    ) -> (SparseMatrix, DVector<f64>) {
        let n = self.num_free_dofs;
        let (triplets, _resistance) = self.assemble_stiffness_triplets();

        let mut k_damp = DVector::<f64>::zeros(n);
        let mut eff_triplets = Vec::with_capacity(triplets.len() + n);
        for t in &triplets {
            k_damp[t.row] += t.val * damp_vector[t.col];
            eff_triplets.push(Triplet::new(t.row, t.col, stiffness_coeff * t.val));
        }
        for i in 0..n {
            eff_triplets.push(Triplet::new(i, i, mass_coeff * mass[i]));
        }

        let k_eff = SparseMatrix::try_new_from_triplets(n, n, &eff_triplets)
            .expect("equation numbers are always in [0, num_free_dofs)");
        (k_eff, k_damp)
    }
}

/// Tie the x-translation DOF of every node in `constrained` to `retained`'s
/// — the standard 2D-frame simplification for a rigid floor diaphragm
/// (Xara/OpenSees's `rigidDiaphragm` is inherently 3D: it ties in-plane
/// translations *and* rotation-about-axis through a lever arm to nodes in a
/// plane perpendicular to a given axis, which doesn't map onto a
/// single-plane 2D model). Equivalent to calling `equal_dof(retained, c,
/// &[0])` for each `c` in `constrained`.
///
/// Planar-only, in its own impl block rather than the shared generic one:
/// naively reusing this identity-tie approach for the spatial profile would
/// silently produce a wrong diaphragm (missing the lever-arm term) instead
/// of refusing to compile — see the spatial-architecture plan's
/// "Constraints" section. A real spatial rigid diaphragm needs a genuine
/// transformation matrix, not exposed yet.
impl Domain {
    pub fn rigid_diaphragm(&mut self, retained: NodeId, constrained: &[NodeId]) {
        for &c in constrained {
            self.equal_dof(retained, c, &[0]);
        }
    }
}

#[cfg(test)]
mod domain3_tests {
    use super::*;
    use crate::analysis::SparseSolver;
    use crate::model::{Material, SpatialDof, Truss3};

    /// A skew 3D truss cantilever: end load projected onto the bar axis
    /// should reproduce simple axial-bar elongation, exercised through the
    /// full `Domain3` assembly/equation-numbering path (not just `Truss3`
    /// in isolation, which `truss.rs`'s own tests already cover).
    #[test]
    fn cantilever_truss3_matches_analytical_axial_stiffness() {
        let mut domain = Domain3::new();
        let fixed = domain.add_node(
            crate::model::Node3::new([0.0, 0.0, 0.0])
                .fix(SpatialDof::Ux as usize)
                .fix(SpatialDof::Uy as usize)
                .fix(SpatialDof::Uz as usize)
                .fix(SpatialDof::Rx as usize)
                .fix(SpatialDof::Ry as usize)
                .fix(SpatialDof::Rz as usize),
        );
        let free = domain.add_node(
            crate::model::Node3::new([3.0, 4.0, 0.0])
                .fix(SpatialDof::Uz as usize)
                .fix(SpatialDof::Rx as usize)
                .fix(SpatialDof::Ry as usize)
                .fix(SpatialDof::Rz as usize),
        );

        let area = 2.0;
        let e = 1000.0;
        let length = 5.0; // 3-4-5 triangle
        domain.add_element(Element3::Truss3(Truss3::new(fixed, free, area, Material::Elastic { e })));

        let force = 100.0;
        domain.load_node(free, SpatialDof::Ux as usize, force * 3.0 / 5.0);
        domain.load_node(free, SpatialDof::Uy as usize, force * 4.0 / 5.0);

        domain.number_dofs();
        assert_eq!(domain.num_free_dofs(), 2);

        let (k, residual) = domain.form_tangent_and_residual(1.0);
        let du = SparseSolver::new().solve(&k, &residual).expect("well-posed system");
        domain.apply_displacement_increment(&du);

        let expected_elongation = force * length / (area * e);
        let node = domain.node(free);
        let axial_disp = node.displacement[SpatialDof::Ux as usize] * 3.0 / 5.0
            + node.displacement[SpatialDof::Uy as usize] * 4.0 / 5.0;
        assert!((axial_disp - expected_elongation).abs() < 1e-9);
    }
}
