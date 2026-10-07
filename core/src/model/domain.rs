use std::collections::HashSet;

use faer::sparse::Triplet;
use nalgebra::{DVector, SMatrix, SVector};
use slotmap::{Key, SlotMap};

use super::constraint::{self, LinearConstraint, ResolveError};
use super::dof_table::{DofEntry, DofTable};
use super::load_pattern::{
    active_element_patterns, effective_element_load, ElementLoadComponents, LoadPattern,
};
use super::{
    Axis3, DofRef, Element, Element3, ElementForce, ElementOps, LoadPatternId, LoadSeries,
    ModelError, Node, Node3Id, NodeId, NodeView, SparseMatrix, TangentSink, VectorSink, NDF,
    PLANAR_NDIM, SPATIAL_NDF, SPATIAL_NDIM,
};

/// Every `(equation, coefficient)` term one node-dof contributes to the
/// assembled free-DOF system: none for a fixed or inactive dof, a single
/// `(eq, 1.0)` for an ordinary free or identity-tied dof, or any number of
/// terms for a constrained dof (a slave's value is a linear combination of
/// free unknowns). Borrowed from the `DofTable`'s term arena, so there is no
/// cap on the number of terms and no allocation in the assembly loop.
#[derive(Debug, Clone, Copy)]
enum DofTerms<'a> {
    None,
    One([(usize, f64); 1]),
    Many(&'a [(usize, f64)]),
}

impl DofTerms<'_> {
    fn as_slice(&self) -> &[(usize, f64)] {
        match self {
            DofTerms::None => &[],
            DofTerms::One(term) => term,
            DofTerms::Many(terms) => terms,
        }
    }

    fn iter(&self) -> std::slice::Iter<'_, (usize, f64)> {
        self.as_slice().iter()
    }
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
    NId = NodeId,
    E = Element,
> where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    nodes: SlotMap<NId, Node<NDIM, NDOF>>,
    elements: SlotMap<E::Id, E>,
    constraints: Vec<LinearConstraint<NId>>,
    /// The first structural problem `number_dofs` found in `constraints`
    /// (duplicate slave, cycle, fixed slave); reported by `validate`.
    constraint_error: Option<ModelError>,
    load_patterns: SlotMap<LoadPatternId, LoadPattern<NDOF, NId, E::Id, E::Load>>,
    default_pattern: LoadPatternId,
    num_free_dofs: usize,
    /// Pseudo-time of the last `commit` — what analyses that run on the
    /// committed state without a time of their own (modal) evaluate element
    /// loads at.
    committed_time: f64,
    /// How every node DOF maps onto the free-DOF system (equation number,
    /// or resolved terms for constrained DOFs),
    /// rebuilt by `number_dofs`.
    dofs: DofTable<NId, NDOF>,
}

/// `Domain`'s spatial instantiation — six-DOF `Node3`/`Element3`. See
/// `Domain`'s doc comment for why this is a type alias over one generic
/// implementation rather than a hand-duplicated struct.
pub type Domain3 = Domain<SPATIAL_NDIM, SPATIAL_NDF, Node3Id, Element3>;

/// Manual rather than `#[derive(Clone)]`: the derive macro only adds `E:
/// Clone`, not `E::Load: Clone` (it can't see through the associated type),
/// so it under-constrains this struct and fails to compile at every call
/// site instead. `Analysis::step`'s snapshot-and-restore-on-failure needs
/// this — see `Domain`'s doc comment.
impl<const NDIM: usize, const NDOF: usize, NId, E> Clone for Domain<NDIM, NDOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    fn clone(&self) -> Self {
        Domain {
            nodes: self.nodes.clone(),
            elements: self.elements.clone(),
            constraints: self.constraints.clone(),
            constraint_error: self.constraint_error,
            load_patterns: self.load_patterns.clone(),
            default_pattern: self.default_pattern,
            num_free_dofs: self.num_free_dofs,
            committed_time: self.committed_time,
            dofs: self.dofs.clone(),
        }
    }
}

impl<const NDIM: usize, const NDOF: usize, NId, E> Default for Domain<NDIM, NDOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<const NDIM: usize, const NDOF: usize, NId, E> Domain<NDIM, NDOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    /// A fresh, empty model — already carrying one `LoadPattern`
    /// (`default_pattern`, `LoadSeries::Linear { slope: 1.0 }`, unscaled)
    /// so simple single-pattern models (the common case) don't need to
    /// name a pattern at all: `domain.load_node(id, dof, value)` reaches
    /// for it implicitly. Reproduces, exactly, the single-pattern behavior
    /// every `Analysis` had before multi-pattern support existed.
    pub fn new() -> Self {
        let mut load_patterns = SlotMap::default();
        let default_pattern =
            load_patterns.insert(LoadPattern::new(LoadSeries::Linear { slope: 1.0 }));
        Domain {
            nodes: SlotMap::default(),
            elements: SlotMap::default(),
            constraints: Vec::new(),
            constraint_error: None,
            load_patterns,
            default_pattern,
            num_free_dofs: 0,
            committed_time: 0.0,
            dofs: DofTable::default(),
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

    fn node_view(&self) -> NodeView<'_, NDIM, NDOF, NId> {
        NodeView::new(&self.nodes)
    }

    /// An element's local nodal force at its current committed state — see
    /// `ElementOps::local_force`'s doc comment. Never called from the
    /// assembly hot loop (only from results recording), so looking the
    /// element and its nodes up here rather than caching anything is fine.
    pub fn element_local_force(&self, id: E::Id) -> ElementForce {
        self.elements[id].local_force(&self.node_view())
    }

    /// Component `component` of the load `id` carries at `pseudo_time` — every
    /// pattern's load on it, scaled by that pattern's factor (frozen patterns
    /// at their frozen value), summed; `0.0` for an unloaded element. Local
    /// axes: `ElementLoad` is `[wx, wy]`, `ElementLoad3` is `[wx, wy, wz]`.
    pub fn element_load_component(&self, id: E::Id, pseudo_time: f64, component: usize) -> f64 {
        let active = active_element_patterns(self.load_patterns.values(), pseudo_time);
        effective_element_load(&active, id).map_or(0.0, |load| load.component(component))
    }

    /// `element_local_force` including the fixed-end effect of element loads:
    /// the member end forces a loaded element actually carries
    /// (`k·d − Σ factor·f_eq`, OpenSees's `K·u + p0`), with every pattern's
    /// factor evaluated at `pseudo_time`. Identical to `element_local_force`
    /// for an element with no load.
    pub fn element_end_force(&self, id: E::Id, pseudo_time: f64) -> ElementForce {
        let element = &self.elements[id];
        let view = self.node_view();
        let mut force = element.local_force(&view);
        for (_, pattern) in self.load_patterns.iter() {
            let factor = pattern.factor(pseudo_time);
            if factor == 0.0 {
                continue;
            }
            if let Some(load) = pattern.element_load(id) {
                force.subtract_scaled(factor, &element.local_load_force(&view, Some(load)));
            }
        }
        force
    }

    /// Every integration point's per-fiber `(strain, stress)` at `id`'s
    /// current committed state — `None` for every element kind that isn't
    /// fiber-discretized (`ElementOps::fiber_responses`'s doc comment).
    /// Outer `Vec` (when `Some`) is one entry per integration point, inner
    /// `Vec` one entry per fiber, in the order each was originally given to
    /// the element's constructor — a fiber recorder resolves its `(section,
    /// fiber)` indices against that same order.
    pub fn element_fiber_responses(&self, id: E::Id) -> Option<Vec<Vec<(f64, f64)>>> {
        self.elements[id].fiber_responses(&self.node_view())
    }

    /// Support/equilibrium reaction at `(node, dof)`: net internal element
    /// resistance there minus the total applied reference load at
    /// `pseudo_time` (Chopra's standard reaction definition) — most useful
    /// at a *fixed* DOF (a free DOF's residual is, by construction, what
    /// `Algorithm::Newton` already drives to ~0 at convergence, so
    /// reading its "reaction" is redundant but harmless). Needs no
    /// equation number — unlike every other accessor here, this bypasses
    /// `dof_terms`/the free-DOF system entirely by summing every element's
    /// and pattern's own contribution to this one `(node, dof)` directly
    /// (the same per-element/per-pattern loops `assemble_stiffness_
    /// triplets`/`assemble_load_with` run, just for one DOF instead of
    /// every free one), so it's well-defined even at a DOF `number_dofs`
    /// never assigned an equation number to.
    ///
    /// Ground-motion excitation (`TransientAnalysis`) never needs a
    /// separate term here: it's an equivalent inertial pseudo-force at free
    /// DOFs only (`GroundMotion`'s own doc comment), never a real applied
    /// load at any DOF — and this crate's dynamics are relative-
    /// displacement (`docs/spatial-architecture.md`), so a fixed DOF's
    /// velocity/acceleration are always exactly zero and contribute nothing
    /// to its reaction regardless of profile or excitation.
    pub fn reaction(&self, node: NId, dof: usize, pseudo_time: f64) -> f64 {
        let view = self.node_view();
        let active = active_element_patterns(self.load_patterns.values(), pseudo_time);
        let touches = |element: &E| element.nodes().contains(&node);

        let mut resistance = PickDof::new(node, dof);
        for (element_id, element) in self.elements.iter() {
            if !touches(element) {
                continue;
            }
            let load = effective_element_load(&active, element_id);
            element.assemble_tangent(&view, load.as_ref(), &mut resistance);
        }

        let mut applied = 0.0;
        for (_, pattern) in self.load_patterns.iter() {
            let factor = pattern.factor(pseudo_time);
            if factor == 0.0 {
                continue;
            }
            if let Some(nodal_load) = pattern.nodal_load(node) {
                applied += factor * nodal_load[dof];
            }
            for (element_id, element) in self.elements.iter() {
                let Some(element_load) = pattern.element_load(element_id) else {
                    continue;
                };
                if !touches(element) {
                    continue;
                }
                let mut load = PickDof::new(node, dof);
                element.assemble_load(&view, Some(element_load), &mut load);
                applied += factor * load.total;
            }
        }

        resistance.total - applied
    }

    /// Imposes `value` as the displacement of a *fixed* (single-point
    /// constrained) DOF — a non-homogeneous SP constraint. A fixed DOF has
    /// no equation number, so the imposed value simply lives in the node's
    /// displacement and every element reads it like any other nodal
    /// displacement; the free DOFs then respond to it on the next
    /// `Analysis::step` (see `Analysis::step_prescribed`, which also makes
    /// the update atomic with respect to a failed step). `reaction` reads
    /// the force needed to hold it there. Returns `false`, changing
    /// nothing, if the DOF isn't fixed or `value` is nonfinite.
    pub fn prescribe_displacement(&mut self, node: NId, dof: usize, value: f64) -> bool {
        if !self.can_prescribe(node, dof, value) {
            return false;
        }
        self.nodes[node].displacement[dof] = value;
        true
    }

    /// Whether `prescribe_displacement` would accept these arguments.
    ///
    /// A DOF that a constraint uses as a master is refused: constraints are
    /// homogeneous and have already eliminated a fixed master's
    /// contribution, so moving it would silently leave its slaves behind.
    pub fn can_prescribe(&self, node: NId, dof: usize, value: f64) -> bool {
        value.is_finite()
            && dof < NDOF
            && self.nodes.get(node).is_some_and(|n| n.fixed[dof])
            && !self
                .constraints
                .iter()
                .any(|c| c.terms.iter().any(|&(n, d, _)| n == node && d == dof))
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
    pub fn add_load_pattern_scaled(
        &mut self,
        series: LoadSeries,
        scale_factor: f64,
    ) -> LoadPatternId {
        self.load_patterns
            .insert(LoadPattern::new(series).with_scale_factor(scale_factor))
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
    /// for the planar profile, `ElementLoad3` for the spatial one (see
    /// `ElementOps`'s doc comment).
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

    /// Constrains `slave`'s DOF to a linear combination of master DOFs on
    /// any nodes: `u_slave = sum(coeff * u_master)`. Repeated masters are
    /// merged and vanishing terms dropped; an empty result means
    /// `u_slave = 0`. Problems that cannot be seen from one constraint (a
    /// slave defined twice, a fixed slave, a dependency cycle, a master that
    /// carries a prescribed value, an inconsistent initial state) are
    /// reported by `validate`. Requires `ConstraintHandler::Transformation`.
    pub fn add_constraint(&mut self, slave: (NId, usize), terms: &[(NId, usize, f64)]) {
        assert!(
            slave.1 < NDOF,
            "constraint slave DOF {} out of range",
            slave.1
        );
        assert!(
            terms.iter().all(|t| t.1 < NDOF),
            "constraint master DOF out of range"
        );
        self.constraints.push(LinearConstraint {
            slave,
            terms: constraint::normalize(terms),
        });
    }

    /// Tie `dofs` of `constrained` exactly to the same DOFs of `retained`
    /// (`u_c = u_r` for each listed dof) — Xara/OpenSees's `equalDOF`.
    /// Requires `ConstraintHandler::Transformation`; see its doc comment
    /// for how this is resolved.
    pub fn equal_dof(&mut self, retained: NId, constrained: NId, dofs: &[usize]) {
        for &dof in dofs {
            self.add_constraint((constrained, dof), &[(retained, dof, 1.0)]);
        }
    }

    /// Equation number of a node's DOF **if it is a free unknown** (a DOF
    /// identity-tied to one shares its equation), or `None` for a fixed DOF,
    /// an inactive DOF, or a constrained DOF that is a combination of
    /// several unknowns. It says where an unknown lives, not what a DOF's
    /// value is: to read a node DOF out of a free-DOF-indexed vector use
    /// `value_at`. Used
    /// internally by `Integrator::DisplacementControl` to locate its
    /// controlled DOF in the free-DOF system, and `pub` (not just
    /// `pub(crate)`) so a caller holding a free-DOF-indexed vector — e.g.
    /// `modal_analysis`'s `Mode::shape` — can map one of its entries back
    /// to a specific node/DOF without `core` needing to expose any of
    /// `gather_displacement`/`scatter_state`'s broader internal-state
    /// surface just for that one lookup.
    pub fn equation_of(&self, node: NId, dof: usize) -> Option<usize> {
        match self.dofs.entry(node, dof) {
            DofEntry::Free(eq) => Some(eq),
            DofEntry::Fixed | DofEntry::Inactive | DofEntry::Affine { .. } => None,
        }
    }

    /// The value of `(node, dof)` in the free-DOF-indexed vector `q`: the
    /// entry for a free unknown, the constraint's combination of entries for
    /// a constrained DOF, and `0.0` for a fixed or inactive DOF. This is how
    /// a mode shape, an influence vector, or any other free-DOF vector is
    /// read back onto nodes; `equation_of` cannot do it for a constrained DOF
    /// that depends on several unknowns (a rigid-diaphragm slave).
    pub fn value_at(&self, q: &DVector<f64>, node: NId, dof: usize) -> f64 {
        match self.dofs.entry(node, dof) {
            DofEntry::Free(eq) => q[eq],
            DofEntry::Affine { start, len } => self
                .dofs
                .terms(start, len)
                .iter()
                .map(|&(eq, coeff)| coeff * q[eq])
                .sum(),
            DofEntry::Fixed | DofEntry::Inactive => 0.0,
        }
    }

    /// `true` for every free-DOF equation that is a rotation (local DOF
    /// index `>= NDIM` — `rz` in the planar profile, `rx`/`ry`/`rz` in the
    /// spatial one). Identity-tied DOFs share their retained equation and
    /// DOF index, and a rigid diaphragm's retained rotation keeps its own,
    /// so the classification is per independent coordinate. Used to pick
    /// moment rather than force tolerances and rotation rather than
    /// translation scales.
    pub(crate) fn rotational_equations(&self) -> Vec<bool> {
        let mut rotational = vec![false; self.num_free_dofs];
        for (id, _) in self.nodes.iter() {
            for dof in NDIM..NDOF {
                if let DofEntry::Free(eq) = self.dofs.entry(id, dof) {
                    rotational[eq] = true;
                }
            }
        }
        rotational
    }

    /// Whether any load pattern still follows an unfrozen
    /// `LoadSeries::Path` (see `LoadPattern::is_unfrozen_path`).
    pub(crate) fn has_unfrozen_path_series(&self) -> bool {
        self.load_patterns
            .iter()
            .any(|(_, pattern)| pattern.is_unfrozen_path())
    }

    /// Whether `(node, dof)` can take part in the analysis: it is fixed,
    /// constrained, stiffened by some element, used as a constraint master,
    /// or carries a nodal mass. A DOF that is none of these has no equation
    /// and nothing to resist a load. Reflects the last numbering
    /// (`validate`/`AnalysisBuilder::build` run it).
    pub fn is_active(&self, node: NId, dof: usize) -> bool {
        self.dofs.entry(node, dof) != DofEntry::Inactive
    }

    /// Numbers the DOFs and checks the model for mistakes that would
    /// otherwise be silently wrong: constraint problems (a slave defined
    /// twice or fixed, a dependency cycle, a master carrying a prescribed
    /// value, mass on a slave that combines several unknowns), an element
    /// rejecting its own geometry, nodal loads on DOFs nothing uses, and a
    /// constrained slave whose initial state contradicts its masters. Also
    /// makes the initial state consistent with the constraints: a slave
    /// whose initial displacement, velocity or acceleration is exactly zero
    /// (the default) takes the value its masters imply; a nonzero value that
    /// disagrees is an error, never overwritten. Called by
    /// `AnalysisBuilder::build`, `modal_analysis` and `TransientAnalysis`.
    pub fn validate(&mut self) -> Result<(), ModelError> {
        self.number_dofs();
        if let Some(error) = self.constraint_error {
            return Err(error);
        }
        let view = self.node_view();
        for (index, (_, element)) in self.elements.iter().enumerate() {
            element
                .validate(&view)
                .map_err(|reason| ModelError::InvalidElement {
                    element: index,
                    reason,
                })?;
        }
        for (index, (id, _)) in self.nodes.iter().enumerate() {
            for dof in 0..NDOF {
                if self.dofs.entry(id, dof) != DofEntry::Inactive {
                    continue;
                }
                let loaded = self
                    .load_patterns
                    .values()
                    .any(|pattern| pattern.nodal_load(id).is_some_and(|load| load[dof] != 0.0));
                if loaded {
                    return Err(ModelError::LoadOnInactiveDof { node: index, dof });
                }
            }
        }
        self.validate_constraint_state()
    }

    fn node_index(&self, id: NId) -> usize {
        self.nodes.keys().position(|key| key == id).unwrap_or(0)
    }

    /// The constraint-related half of `validate` that needs node state.
    fn validate_constraint_state(&mut self) -> Result<(), ModelError> {
        for constraint in &self.constraints {
            for &(node, dof, _) in &constraint.terms {
                // A fixed master's contribution has been eliminated, which is
                // only right while it stays zero.
                if self.nodes[node].fixed[dof] && self.nodes[node].displacement[dof] != 0.0 {
                    return Err(ModelError::ConstraintOnPrescribedDof {
                        node: self.node_index(node),
                        dof,
                    });
                }
            }
            let (slave, dof) = constraint.slave;
            let several_unknowns = matches!(
                self.dofs.entry(slave, dof),
                DofEntry::Affine { len, .. } if len > 1
            );
            if several_unknowns && self.nodes[slave].mass[dof] != 0.0 {
                return Err(ModelError::MassOnConstrainedDof {
                    node: self.node_index(slave),
                    dof,
                });
            }
        }

        // Project zero-valued slave state from the masters, in dependency
        // order (repeat until nothing changes; a chain settles in at most
        // one pass per link).
        for _ in 0..=self.constraints.len() {
            let mut changed = false;
            for index in 0..self.constraints.len() {
                let (slave, dof) = self.constraints[index].slave;
                let implied = self.implied_state(index);
                let node = &mut self.nodes[slave];
                let targets = [
                    &mut node.displacement[dof],
                    &mut node.velocity[dof],
                    &mut node.acceleration[dof],
                ];
                for (target, value) in targets.into_iter().zip(implied) {
                    if *target == 0.0 && value != 0.0 {
                        *target = value;
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }
        for index in 0..self.constraints.len() {
            let (slave, dof) = self.constraints[index].slave;
            let implied = self.implied_state(index);
            let node = &self.nodes[slave];
            let actual = [
                node.displacement[dof],
                node.velocity[dof],
                node.acceleration[dof],
            ];
            for (actual, implied) in actual.into_iter().zip(implied) {
                if (actual - implied).abs() > 1e-12 * (1.0 + actual.abs() + implied.abs()) {
                    return Err(ModelError::InconsistentInitialState {
                        node: self.node_index(slave),
                        dof,
                    });
                }
            }
        }
        Ok(())
    }

    /// The displacement, velocity and acceleration constraint `index`'s
    /// masters imply for its slave.
    fn implied_state(&self, index: usize) -> [f64; 3] {
        let mut implied = [0.0; 3];
        for &(node, dof, coeff) in &self.constraints[index].terms {
            let node = &self.nodes[node];
            implied[0] += coeff * node.displacement[dof];
            implied[1] += coeff * node.velocity[dof];
            implied[2] += coeff * node.acceleration[dof];
        }
        implied
    }

    /// Assign a sequential equation number to every free, unconstrained
    /// DOF, in node insertion order, then alias every multi-point-
    /// constrained DOF to its retained node's equation number for that DOF
    /// (`ConstraintHandler::Transformation`). Fixed
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
        let slaves: HashSet<(NId, usize)> = self.constraints.iter().map(|c| c.slave).collect();
        self.constraint_error = None;

        // A DOF is active if an element stiffens it, a constraint uses it as
        // a master (a diaphragm's retained node typically has no element at
        // all), or it carries a nodal mass. Everything else is not an
        // equation.
        let mut active: HashSet<(NId, usize)> = HashSet::new();
        for (_, element) in self.elements.iter() {
            let mask = element.dof_mask();
            for node in element.nodes() {
                for dof in (0..NDOF).filter(|&dof| mask.contains(dof)) {
                    active.insert((node, dof));
                }
            }
        }
        // A nodal mass gives its DOF inertia even with no stiffness (a free
        // mass integrates fine in a transient analysis), so it is active.
        for (id, node) in self.nodes.iter() {
            for dof in (0..NDOF).filter(|&dof| node.mass[dof] != 0.0) {
                active.insert((id, dof));
            }
        }
        for constraint in &self.constraints {
            for &(node, dof, _) in &constraint.terms {
                active.insert((node, dof));
            }
        }

        let mut table = DofTable::<NId, NDOF>::for_nodes(self.nodes.keys());
        let mut next = 0;
        for (index, (id, node)) in self.nodes.iter().enumerate() {
            for dof in 0..NDOF {
                if slaves.contains(&(id, dof)) {
                    if node.fixed[dof] {
                        self.constraint_error
                            .get_or_insert(ModelError::SlaveIsFixed { node: index, dof });
                    }
                    continue;
                }
                if node.fixed[dof] {
                    continue;
                }
                if !active.contains(&(id, dof)) {
                    table.set(id, dof, DofEntry::Inactive);
                    continue;
                }
                table.set(id, dof, DofEntry::Free(next));
                next += 1;
            }
        }

        if let Err(error) = constraint::resolve(&self.constraints, &mut table) {
            let (ResolveError::DuplicateSlave((id, dof)) | ResolveError::Cycle((id, dof))) = error;
            let node = self.nodes.keys().position(|k| k == id).unwrap_or(0);
            self.constraint_error.get_or_insert(match error {
                ResolveError::DuplicateSlave(_) => ModelError::DuplicateSlave { node, dof },
                ResolveError::Cycle(_) => ModelError::ConstraintCycle { node, dof },
            });
        }

        self.dofs = table;
        self.num_free_dofs = next;
        next
    }

    /// Pseudo-time of the last `commit` (zero before any).
    pub(crate) fn committed_time(&self) -> f64 {
        self.committed_time
    }

    pub fn num_free_dofs(&self) -> usize {
        self.num_free_dofs
    }

    /// Whether any `equal_dof`/`rigid_diaphragm` constraint has been added —
    /// used by `AnalysisBuilder<Ready>::build` to reject `ConstraintHandler
    /// ::Plain` (which can't resolve them) early, at model-construction
    /// time rather than as a silently-wrong solve. Also used by callers
    /// (e.g. `carapace-wasm`'s decoder) to pick `ConstraintHandler::
    /// Transformation` automatically whenever a domain actually has any.
    pub fn has_mp_constraints(&self) -> bool {
        !self.constraints.is_empty()
    }

    /// Every `(equation, coefficient)` term `node_id`'s `dof` contributes
    /// to the assembled free-DOF system — the single source every assembly
    /// hot loop (`assemble_stiffness_triplets`, `assemble_mass_diagonal`,
    /// `assemble_load_with`) reads instead of `node.equation[dof]`
    /// directly, so each of them handles ordinary free dofs, identity-tied
    /// dofs, and affine-tied (rigid-diaphragm) dofs uniformly rather than
    /// three separate code paths.
    fn dof_terms(&self, node_id: NId, dof: usize) -> DofTerms<'_> {
        match self.dofs.entry(node_id, dof) {
            DofEntry::Free(eq) => DofTerms::One([(eq, 1.0)]),
            DofEntry::Affine { start, len } => DofTerms::Many(self.dofs.terms(start, len)),
            DofEntry::Fixed | DofEntry::Inactive => DofTerms::None,
        }
    }

    /// Debug-build guard behind `ElementOps::dof_mask`: an element that
    /// produces a stiffness, resistance, load or mass on a DOF that is
    /// inactive (not fixed, not constrained, not stiffened by any element,
    /// not a constraint master, no nodal mass) has under-declared its mask, and the
    /// contribution would be silently dropped. Contributions below a
    /// relative 1e-10 of the element's largest are treated as zero (an
    /// orientation frame leaves round-off in entries that are zero in exact
    /// arithmetic).
    #[cfg(debug_assertions)]
    fn assert_declared<const N: usize>(
        &self,
        dofs: &[DofRef<NId>; N],
        k: Option<&SMatrix<f64, N, N>>,
        values: &[f64],
    ) {
        let scale = |it: &mut dyn Iterator<Item = f64>| it.fold(0.0_f64, |m, v| m.max(v.abs()));
        let vector_scale = scale(&mut values.iter().copied());
        let diagonal_scale = k.map_or(0.0, |k| scale(&mut (0..N).map(|i| k[(i, i)])));
        for (a, &(node, slot)) in dofs.iter().enumerate() {
            if self.dofs.entry(node, slot as usize) != DofEntry::Inactive {
                continue;
            }
            let stray_value = values[a].abs() > 1e-10 * vector_scale;
            let stray_stiffness = k.is_some_and(|k| k[(a, a)].abs() > 1e-10 * diagonal_scale);
            assert!(
                !(stray_value || stray_stiffness),
                "an element contributes to DOF {slot} of a node on which it did not declare a DOF in \
                 `ElementOps::dof_mask` (and nothing else uses that DOF), so the contribution would be \
                 silently dropped"
            );
        }
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

        for (id, node) in self.nodes.iter() {
            for dof in 0..NDOF {
                accumulate_mass(&mut mass, self.dof_terms(id, dof), node.mass[dof]);
            }
        }

        let view = self.node_view();
        let mut sink = MassSink {
            domain: self,
            mass: &mut mass,
        };
        for (_, element) in self.elements.iter() {
            element.assemble_mass(&view, &mut sink);
        }

        mass
    }

    /// Per-element tangent-stiffness triplets and the assembled internal
    /// resisting force, over free DOFs, at the current nodal displacement
    /// state. Duplicate `(row, col)` triplets (every DOF shared by more
    /// than one element) are summed by whoever consumes them — `faer`'s
    /// triplet constructor does this automatically.
    ///
    /// `pseudo_time` fixes the element loads that state-dependent elements
    /// (`ForceBeamColumn`) fold into their own resistance — see
    /// `ElementOps::form_tangent_and_resistance`.
    fn assemble_stiffness_triplets(
        &self,
        pseudo_time: f64,
    ) -> (Vec<Triplet<usize, usize, f64>>, DVector<f64>) {
        let n = self.num_free_dofs;
        let active = active_element_patterns(self.load_patterns.values(), pseudo_time);
        let view = self.node_view();
        let mut sink = TripletSink {
            domain: self,
            triplets: Vec::new(),
            resistance: DVector::<f64>::zeros(n),
        };
        for (element_id, element) in self.elements.iter() {
            let load = effective_element_load(&active, element_id);
            element.assemble_tangent(&view, load.as_ref(), &mut sink);
        }
        (sink.triplets, sink.resistance)
    }

    /// Assemble the global tangent stiffness (sparse) and internal
    /// resisting force over free DOFs only, at the current nodal
    /// displacement state. No load-pattern contribution — see
    /// `assemble_reference_load`.
    pub(crate) fn assemble_tangent_and_resistance(
        &self,
        pseudo_time: f64,
    ) -> (SparseMatrix, DVector<f64>) {
        let n = self.num_free_dofs;
        let (triplets, resistance) = self.assemble_stiffness_triplets(pseudo_time);
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
    fn assemble_load_with(
        &self,
        pattern_factor: impl Fn(&LoadPattern<NDOF, NId, E::Id, E::Load>) -> f64,
    ) -> DVector<f64> {
        let n = self.num_free_dofs;
        let mut load = DVector::<f64>::zeros(n);

        for (_, pattern) in self.load_patterns.iter() {
            let factor = pattern_factor(pattern);
            if factor == 0.0 {
                continue;
            }

            for (node_id, _node) in self.nodes.iter() {
                let Some(nodal_load) = pattern.nodal_load(node_id) else {
                    continue;
                };
                for (dof, &value) in nodal_load.iter().enumerate() {
                    for &(eq, coeff) in self.dof_terms(node_id, dof).iter() {
                        load[eq] += coeff * factor * value;
                    }
                }
            }

            let view = self.node_view();
            for (element_id, element) in self.elements.iter() {
                let Some(element_load) = pattern.element_load(element_id) else {
                    continue;
                };
                let mut sink = LoadSink {
                    domain: self,
                    factor,
                    load: &mut load,
                };
                element.assemble_load(&view, Some(element_load), &mut sink);
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
    pub(crate) fn form_tangent_and_residual(
        &self,
        pseudo_time: f64,
    ) -> (SparseMatrix, DVector<f64>) {
        let (k, resistance) = self.assemble_tangent_and_resistance(pseudo_time);
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
    pub fn direction_incidence(&self, dof_direction: usize) -> DVector<f64> {
        let mut incidence = DVector::<f64>::zeros(self.num_free_dofs);
        for (id, _) in self.nodes.iter() {
            if let DofEntry::Free(eq) = self.dofs.entry(id, dof_direction) {
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
    ///
    /// `pseudo_time` is the converged step's pseudo-time: a state-dependent
    /// element (`ForceBeamColumn`) commits the state that is consistent with
    /// the element load at that time.
    pub(crate) fn commit(&mut self, pseudo_time: f64) {
        self.committed_time = pseudo_time;
        let active = active_element_patterns(self.load_patterns.values(), pseudo_time);
        // `self.nodes`, `self.elements` and `self.load_patterns` are
        // disjoint fields, so the node view can borrow `nodes` while
        // `elements` is iterated mutably.
        let view = NodeView::new(&self.nodes);
        for (element_id, element) in self.elements.iter_mut() {
            let load = effective_element_load(&active, element_id);
            element.commit(&view, load.as_ref());
        }
    }

    /// Scatter a global free-DOF displacement increment back onto nodes.
    /// Every ordinarily-numbered dof (free, or identity-tied via
    /// an identity tie) reads its increment directly out of `du`; every
    /// constrained dof (a rigid diaphragm's lever-arm
    /// coupling) then gets its increment *derived* from the same `du`,
    /// applying the constraint's linear combination directly to the
    /// increment rather than to an absolute value — valid because the
    /// relation is linear and homogeneous (constant geometry-derived
    /// coefficients, no offset term), so it holds identically for an
    /// increment as it does for an absolute displacement.
    pub(crate) fn apply_displacement_increment(&mut self, du: &DVector<f64>) {
        for (id, node) in self.nodes.iter_mut() {
            for dof in 0..NDOF {
                match self.dofs.entry(id, dof) {
                    DofEntry::Free(eq) => node.displacement[dof] += du[eq],
                    DofEntry::Affine { start, len } => {
                        let delta: f64 = self
                            .dofs
                            .terms(start, len)
                            .iter()
                            .map(|&(eq, coeff)| coeff * du[eq])
                            .sum();
                        node.displacement[dof] += delta;
                    }
                    DofEntry::Fixed | DofEntry::Inactive => {}
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
        for (id, node) in self.nodes.iter() {
            for dof in 0..NDOF {
                if let DofEntry::Free(eq) = self.dofs.entry(id, dof) {
                    v[eq] = get(node, dof);
                }
            }
        }
        v
    }

    /// Replace (not increment — every `TransientAnalysis` step computes a
    /// full new state, not a correction) each free DOF's displacement,
    /// velocity, and acceleration. Affine-tied dofs (see
    /// `apply_displacement_increment`'s doc comment) are derived from `u`/
    /// `v`/`a` the same way — the linear relation holds identically for
    /// displacement, velocity, and acceleration.
    pub(crate) fn scatter_state(&mut self, u: &DVector<f64>, v: &DVector<f64>, a: &DVector<f64>) {
        for (id, node) in self.nodes.iter_mut() {
            for dof in 0..NDOF {
                match self.dofs.entry(id, dof) {
                    DofEntry::Free(eq) => {
                        node.displacement[dof] = u[eq];
                        node.velocity[dof] = v[eq];
                        node.acceleration[dof] = a[eq];
                    }
                    DofEntry::Affine { start, len } => {
                        let (mut du, mut dv, mut da) = (0.0, 0.0, 0.0);
                        for &(eq, coeff) in self.dofs.terms(start, len) {
                            du += coeff * u[eq];
                            dv += coeff * v[eq];
                            da += coeff * a[eq];
                        }
                        node.displacement[dof] = du;
                        node.velocity[dof] = dv;
                        node.acceleration[dof] = da;
                    }
                    DofEntry::Fixed | DofEntry::Inactive => {}
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
    pub(crate) fn multiply_stiffness(&self, v: &DVector<f64>, pseudo_time: f64) -> DVector<f64> {
        let (triplets, _resistance) = self.assemble_stiffness_triplets(pseudo_time);
        let mut result = DVector::<f64>::zeros(self.num_free_dofs);
        for t in &triplets {
            result[t.row] += t.val * v[t.col];
        }
        result
    }

    /// Newmark's effective dynamic operator and internal resistance, at the
    /// domain's *current* state (whatever trial displacement is on it when
    /// called): `K_eff = stiffness_coeff*K(u) + mass_coeff*diag(mass)`,
    /// `resistance = F_int(u)` — the `Op`/`F_int` half of `docs/
    /// algorithms.md` §6.1's Newton corrector residual (`TransientAnalysis::
    /// try_step` composes the rest — the damping/inertial terms, which need
    /// `v_trial`/`a_trial`, not just `u` — around this). `stiffness_coeff`/
    /// `mass_coeff` are `docs/algorithms.md` §6.1's `c1`/`c3` (Newmark's
    /// `a1 + a4*alpha_m`/`1 + a4*beta_k`).
    pub(crate) fn assemble_newmark_system(
        &self,
        mass: &DVector<f64>,
        mass_coeff: f64,
        stiffness_coeff: f64,
        pseudo_time: f64,
    ) -> (SparseMatrix, DVector<f64>) {
        let n = self.num_free_dofs;
        let (triplets, resistance) = self.assemble_stiffness_triplets(pseudo_time);

        let mut eff_triplets = Vec::with_capacity(triplets.len() + n);
        for t in &triplets {
            eff_triplets.push(Triplet::new(t.row, t.col, stiffness_coeff * t.val));
        }
        for i in 0..n {
            eff_triplets.push(Triplet::new(i, i, mass_coeff * mass[i]));
        }

        let k_eff = SparseMatrix::try_new_from_triplets(n, n, &eff_triplets)
            .expect("equation numbers are always in [0, num_free_dofs)");
        (k_eff, resistance)
    }
}

/// Adds a lumped-mass value on one node DOF to the diagonal mass vector.
/// A DOF that is a single scaled unknown (`u = c * q`) contributes
/// `c^2 * m`; one that combines several unknowns would need the full
/// `T^T M T`, which a diagonal mass vector cannot hold.
fn accumulate_mass(mass: &mut DVector<f64>, terms: DofTerms<'_>, value: f64) {
    if value == 0.0 {
        return;
    }
    match terms.as_slice() {
        [] => {}
        [(eq, coeff)] => mass[*eq] += coeff * coeff * value,
        _ => panic!(
            "assemble_mass_diagonal: mass on a DOF constrained to a combination of several unknowns (a \
             rigid diaphragm or rigid link slave) isn't supported — the lumped mass diagonal can't \
             represent the rotational inertia this would induce at the master; assign the mass at the \
             master node (or an unconstrained node) instead"
        ),
    }
}

/// Scatters element tangents into triplets and element resistances into the
/// free-DOF resistance vector, through each DOF's equation terms.
struct TripletSink<'a, const NDIM: usize, const NDOF: usize, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    domain: &'a Domain<NDIM, NDOF, NId, E>,
    triplets: Vec<Triplet<usize, usize, f64>>,
    resistance: DVector<f64>,
}

impl<const NDIM: usize, const NDOF: usize, NId, E> TangentSink<NId>
    for TripletSink<'_, NDIM, NDOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    fn add<const N: usize>(
        &mut self,
        dofs: &[DofRef<NId>; N],
        k: &SMatrix<f64, N, N>,
        r: &SVector<f64, N>,
    ) {
        #[cfg(debug_assertions)]
        self.domain.assert_declared(dofs, Some(k), r.as_slice());
        let dof_terms: [DofTerms<'_>; N] =
            std::array::from_fn(|a| self.domain.dof_terms(dofs[a].0, dofs[a].1 as usize));
        for (a, terms_a) in dof_terms.iter().enumerate() {
            for &(eq_a, coeff_a) in terms_a.iter() {
                self.resistance[eq_a] += coeff_a * r[a];
                for (b, terms_b) in dof_terms.iter().enumerate() {
                    for &(eq_b, coeff_b) in terms_b.iter() {
                        self.triplets
                            .push(Triplet::new(eq_a, eq_b, coeff_a * coeff_b * k[(a, b)]));
                    }
                }
            }
        }
    }
}

/// Scatters scaled element load vectors into the free-DOF load vector.
struct LoadSink<'a, const NDIM: usize, const NDOF: usize, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    domain: &'a Domain<NDIM, NDOF, NId, E>,
    factor: f64,
    load: &'a mut DVector<f64>,
}

impl<const NDIM: usize, const NDOF: usize, NId, E> VectorSink<NId>
    for LoadSink<'_, NDIM, NDOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    fn add<const N: usize>(&mut self, dofs: &[DofRef<NId>; N], v: &SVector<f64, N>) {
        #[cfg(debug_assertions)]
        self.domain.assert_declared(dofs, None, v.as_slice());
        for (a, &(node, slot)) in dofs.iter().enumerate() {
            for &(eq, coeff) in self.domain.dof_terms(node, slot as usize).iter() {
                self.load[eq] += coeff * self.factor * v[a];
            }
        }
    }
}

/// Scatters element lumped masses into the free-DOF mass diagonal.
struct MassSink<'a, const NDIM: usize, const NDOF: usize, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    domain: &'a Domain<NDIM, NDOF, NId, E>,
    mass: &'a mut DVector<f64>,
}

impl<const NDIM: usize, const NDOF: usize, NId, E> VectorSink<NId>
    for MassSink<'_, NDIM, NDOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    fn add<const N: usize>(&mut self, dofs: &[DofRef<NId>; N], v: &SVector<f64, N>) {
        #[cfg(debug_assertions)]
        self.domain.assert_declared(dofs, None, v.as_slice());
        for (a, &(node, slot)) in dofs.iter().enumerate() {
            accumulate_mass(self.mass, self.domain.dof_terms(node, slot as usize), v[a]);
        }
    }
}

/// Picks out the entry of an element's resistance or load vector that
/// belongs to one `(node, dof)` — what `Domain::reaction` needs, without
/// involving equation numbers at all.
struct PickDof<NId> {
    node: NId,
    dof: u8,
    total: f64,
}

impl<NId> PickDof<NId> {
    fn new(node: NId, dof: usize) -> Self {
        PickDof {
            node,
            dof: dof as u8,
            total: 0.0,
        }
    }
}

impl<NId: Key> TangentSink<NId> for PickDof<NId> {
    fn add<const N: usize>(
        &mut self,
        dofs: &[DofRef<NId>; N],
        _k: &SMatrix<f64, N, N>,
        r: &SVector<f64, N>,
    ) {
        for (a, &(node, slot)) in dofs.iter().enumerate() {
            if node == self.node && slot == self.dof {
                self.total += r[a];
            }
        }
    }
}

impl<NId: Key> VectorSink<NId> for PickDof<NId> {
    fn add<const N: usize>(&mut self, dofs: &[DofRef<NId>; N], v: &SVector<f64, N>) {
        for (a, &(node, slot)) in dofs.iter().enumerate() {
            if node == self.node && slot == self.dof {
                self.total += v[a];
            }
        }
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
/// of refusing to compile — see `Domain3::rigid_diaphragm` for the real
/// (affine, lever-arm-aware) spatial version below.
impl Domain {
    pub fn rigid_diaphragm(&mut self, retained: NodeId, constrained: &[NodeId]) {
        for &c in constrained {
            self.equal_dof(retained, c, &[0]);
        }
    }

    /// Rigid link: `slave` moves with `master` as a rigid body. Its rotation
    /// equals the master's and its translations equal the master's plus the
    /// small-rotation lever arm, `u_s = u_m + theta_m x d` with
    /// `d = pos(slave) - pos(master)` from the reference coordinates at the
    /// time of the call (`ux_s = ux_m - theta*dy`, `uy_s = uy_m + theta*dx`).
    /// Requires `ConstraintHandler::Transformation`.
    pub fn rigid_link(&mut self, master: NodeId, slave: NodeId) {
        let m = self.nodes[master].coords;
        let s = self.nodes[slave].coords;
        let (dx, dy) = (s[0] - m[0], s[1] - m[1]);
        self.add_constraint((slave, 0), &[(master, 0, 1.0), (master, 2, -dy)]);
        self.add_constraint((slave, 1), &[(master, 1, 1.0), (master, 2, dx)]);
        self.add_constraint((slave, 2), &[(master, 2, 1.0)]);
    }
}

/// Spatial rigid diaphragm: for each node in `constrained`, ties its two
/// in-plane translational dofs (the two coordinate axes other than
/// `normal`) to `retained`'s same in-plane translations plus a lever-arm
/// term from `retained`'s rotation about `normal` — a genuine
/// general linear constraint (`ConstraintHandler::Transformation`), not a further
/// identity alias like planar `Domain::rigid_diaphragm` (see its doc
/// comment for why identity aliasing can't express a real diaphragm's
/// lever-arm kinematics).
///
/// Deliberately narrower than Xara/OpenSees's `rigidDiaphragm`: only the
/// two in-plane translational dofs of each constrained node are tied.
/// Every rotational dof of the constrained node — including its rotation
/// about `normal` itself — and its out-of-plane translation are left
/// completely free, unconstrained by this call. A floor diaphragm is rigid
/// in its own plane; it doesn't impose a shared rotation on whatever's
/// attached to it, and a linear frame model gives it no out-of-plane
/// (vertical) stiffness of its own. This matches planar `rigid_diaphragm`'s
/// existing scope (which likewise only ties translation, never `rz`), just
/// generalized from "the one translational in-plane dof a 2D model has" to
/// "the two a 3D diaphragm plane has, with the lever-arm term a real 3D
/// diaphragm needs".
///
/// Kinematics (small-rotation rigid-body motion in the plane perpendicular
/// to `normal`, at retained node `r`, evaluated at constrained node `c`):
/// with `n̂` the unit normal and `d = pos(c) - pos(r)`, a rotation
/// `theta_n` about `n̂` displaces `c` by `theta_n * (n̂ × d)`. Writing `a =
/// (normal+1) % 3`, `b = (normal+2) % 3` for the other two axes (so `(n, a,
/// b)` is a right-handed triad):
///
/// ```text
/// u_c[a] = u_r[a] - theta_r[normal] * (coord[b]_c - coord[b]_r)
/// u_c[b] = u_r[b] + theta_r[normal] * (coord[a]_c - coord[a]_r)
/// ```
///
/// Cross-checked against Xara/OpenSees's `RigidDiaphragm` constraint
/// equations for the z-normal case (`perpDirn = 3`, their default: `Ux_c =
/// Ux_r - Rz_r*(Yc-Yr)`, `Uy_c = Uy_r + Rz_r*(Xc-Xr)`) — substituting
/// `normal = Z` (`a = X, b = Y`) into the formula above reproduces exactly
/// those two equations. Verified through the full `Domain3`/`Analysis3`
/// stack for rigid-rotation invariance and lever-arm displacement magnitude
/// (`core/tests/m20_rigid_diaphragm3.rs`), not just asserted by that
/// cross-check.
///
/// Lever arms are computed once here, from each node's fixed reference
/// `coords` at the time this is called — not re-evaluated from current
/// displacement, matching every other `Linear`/`Linear3` (small-
/// displacement) geometric transform's scope in this crate.
///
/// Planar-only, in its own impl block: `normal` requires a full 3D
/// coordinate/rotation set (`Axis3`, and a rotation-about-`normal` dof)
/// that only the spatial profile's `Node3`/`SpatialDof` provide.
impl Domain3 {
    /// `rigid_diaphragm_about` with `normal = Axis3::Y` — the "up" axis in
    /// this crate's y-up global convention (see spatial-architecture.md),
    /// so the diaphragm plane is the horizontal `x`-`z` plane: the common
    /// case (an ordinary floor diaphragm).
    pub fn rigid_diaphragm(&mut self, retained: Node3Id, constrained: &[Node3Id]) {
        self.rigid_diaphragm_about(retained, constrained, Axis3::Y);
    }

    /// `rigid_diaphragm` with an explicit diaphragm-plane normal — see this
    /// impl block's doc comment for the kinematics and scope.
    pub fn rigid_diaphragm_about(
        &mut self,
        retained: Node3Id,
        constrained: &[Node3Id],
        normal: Axis3,
    ) {
        let n = normal as usize;
        let a = (n + 1) % SPATIAL_NDIM;
        let b = (n + 2) % SPATIAL_NDIM;
        let rot_normal = SPATIAL_NDIM + n;

        let retained_coords = self.nodes[retained].coords;
        for &c in constrained {
            let coords = self.nodes[c].coords;
            let lever_a = coords[a] - retained_coords[a];
            let lever_b = coords[b] - retained_coords[b];
            self.add_constraint(
                (c, a),
                &[(retained, a, 1.0), (retained, rot_normal, -lever_b)],
            );
            self.add_constraint(
                (c, b),
                &[(retained, b, 1.0), (retained, rot_normal, lever_a)],
            );
        }
    }

    /// Rigid link: `slave` moves with `master` as a rigid body. All six slave
    /// DOFs are tied: its rotations equal the master's, and its translations
    /// equal the master's plus the small-rotation lever arm,
    /// `u_s = u_m + theta_m x d` with `d = pos(slave) - pos(master)` taken
    /// from the reference coordinates at the time of the call. Joins, for
    /// example, a beam end to the edge of a shell, or an eccentric load point
    /// to the member it loads. Requires `ConstraintHandler::Transformation`.
    pub fn rigid_link(&mut self, master: Node3Id, slave: Node3Id) {
        let m = self.nodes[master].coords;
        let s = self.nodes[slave].coords;
        let (dx, dy, dz) = (s[0] - m[0], s[1] - m[1], s[2] - m[2]);
        // Rotation dofs are 3 (rx), 4 (ry), 5 (rz).
        self.add_constraint(
            (slave, 0),
            &[(master, 0, 1.0), (master, 4, dz), (master, 5, -dy)],
        );
        self.add_constraint(
            (slave, 1),
            &[(master, 1, 1.0), (master, 5, dx), (master, 3, -dz)],
        );
        self.add_constraint(
            (slave, 2),
            &[(master, 2, 1.0), (master, 3, dy), (master, 4, -dx)],
        );
        for rotation in 3..6 {
            self.add_constraint((slave, rotation), &[(master, rotation, 1.0)]);
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
        domain.add_element(Element3::Truss3(Truss3::new(
            fixed,
            free,
            area,
            Material::Elastic { e },
        )));

        let force = 100.0;
        domain.load_node(free, SpatialDof::Ux as usize, force * 3.0 / 5.0);
        domain.load_node(free, SpatialDof::Uy as usize, force * 4.0 / 5.0);

        domain.number_dofs();
        assert_eq!(domain.num_free_dofs(), 2);

        let (k, residual) = domain.form_tangent_and_residual(1.0);
        let du = SparseSolver::new()
            .solve(&k, &residual)
            .expect("well-posed system");
        domain.apply_displacement_increment(&du);

        let expected_elongation = force * length / (area * e);
        let node = domain.node(free);
        let axial_disp = node.displacement[SpatialDof::Ux as usize] * 3.0 / 5.0
            + node.displacement[SpatialDof::Uy as usize] * 4.0 / 5.0;
        assert!((axial_disp - expected_elongation).abs() < 1e-9);
    }
}

#[cfg(test)]
mod numbering_tests {
    use super::*;
    use crate::model::{Material, Node3, SpatialDof, Truss, ZeroLength3};

    /// Equation numbers are assigned to free, unconstrained DOFs in node
    /// insertion order; an identity-tied DOF shares its retained DOF's
    /// equation; a fixed retained DOF leaves its tied DOF without one.
    #[test]
    fn identity_ties_alias_equations_and_fixed_retained_dofs_drop_out() {
        let mut domain = Domain::new();
        let a = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let b = domain.add_node(Node::new([1.0, 0.0]).fix(2));
        let c = domain.add_node(Node::new([2.0, 0.0]).fix(1).fix(2));
        let d = domain.add_node(Node::new([3.0, 0.0]).fix(2));
        domain.add_element(Element::Truss(Truss::new(
            a,
            b,
            1.0,
            Material::Elastic { e: 1.0 },
        )));
        domain.equal_dof(b, c, &[0]);
        // `a` is fully fixed, so `d`'s tied DOF has no equation.
        domain.equal_dof(a, d, &[1]);

        assert_eq!(domain.number_dofs(), 2);
        assert_eq!(domain.equation_of(a, 0), None);
        assert_eq!(domain.equation_of(b, 0), Some(0));
        assert_eq!(domain.equation_of(b, 1), Some(1));
        assert_eq!(
            domain.equation_of(c, 0),
            Some(0),
            "tied dof shares the retained equation"
        );
        assert_eq!(domain.equation_of(c, 1), None);
        // No element stiffens `d`'s ux, so it is not an equation at all.
        assert_eq!(domain.equation_of(d, 0), None);
        assert!(!domain.is_active(d, 0));
        assert_eq!(domain.equation_of(d, 1), None);
        assert_eq!(
            domain.dof_terms(c, 0).iter().copied().collect::<Vec<_>>(),
            vec![(0, 1.0)]
        );
        assert!(domain.dof_terms(d, 1).iter().next().is_none());
    }

    /// A 3D rigid diaphragm: the slave's in-plane translations are affine
    /// in the master's translation and rotation about the normal, with the
    /// lever-arm coefficients, and the slave has no equation of its own.
    #[test]
    fn rigid_diaphragm_slave_dofs_resolve_to_affine_terms() {
        let fixed_but = |free: &[SpatialDof]| {
            let mut node = Node3::new([0.0, 0.0, 0.0]);
            for dof in [
                SpatialDof::Ux,
                SpatialDof::Uy,
                SpatialDof::Uz,
                SpatialDof::Rx,
                SpatialDof::Ry,
                SpatialDof::Rz,
            ] {
                if !free.iter().any(|f| *f as usize == dof as usize) {
                    node = node.fix(dof as usize);
                }
            }
            node
        };
        let mut domain = Domain3::new();
        let master = domain.add_node(fixed_but(&[SpatialDof::Ux, SpatialDof::Uz, SpatialDof::Ry]));
        let mut slave_node = fixed_but(&[SpatialDof::Ux, SpatialDof::Uz]);
        slave_node.coords = [2.0, 0.0, 3.0];
        let slave = domain.add_node(slave_node);
        let ground = domain.add_node(fixed_but(&[]));
        domain.add_element(Element3::ZeroLength3(
            ZeroLength3::new(ground, slave).with_material(0, Material::Elastic { e: 1.0 }),
        ));
        domain.rigid_diaphragm(master, &[slave]);

        assert_eq!(domain.number_dofs(), 3);
        let (ux, uz, ry) = (
            domain.equation_of(master, SpatialDof::Ux as usize).unwrap(),
            domain.equation_of(master, SpatialDof::Uz as usize).unwrap(),
            domain.equation_of(master, SpatialDof::Ry as usize).unwrap(),
        );
        assert_eq!(domain.equation_of(slave, SpatialDof::Ux as usize), None);
        // Normal Y: a = Z, b = X, so u_c[Z] = u_r[Z] - theta_Y * (x_c - x_r)
        // and u_c[X] = u_r[X] + theta_Y * (z_c - z_r).
        let terms = |dof: SpatialDof| {
            domain
                .dof_terms(slave, dof as usize)
                .iter()
                .copied()
                .collect::<Vec<_>>()
        };
        assert_eq!(terms(SpatialDof::Ux), vec![(ux, 1.0), (ry, 3.0)]);
        assert_eq!(terms(SpatialDof::Uz), vec![(uz, 1.0), (ry, -2.0)]);
    }
}
