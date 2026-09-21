use faer::sparse::Triplet;
use nalgebra::DVector;
use slotmap::SlotMap;

use super::{Element, ElementId, Node, NodeId, SparseMatrix, ELEMENT_DOF, NDF};

/// Owns all nodes and elements. No serialization/broker machinery (§3.3) —
/// this is the whole model, in memory, for one worker.
#[derive(Debug, Default)]
pub struct Domain {
    nodes: SlotMap<NodeId, Node>,
    elements: SlotMap<ElementId, Element>,
    num_free_dofs: usize,
}

impl Domain {
    pub fn new() -> Self {
        Domain::default()
    }

    pub fn add_node(&mut self, node: Node) -> NodeId {
        self.nodes.insert(node)
    }

    pub fn add_element(&mut self, element: Element) -> ElementId {
        self.elements.insert(element)
    }

    pub fn node(&self, id: NodeId) -> &Node {
        &self.nodes[id]
    }

    /// Equation number of a node's DOF, or `None` if it's fixed. Used by
    /// `Integrator::DisplacementControl` to locate its controlled DOF in
    /// the free-DOF system.
    pub(crate) fn equation_of(&self, node: NodeId, dof: usize) -> Option<usize> {
        self.nodes[node].equation[dof]
    }

    /// Assign a sequential equation number to every free DOF, in node
    /// insertion order. Fixed DOFs get no equation number. This is the
    /// entire numbering pass — done once, not re-checked every step (§5.4:
    /// no live re-solve means no `hasDomainChanged()`-style machinery).
    ///
    /// Returns the number of free DOFs (the size of the global system).
    pub(crate) fn number_dofs(&mut self) -> usize {
        let mut next = 0;
        for (_, node) in self.nodes.iter_mut() {
            for dof in 0..NDF {
                node.equation[dof] = if node.fixed[dof] {
                    None
                } else {
                    let eq = next;
                    next += 1;
                    Some(eq)
                };
            }
        }
        self.num_free_dofs = next;
        next
    }

    pub fn num_free_dofs(&self) -> usize {
        self.num_free_dofs
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
            for dof in 0..NDF {
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
            equations[..NDF].copy_from_slice(&node_i.equation);
            equations[NDF..].copy_from_slice(&node_j.equation);

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
            equations[..NDF].copy_from_slice(&node_i.equation);
            equations[NDF..].copy_from_slice(&node_j.equation);

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

    /// Assemble the reference load pattern (nodal loads + element
    /// equivalent loads) over free DOFs, unscaled by any load factor —
    /// what `Integrator`s scale by `load_factor` to get the external load
    /// side of equilibrium.
    pub(crate) fn assemble_reference_load(&self) -> DVector<f64> {
        let n = self.num_free_dofs;
        let mut load = DVector::<f64>::zeros(n);

        for (_, node) in self.nodes.iter() {
            for dof in 0..NDF {
                if let Some(eq) = node.equation[dof] {
                    load[eq] += node.load[dof];
                }
            }
        }

        for (_, element) in self.elements.iter() {
            let [id_i, id_j] = element.nodes();
            let node_i = &self.nodes[id_i];
            let node_j = &self.nodes[id_j];
            let load_local = element.form_load_vector(node_i, node_j);

            let mut equations = [None; ELEMENT_DOF];
            equations[..NDF].copy_from_slice(&node_i.equation);
            equations[NDF..].copy_from_slice(&node_j.equation);

            for (a, eq_a) in equations.iter().enumerate() {
                let Some(eq_a) = eq_a else { continue };
                load[*eq_a] += load_local[a];
            }
        }

        load
    }

    /// Assemble the global tangent stiffness and unbalanced-force residual
    /// (`load_factor * reference_load - internal_resistance`) over free
    /// DOFs only. Called once per Newton iteration by `Analysis::step`.
    pub(crate) fn form_tangent_and_residual(&self, load_factor: f64) -> (SparseMatrix, DVector<f64>) {
        let (k, resistance) = self.assemble_tangent_and_resistance();
        let residual = load_factor * self.assemble_reference_load() - resistance;
        (k, residual)
    }

    /// Scatter a global free-DOF displacement increment back onto nodes.
    pub(crate) fn apply_displacement_increment(&mut self, du: &DVector<f64>) {
        for (_, node) in self.nodes.iter_mut() {
            for dof in 0..NDF {
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
    fn gather(&self, get: impl Fn(&Node, usize) -> f64) -> DVector<f64> {
        let mut v = DVector::<f64>::zeros(self.num_free_dofs);
        for (_, node) in self.nodes.iter() {
            for dof in 0..NDF {
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
            for dof in 0..NDF {
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
