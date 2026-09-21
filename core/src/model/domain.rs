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
    /// lumped-mass matrix's diagonal over free DOFs — needed by modal
    /// analysis (M5). No element-consistent mass matrices yet (M6); entirely
    /// `Node::mass`. Just the diagonal, not a full (dense or sparse) N×N
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
}
