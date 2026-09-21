use nalgebra::{DMatrix, DVector};
use slotmap::SlotMap;

use super::{Element, ElementId, Node, NodeId, ELEMENT_DOF, NDF};

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

    /// Assemble the global tangent stiffness and internal resisting force
    /// over free DOFs only, at the current nodal displacement state. No
    /// load-pattern contribution — see `assemble_reference_load`.
    pub(crate) fn assemble_tangent_and_resistance(&self) -> (DMatrix<f64>, DVector<f64>) {
        let n = self.num_free_dofs;
        let mut k = DMatrix::<f64>::zeros(n, n);
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
                    k[(*eq_a, *eq_b)] += k_local[(a, b)];
                }
            }
        }

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
    pub(crate) fn form_tangent_and_residual(&self, load_factor: f64) -> (DMatrix<f64>, DVector<f64>) {
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
