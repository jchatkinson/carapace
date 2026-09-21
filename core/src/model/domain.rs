use nalgebra::{DMatrix, DVector};
use slotmap::SlotMap;

use super::{Element, ElementId, Node, NodeId, NDF};

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

    /// Assemble the global tangent stiffness and unbalanced-force residual
    /// (`load_factor * external_load - internal_resistance`) over free DOFs
    /// only. Called once per Newton iteration by `Analysis::step`.
    pub(crate) fn form_tangent_and_residual(&self, load_factor: f64) -> (DMatrix<f64>, DVector<f64>) {
        let n = self.num_free_dofs;
        let mut k = DMatrix::<f64>::zeros(n, n);
        let mut residual = DVector::<f64>::zeros(n);

        for (_, node) in self.nodes.iter() {
            for dof in 0..NDF {
                if let Some(eq) = node.equation[dof] {
                    residual[eq] += load_factor * node.load[dof];
                }
            }
        }

        for (_, element) in self.elements.iter() {
            let [id_i, id_j] = element.nodes();
            let node_i = &self.nodes[id_i];
            let node_j = &self.nodes[id_j];
            let (k_local, r_local) = element.form_tangent_and_resistance(node_i, node_j);

            let equations: [Option<usize>; 4] = [
                node_i.equation[0],
                node_i.equation[1],
                node_j.equation[0],
                node_j.equation[1],
            ];

            for (a, eq_a) in equations.iter().enumerate() {
                let Some(eq_a) = eq_a else { continue };
                residual[*eq_a] -= r_local[a];
                for (b, eq_b) in equations.iter().enumerate() {
                    let Some(eq_b) = eq_b else { continue };
                    k[(*eq_a, *eq_b)] += k_local[(a, b)];
                }
            }
        }

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
