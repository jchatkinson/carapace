use std::fmt;

/// A model that cannot be analyzed as built, found by `Domain::validate`
/// before any equation is assembled. Node and element indices are insertion
/// order (the order the objects were added to the `Domain`), so a caller
/// that created them from a table can map them straight back to a row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ModelError {
    /// A nodal load is applied to a DOF that no element stiffens, that is
    /// not fixed, that carries no mass and that no constraint uses: the load
    /// would have nothing to push against and would be silently dropped.
    LoadOnInactiveDof { node: usize, dof: usize },
    /// An element rejected its own geometry or parameters.
    InvalidElement {
        element: usize,
        reason: &'static str,
    },
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModelError::LoadOnInactiveDof { node, dof } => write!(
                f,
                "a nodal load is applied to DOF {dof} of node {node}, which no element stiffens \
                 (and which is not fixed or used by a constraint)"
            ),
            ModelError::InvalidElement { element, reason } => {
                write!(f, "element {element} is invalid: {reason}")
            }
        }
    }
}

impl std::error::Error for ModelError {}
