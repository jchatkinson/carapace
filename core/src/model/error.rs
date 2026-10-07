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
    /// Two constraints define the same slave DOF.
    DuplicateSlave { node: usize, dof: usize },
    /// A constraint's slave DOF is also fixed by the user.
    SlaveIsFixed { node: usize, dof: usize },
    /// A slave depends on itself, directly or through other slaves.
    ConstraintCycle { node: usize, dof: usize },
    /// A constraint master is fixed with a nonzero value: constraints are
    /// homogeneous and have already eliminated its contribution, so its
    /// slaves would not follow it. (Prescribing such a DOF later is refused
    /// by `Domain::can_prescribe` for the same reason.)
    ConstraintOnPrescribedDof { node: usize, dof: usize },
    /// A slave's nonzero initial displacement, velocity or acceleration
    /// disagrees with what its masters imply.
    InconsistentInitialState { node: usize, dof: usize },
    /// Nodal mass on a slave that combines several unknowns, which a
    /// diagonal mass matrix cannot represent.
    MassOnConstrainedDof { node: usize, dof: usize },
    /// An element load assigned to an element that cannot carry that kind of
    /// load (a body force on a beam, a beam load on a continuum element, an
    /// edge load past the element's last edge).
    IncompatibleElementLoad { element: usize },
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
            ModelError::DuplicateSlave { node, dof } => {
                write!(f, "DOF {dof} of node {node} is the slave of more than one constraint")
            }
            ModelError::SlaveIsFixed { node, dof } => write!(
                f,
                "DOF {dof} of node {node} is a constraint slave and is also fixed"
            ),
            ModelError::ConstraintCycle { node, dof } => write!(
                f,
                "DOF {dof} of node {node} depends on itself through its constraints"
            ),
            ModelError::ConstraintOnPrescribedDof { node, dof } => write!(
                f,
                "DOF {dof} of node {node} is a constraint master and is fixed at a nonzero value; \
                 constraints cannot follow a prescribed master"
            ),
            ModelError::InconsistentInitialState { node, dof } => write!(
                f,
                "DOF {dof} of node {node} is a constraint slave whose initial state disagrees with its masters"
            ),
            ModelError::MassOnConstrainedDof { node, dof } => write!(
                f,
                "DOF {dof} of node {node} is constrained to several unknowns and carries a nodal mass; \
                 assign the mass at a master node"
            ),
            ModelError::IncompatibleElementLoad { element } => write!(
                f,
                "element {element} cannot carry the kind of element load assigned to it"
            ),
            ModelError::InvalidElement { element, reason } => {
                write!(f, "element {element} is invalid: {reason}")
            }
        }
    }
}

impl std::error::Error for ModelError {}
