mod domain;
mod element;
mod material;
mod node;

pub use domain::Domain;
pub use element::{Element, ElementId, Truss, ZeroLength};
pub use material::Material;
pub use node::{Node, NodeId};

/// Degrees of freedom per node. Fixed at 2 (2D: x, y translation) for M1's
/// scope (a 2D truss). Revisit if/when a 3D element lands.
pub const NDF: usize = 2;
