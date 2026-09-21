mod beam;
mod domain;
mod element;
mod material;
mod node;
mod transform;

pub use beam::ElasticBeamColumn;
pub use domain::Domain;
pub use element::{Element, ElementId, Truss, ZeroLength};
pub use material::Material;
pub use node::{Node, NodeId};
pub use transform::GeomTransf;

/// Degrees of freedom per node: 2D translation (x, y) plus in-plane
/// rotation (z). Bumped from 2 to 3 at M3 to carry `ElasticBeamColumn`'s
/// bending DOF — `Truss`/`ZeroLength` only ever populate the translational
/// entries, so existing models need the rotation DOF fixed at any node not
/// otherwise connected to a beam-column element (see M3 note in the
/// implementation plan), or the global system is singular in that DOF.
pub const NDF: usize = 3;

/// `2 * NDF` — the size of a 2-node element's local/global stiffness and
/// force vectors. A named constant instead of a `const` expression because
/// stable Rust's const generics can't yet compute `2 * NDF` inline at each
/// `SMatrix`/`SVector` use site.
pub const ELEMENT_DOF: usize = 2 * NDF;
