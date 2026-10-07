mod constraint;
pub mod continuum;
mod dof_table;
mod domain;
mod elements;
mod error;
mod fiber_section;
mod integration;
mod load_pattern;
mod materials;
mod node;
mod transform;

pub use domain::{Domain, Domain3};
pub use elements::{
    two_node_dofs, DispBeamColumn, DispBeamColumn3, DofMask, DofRef, ElasticBeamColumn,
    ElasticBeamColumn3, Element, Element3, Element3Id, ElementForce, ElementId, ElementOps,
    ForceBeamColumn, ForceBeamColumn3, Friction, Friction3, NodeList, NodeView, Orientation,
    OrientationError, SpatialElementMatrix, SpatialElementVector, TangentSink, Truss, Truss3,
    VectorSink, ZeroLength, ZeroLength3, ZeroLengthSection, ZeroLengthSection3, MAX_ELEMENT_NODES,
};
pub use error::ModelError;
pub use fiber_section::{Fiber, Fiber3, FiberSection, FiberSection3};
pub use integration::BeamIntegration;
pub use load_pattern::{
    ElementLoad, ElementLoad3, ElementLoadComponents, LoadPatternId, LoadSeries,
};
pub use materials::{Material, Pinching4DmgCyc, Pinching4State};
pub use node::{
    Axis3, Node, Node2, Node3, Node3Id, NodeId, SpatialDof, PLANAR_NDIM, SPATIAL_ELEMENT_DOF,
    SPATIAL_NDF, SPATIAL_NDIM,
};
pub use transform::{GeomTransf, GeomTransf3};

/// Degrees of freedom per node: 2D translation (x, y) plus in-plane
/// rotation (z). Bumped from 2 to 3 at M3 to carry `ElasticBeamColumn`'s
/// bending DOF. A node's DOF becomes an equation only if something uses it
/// (an element stiffens it, a constraint ties to it, it carries a mass, or
/// it is fixed), so a node connected only to trusses has no rotation
/// unknown and needs none fixed by hand (see `Domain::is_active`).
pub const NDF: usize = 3;

/// `2 * NDF` — the size of a 2-node element's local/global stiffness and
/// force vectors. A named constant instead of a `const` expression because
/// stable Rust's const generics can't yet compute `2 * NDF` inline at each
/// `SMatrix`/`SVector` use site.
pub const ELEMENT_DOF: usize = 2 * NDF;

/// The global tangent stiffness matrix's representation: sparse, not a
/// dense `nalgebra::DMatrix` — resolves implementation-plan §5 decision #1
/// for real. `faer` over `nalgebra-sparse`: a built-in sparse LU with
/// COLAMD/AMD fill-reducing ordering (so no separate `DOF_Numberer`
/// abstraction is needed, per the original §5 question), pure Rust (no C
/// dependency to fight through `wasm32-unknown-unknown`), and confirmed to
/// build and solve correctly under wasm32 + Node with a minimal feature set
/// (`std` + `sparse-linalg`, dropping `rand`/`rayon`/`npy`, which otherwise
/// pull in a `getrandom` wasm build failure).
pub(crate) type SparseMatrix = faer::sparse::SparseColMat<usize, f64>;
