//! Isoparametric machinery shared by the continuum elements: shape
//! functions and their reference derivatives (`shape`), Gauss rules
//! (`quadrature`), the Jacobian and physical derivatives (`jacobian`), and
//! scale-aware geometry validation (`validate`). Everything here is
//! dimension-agnostic bookkeeping; what is specific to an element (the `B`
//! matrix shape, the constitutive reduction) stays in the element.

pub mod embed;
pub mod jacobian;
pub mod natural_coords;
pub mod quadrature;
pub mod shape;
pub mod validate;

pub use embed::embed_polyline;
pub use jacobian::{physical_derivatives, Jacobian};
pub use natural_coords::{
    quad4_contains, quad4_natural_coords, tri3_contains, tri3_natural_coords,
};
pub use quadrature::{EDGE_GAUSS_2, QUAD4_GAUSS_2X2, TRI3_GAUSS_1, TRI3_GAUSS_3};
pub use shape::{quad4_dshape, quad4_shape, tri3_dshape, tri3_shape, QUAD4_CORNERS};
pub use validate::{validate_quad, validate_triangle};
