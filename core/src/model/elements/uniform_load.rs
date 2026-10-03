//! Consistent (virtual-work) equivalent nodal loads for a uniform local-axis
//! load on a prismatic 2-node beam: cubic-Hermite transverse and linear axial
//! shape functions integrated against the load. This is purely a function of
//! the load and the member length — it does not depend on the element's
//! constitutive model or integration rule — so every beam element
//! (`ElasticBeamColumn`, `DispBeamColumn`, ...) shares it, as OpenSees's
//! fixed-end forces `p0` are shared across them.

use nalgebra::SVector;

use super::truss::SpatialElementVector;

/// Planar: local DOF order `[u, v, rz]` per node. `wx` axial, `wy` transverse
/// (force/length).
pub(super) fn planar_local(l: f64, wx: f64, wy: f64) -> SVector<f64, 6> {
    SVector::<f64, 6>::from_row_slice(&[
        wx * l / 2.0,
        wy * l / 2.0,
        wy * l * l / 12.0,
        wx * l / 2.0,
        wy * l / 2.0,
        -wy * l * l / 12.0,
    ])
}

/// Spatial: local DOF order `[u, v, w, rx, ry, rz]` per node. The `wy`-induced
/// terms (`[v, rz]`, indices `[1,5,7,11]`) equal the planar formula; the
/// `wz`-induced terms (`[w, ry]`, indices `[2,4,8,10]`) carry the rotation
/// sign flip from `ry = -dw/dx` (vs `rz = +dv/dx`) — verified against
/// Xara/OpenSees's `Beam3dUniformLoad` fixed-end forces.
pub(super) fn spatial_local(l: f64, wx: f64, wy: f64, wz: f64) -> SpatialElementVector {
    let mut local = SpatialElementVector::zeros();
    local[0] = wx * l / 2.0;
    local[6] = wx * l / 2.0;

    local[1] = wy * l / 2.0;
    local[5] = wy * l * l / 12.0;
    local[7] = wy * l / 2.0;
    local[11] = -wy * l * l / 12.0;

    local[2] = wz * l / 2.0;
    local[4] = -wz * l * l / 12.0;
    local[8] = wz * l / 2.0;
    local[10] = wz * l * l / 12.0;
    local
}
