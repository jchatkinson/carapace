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

/// Planar, for the *linear part* of a trapezoid load that is zero at node i and `dwx`, `dwy` at node j
/// (force/length): cubic-Hermite and linear shape functions integrated against `w(s) = dw * s / l`.
/// Transverse: shears `3 dw l / 20` and `7 dw l / 20`, end moments `dw l^2 / 30` and `-dw l^2 / 20`
/// (same sign convention as `planar_local`); axial: `dwx l / 6` and `dwx l / 3`.
pub(super) fn planar_local_taper(l: f64, dwx: f64, dwy: f64) -> SVector<f64, 6> {
    SVector::<f64, 6>::from_row_slice(&[
        dwx * l / 6.0,
        3.0 * dwy * l / 20.0,
        dwy * l * l / 30.0,
        dwx * l / 3.0,
        7.0 * dwy * l / 20.0,
        -dwy * l * l / 20.0,
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The taper's equivalent nodal loads are statically equivalent to the load `dw * s / l`: the shears
    /// sum to its resultant `dw l / 2` and, with the end moments, reproduce its moment `dw l^2 / 3` about node i.
    #[test]
    fn planar_taper_is_statically_equivalent_to_the_load() {
        let (l, dwx, dwy) = (7.5, 2.0, -3.0);
        let f = planar_local_taper(l, dwx, dwy);
        assert!((f[0] + f[3] - dwx * l / 2.0).abs() < 1e-12);
        assert!((f[1] + f[4] - dwy * l / 2.0).abs() < 1e-12);
        assert!((f[2] + f[5] + f[4] * l - dwy * l * l / 3.0).abs() < 1e-12);
        // Equal ends reduce to the uniform formula: uniform + 0 taper, and a taper of a flat load is the difference.
        assert_eq!(planar_local_taper(l, 0.0, 0.0), SVector::<f64, 6>::zeros());
    }
}
