//! Helpers shared by the 2D continuum elements (`Tri3`, `Quad4`): DOF
//! layout, nodal displacement and coordinate gathering, and the strain-
//! displacement matrix.

use nalgebra::{SMatrix, SVector};

use super::{DofRef, NodeView};
use crate::model::continuum::EDGE_GAUSS_2;
use crate::model::{ElementLoad, NodeId, NDF, PLANAR_NDIM};

pub(super) type PlaneView<'a> = NodeView<'a, PLANAR_NDIM, NDF, NodeId>;

/// `[(n0, ux), (n0, uy), (n1, ux), ...]`; `N` must be `2 * NN`.
pub(super) fn plane_dofs<const NN: usize, const N: usize>(
    nodes: &[NodeId; NN],
) -> [DofRef<NodeId>; N] {
    debug_assert_eq!(N, 2 * NN);
    std::array::from_fn(|a| (nodes[a / 2], (a % 2) as u8))
}

pub(super) fn plane_coords<const NN: usize>(
    nodes: &[NodeId; NN],
    view: &PlaneView<'_>,
) -> [[f64; 2]; NN] {
    std::array::from_fn(|i| view.get(nodes[i]).coords)
}

/// The element's nodal translations, `[ux0, uy0, ux1, ...]`.
pub(super) fn plane_displacement<const NN: usize, const N: usize>(
    nodes: &[NodeId; NN],
    view: &PlaneView<'_>,
) -> SVector<f64, N> {
    SVector::from_fn(|a, _| view.get(nodes[a / 2]).displacement[a % 2])
}

/// `B` (strain `[eps_x, eps_y, gamma_xy]` per nodal translation) from the
/// physical shape-function derivatives `[dN/dx, dN/dy]`.
pub(super) fn b_matrix<const NN: usize, const N: usize>(dx: &[[f64; 2]; NN]) -> SMatrix<f64, 3, N> {
    debug_assert_eq!(N, 2 * NN);
    let mut b = SMatrix::<f64, 3, N>::zeros();
    for (i, d) in dx.iter().enumerate() {
        b[(0, 2 * i)] = d[0];
        b[(1, 2 * i + 1)] = d[1];
        b[(2, 2 * i)] = d[1];
        b[(2, 2 * i + 1)] = d[0];
    }
    b
}

pub(super) fn check_section(
    thickness: f64,
    density: f64,
    material: &crate::model::PlaneMaterial,
) -> Result<(), &'static str> {
    if !thickness.is_finite() || thickness <= 0.0 {
        return Err("thickness must be positive");
    }
    if !density.is_finite() || density < 0.0 {
        return Err("density must be non-negative");
    }
    material.validate().map_err(|error| error.message())
}

/// Consistent nodal forces of a continuum element's load, in global DOF order.
/// `weights[i]` is `t * integral(N_i dA)` for node `i` (the body-force share, which
/// is also the lumped-mass share per unit density); edge loads are integrated on the
/// edge with the 2-point Gauss rule and the edge's length scale, edge `k` joining
/// local node `k` to `k + 1`. A pressure acts along the inward normal, which for
/// counter-clockwise nodes is the left normal of the edge direction.
pub(super) fn load_vector<const NN: usize, const N: usize>(
    coords: &[[f64; 2]; NN],
    weights: &[f64; NN],
    thickness: f64,
    load: &ElementLoad,
) -> SVector<f64, N> {
    debug_assert_eq!(N, 2 * NN);
    let mut f = SVector::<f64, N>::zeros();
    for (i, w) in weights.iter().enumerate() {
        f[2 * i] += w * load.body[0];
        f[2 * i + 1] += w * load.body[1];
    }
    for edge in 0..NN {
        let (a, b) = (edge, (edge + 1) % NN);
        let (dx, dy) = (coords[b][0] - coords[a][0], coords[b][1] - coords[a][1]);
        let length = dx.hypot(dy);
        let p = load.edge_pressure[edge];
        let traction = [
            load.edge_traction[edge][0] - p * dy / length,
            load.edge_traction[edge][1] + p * dx / length,
        ];
        if traction == [0.0, 0.0] {
            continue;
        }
        for &(s, w) in &EDGE_GAUSS_2 {
            let scale = w * 0.5 * length * thickness;
            for (node, shape) in [(a, 0.5 * (1.0 - s)), (b, 0.5 * (1.0 + s))] {
                f[2 * node] += scale * shape * traction[0];
                f[2 * node + 1] += scale * shape * traction[1];
            }
        }
    }
    f
}
