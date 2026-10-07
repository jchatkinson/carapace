//! Shape functions and reference derivatives of the 3-node triangle (natural
//! coordinates `xi, eta` on the unit triangle) and the 4-node quadrilateral
//! (`xi, eta` in `[-1, 1]`). Node order is counter-clockwise; quad node `k`
//! sits at the reference corner `QUAD4_CORNERS[k]`.

/// Reference corners of the bilinear quadrilateral, in node order.
pub const QUAD4_CORNERS: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];

/// `N_i(xi, eta)` for the triangle: `1 - xi - eta`, `xi`, `eta`.
pub fn tri3_shape(xi: f64, eta: f64) -> [f64; 3] {
    [1.0 - xi - eta, xi, eta]
}

/// `[dN_i/dxi, dN_i/deta]` for the triangle (constant).
pub fn tri3_dshape() -> [[f64; 2]; 3] {
    [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]
}

/// `N_i(xi, eta) = (1 + xi_i xi)(1 + eta_i eta) / 4`.
pub fn quad4_shape(xi: f64, eta: f64) -> [f64; 4] {
    std::array::from_fn(|i| {
        let [xi_i, eta_i] = QUAD4_CORNERS[i];
        0.25 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta)
    })
}

/// `[dN_i/dxi, dN_i/deta]` for the quadrilateral.
pub fn quad4_dshape(xi: f64, eta: f64) -> [[f64; 2]; 4] {
    std::array::from_fn(|i| {
        let [xi_i, eta_i] = QUAD4_CORNERS[i];
        [
            0.25 * xi_i * (1.0 + eta_i * eta),
            0.25 * eta_i * (1.0 + xi_i * xi),
        ]
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const POINTS: [[f64; 2]; 4] = [[0.1, 0.2], [0.3, 0.3], [-0.7, 0.4], [0.5, -0.9]];

    #[test]
    fn shape_functions_partition_unity_and_derivatives_sum_to_zero() {
        for [xi, eta] in POINTS {
            let tri = [xi.abs() * 0.4, eta.abs() * 0.4];
            assert!((tri3_shape(tri[0], tri[1]).iter().sum::<f64>() - 1.0).abs() < 1e-15);
            assert!((quad4_shape(xi, eta).iter().sum::<f64>() - 1.0).abs() < 1e-15);
            for axis in 0..2 {
                assert!(tri3_dshape().iter().map(|d| d[axis]).sum::<f64>().abs() < 1e-15);
                assert!(
                    quad4_dshape(xi, eta)
                        .iter()
                        .map(|d| d[axis])
                        .sum::<f64>()
                        .abs()
                        < 1e-15
                );
            }
        }
    }

    #[test]
    fn shape_functions_are_one_at_their_own_node_and_zero_elsewhere() {
        let tri_nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        for (k, [xi, eta]) in tri_nodes.into_iter().enumerate() {
            for (i, n) in tri3_shape(xi, eta).into_iter().enumerate() {
                assert_eq!(n, if i == k { 1.0 } else { 0.0 });
            }
        }
        for (k, [xi, eta]) in QUAD4_CORNERS.into_iter().enumerate() {
            for (i, n) in quad4_shape(xi, eta).into_iter().enumerate() {
                assert_eq!(n, if i == k { 1.0 } else { 0.0 });
            }
        }
    }

    #[test]
    fn quad_derivatives_match_central_finite_differences() {
        let h = 1e-6;
        for [xi, eta] in POINTS {
            let d = quad4_dshape(xi, eta);
            let (dxi_p, dxi_m) = (quad4_shape(xi + h, eta), quad4_shape(xi - h, eta));
            let (deta_p, deta_m) = (quad4_shape(xi, eta + h), quad4_shape(xi, eta - h));
            for i in 0..4 {
                assert!((d[i][0] - (dxi_p[i] - dxi_m[i]) / (2.0 * h)).abs() < 1e-9);
                assert!((d[i][1] - (deta_p[i] - deta_m[i]) / (2.0 * h)).abs() < 1e-9);
            }
        }
        let d = tri3_dshape();
        let (a, b) = (tri3_shape(0.3 + h, 0.2), tri3_shape(0.3 - h, 0.2));
        let (c, e) = (tri3_shape(0.3, 0.2 + h), tri3_shape(0.3, 0.2 - h));
        for i in 0..3 {
            assert!((d[i][0] - (a[i] - b[i]) / (2.0 * h)).abs() < 1e-9);
            assert!((d[i][1] - (c[i] - e[i]) / (2.0 * h)).abs() < 1e-9);
        }
    }
}
