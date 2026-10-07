//! The 2D isoparametric Jacobian: `J = [dx/dxi dy/dxi; dx/deta dy/deta]`
//! (row `a` is the derivative with respect to natural coordinate `a`), its
//! determinant and inverse, and physical shape-function derivatives.

/// `J`, `det J` and `J^{-1}` at one point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Jacobian {
    pub matrix: [[f64; 2]; 2],
    pub det: f64,
    /// `J^{-1}`; meaningless (non-finite entries) when `det` is zero.
    pub inverse: [[f64; 2]; 2],
}

impl Jacobian {
    /// `J[a][c] = sum_i dN_i/dxi_a * x_i[c]` for `NN` nodes with reference
    /// derivatives `dn` (`[dN/dxi, dN/deta]` per node).
    pub fn new<const NN: usize>(coords: &[[f64; 2]; NN], dn: &[[f64; 2]; NN]) -> Self {
        let mut j = [[0.0; 2]; 2];
        for (x, d) in coords.iter().zip(dn) {
            for a in 0..2 {
                for c in 0..2 {
                    j[a][c] += d[a] * x[c];
                }
            }
        }
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let inverse = [
            [j[1][1] / det, -j[0][1] / det],
            [-j[1][0] / det, j[0][0] / det],
        ];
        Jacobian {
            matrix: j,
            det,
            inverse,
        }
    }
}

/// `dN_i/dx_c = sum_a J^{-1}[c][a] dN_i/dxi_a`, one `[dN/dx, dN/dy]` per node.
pub fn physical_derivatives<const NN: usize>(
    jacobian: &Jacobian,
    dn: &[[f64; 2]; NN],
) -> [[f64; 2]; NN] {
    let inv = &jacobian.inverse;
    std::array::from_fn(|i| {
        [
            inv[0][0] * dn[i][0] + inv[0][1] * dn[i][1],
            inv[1][0] * dn[i][0] + inv[1][1] * dn[i][1],
        ]
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::continuum::shape::{quad4_dshape, tri3_dshape};

    #[test]
    fn affine_maps_have_constant_determinant_equal_to_the_area_scale() {
        // Triangle (0,0), (4,0), (1,3): area 6, reference area 1/2 -> detJ = 12.
        let tri = [[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]];
        let j = Jacobian::new(&tri, &tri3_dshape());
        assert!((j.det - 12.0).abs() < 1e-14);
        // Parallelogram (0,0), (2,0), (3,1.5), (1,1.5): area 3, reference area 4 -> 0.75 everywhere.
        let par = [[0.0, 0.0], [2.0, 0.0], [3.0, 1.5], [1.0, 1.5]];
        for [xi, eta] in [[0.0, 0.0], [0.6, -0.3], [-1.0, 1.0]] {
            let j = Jacobian::new(&par, &quad4_dshape(xi, eta));
            assert!((j.det - 0.75).abs() < 1e-14);
        }
    }

    #[test]
    fn physical_derivatives_reproduce_a_linear_field_gradient() {
        // u(x, y) = 3x - 2y + 1 sampled at the nodes of a skewed quad has gradient (3, -2) everywhere.
        let quad = [[0.0, 0.0], [2.0, 0.3], [2.4, 1.9], [-0.2, 1.5]];
        let u: Vec<f64> = quad.iter().map(|p| 3.0 * p[0] - 2.0 * p[1] + 1.0).collect();
        for [xi, eta] in [[0.0, 0.0], [0.57, -0.57], [-0.8, 0.2]] {
            let dn = quad4_dshape(xi, eta);
            let dx = physical_derivatives(&Jacobian::new(&quad, &dn), &dn);
            let grad: [f64; 2] =
                std::array::from_fn(|c| dx.iter().zip(&u).map(|(d, u)| d[c] * u).sum());
            assert!((grad[0] - 3.0).abs() < 1e-13 && (grad[1] + 2.0).abs() < 1e-13);
        }
    }

    #[test]
    fn inverse_times_matrix_is_identity() {
        let quad = [[0.0, 0.0], [2.0, 0.3], [2.4, 1.9], [-0.2, 1.5]];
        let j = Jacobian::new(&quad, &quad4_dshape(0.3, -0.4));
        for a in 0..2 {
            for b in 0..2 {
                let v: f64 = (0..2).map(|c| j.matrix[a][c] * j.inverse[c][b]).sum();
                assert!((v - if a == b { 1.0 } else { 0.0 }).abs() < 1e-14);
            }
        }
    }
}
