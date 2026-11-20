//! Natural coordinates of a physical point inside a T3 or Q4 element: a
//! closed-form inverse for the (affine) triangle, and a Newton inverse for
//! the (generally non-affine) bilinear quadrilateral — used by
//! `embed_polyline` to bind an embedded-bar segment to its host.
//!
//! This is deliberately not a port of `postprocess.ts`'s
//! `naturalCoordinatesAt`: that function returns after 10 iterations
//! without checking that the Newton loop actually converged or that the
//! coordinates it returns are finite (plan `docs/csfm-plan.md` §2.3).
//! Both inverses here require finite output and a physical residual
//! (`|x(xi, eta) - point|`) below `physical_tol`, and report
//! non-convergence as an error instead of returning a stale guess.

use super::jacobian::Jacobian;
use super::shape::{quad4_dshape, quad4_shape};

/// `(xi, eta)` of `point` inside the (affine) triangle `coords`, solved in
/// closed form from `point = coords[0] + xi * (coords[1] - coords[0]) +
/// eta * (coords[2] - coords[0])` by Cramer's rule. Errs if the triangle
/// is degenerate (zero or non-finite area) or if `point` is not finite.
pub fn tri3_natural_coords(
    coords: &[[f64; 2]; 3],
    point: [f64; 2],
) -> Result<[f64; 2], &'static str> {
    if !point[0].is_finite() || !point[1].is_finite() {
        return Err("point to invert is not finite");
    }
    let d1 = [coords[1][0] - coords[0][0], coords[1][1] - coords[0][1]];
    let d2 = [coords[2][0] - coords[0][0], coords[2][1] - coords[0][1]];
    let rhs = [point[0] - coords[0][0], point[1] - coords[0][1]];
    let det = d1[0] * d2[1] - d1[1] * d2[0];
    if !det.is_finite() || det == 0.0 {
        return Err("triangle is degenerate");
    }
    let xi = (d2[1] * rhs[0] - d2[0] * rhs[1]) / det;
    let eta = (d1[0] * rhs[1] - d1[1] * rhs[0]) / det;
    if !xi.is_finite() || !eta.is_finite() {
        return Err("triangle natural coordinates are not finite");
    }
    Ok([xi, eta])
}

/// `(xi, eta)` of `point` inside the (generally non-affine) quadrilateral
/// `coords`, by Newton iteration from the element center. Errs if any
/// iterate is non-finite (a singular or near-singular Jacobian), or if
/// `max_iterations` is exhausted without the physical residual `|x(xi,
/// eta) - point|` dropping below `physical_tol`.
pub fn quad4_natural_coords(
    coords: &[[f64; 2]; 4],
    point: [f64; 2],
    physical_tol: f64,
    max_iterations: usize,
) -> Result<[f64; 2], &'static str> {
    if !point[0].is_finite() || !point[1].is_finite() {
        return Err("point to invert is not finite");
    }
    let mut natural = [0.0, 0.0];
    for _ in 0..max_iterations {
        let [xi, eta] = natural;
        let n = quad4_shape(xi, eta);
        let physical: [f64; 2] =
            std::array::from_fn(|c| n.iter().zip(coords).map(|(ni, x)| ni * x[c]).sum());
        let residual = [point[0] - physical[0], point[1] - physical[1]];
        let residual_norm = residual[0].hypot(residual[1]);
        if !residual_norm.is_finite() {
            return Err("quad4 Newton inverse produced a non-finite iterate");
        }
        if residual_norm < physical_tol {
            return Ok(natural);
        }
        let dn = quad4_dshape(xi, eta);
        let jacobian = Jacobian::new(coords, &dn);
        if !jacobian.det.is_finite() || jacobian.det == 0.0 {
            return Err("quad4 Newton inverse hit a singular Jacobian");
        }
        // d(physical) = J^T * d(natural), so d(natural) = (J^{-1})^T * d(physical).
        let inv = &jacobian.inverse;
        let delta = [
            inv[0][0] * residual[0] + inv[1][0] * residual[1],
            inv[0][1] * residual[0] + inv[1][1] * residual[1],
        ];
        if !delta[0].is_finite() || !delta[1].is_finite() {
            return Err("quad4 Newton inverse produced a non-finite step");
        }
        natural = [xi + delta[0], eta + delta[1]];
    }
    Err("quad4 Newton inverse did not converge")
}

/// Whether `[xi, eta]` lies inside the reference triangle (barycentric
/// coordinates `1 - xi - eta`, `xi`, `eta` all `>= -tol`).
pub fn tri3_contains(natural: [f64; 2], tol: f64) -> bool {
    let [xi, eta] = natural;
    xi >= -tol && eta >= -tol && (1.0 - xi - eta) >= -tol
}

/// Whether `[xi, eta]` lies inside the reference square `[-1, 1]^2`
/// (expanded by `tol`).
pub fn quad4_contains(natural: [f64; 2], tol: f64) -> bool {
    let [xi, eta] = natural;
    (-1.0 - tol..=1.0 + tol).contains(&xi) && (-1.0 - tol..=1.0 + tol).contains(&eta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::continuum::shape::tri3_shape;

    #[test]
    fn tri3_inverse_round_trips_through_the_shape_functions() {
        let coords = [[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]];
        for natural in [[0.2, 0.3], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.6, 0.1]] {
            let n = tri3_shape(natural[0], natural[1]);
            let point: [f64; 2] =
                std::array::from_fn(|c| n.iter().zip(&coords).map(|(ni, x)| ni * x[c]).sum());
            let recovered = tri3_natural_coords(&coords, point).unwrap();
            assert!((recovered[0] - natural[0]).abs() < 1e-12);
            assert!((recovered[1] - natural[1]).abs() < 1e-12);
        }
    }

    #[test]
    fn tri3_degenerate_triangle_is_an_error() {
        let coords = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        assert!(tri3_natural_coords(&coords, [0.5, 0.0]).is_err());
    }

    #[test]
    fn quad4_newton_inverse_round_trips_on_a_distorted_element() {
        let coords = [[0.0, 0.0], [2.0, 0.2], [2.4, 1.9], [-0.2, 1.5]];
        for natural in [[0.0, 0.0], [0.5, -0.3], [-0.8, 0.7], [0.9, 0.9]] {
            let n = quad4_shape(natural[0], natural[1]);
            let point: [f64; 2] =
                std::array::from_fn(|c| n.iter().zip(&coords).map(|(ni, x)| ni * x[c]).sum());
            let recovered = quad4_natural_coords(&coords, point, 1e-9 * 2.4, 50).unwrap();
            assert!((recovered[0] - natural[0]).abs() < 1e-6, "{recovered:?}");
            assert!((recovered[1] - natural[1]).abs() < 1e-6, "{recovered:?}");
        }
    }

    #[test]
    fn quad4_inverse_of_a_point_outside_the_element_does_not_panic() {
        let coords = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let result = quad4_natural_coords(&coords, [50.0, 50.0], 1e-9, 50);
        // Either converges to a (far outside [-1,1]) natural point, or
        // reports non-convergence - either way, no panic and no silent
        // garbage.
        if let Ok(natural) = result {
            assert!(natural[0].is_finite() && natural[1].is_finite());
        }
    }

    #[test]
    fn contains_predicates_accept_the_boundary_and_reject_clearly_outside() {
        assert!(tri3_contains([0.0, 0.0], 1e-9));
        assert!(tri3_contains([1.0, 0.0], 1e-9));
        assert!(!tri3_contains([0.6, 0.6], 1e-9));
        assert!(quad4_contains([1.0, 1.0], 1e-9));
        assert!(!quad4_contains([1.5, 0.0], 1e-9));
    }
}
