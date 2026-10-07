//! Geometry validation with a scale-aware tolerance: every threshold is
//! relative to the element's characteristic length squared (its longest
//! node-to-node distance squared), never an absolute number, so a unit
//! element and a 1e6-scaled copy are judged identically.

/// A distance below `COINCIDENT * L` is "the same node".
const COINCIDENT: f64 = 1e-8;
/// A determinant (or twice the triangle area) below `DEGENERATE * L^2` is zero.
const DEGENERATE: f64 = 1e-10;

fn length_squared<const N: usize>(p: &[[f64; 2]; N]) -> f64 {
    let mut longest = 0.0_f64;
    for i in 0..N {
        for j in i + 1..N {
            let (dx, dy) = (p[j][0] - p[i][0], p[j][1] - p[i][1]);
            longest = longest.max(dx * dx + dy * dy);
        }
    }
    longest
}

fn check_distinct<const N: usize>(p: &[[f64; 2]; N], l2: f64) -> Result<(), &'static str> {
    for i in 0..N {
        for j in i + 1..N {
            let (dx, dy) = (p[j][0] - p[i][0], p[j][1] - p[i][1]);
            if dx * dx + dy * dy <= COINCIDENT * COINCIDENT * l2 {
                return Err("repeated or coincident nodes");
            }
        }
    }
    Ok(())
}

/// A triangle needs distinct nodes, area above tolerance, and counter-clockwise order.
pub fn validate_triangle(p: &[[f64; 2]; 3]) -> Result<(), &'static str> {
    let l2 = length_squared(p);
    if !l2.is_finite() || l2 == 0.0 {
        return Err("repeated or coincident nodes");
    }
    check_distinct(p, l2)?;
    let twice_area =
        (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
    if twice_area.abs() <= DEGENERATE * l2 {
        return Err("triangle area is zero (collinear nodes)");
    }
    if twice_area < 0.0 {
        return Err("triangle nodes are ordered clockwise");
    }
    Ok(())
}

/// `det J` at reference corner `k` of a bilinear quad, from the two edge
/// vectors leaving the corner (`J = [e1; e2] / 2` there).
fn corner_det(p: &[[f64; 2]; 4], k: usize) -> f64 {
    let (next, prev) = (p[(k + 1) % 4], p[(k + 3) % 4]);
    let (e1, e2) = (
        [next[0] - p[k][0], next[1] - p[k][1]],
        [prev[0] - p[k][0], prev[1] - p[k][1]],
    );
    // Counter-clockwise order: the edge to the next node cross the edge to the previous is positive.
    0.25 * (e1[0] * e2[1] - e1[1] * e2[0])
}

/// A quad needs distinct nodes and a positive `det J` at all four reference corners.
///
/// `det J` of a bilinear map is linear in `(xi, eta)`, so positivity at the
/// corners is necessary and sufficient for positivity everywhere; positivity
/// at the Gauss points alone is not (a quad with one reflex corner can have
/// all Gauss-point determinants positive).
pub fn validate_quad(p: &[[f64; 2]; 4]) -> Result<(), &'static str> {
    let l2 = length_squared(p);
    if !l2.is_finite() || l2 == 0.0 {
        return Err("repeated or coincident nodes");
    }
    check_distinct(p, l2)?;
    let tol = DEGENERATE * l2;
    let dets: [f64; 4] = std::array::from_fn(|k| corner_det(p, k));
    if dets.iter().all(|&d| d < -tol) {
        return Err("quad nodes are ordered clockwise");
    }
    const MESSAGES: [&str; 4] = [
        "quad corner 0 has a non-positive Jacobian determinant (reflex or degenerate corner)",
        "quad corner 1 has a non-positive Jacobian determinant (reflex or degenerate corner)",
        "quad corner 2 has a non-positive Jacobian determinant (reflex or degenerate corner)",
        "quad corner 3 has a non-positive Jacobian determinant (reflex or degenerate corner)",
    ];
    for (k, &d) in dets.iter().enumerate() {
        if d <= tol {
            return Err(MESSAGES[k]);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::continuum::{quad4_dshape, Jacobian, QUAD4_GAUSS_2X2};

    fn scaled<const N: usize>(p: [[f64; 2]; N], s: f64) -> [[f64; 2]; N] {
        p.map(|[x, y]| [x * s, y * s])
    }

    #[test]
    fn the_reflex_corner_quad_is_rejected_though_its_gauss_point_determinants_are_positive() {
        let reflex = [[0.0, 0.0], [1.0, 0.0], [0.4, 0.4], [0.0, 1.0]];
        let mut gauss: Vec<f64> = QUAD4_GAUSS_2X2
            .iter()
            .map(|&(xi, eta, _)| Jacobian::new(&reflex, &quad4_dshape(xi, eta)).det)
            .collect();
        gauss.sort_by(f64::total_cmp);
        assert!(gauss.iter().all(|&d| d > 0.0), "{gauss:?}");
        for (got, expected) in gauss.iter().zip([0.013, 0.1, 0.1, 0.187]) {
            assert!((got - expected).abs() < 1e-3, "{gauss:?}");
        }
        let corner = Jacobian::new(&reflex, &quad4_dshape(1.0, 1.0)).det;
        assert!((corner + 0.05).abs() < 1e-12);
        let error = validate_quad(&reflex).unwrap_err();
        assert!(error.contains("corner 2"), "{error}");
    }

    #[test]
    fn convex_trapezoid_and_parallelogram_are_accepted() {
        assert!(validate_quad(&[[0.0, 0.0], [2.0, 0.0], [1.5, 1.0], [0.5, 1.0]]).is_ok());
        assert!(validate_quad(&[[0.0, 0.0], [2.0, 0.0], [3.0, 1.5], [1.0, 1.5]]).is_ok());
        assert!(validate_triangle(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]).is_ok());
    }

    #[test]
    fn repeated_clockwise_and_degenerate_inputs_are_rejected() {
        let quad = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let clockwise = [quad[0], quad[3], quad[2], quad[1]];
        assert!(validate_quad(&clockwise).unwrap_err().contains("clockwise"));
        let repeated = [quad[0], quad[1], quad[1], quad[3]];
        assert!(validate_quad(&repeated).unwrap_err().contains("coincident"));
        assert!(validate_quad(&[[0.0, 0.0]; 4]).is_err());
        // A bow-tie (self-intersecting) quad.
        assert!(validate_quad(&[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0]]).is_err());

        let tri = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        assert!(validate_triangle(&[tri[0], tri[2], tri[1]])
            .unwrap_err()
            .contains("clockwise"));
        assert!(validate_triangle(&[tri[0], tri[1], tri[1]])
            .unwrap_err()
            .contains("coincident"));
        assert!(validate_triangle(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
            .unwrap_err()
            .contains("collinear"));
    }

    #[test]
    fn the_tolerance_is_scale_invariant() {
        let near_degenerate = [[0.0, 0.0], [1.0, 0.0], [0.5, 1e-11]];
        let sliver_quad = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        for s in [1.0, 1e-6, 1e6] {
            assert!(
                validate_triangle(&scaled(near_degenerate, s)).is_err(),
                "scale {s}"
            );
            assert!(validate_quad(&scaled(sliver_quad, s)).is_ok(), "scale {s}");
            assert!(
                validate_triangle(&scaled([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], s)).is_ok(),
                "scale {s}"
            );
        }
    }
}
