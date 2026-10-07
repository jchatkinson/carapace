//! Gauss rules as `(xi, eta, weight)` triples (`(t, weight)` for the edge
//! rule). Triangle weights sum to the reference area `1/2`, quadrilateral
//! weights to `4`, edge weights to `2`.

const G: f64 = 0.577_350_269_189_625_8; // 1 / sqrt(3)

/// One-point triangle rule (centroid), exact for degree 1.
pub const TRI3_GAUSS_1: [(f64, f64, f64); 1] = [(1.0 / 3.0, 1.0 / 3.0, 0.5)];

/// Three-point triangle rule (interior points), exact for degree 2.
pub const TRI3_GAUSS_3: [(f64, f64, f64); 3] = [
    (1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0),
    (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
    (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
];

/// 2x2 rule on `[-1, 1]^2`, exact for degree 3 per axis. Point order is
/// `(-,-), (+,-), (+,+), (-,+)`, matching the quad's node order.
pub const QUAD4_GAUSS_2X2: [(f64, f64, f64); 4] =
    [(-G, -G, 1.0), (G, -G, 1.0), (G, G, 1.0), (-G, G, 1.0)];

/// Two-point rule on `[-1, 1]` for element edges, exact for degree 3.
pub const EDGE_GAUSS_2: [(f64, f64); 2] = [(-G, 1.0), (G, 1.0)];

#[cfg(test)]
mod tests {
    use super::*;

    fn integrate(rule: &[(f64, f64, f64)], f: impl Fn(f64, f64) -> f64) -> f64 {
        rule.iter().map(|&(x, y, w)| w * f(x, y)).sum()
    }

    #[test]
    fn weights_give_the_reference_areas() {
        assert!((integrate(&TRI3_GAUSS_1, |_, _| 1.0) - 0.5).abs() < 1e-15);
        assert!((integrate(&TRI3_GAUSS_3, |_, _| 1.0) - 0.5).abs() < 1e-15);
        assert!((integrate(&QUAD4_GAUSS_2X2, |_, _| 1.0) - 4.0).abs() < 1e-15);
        assert!((EDGE_GAUSS_2.iter().map(|p| p.1).sum::<f64>() - 2.0).abs() < 1e-15);
    }

    #[test]
    fn triangle_rules_are_exact_to_their_degree() {
        // Exact: int x^a y^b over the unit triangle = a! b! / (a + b + 2)!.
        assert!((integrate(&TRI3_GAUSS_1, |x, y| x + 2.0 * y) - 3.0 / 6.0).abs() < 1e-15);
        for (a, b, exact) in [(2, 0, 1.0 / 12.0), (1, 1, 1.0 / 24.0), (0, 2, 1.0 / 12.0)] {
            let got = integrate(&TRI3_GAUSS_3, |x, y| x.powi(a) * y.powi(b));
            assert!((got - exact).abs() < 1e-15, "x^{a} y^{b}: {got} vs {exact}");
        }
    }

    #[test]
    fn quad_and_edge_rules_integrate_cubics_exactly() {
        // Exact over [-1,1]^2: x^2 y^2 -> 4/9; x^3 y -> 0; x^2 + y^3 + 1 -> 4/3 + 0 + 4.
        assert!((integrate(&QUAD4_GAUSS_2X2, |x, y| x * x * y * y) - 4.0 / 9.0).abs() < 1e-15);
        assert!(integrate(&QUAD4_GAUSS_2X2, |x, y| x.powi(3) * y).abs() < 1e-15);
        assert!(
            (integrate(&QUAD4_GAUSS_2X2, |x, y| x * x + y.powi(3) + 1.0) - (4.0 / 3.0 + 4.0)).abs()
                < 1e-14
        );
        let edge = |f: fn(f64) -> f64| EDGE_GAUSS_2.iter().map(|&(t, w)| w * f(t)).sum::<f64>();
        assert!((edge(|t| t * t) - 2.0 / 3.0).abs() < 1e-15);
        assert!(edge(|t| t.powi(3)).abs() < 1e-15);
    }
}
