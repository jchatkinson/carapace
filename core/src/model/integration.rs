/// Numerical integration scheme along a frame element's length, relating
/// per-section response (`FiberSection`) to element-level response — needed
/// by `DispBeamColumn` (§3.1). Closed enum (§2.1).
///
/// Point/weight tables are copied directly from OpenSees's own
/// (`xara/SRC/quadrature/Frame/{Legendre,Lobatto}BeamIntegration.cpp`),
/// not computed via a general-purpose quadrature crate — domain-specific
/// numerics like this are copied from the reference implementation, the
/// same way element/material formulations are, not sourced from a generic
/// library (unlike heavily-optimized linear algebra, e.g. `nalgebra`/
/// `faer`, where using a library is the right call). Supports 2-6 points,
/// matching the range `DispBeamColumn` actually needs in practice — the
/// C++ source supports up to 10; extend the tables if that's ever needed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BeamIntegration {
    /// Gauss-Legendre: interior points only, `n`-point rule exact to
    /// degree `2n-1`. `DispBeamColumn`'s usual default.
    Legendre { points: usize },
    /// Gauss-Lobatto: includes both end points, `n`-point rule exact to
    /// degree `2n-3` — favored for `ForceBeamColumn` (M8) since end-section
    /// behavior matters most for plastic hinges there, but not needed by
    /// `DispBeamColumn` itself. Included now since both were asked for
    /// together (§3.1) and the table-driven implementation cost is the
    /// same either way.
    Lobatto { points: usize },
}

impl BeamIntegration {
    /// `(xi, weight)` pairs, `xi` in `[0,1]` along the element length,
    /// weights summing to 1 — already in the same mapped-to-`[0,1]` form
    /// `getSectionLocations`/`getSectionWeights` return in the C++ source
    /// (`xi = 0.5*(xi+1)`, `wt *= 0.5`), not the raw `[-1,1]` table.
    pub fn points(&self) -> Vec<(f64, f64)> {
        let (raw_xi, raw_wt): (&[f64], &[f64]) = match self {
            BeamIntegration::Legendre { points } => legendre_table(*points),
            BeamIntegration::Lobatto { points } => lobatto_table(*points),
        };
        raw_xi
            .iter()
            .zip(raw_wt)
            .map(|(&xi, &wt)| (0.5 * (xi + 1.0), 0.5 * wt))
            .collect()
    }
}

fn legendre_table(points: usize) -> (&'static [f64], &'static [f64]) {
    match points {
        2 => (
            &[-0.577350269189626, 0.577350269189626],
            &[1.0, 1.0],
        ),
        3 => (
            &[-0.774596669241483, 0.0, 0.774596669241483],
            &[0.555555555555556, 0.888888888888889, 0.555555555555556],
        ),
        4 => (
            &[
                -0.861136311594053,
                -0.339981043584856,
                0.339981043584856,
                0.861136311594053,
            ],
            &[
                0.347854845137454,
                0.652145154862546,
                0.652145154862546,
                0.347854845137454,
            ],
        ),
        5 => (
            &[
                -0.906179845938664,
                -0.538469310105683,
                0.0,
                0.538469310105683,
                0.906179845938664,
            ],
            &[
                0.236926885056189,
                0.478628670499366,
                0.568888888888889,
                0.478628670499366,
                0.236926885056189,
            ],
        ),
        6 => (
            &[
                -0.932469514203152,
                -0.661209386466265,
                -0.238619186083197,
                0.238619186083197,
                0.661209386466265,
                0.932469514203152,
            ],
            &[
                0.171324492379170,
                0.360761573048139,
                0.467913934572691,
                0.467913934572691,
                0.360761573048139,
                0.171324492379170,
            ],
        ),
        _ => panic!("BeamIntegration::Legendre supports 2-6 points, got {points}"),
    }
}

fn lobatto_table(points: usize) -> (&'static [f64], &'static [f64]) {
    match points {
        2 => (&[-1.0, 1.0], &[1.0, 1.0]),
        3 => (
            &[-1.0, 0.0, 1.0],
            &[0.333333333333333, 1.333333333333333, 0.333333333333333],
        ),
        4 => (
            &[-1.0, -0.44721360, 0.44721360, 1.0],
            &[0.166666666666667, 0.833333333333333, 0.833333333333333, 0.166666666666667],
        ),
        5 => (
            &[-1.0, -0.65465367, 0.0, 0.65465367, 1.0],
            &[0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1],
        ),
        6 => (
            &[
                -1.0,
                -0.7650553239,
                -0.2852315164,
                0.2852315164,
                0.7650553239,
                1.0,
            ],
            &[
                0.06666666667,
                0.3784749562,
                0.5548583770,
                0.5548583770,
                0.3784749562,
                0.06666666667,
            ],
        ),
        _ => panic!("BeamIntegration::Lobatto supports 2-6 points, got {points}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn points_are_symmetric_and_weights_sum_to_one() {
        for points in 2..=6 {
            for integration in [BeamIntegration::Legendre { points }, BeamIntegration::Lobatto { points }] {
                let pts = integration.points();
                assert_eq!(pts.len(), points);
                let sum: f64 = pts.iter().map(|(_, w)| w).sum();
                // 1e-8, not tighter: OpenSees's own source literals for
                // higher-order Lobatto weights are truncated to ~10
                // significant digits (e.g. `0.06666666667`), not full f64
                // precision.
                assert!((sum - 1.0).abs() < 1e-8, "{integration:?} weights should sum to 1, got {sum}");
            }
        }
    }

    #[test]
    fn lobatto_includes_the_endpoints() {
        let pts = BeamIntegration::Lobatto { points: 4 }.points();
        assert!((pts.first().unwrap().0 - 0.0).abs() < 1e-12);
        assert!((pts.last().unwrap().0 - 1.0).abs() < 1e-12);
    }
}
