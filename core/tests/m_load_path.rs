use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, LoadSeries, Material, Node, Truss};

/// A `LoadSeries::Path` pattern driven by `Integrator::LoadControl` walking
/// pseudo-time along the path's own control points — the load-controlled
/// equivalent of Xara's dedicated `Integrator::LoadPath` (a `StaticIntegrator`
/// that walks a prescribed load-factor vector one entry per step, see
/// `SRC/analysis/integrator/Static/LoadPath.h`): no separate `Integrator`
/// variant needed here, since a `Path` series already *is* a prescribed
/// load-factor-vs-pseudo-time path, and `LoadControl` already knows how to
/// step pseudo-time forward by a fixed amount each step.
///
/// A single elastic `Truss` makes this an exact, closed-form check at every
/// step (`u = F/k`, `k = E*A/L`), not just a "did it run" smoke test — the
/// load factor at pseudo-time `t` is exactly `path.factor(t)` (piecewise-
/// linear interpolation over the control points), so the expected
/// displacement at each step is computable directly from the same control
/// points the test hands to the pattern.
#[test]
fn load_control_walks_a_path_series_pattern_through_its_prescribed_load_factors() {
    let (e, area, length) = (30000.0, 2.0, 100.0);
    let k = e * area / length;

    let times = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let factors = vec![0.0, 50.0, -20.0, -20.0, 30.0];

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.add_element(Element::Truss(Truss::new(node_i, node_j, area, Material::Elastic { e })));

    let pattern = domain.add_load_pattern(LoadSeries::Path {
        times: times.clone(),
        factors: factors.clone(),
    });
    domain.add_nodal_load(pattern, node_j, 0, 1.0);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);

    // Step through every control point exactly (pseudo-time advances by
    // 1.0 each step, matching `times`' spacing) and check the closed form
    // at each one — including the flat segment (2.0 -> 3.0, factor held at
    // -20.0) and the final reversal back up to 30.0.
    for (i, &expected_factor) in factors.iter().enumerate().skip(1) {
        let result = analysis.step().expect("path-controlled truss should solve every step");
        assert_eq!(result.load_factor, times[i]);

        let expected_u = expected_factor / k;
        let u = analysis.domain().node(node_j).displacement[0];
        assert!(
            (u - expected_u).abs() < 1e-9,
            "step {i} (t={}): expected u={expected_u} (F={expected_factor}), got u={u}",
            times[i]
        );
    }
}
