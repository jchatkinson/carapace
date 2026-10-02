use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{Domain, Element, Material, Node, ZeroLength};

// A yielding ground spring in series with an elastic spring exercises
// both the controlled DOF and an unconstrained equilibrium correction.
// The post-yield load sensitivity differs from the initial one.
#[test]
fn displacement_control_converges_across_yield_with_each_tangent_strategy() {
    for algorithm in [
        Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        },
        Algorithm::Newton {
            tangent: TangentStrategy::ReuseAtStepStart,
            line_search: None,
        },
        Algorithm::Newton {
            tangent: TangentStrategy::Initial,
            line_search: None,
        },
        Algorithm::KrylovNewton {
            tangent: TangentStrategy::Initial,
            max_dimension: 2,
        },
    ] {
        let mut domain = Domain::new();
        let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let middle = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
        let tip = domain.add_node(Node::new([2.0, 0.0]).fix(1).fix(2));
        for material in [
            Material::Elastic { e: 50.0 },
            Material::elastic_pp(100.0, 0.01),
        ] {
            domain.add_element(Element::ZeroLength(
                ZeroLength::new(ground, middle).with_material(0, material),
            ));
        }
        domain.add_element(Element::ZeroLength(
            ZeroLength::new(middle, tip).with_material(0, Material::Elastic { e: 100.0 }),
        ));
        domain.load_node(tip, 0, 1.0);
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::DisplacementControl {
                node: tip,
                dof: 0,
                increment: 0.03,
            })
            .algorithm(algorithm)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-10,
                max_iter: 100,
            })
            .build(domain);
        let result = analysis
            .step()
            .expect("yielding displacement-controlled step should converge");
        // lambda = 50*u_middle + 1 = 100*(u_tip-u_middle).
        assert!((analysis.domain().node(tip).displacement[0] - 0.03).abs() < 1e-10);
        assert!((analysis.domain().node(middle).displacement[0] - 2.0 / 150.0).abs() < 1e-10);
        assert!((result.load_factor - 5.0 / 3.0).abs() < 1e-9);
        assert!(result.iterations > 1);
    }
}
