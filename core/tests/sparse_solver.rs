use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, ElasticBeamColumn, Element, GeomTransf, Node};

/// Not tied to a specific milestone — a correctness + scale check for the
/// `SparseSolver` swap to `faer` (see implementation-plan §6 decision #1).
/// A single 2-node `ElasticBeamColumn` is exact for M3's simply-supported
/// UDL case (§7); this does the equivalent check for a *chain* of many
/// elements, which exercises the sparse assembly/solve path at a
/// non-trivial system size (150 free DOFs) with a real banded/sparse
/// pattern, rather than the single- or few-element systems every other
/// test uses. Chaining `n` identical prismatic Euler-Bernoulli elements
/// under only a tip point load (no interior element loads) is exact
/// against the classical cantilever closed form `delta = P*L^3/(3*E*I)`
/// regardless of `n` — cubic Hermite elements reproduce a cubic
/// deflection field exactly, so subdividing adds no discretization error.
#[test]
fn long_cantilever_chain_matches_closed_form_under_sparse_solve() {
    let (e, iz, area) = (30000.0, 1000.0, 100.0);
    let total_length = 1000.0;
    let n = 50;
    let element_length = total_length / n as f64;
    let tip_load = -10.0;

    let mut domain = Domain::new();
    let mut nodes = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let x = element_length * i as f64;
        let mut node = Node::new([x, 0.0]);
        if i == 0 {
            node = node.fix(0).fix(1).fix(2); // fixed support
        } else if i == n {
            node = node.with_load(1, tip_load);
        }
        nodes.push(domain.add_node(node));
    }
    for i in 0..n {
        domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            nodes[i],
            nodes[i + 1],
            e,
            area,
            iz,
            GeomTransf::Linear,
        )));
    }

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    analysis.step().expect("long cantilever chain should solve");

    let v_tip = analysis.domain().node(*nodes.last().unwrap()).displacement[1];
    let expected = tip_load * total_length.powi(3) / (3.0 * e * iz);
    assert!(
        (v_tip - expected).abs() < 1e-6,
        "expected tip deflection={expected}, got {v_tip}"
    );
}
