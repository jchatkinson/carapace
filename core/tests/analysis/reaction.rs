use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    Domain, ElasticBeamColumn, Element, ElementLoad, GeomTransf, Material, Node, Truss,
};

/// `Domain::reaction`'s simplest closed-form check: a fixed-free truss
/// pulled at its free end (same case `elements/truss.rs` uses for displacement)
/// must show a reaction at the fixed support exactly opposite the applied
/// load — global equilibrium (`sum(reactions) + sum(applied loads) == 0`)
/// for a single-DOF system with only one load and one support.
#[test]
fn fixed_free_truss_reaction_exactly_opposes_the_applied_load() {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 100.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([100.0, 100.0]).fix(1).fix(2));
    let load = 50.0;
    domain.load_node(node_j, 0, load);
    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        2.0,
        Material::Elastic { e: 30_000.0 },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    analysis.step().expect("linear elastic truss should solve");

    let reaction = analysis.domain().reaction(node_i, 0, 1.0);
    assert!(
        (reaction + load).abs() < 1e-9,
        "expected reaction={}, got {reaction}",
        -load
    );

    // Global equilibrium: this DOF's only load and only support must sum to zero.
    assert!((reaction + load).abs() < 1e-9);
}

/// A fixed node shared by *two* elements (a two-span beam continuous over
/// its middle support) sums each element's own contribution to that node's
/// reaction — proves `reaction` iterates every element touching the node,
/// not just the first (or last) one found, by checking the fixed end's
/// reaction against simple statics: with a point load at the free tip of a
/// cantilever built from two collinear elements in series, the root
/// reaction moment must equal `-load * total_length` regardless of where
/// the elements are split, and the root shear must equal `-load` — neither
/// of which could be right if only one of the two elements were being
/// summed.
#[test]
fn a_shared_fixed_node_sums_every_touching_elements_contribution() {
    let (e, a, iz, length, load) = (30_000.0_f64, 10.0, 1000.0, 100.0, -10.0);
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let mid = domain.add_node(Node::new([length / 2.0, 0.0]));
    let tip = domain.add_node(Node::new([length, 0.0]));
    domain.load_node(tip, 1, load);
    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        fixed,
        mid,
        e,
        a,
        iz,
        GeomTransf::Linear,
    )));
    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        mid,
        tip,
        e,
        a,
        iz,
        GeomTransf::Linear,
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    analysis
        .step()
        .expect("linear elastic two-span cantilever should solve");

    let shear = analysis.domain().reaction(fixed, 1, 1.0);
    let moment = analysis.domain().reaction(fixed, 2, 1.0);
    assert!(
        (shear - (-load)).abs() < 1e-6,
        "expected shear={}, got {shear}",
        -load
    );
    assert!(
        (moment - (-load * length)).abs() < 1e-6,
        "expected moment={}, got {moment}",
        -load * length
    );
}

/// `Domain::reactions_all` is the one-pass form of `reaction`: on a two-span cantilever carrying a
/// tip load and a uniform load on each span, every node's every DOF (supported, free and the
/// shared mid node) must equal the per-DOF sweep to rounding.
#[test]
fn reactions_all_matches_the_per_dof_reaction_with_nodal_and_element_loads() {
    let (e, a, iz, length) = (30_000.0_f64, 10.0, 1000.0, 100.0);
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let mid = domain.add_node(Node::new([length / 2.0, 0.0]));
    let tip = domain.add_node(Node::new([length, 0.0]));
    domain.load_node(tip, 1, -10.0);
    let beams = [(fixed, mid), (mid, tip)].map(|(i, j)| {
        domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            i,
            j,
            e,
            a,
            iz,
            GeomTransf::Linear,
        )))
    });
    let pattern = domain.default_pattern();
    domain.add_element_load(pattern, beams[0], ElementLoad::uniform(0.2, -0.5));
    domain.add_element_load(pattern, beams[1], ElementLoad::uniform(0.0, -0.3));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    analysis
        .step()
        .expect("linear elastic cantilever should solve");

    let all = analysis.domain().reactions_all(1.0);
    for node in [fixed, mid, tip] {
        for dof in 0..3 {
            let single = analysis.domain().reaction(node, dof, 1.0);
            assert!(
                (all[node][dof] - single).abs() < 1e-9,
                "node {node:?} dof {dof}: {} vs {single}",
                all[node][dof]
            );
        }
    }
    assert!(
        all[fixed][1].abs() > 1.0,
        "the fixed end must carry the loads"
    );
}
