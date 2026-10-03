use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    Domain, Domain3, ElasticBeamColumn, ElasticBeamColumn3, Element, Element3, ElementLoad,
    ElementLoad3, GeomTransf, GeomTransf3, Node, Node3, SpatialDof,
};

const E: f64 = 30_000.0;
const A: f64 = 10.0;
const IZ: f64 = 1000.0;
const L: f64 = 100.0;

fn solve_planar_cantilever(
    loads: &[(f64, f64)],
) -> (
    Domain,
    carapace_core::model::ElementId,
    carapace_core::model::NodeId,
) {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([L, 0.0]));
    let beam = domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        node_i,
        node_j,
        E,
        A,
        IZ,
        GeomTransf::Linear,
    )));
    let pattern = domain.default_pattern();
    for &(wx, wy) in loads {
        domain.add_element_load(pattern, beam, ElementLoad::Uniform { wx, wy });
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
    analysis.step().expect("cantilever should solve");
    (analysis.domain().clone(), beam, node_j)
}

/// Axial + transverse UDL on a cantilever: tip `ux = wx L² / 2EA`, tip
/// `uy = wy L⁴ / 8EI`, and the member end forces equal the statics result —
/// the free end carries nothing, the fixed end carries the whole load.
#[test]
fn planar_cantilever_axial_and_transverse_udl_match_closed_form() {
    let (wx, wy) = (2.0, -1.0);
    let (domain, beam, node_j) = solve_planar_cantilever(&[(wx, wy)]);

    let ux_tip = domain.node(node_j).displacement[0];
    assert!((ux_tip - wx * L * L / (2.0 * E * A)).abs() < 1e-12);
    let uy_tip = domain.node(node_j).displacement[1];
    assert!((uy_tip - wy * L.powi(4) / (8.0 * E * IZ)).abs() < 1e-9);

    let f = domain.element_end_force(beam, 1.0);
    let expected = [-wx * L, -wy * L, -wy * L * L / 2.0, 0.0, 0.0, 0.0];
    for i in 0..6 {
        assert!(
            (f[i] - expected[i]).abs() < 1e-6,
            "component {i}: got {}, expected {}",
            f[i],
            expected[i]
        );
    }
}

/// Several loads on one element in one pattern sum rather than the last one
/// winning (OpenSees semantics).
#[test]
fn loads_on_the_same_element_accumulate() {
    let (d, _, n) = solve_planar_cantilever(&[(1.0, 0.0), (1.0, 0.0)]);
    let ux_split = d.node(n).displacement[0];
    let (d, _, n) = solve_planar_cantilever(&[(2.0, 0.0)]);
    let ux_whole = d.node(n).displacement[0];
    assert!(ux_whole.abs() > 0.0);
    assert!((ux_split - ux_whole).abs() < 1e-15);
}

/// Spatial counterpart: a member along global x with `vec_xz = [0,0,1]` has
/// local axes equal to global, so `wx` is plain axial.
#[test]
fn spatial_cantilever_axial_udl_matches_closed_form() {
    let (g, j, iy, iz) = (12_000.0, 500.0, 1500.0, 1000.0);
    let wx = 3.0;
    let mut domain = Domain3::new();
    let mut fixed = Node3::new([0.0, 0.0, 0.0]);
    for dof in [
        SpatialDof::Ux,
        SpatialDof::Uy,
        SpatialDof::Uz,
        SpatialDof::Rx,
        SpatialDof::Ry,
        SpatialDof::Rz,
    ] {
        fixed = fixed.fix(dof as usize);
    }
    let node_i = domain.add_node(fixed);
    let node_j = domain.add_node(Node3::new([L, 0.0, 0.0]));
    let beam = domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
        node_i,
        node_j,
        E,
        g,
        A,
        j,
        iy,
        iz,
        GeomTransf3::linear([0.0, 0.0, 1.0]),
    )));
    let pattern = domain.default_pattern();
    domain.add_element_load(
        pattern,
        beam,
        ElementLoad3::Uniform {
            wx,
            wy: 0.0,
            wz: 0.0,
        },
    );
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    analysis.step().expect("spatial cantilever should solve");

    let ux = analysis.domain().node(node_j).displacement[SpatialDof::Ux as usize];
    assert!((ux - wx * L * L / (2.0 * E * A)).abs() < 1e-12);
    let f = analysis.domain().element_end_force(beam, 1.0);
    assert!((f[0] + wx * L).abs() < 1e-6, "fixed-end axial {}", f[0]);
    assert!(f[6].abs() < 1e-6, "free-end axial {}", f[6]);
}
