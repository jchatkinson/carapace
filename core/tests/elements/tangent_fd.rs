//! The generic element-level tangent-vs-central-difference checker
//! (`carapace_core::testkit::assert_element_tangent_matches_fd`,
//! `docs/csfm-plan.md` phase 0), driven through `ElementOps` on `Element`
//! exactly as `Domain` assembles it — the exit criterion "the checker
//! passes on existing elements". `EmbeddedBar` (phase 2) and any nonlinear
//! `PlaneMaterial` (phase 1) will reuse this same checker.

use carapace_core::model::{Element, Material, Node, PlaneMaterial, Quad4};
use carapace_core::testkit::{assert_element_tangent_matches_fd, FdConfig};
use slotmap::SlotMap;

/// `Truss` with a nonlinear `Steel01` material, at a displacement state
/// whose axial strain sits on the elastic branch (well inside `fy / e0`).
/// Reference: central difference of the resisting force, `rel` 1e-7.
#[test]
fn truss_with_steel01_matches_central_difference_on_the_elastic_branch() {
    let (length, area) = (100.0, 2.0);
    let (fy, e0, b) = (250.0, 200_000.0, 0.02);

    let mut nodes = SlotMap::with_key();
    let node_i = nodes.insert(Node::new([0.0, 0.0]));
    let node_j = nodes.insert(Node::new([length, 0.0]));
    // Axial strain = 0.0005, well under fy/e0 = 0.00125.
    nodes[node_j].displacement = [0.05, 0.02, 0.0];

    let mut element = Element::Truss(carapace_core::model::Truss::new(
        node_i,
        node_j,
        area,
        Material::steel01(fy, e0, b, 0.0, 1.0, 0.0, 1.0),
    ));

    assert_element_tangent_matches_fd(&mut element, &mut nodes, None, &FdConfig::default());
}

/// `ZeroLength` with a nonlinear `ElasticPP` material on `dof 0`, at a
/// trial deformation on the elastic branch (strictly under the yield
/// deformation `eyp`). Reference: central difference of the resisting
/// force, `rel` 1e-7.
#[test]
fn zero_length_with_elastic_pp_matches_central_difference_on_the_elastic_branch() {
    let (e, eyp) = (100.0, 0.05);

    let mut nodes = SlotMap::with_key();
    let node_i = nodes.insert(Node::new([0.0, 0.0]));
    let node_j = nodes.insert(Node::new([0.0, 0.0]));
    nodes[node_j].displacement = [0.02, 0.0, 0.0];

    let mut element = Element::ZeroLength(
        carapace_core::model::ZeroLength::new(node_i, node_j)
            .with_material(0, Material::elastic_pp(e, eyp)),
    );

    assert_element_tangent_matches_fd(&mut element, &mut nodes, None, &FdConfig::default());
}

/// `Quad4` (elastic `PlaneStress`) on a skewed quadrilateral under a
/// generic displacement state. The material is linear, so this checks the
/// element's own `B`/`assemble_tangent` wiring (shape-function derivatives,
/// Gauss integration, the sink) against an independent finite difference
/// of the resisting force — not the material law. Reference: central
/// difference, `rel` 1e-7.
#[test]
fn quad4_matches_central_difference_on_a_skewed_element() {
    let coords = [[0.0, 0.0], [1.2, -0.1], [1.4, 1.3], [-0.1, 1.1]];
    let material = PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap();

    let mut nodes = SlotMap::with_key();
    let ids = coords.map(|c| nodes.insert(Node::new(c)));
    for (k, id) in ids.into_iter().enumerate() {
        let [x, y] = nodes[id].coords;
        nodes[id].displacement = [1e-3 * x * y + 1e-4 * k as f64, -5e-4 * x * x, 0.0];
    }

    let mut element = Element::Quad4(Quad4::new(ids, 0.4, material));

    assert_element_tangent_matches_fd(&mut element, &mut nodes, None, &FdConfig::default());
}
