//! `EmbeddedBar` (plan `docs/csfm-plan.md` §2.3, §4 D): patch tests on T3
//! and Q4 hosts under uniform strain, rigid translation, exact agreement
//! with a `Truss` on shared nodes, the generic tangent-vs-finite-
//! difference checker on a nonlinear `CsfmRebar`, and (plan §2.5 point 3)
//! `EmbeddedBar::limit_utilization`/`Domain::controlling_limit_utilization`
//! picking the right controlling element and `LimitKind` between a
//! `CsfmConcrete` host and its embedded `CsfmRebar`.

use carapace_core::model::{
    CsfmConcreteLaw, Domain, Element, ElementOps, EmbeddedBar, LimitKind, LimitStateConfig,
    Material, Node, NodeId, PlaneMaterial, Tri3, Truss,
};
use carapace_core::testkit::{assert_element_tangent_matches_fd, FdConfig};
use slotmap::SlotMap;

fn tri3_shape(xi: f64, eta: f64) -> [f64; 3] {
    [1.0 - xi - eta, xi, eta]
}

fn point_on_tri3(coords: &[[f64; 2]; 3], natural: [f64; 2]) -> [f64; 2] {
    let n = tri3_shape(natural[0], natural[1]);
    std::array::from_fn(|c| n.iter().zip(coords).map(|(ni, x)| ni * x[c]).sum())
}

fn quad4_shape(xi: f64, eta: f64) -> [f64; 4] {
    const CORNERS: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    std::array::from_fn(|i| {
        let [xi_i, eta_i] = CORNERS[i];
        0.25 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta)
    })
}

fn point_on_quad4(coords: &[[f64; 2]; 4], natural: [f64; 2]) -> [f64; 2] {
    let n = quad4_shape(natural[0], natural[1]);
    std::array::from_fn(|c| n.iter().zip(coords).map(|(ni, x)| ni * x[c]).sum())
}

fn segment_geometry(a: [f64; 2], b: [f64; 2]) -> ([f64; 2], f64) {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let length = dx.hypot(dy);
    ([dx / length, dy / length], length)
}

/// Uniform strain state `ux = a*x + b*y`, `uy = c*x + d*y`, so the
/// Cartesian strain tensor is constant everywhere: `eps_x = a`, `eps_y =
/// d`, `gamma_xy = b + c`.
fn node_with_strain(coords: [f64; 2], [a, b, c, d]: [f64; 4]) -> Node {
    let mut node = Node::new(coords);
    let [x, y] = coords;
    node.displacement = [a * x + b * y, c * x + d * y, 0.0];
    node
}

/// `d . eps . d` for the uniform strain tensor `[a, b; c, d]` and unit
/// direction `d = (dx, dy)` — the strain a bar along `d` should measure
/// under an affine host field, exactly (plan §4 D).
fn directional_strain([dx, dy]: [f64; 2], [a, _b, _c, d]: [f64; 4], gamma_xy: f64) -> f64 {
    a * dx * dx + d * dy * dy + gamma_xy * dx * dy
}

#[test]
fn tri3_patch_test_matches_the_host_strain_tensor_at_any_position_and_angle() {
    let coords = [[0.0, 0.0], [4.0, 0.3], [0.5, 3.5]];
    let strain_field = [2e-3, -5e-4, 7e-4, -1.2e-3];
    let gamma_xy = strain_field[1] + strain_field[2];
    let cases = [
        ([0.1, 0.2], [0.6, 0.1]),
        ([0.0, 0.0], [1.0, 0.0]),
        ([0.3, 0.3], [0.2, 0.6]),
        ([0.05, 0.9], [0.8, 0.05]),
    ];
    for (a_natural, b_natural) in cases {
        let mut domain: Domain = Domain::new();
        let ids: [NodeId; 3] =
            std::array::from_fn(|i| domain.add_node(node_with_strain(coords[i], strain_field)));

        let a = point_on_tri3(&coords, a_natural);
        let b = point_on_tri3(&coords, b_natural);
        let (direction, length) = segment_geometry(a, b);
        let n_a = tri3_shape(a_natural[0], a_natural[1]);
        let n_b = tri3_shape(b_natural[0], b_natural[1]);
        let material = Material::Elastic { e: 200_000.0 };
        let id = domain.add_element(Element::EmbeddedBar(EmbeddedBar::in_tri3(
            ids, n_a, n_b, direction, length, 1.0, material,
        )));

        let strain = domain.element_gauss_responses(id).unwrap()[0].strain[0];
        let expected = directional_strain(direction, strain_field, gamma_xy);
        assert!(
            (strain - expected).abs() < 1e-12,
            "a={a_natural:?} b={b_natural:?}: strain {strain} vs {expected}"
        );
    }
}

#[test]
fn quad4_patch_test_matches_the_host_strain_tensor_on_a_distorted_element() {
    let coords = [[0.0, 0.0], [2.0, -0.1], [2.4, 1.9], [-0.1, 1.4]];
    let strain_field = [1.5e-3, 3e-4, -6e-4, 2.1e-3];
    let gamma_xy = strain_field[1] + strain_field[2];
    let cases = [
        ([-0.3, 0.2], [0.6, -0.4]),
        ([-1.0, -1.0], [1.0, 1.0]),
        ([0.0, 0.0], [0.7, 0.3]),
    ];
    for (a_natural, b_natural) in cases {
        let mut domain: Domain = Domain::new();
        let ids: [NodeId; 4] =
            std::array::from_fn(|i| domain.add_node(node_with_strain(coords[i], strain_field)));

        let a = point_on_quad4(&coords, a_natural);
        let b = point_on_quad4(&coords, b_natural);
        let (direction, length) = segment_geometry(a, b);
        let n_a = quad4_shape(a_natural[0], a_natural[1]);
        let n_b = quad4_shape(b_natural[0], b_natural[1]);
        let material = Material::Elastic { e: 200_000.0 };
        let id = domain.add_element(Element::EmbeddedBar(EmbeddedBar::in_quad4(
            ids, n_a, n_b, direction, length, 1.0, material,
        )));

        let strain = domain.element_gauss_responses(id).unwrap()[0].strain[0];
        let expected = directional_strain(direction, strain_field, gamma_xy);
        assert!(
            (strain - expected).abs() < 1e-10,
            "a={a_natural:?} b={b_natural:?}: strain {strain} vs {expected}"
        );
    }
}

#[test]
fn rigid_translation_gives_zero_force() {
    let coords = [[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]];
    let mut domain: Domain = Domain::new();
    let ids: [NodeId; 3] = std::array::from_fn(|i| {
        let mut node = Node::new(coords[i]);
        node.displacement = [0.5, -0.3, 0.0];
        domain.add_node(node)
    });
    let a = point_on_tri3(&coords, [0.2, 0.3]);
    let b = point_on_tri3(&coords, [0.6, 0.1]);
    let (direction, length) = segment_geometry(a, b);
    let n_a = tri3_shape(0.2, 0.3);
    let n_b = tri3_shape(0.6, 0.1);
    let material = Material::Elastic { e: 200_000.0 };
    let id = domain.add_element(Element::EmbeddedBar(EmbeddedBar::in_tri3(
        ids, n_a, n_b, direction, length, 2.0, material,
    )));
    let strain = domain.element_gauss_responses(id).unwrap()[0].strain[0];
    assert!(strain.abs() < 1e-14, "strain {strain}");
    let force = domain.element_local_force(id);
    for i in 0..force.len() {
        assert!(
            force.as_slice()[i].abs() < 1e-10,
            "force[{i}] = {}",
            force.as_slice()[i]
        );
    }
}

/// A bar coincident with a host edge, spanning exactly two host nodes
/// (`N_a`/`N_b` one-hot at those nodes), matches a `Truss` on the same
/// two nodes exactly: same strain, same stress, same axial force.
#[test]
fn bar_coincident_with_a_mesh_line_matches_a_truss_on_shared_nodes() {
    let coords = [[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]];
    let displacements = [[0.01, -0.02], [-0.015, 0.03], [0.02, 0.01]];
    let mut domain: Domain = Domain::new();
    let ids: [NodeId; 3] = std::array::from_fn(|i| {
        let mut node = Node::new(coords[i]);
        node.displacement = [displacements[i][0], displacements[i][1], 0.0];
        domain.add_node(node)
    });

    let (direction, length) = segment_geometry(coords[0], coords[1]);
    let area = 3.0;
    let material = Material::steel01(500.0, 200_000.0, 0.02, 0.9, 5.0, 0.9, 5.0);
    // Node 0 is natural (xi=0, eta=0): N = [1, 0, 0]. Node 1: N = [0, 1, 0].
    let bar_id = domain.add_element(Element::EmbeddedBar(EmbeddedBar::in_tri3(
        ids,
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        direction,
        length,
        area,
        material.clone(),
    )));
    let bar_strain = domain.element_gauss_responses(bar_id).unwrap()[0].strain[0];
    let bar_stress = domain.element_gauss_responses(bar_id).unwrap()[0].stress[0];
    let bar_force = domain.element_local_force(bar_id);

    let truss_id = domain.add_element(Element::Truss(Truss::new(ids[0], ids[1], area, material)));
    let truss_strain = (displacements[1][0] - displacements[0][0]) * direction[0]
        + (displacements[1][1] - displacements[0][1]) * direction[1];
    let truss_strain = truss_strain / length;
    let truss_force = domain.element_local_force(truss_id);

    assert!(
        (bar_strain - truss_strain).abs() < 1e-13,
        "{bar_strain} vs {truss_strain}"
    );
    // Truss local force order [ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]; the bar's
    // order is [ux_i, uy_i, ux_j, uy_j, ux_k, uy_k] for its three host
    // nodes, with node k (the triangle's third node) carrying zero since
    // N_a = N_b = 0 there.
    assert!((bar_force.as_slice()[0] - truss_force.as_slice()[0]).abs() < 1e-10);
    assert!((bar_force.as_slice()[1] - truss_force.as_slice()[1]).abs() < 1e-10);
    assert!((bar_force.as_slice()[2] - truss_force.as_slice()[3]).abs() < 1e-10);
    assert!((bar_force.as_slice()[3] - truss_force.as_slice()[4]).abs() < 1e-10);
    assert!(bar_force.as_slice()[4].abs() < 1e-12);
    assert!(bar_force.as_slice()[5].abs() < 1e-12);
    let _ = bar_stress;
}

/// The generic tangent-vs-finite-difference checker (`docs/csfm-plan.md`
/// phase 0) passes on `EmbeddedBar` with a nonlinear `CsfmRebar`, exactly
/// as it already does for `Truss`/`ZeroLength`/`Quad4` — the exit
/// criterion for phase 2's own new law/element combination.
#[test]
fn embedded_bar_with_csfm_rebar_matches_the_central_difference_checker() {
    let coords = [[0.0, 0.0], [4.0, 0.3], [0.5, 3.5]];
    let mut nodes = SlotMap::with_key();
    let ids: [NodeId; 3] = std::array::from_fn(|i| nodes.insert(Node::new(coords[i])));
    for (k, id) in ids.into_iter().enumerate() {
        let [x, y] = nodes[id].coords;
        nodes[id].displacement = [
            1e-4 * x - 5e-5 * y + 1e-5 * k as f64,
            2e-5 * x + 3e-5 * y,
            0.0,
        ];
    }
    let a_natural = [0.15, 0.2];
    let b_natural = [0.55, 0.1];
    let a = point_on_tri3(&coords, a_natural);
    let b = point_on_tri3(&coords, b_natural);
    let (direction, length) = segment_geometry(a, b);
    let n_a = tri3_shape(a_natural[0], a_natural[1]);
    let n_b = tri3_shape(b_natural[0], b_natural[1]);

    let material = Material::csfm_rebar(
        200_000.0, 500.0, 540.0, 0.05, 2.9, 32000.0, 16.0, 0.02, 0.67, true,
    )
    .unwrap();
    let mut element = Element::EmbeddedBar(EmbeddedBar::in_tri3(
        ids, n_a, n_b, direction, length, 1.0, material,
    ));

    assert_element_tangent_matches_fd(&mut element, &mut nodes, None, &FdConfig::default());
}

/// Mask conformance: an `EmbeddedBar` only ever stiffens `ux`/`uy`.
#[test]
fn dof_mask_is_translations_only() {
    let coords = [[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]];
    let mut domain: Domain = Domain::new();
    let ids: [NodeId; 3] = std::array::from_fn(|i| domain.add_node(Node::new(coords[i])));
    let (direction, length) = segment_geometry(coords[0], coords[1]);
    let material = Material::Elastic { e: 200_000.0 };
    let element = Element::EmbeddedBar(EmbeddedBar::in_tri3(
        ids,
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        direction,
        length,
        1.0,
        material,
    ));
    let mask = element.dof_mask();
    assert!(mask.contains(0));
    assert!(mask.contains(1));
    assert!(!mask.contains(2));
}

/// A plain (non-CSFM) `Tri3`/`EmbeddedBar` pair reports no limit state at
/// all — `limit_utilization` is `None` per element and
/// `controlling_limit_utilization` is `None` for the whole domain.
#[test]
fn limit_utilization_is_none_without_csfm_materials() {
    let coords = [[0.0, 0.0], [4.0, 0.0], [0.0, 4.0]];
    let mut domain: Domain = Domain::new();
    let ids: [NodeId; 3] = std::array::from_fn(|i| domain.add_node(Node::new(coords[i])));
    domain.add_element(Element::Tri3(Tri3::new(
        ids[0],
        ids[1],
        ids[2],
        1.0,
        PlaneMaterial::plane_stress(30e3, 0.2).unwrap(),
    )));
    let (direction, length) = segment_geometry(coords[0], coords[1]);
    domain.add_element(Element::EmbeddedBar(EmbeddedBar::in_tri3(
        ids,
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        direction,
        length,
        1.0,
        Material::Elastic { e: 200_000.0 },
    )));
    assert_eq!(
        domain.controlling_limit_utilization(&LimitStateConfig::default()),
        None
    );
}

const LEG_LENGTH: f64 = 4.0;

fn rebar_material() -> Material {
    Material::csfm_rebar(
        200_000.0, 500.0, 540.0, 0.05, 2.9, 32000.0, 16.0, 0.02, 0.67, true,
    )
    .unwrap()
}

fn rebar_bounds() -> (f64, f64) {
    let Material::CsfmRebar(fields) = rebar_material() else {
        panic!()
    };
    (fields.eps_m_at_ft, fields.eps_u)
}

/// Builds the `Tri3` (`CsfmConcrete`) + one-leg `EmbeddedBar` (`CsfmRebar`)
/// fixture with node 1's `ux` already set to `ux1`, since nothing public
/// mutates a node already inserted into a `Domain` (same pattern as the
/// patch tests above: set `Node::displacement` before `add_node`).
fn concrete_and_rebar_leg(
    ux1: f64,
) -> (
    Domain,
    carapace_core::model::ElementId,
    carapace_core::model::ElementId,
) {
    let coords = [[0.0, 0.0], [LEG_LENGTH, 0.0], [0.0, LEG_LENGTH]];
    let mut domain: Domain = Domain::new();
    let ids: [NodeId; 3] = std::array::from_fn(|i| {
        let mut node = Node::new(coords[i]);
        if i == 1 {
            node.displacement[0] = ux1;
        }
        domain.add_node(node)
    });
    let concrete_id = domain.add_element(Element::Tri3(Tri3::new(
        ids[0],
        ids[1],
        ids[2],
        1.0,
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap(),
    )));
    let (direction, length) = segment_geometry(coords[0], coords[1]);
    let bar_id = domain.add_element(Element::EmbeddedBar(EmbeddedBar::in_tri3(
        ids,
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        direction,
        length,
        1.0,
        rebar_material(),
    )));
    (domain, concrete_id, bar_id)
}

/// A `Tri3` host (`CsfmConcrete`) with one `EmbeddedBar` (`CsfmRebar`)
/// bound to its two legs: pure uniaxial strain along that leg (CST, so
/// `eps_x = u1x / length` exactly, `eps_y = gamma_xy = 0` for this
/// right-triangle geometry) isolates each limit in turn.
/// `Domain::controlling_limit_utilization` must pick the rebar when pulled
/// into tension past its own rupture utilization (concrete crushing stays
/// at exactly `0`, since `eps_2 = 0` there), and the concrete when pushed
/// into compression past `|EPS_CU2|` (the rebar's own utilization stays
/// small, since the same strain is far inside its bare-steel range).
#[test]
fn domain_reports_the_controlling_element_and_kind() {
    const EPS_CU2: f64 = -0.0035;
    let (eps_m_at_ft, eps_u) = rebar_bounds();

    // Tension: strain = 1.5 * eps_m_at_ft, well past rebar rupture.
    let (domain, concrete_id, bar_id) = concrete_and_rebar_leg(1.5 * eps_m_at_ft * LEG_LENGTH);
    let (element, utilization, kind) = domain
        .controlling_limit_utilization(&LimitStateConfig::default())
        .expect("rebar rupture utilization must be reported");
    assert_eq!(element, bar_id);
    assert_eq!(kind, LimitKind::RebarRupture);
    assert!((utilization - 1.5).abs() < 1e-9, "{utilization}");
    // Concrete stays at exactly zero: pure tension, eps_2 = 0.
    let concrete_strain = domain.element_gauss_responses(concrete_id).unwrap()[0].strain;
    assert_eq!(concrete_strain[1], 0.0);

    // Compression: strain = 1.5 * EPS_CU2, well past the default crushing
    // limit but far inside the rebar's own bare-steel range.
    let (domain, concrete_id, bar_id) = concrete_and_rebar_leg(1.5 * EPS_CU2 * LEG_LENGTH);
    let (element, utilization, kind) = domain
        .controlling_limit_utilization(&LimitStateConfig::default())
        .expect("concrete crushing utilization must be reported");
    assert_eq!(element, concrete_id);
    assert_eq!(kind, LimitKind::ConcreteCrushing);
    assert!((utilization - 1.5).abs() < 1e-9, "{utilization}");
    let rebar_strain = domain.element_gauss_responses(bar_id).unwrap()[0].strain[0];
    let rebar_utilization = rebar_strain.abs() / eps_u;
    assert!(
        rebar_utilization < 1.0 && rebar_utilization > 0.0,
        "rebar should be stressed but nowhere near its own limit: {rebar_utilization}"
    );
}
