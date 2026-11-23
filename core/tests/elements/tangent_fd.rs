//! The generic element-level tangent-vs-central-difference checker
//! (`carapace_core::testkit::assert_element_tangent_matches_fd`,
//! `docs/csfm-plan.md` phase 0), driven through `ElementOps` on `Element`
//! exactly as `Domain` assembles it — the exit criterion "the checker
//! passes on existing elements". `EmbeddedBar` (phase 2) and any nonlinear
//! `PlaneMaterial` (phase 1) will reuse this same checker. The nonlinear
//! stud-wall cases include a door, slot loads, O86 and hysteretic sheathing;
//! central differences use a 1e-6 relative tolerance with a small absolute
//! floor for nearly zero tangent entries.

use carapace_core::model::{
    ChordDistribution, Element, Material, Node, PlaneMaterial, Quad4, StudWall, StudWallLayout,
    WallChord, WallShearLaw, WallSheathing, WallStud,
};
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

/// O86 pier condensation uses the implicit local-law derivative in the exterior tangent.
/// Reference: central difference of the resisting force; relative tolerance 1e-6.
#[test]
fn o86_stud_wall_matches_central_difference() {
    let mut nodes = SlotMap::with_key();
    let ids = [
        nodes.insert(Node::new([0., 0.])),
        nodes.insert(Node::new([10., 0.])),
        nodes.insert(Node::new([10., 8.])),
        nodes.insert(Node::new([0., 8.])),
    ];
    for (i, id) in ids.iter().enumerate() {
        nodes[*id].displacement = [i as f64 * 0.001, (3 - i) as f64 * -0.0007, 0.];
    }
    let layout = StudWallLayout {
        length: 10.,
        height: 8.,
        studs: vec![
            WallStud { x: 0., ea: 30_000. },
            WallStud { x: 5., ea: 10_000. },
            WallStud {
                x: 10.,
                ea: 30_000.,
            },
        ],
        top_ea: 1e8,
        top_ei: Some(1e8),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: 1e8,
        bottom_ei: 1e8,
        openings: vec![],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::O86 {
                bv: 100.,
                nail_spacing: 100.,
                nail_diameter: 3.25,
            },
        }),
        foundation: false,
        floor: None,
    };
    let wall = StudWall::new(ids, vec![], layout, std::array::from_fn(|_| vec![])).unwrap();
    let mut element = Element::StudWall(wall);
    assert_element_tangent_matches_fd(
        &mut element,
        &mut nodes,
        None,
        &FdConfig {
            h_rel: 1e-4,
            h_min: 1e-7,
            tol_abs: 1e-9,
            tol_rel: 1e-6,
        },
    );
}

/// A two-pier door wall with a slot load has a consistent tangent for O86 and hysteretic sheathing.
/// Reference: central difference of exterior resisting force; relative tolerance 1e-6.
#[test]
fn nonlinear_door_wall_with_slot_load_matches_central_difference() {
    for law in [
        WallShearLaw::O86 {
            bv: 100.,
            nail_spacing: 100.,
            nail_diameter: 3.25,
        },
        WallShearLaw::Material(Material::hysteretic(
            10., 0.01, 15., 0.03, 16., 0.08, -10., -0.01, -15., -0.03, -16., -0.08, 0.2, 0.2, 0.,
            0., 0.,
        )),
    ] {
        let mut nodes = SlotMap::with_key();
        let ids = [
            nodes.insert(Node::new([0., 0.])),
            nodes.insert(Node::new([12., 0.])),
            nodes.insert(Node::new([12., 8.])),
            nodes.insert(Node::new([0., 8.])),
        ];
        for (i, id) in ids.iter().enumerate() {
            nodes[*id].displacement = [i as f64 * 0.2, (3 - i) as f64 * -0.1, 0.];
        }
        let stations = [
            (WallChord::Top, 2., nodes.insert(Node::new([2., 8.]))),
            (WallChord::Bottom, 10., nodes.insert(Node::new([10., 0.]))),
        ]
        .map(|(chord, x, node)| carapace_core::model::WallStation { chord, x, node });
        let layout = StudWallLayout {
            length: 12.,
            height: 8.,
            studs: [0., 2., 3., 5., 7., 9., 10., 12.]
                .into_iter()
                .map(|x| WallStud { x, ea: 30_000. })
                .collect(),
            top_ea: 1e8,
            top_ei: Some(1e8),
            top_distribution: ChordDistribution::Beam,
            bottom_ea: 1e8,
            bottom_ei: 1e8,
            openings: vec![carapace_core::model::WallOpening {
                x0: 4.,
                x1: 8.,
                head_height: 4.,
                header_ei: 1e5,
                jack_ea: 1e4,
                door: true,
            }],
            sheathing: Some(WallSheathing { law }),
            foundation: false,
            floor: None,
        };
        let loads = std::array::from_fn(|slot| {
            if slot == 0 {
                vec![carapace_core::model::WallLoad::TopLine {
                    x0: 0.,
                    x1: 12.,
                    wx0: 1.,
                    wx1: 1.,
                    wy0: 0.,
                    wy1: 0.,
                }]
            } else {
                vec![]
            }
        });
        let wall = StudWall::new(ids, stations.to_vec(), layout, loads).unwrap();
        let mut element = Element::StudWall(wall);
        let load = carapace_core::model::ElementLoad {
            wall: [1., 0., 0., 0., 0., 0., 0., 0.],
            ..Default::default()
        };
        assert_element_tangent_matches_fd(
            &mut element,
            &mut nodes,
            Some(&load),
            &FdConfig {
                h_rel: 1e-4,
                h_min: 2e-3,
                tol_abs: 2e-6,
                tol_rel: 1e-6,
            },
        );
    }
}
