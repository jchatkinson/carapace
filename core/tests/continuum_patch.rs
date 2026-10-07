//! Patch tests for the continuum elements: a mesh with irregular interior
//! nodes, boundary displacements from a linear field, interior nodes free.
//! A converged element must reproduce the field exactly at the interior
//! nodes and the constant stress `D * strain` at every Gauss point.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{Domain, Element, ElementId, Node, NodeId, PlaneMaterial, Quad4, Tri3};

const NX: usize = 4;
const NY: usize = 3;

/// Grid node coordinates (row-major, `NX + 1` by `NY + 1`), interior nodes perturbed.
fn grid() -> Vec<[f64; 2]> {
    let mut coords = Vec::new();
    for j in 0..=NY {
        for i in 0..=NX {
            let (mut x, mut y) = (2.0 * i as f64, 1.5 * j as f64);
            if i != 0 && i != NX && j != 0 && j != NY {
                x += 0.45 * ((3 * i + 5 * j) as f64).sin();
                y += 0.35 * ((7 * i + 2 * j) as f64).cos();
            }
            coords.push([x, y]);
        }
    }
    coords
}

fn is_boundary(index: usize) -> bool {
    let (i, j) = (index % (NX + 1), index / (NX + 1));
    i == 0 || i == NX || j == 0 || j == NY
}

/// `u = (a x + b y + c, d x + e y + f)`.
const FIELD: [f64; 6] = [1.2e-3, 4e-4, 3e-4, -7e-4, 9e-4, -2e-4];

fn field(p: [f64; 2]) -> [f64; 2] {
    [
        FIELD[0] * p[0] + FIELD[1] * p[1] + FIELD[2],
        FIELD[3] * p[0] + FIELD[4] * p[1] + FIELD[5],
    ]
}

fn build(domain: Domain) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-12,
            max_iter: 5,
        })
        .build(domain)
}

/// Adds the nodes (boundary translations fixed) and returns their ids.
fn add_nodes(domain: &mut Domain) -> Vec<NodeId> {
    grid()
        .into_iter()
        .enumerate()
        .map(|(index, c)| {
            let mut node = Node::new(c);
            if is_boundary(index) {
                node = node.fix(0).fix(1);
            }
            domain.add_node(node)
        })
        .collect()
}

fn run_patch(domain: Domain, ids: &[NodeId]) -> Analysis {
    let mut analysis = build(domain);
    let prescribed: Vec<(NodeId, usize, f64)> = grid()
        .into_iter()
        .enumerate()
        .filter(|&(index, _)| is_boundary(index))
        .flat_map(|(index, c)| {
            let u = field(c);
            [(ids[index], 0, u[0]), (ids[index], 1, u[1])]
        })
        .collect();
    analysis.step_prescribed(&prescribed).expect("patch solves");
    analysis
}

fn assert_field_reproduced(analysis: &Analysis, ids: &[NodeId], tol: f64) {
    for (index, c) in grid().into_iter().enumerate() {
        let u = field(c);
        let got = analysis.domain().node(ids[index]).displacement;
        for axis in 0..2 {
            assert!(
                (got[axis] - u[axis]).abs() < tol,
                "node {index} axis {axis}: {} vs {}",
                got[axis],
                u[axis]
            );
        }
    }
}

fn assert_constant_stress(
    analysis: &Analysis,
    elements: &[ElementId],
    material: &PlaneMaterial,
    tol: f64,
) {
    let strain = [FIELD[0], FIELD[4], FIELD[1] + FIELD[3]];
    let stress = material.d_matrix() * nalgebra::Vector3::from(strain);
    for &id in elements {
        for point in analysis
            .domain()
            .element_gauss_responses(id)
            .expect("continuum element")
        {
            for c in 0..3 {
                assert!((point.strain[c] - strain[c]).abs() < tol, "{point:?}");
                assert!((point.stress[c] - stress[c]).abs() < tol * 1e5, "{point:?}");
            }
        }
    }
}

fn tri3_mesh(material: &PlaneMaterial) -> (Domain, Vec<NodeId>, Vec<ElementId>) {
    let mut domain = Domain::new();
    let ids = add_nodes(&mut domain);
    let at = |i: usize, j: usize| ids[j * (NX + 1) + i];
    let mut elements = Vec::new();
    for j in 0..NY {
        for i in 0..NX {
            // Alternate the diagonal so the mesh is not a uniform pattern.
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            let cells: [[NodeId; 3]; 2] = if (i + j) % 2 == 0 {
                [[a, b, c], [a, c, d]]
            } else {
                [[a, b, d], [b, c, d]]
            };
            for [n1, n2, n3] in cells {
                elements.push(domain.add_element(Element::Tri3(Tri3::new(
                    n1,
                    n2,
                    n3,
                    0.3,
                    material.clone(),
                ))));
            }
        }
    }
    (domain, ids, elements)
}

#[test]
fn tri3_constant_strain_patch_test_is_exact_on_an_irregular_mesh() {
    for material in [
        PlaneMaterial::plane_stress(30e3, 0.2).unwrap(),
        PlaneMaterial::plane_strain(30e3, 0.3).unwrap(),
        PlaneMaterial::orthotropic(30e3, 10e3, 0.2, 4e3, 0.5).unwrap(),
    ] {
        let (domain, ids, elements) = tri3_mesh(&material);
        let analysis = run_patch(domain, &ids);
        assert_field_reproduced(&analysis, &ids, 1e-12);
        assert_constant_stress(&analysis, &elements, &material, 1e-12);
    }
}

#[test]
fn a_continuum_model_with_no_rotation_fixed_solves() {
    // `rz` is not stiffened by a Tri3, so it is not an equation: no rotational fixity is needed
    // (the patch above fixes only translations), and loading `rz` is a structured error.
    let material = PlaneMaterial::plane_stress(30e3, 0.2).unwrap();
    let (mut domain, ids, _) = tri3_mesh(&material);
    domain.load_node(ids[NX + 2], 2, 1.0);
    assert!(matches!(
        domain.validate(),
        Err(carapace_core::model::ModelError::LoadOnInactiveDof { dof: 2, .. })
    ));
}

fn quad4_mesh(material: &PlaneMaterial) -> (Domain, Vec<NodeId>, Vec<ElementId>) {
    let mut domain = Domain::new();
    let ids = add_nodes(&mut domain);
    let at = |i: usize, j: usize| ids[j * (NX + 1) + i];
    let mut elements = Vec::new();
    for j in 0..NY {
        for i in 0..NX {
            let nodes = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
            elements.push(domain.add_element(Element::Quad4(Quad4::new(
                nodes,
                0.3,
                material.clone(),
            ))));
        }
    }
    (domain, ids, elements)
}

#[test]
fn quad4_patch_test_is_exact_on_a_distorted_mesh() {
    for material in [
        PlaneMaterial::plane_stress(30e3, 0.2).unwrap(),
        PlaneMaterial::plane_strain(30e3, 0.3).unwrap(),
        PlaneMaterial::orthotropic(30e3, 10e3, 0.2, 4e3, 0.5).unwrap(),
    ] {
        let (domain, ids, elements) = quad4_mesh(&material);
        let analysis = run_patch(domain, &ids);
        assert_field_reproduced(&analysis, &ids, 1e-12);
        assert_constant_stress(&analysis, &elements, &material, 1e-12);
    }
}
