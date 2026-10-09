//! Element loads on continuum elements: accumulation, compatibility checks,
//! consistent nodal forces (total force and pressure sign), self-weight, and a
//! patch test driven entirely by boundary tractions.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    Domain, ElasticBeamColumn, Element, ElementId, ElementLoad, ElementLoadComponents, GeomTransf,
    LoadSeries, Material, ModelError, Node, NodeId, PlaneMaterial, Quad4, Quad4Formulation, Tri3,
    Truss,
};

fn analysis(domain: Domain) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain)
}

fn material() -> PlaneMaterial {
    PlaneMaterial::plane_stress(1000.0, 0.25).unwrap()
}

/// `ElementLoad` values add and scale field by field (beam uniform, body force, edge traction, edge
/// pressure), and expose a fixed 16-component layout via `component(i)`, which is what element-load
/// recorders index into.
#[test]
fn loads_accumulate_field_by_field_and_expose_a_fixed_component_layout() {
    let both = ElementLoad::uniform(1.0, 2.0) + ElementLoad::body(3.0, 4.0);
    assert_eq!(both.beam_uniform, [1.0, 2.0]);
    assert_eq!(both.body, [3.0, 4.0]);
    let summed = ElementLoad::edge_traction(1, 5.0, 6.0)
        + ElementLoad::edge_pressure(1, 7.0)
        + ElementLoad::edge_pressure(1, 1.0);
    assert_eq!(summed.edge_traction[1], [5.0, 6.0]);
    assert_eq!(summed.edge_pressure[1], 8.0);
    let scaled = both * 2.0;
    assert_eq!((scaled.beam_uniform, scaled.body), ([2.0, 4.0], [6.0, 8.0]));

    // [wx, wy, bx, by, t0x, t0y, p0, t1x, t1y, p1, ...]
    assert_eq!(ElementLoad::COUNT, 16);
    let load = both + summed;
    let expect = [1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 5.0, 6.0, 8.0];
    for (i, e) in expect.into_iter().enumerate() {
        assert_eq!(load.component(i), e, "component {i}");
    }
    assert_eq!(load.component(16), 0.0);
}

/// Body-force loads from several patterns on the same Quad4 add up with each pattern's load factor
/// and scale: the displacement from split patterns equals that from one pattern carrying the summed
/// load.
#[test]
fn several_body_force_patterns_sum_with_their_factors() {
    let build = |patterns: &[(LoadSeries, f64, [f64; 2])]| {
        let mut domain = Domain::new();
        let n: Vec<NodeId> = [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]
            .into_iter()
            .map(|c| {
                let mut node = Node::new(c);
                if c[0] == 0.0 {
                    node = node.fix(0).fix(1);
                }
                domain.add_node(node)
            })
            .collect();
        let quad = domain.add_element(Element::Quad4(Quad4::new(
            [n[0], n[1], n[2], n[3]],
            0.5,
            material(),
        )));
        for (series, scale, body) in patterns {
            let pattern = domain.add_load_pattern_scaled(series.clone(), *scale);
            domain.add_element_load(pattern, quad, ElementLoad::body(body[0], body[1]));
        }
        let mut analysis = analysis(domain);
        analysis.step().unwrap();
        analysis.domain().node(n[2]).displacement
    };
    // At pseudo-time 1: (0, -1) * 1 + (0, -2) * 0.5 + (3, 0) * 2 = (6, -2).
    let summed = build(&[
        (LoadSeries::Linear { slope: 1.0 }, 1.0, [0.0, -1.0]),
        (LoadSeries::Constant, 0.5, [0.0, -2.0]),
        (LoadSeries::Linear { slope: 1.0 }, 2.0, [3.0, 0.0]),
    ]);
    let single = build(&[(LoadSeries::Linear { slope: 1.0 }, 1.0, [6.0, -2.0])]);
    for i in 0..2 {
        assert!(
            (summed[i] - single[i]).abs() < 1e-14,
            "{summed:?} vs {single:?}"
        );
    }
    assert!(single[0] != 0.0 && single[1] != 0.0);
}

/// Loads that do not fit their element (beam loads on continuum elements, continuum loads on a
/// beam, edges out of range) are reported as model errors by `validate`, before any load is
/// accumulated.
#[test]
fn incompatible_loads_are_a_model_error_before_anything_is_accumulated() {
    let incompatible = |element_of: &dyn Fn(&mut Domain) -> ElementId, load: ElementLoad| {
        let mut domain = Domain::new();
        let id = element_of(&mut domain);
        let pattern = domain.add_load_pattern_scaled(LoadSeries::Constant, 1.0);
        domain.add_element_load(pattern, id, load);
        domain.validate()
    };
    let square = |domain: &mut Domain| -> [NodeId; 4] {
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]].map(|c| domain.add_node(Node::new(c)))
    };
    let quad = |domain: &mut Domain| {
        let n = square(domain);
        domain.add_element(Element::Quad4(Quad4::new(n, 1.0, material())))
    };
    let tri = |domain: &mut Domain| {
        let n = square(domain);
        domain.add_element(Element::Tri3(Tri3::new(n[0], n[1], n[2], 1.0, material())))
    };
    let beam = |domain: &mut Domain| {
        let a = domain.add_node(Node::new([0.0, 0.0]));
        let b = domain.add_node(Node::new([3.0, 0.0]));
        domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            a,
            b,
            1000.0,
            1.0,
            1.0,
            GeomTransf::Linear,
        )))
    };
    let truss = |domain: &mut Domain| {
        let a = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1));
        let b = domain.add_node(Node::new([3.0, 0.0]).fix(1));
        domain.add_element(Element::Truss(Truss::new(
            a,
            b,
            1.0,
            Material::Elastic { e: 1000.0 },
        )))
    };
    let rejected = Err(ModelError::IncompatibleElementLoad { element: 0 });

    // A continuum element rejects a beam load; a beam rejects body and edge loads; a
    // truss carries none; a triangle has no fourth edge.
    assert_eq!(
        incompatible(&quad, ElementLoad::uniform(1.0, 0.0)),
        rejected
    );
    assert_eq!(incompatible(&tri, ElementLoad::uniform(1.0, 0.0)), rejected);
    assert_eq!(incompatible(&beam, ElementLoad::body(0.0, -1.0)), rejected);
    assert_eq!(
        incompatible(&beam, ElementLoad::edge_pressure(0, 1.0)),
        rejected
    );
    assert_eq!(
        incompatible(&truss, ElementLoad::uniform(0.0, 1.0)),
        rejected
    );
    assert_eq!(
        incompatible(&tri, ElementLoad::edge_traction(3, 1.0, 0.0)),
        rejected
    );
    // And what each does carry is accepted.
    assert_eq!(
        incompatible(
            &quad,
            ElementLoad::body(0.0, -1.0) + ElementLoad::edge_pressure(3, 1.0)
        ),
        Ok(())
    );
    assert_eq!(
        incompatible(&tri, ElementLoad::edge_traction(2, 1.0, 0.0)),
        Ok(())
    );
    assert_eq!(incompatible(&beam, ElementLoad::uniform(0.0, -1.0)), Ok(()));
}

/// A column of `n` stacked quads, base fixed, self-weight as a body force: the base
/// reaction is the total weight `b t A`, and the column compresses.
#[test]
fn self_weight_column_reaction_equals_the_total_weight() {
    let (n, w, h, t, b) = (6usize, 0.5, 4.0, 0.3, -2.5);
    let mut domain = Domain::new();
    let nodes: Vec<[NodeId; 2]> = (0..=n)
        .map(|j| {
            let y = h * j as f64 / n as f64;
            [0.0, w].map(|x| {
                let mut node = Node::new([x, y]);
                if j == 0 {
                    node = node.fix(0).fix(1);
                }
                domain.add_node(node)
            })
        })
        .collect();
    let pattern = domain.add_load_pattern_scaled(LoadSeries::Linear { slope: 1.0 }, 1.0);
    for j in 0..n {
        let quad = domain.add_element(Element::Quad4(Quad4::new(
            [nodes[j][0], nodes[j][1], nodes[j + 1][1], nodes[j + 1][0]],
            t,
            material(),
        )));
        domain.add_element_load(pattern, quad, ElementLoad::body(0.0, b));
    }
    let mut analysis = analysis(domain);
    analysis.step().unwrap();
    let total_weight = b * t * w * h;
    let reaction: f64 = nodes[0]
        .iter()
        .map(|&id| analysis.domain().reaction(id, 1, 1.0))
        .sum();
    assert!(
        (reaction + total_weight).abs() < 1e-9,
        "{reaction} vs {total_weight}"
    );
    assert!(analysis.domain().node(nodes[n][0]).displacement[1] < 0.0);
}

/// Boundary of the patch-test mesh: a rectangle with perturbed interior nodes.
const NX: usize = 4;
const NY: usize = 3;

fn coords(i: usize, j: usize) -> [f64; 2] {
    let (mut x, mut y) = (2.0 * i as f64, 1.5 * j as f64);
    if i != 0 && i != NX && j != 0 && j != NY {
        x += 0.45 * ((3 * i + 5 * j) as f64).sin();
        y += 0.35 * ((7 * i + 2 * j) as f64).cos();
    }
    [x, y]
}

/// Constant stress `[sx, sy, txy]` applied purely as boundary tractions `sigma . n`; the
/// interior must then carry exactly that stress at every Gauss point.
fn edge_load_patch(quads: bool, formulation: Quad4Formulation) {
    let sigma = [2.0, -1.0, 0.5];
    let mut domain = Domain::new();
    let mut ids = Vec::new();
    for j in 0..=NY {
        for i in 0..=NX {
            let mut node = Node::new(coords(i, j));
            // Supports that remove only the rigid-body modes.
            if (i, j) == (0, 0) {
                node = node.fix(0).fix(1);
            }
            if (i, j) == (NX, 0) {
                node = node.fix(1);
            }
            ids.push(domain.add_node(node));
        }
    }
    let at = |i: usize, j: usize| ids[j * (NX + 1) + i];
    let pattern = domain.add_load_pattern_scaled(LoadSeries::Linear { slope: 1.0 }, 1.0);
    let mat = material();
    // The traction on a boundary edge from grid indices `p -> q`, if both lie on one side.
    let traction = |p: (usize, usize), q: (usize, usize)| -> Option<(f64, f64)> {
        let normal = if p.1 == 0 && q.1 == 0 {
            (0.0, -1.0)
        } else if p.1 == NY && q.1 == NY {
            (0.0, 1.0)
        } else if p.0 == 0 && q.0 == 0 {
            (-1.0, 0.0)
        } else if p.0 == NX && q.0 == NX {
            (1.0, 0.0)
        } else {
            return None;
        };
        Some((
            sigma[0] * normal.0 + sigma[2] * normal.1,
            sigma[2] * normal.0 + sigma[1] * normal.1,
        ))
    };
    let mut elements: Vec<ElementId> = Vec::new();
    for j in 0..NY {
        for i in 0..NX {
            let corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)];
            let cells: Vec<Vec<(usize, usize)>> = if quads {
                vec![corners.to_vec()]
            } else if (i + j) % 2 == 0 {
                vec![
                    vec![corners[0], corners[1], corners[2]],
                    vec![corners[0], corners[2], corners[3]],
                ]
            } else {
                vec![
                    vec![corners[0], corners[1], corners[3]],
                    vec![corners[1], corners[2], corners[3]],
                ]
            };
            for cell in cells {
                let n: Vec<NodeId> = cell.iter().map(|&(a, b)| at(a, b)).collect();
                let id = if quads {
                    domain.add_element(Element::Quad4(
                        Quad4::new([n[0], n[1], n[2], n[3]], 0.3, mat.clone())
                            .with_formulation(formulation),
                    ))
                } else {
                    domain.add_element(Element::Tri3(Tri3::new(n[0], n[1], n[2], 0.3, mat.clone())))
                };
                elements.push(id);
                for k in 0..cell.len() {
                    if let Some((tx, ty)) = traction(cell[k], cell[(k + 1) % cell.len()]) {
                        domain.add_element_load(pattern, id, ElementLoad::edge_traction(k, tx, ty));
                    }
                }
            }
        }
    }
    let mut analysis = analysis(domain);
    analysis.step().unwrap();
    for id in elements {
        for point in analysis.domain().element_gauss_responses(id).unwrap() {
            for c in 0..3 {
                assert!(
                    (point.stress[c] - sigma[c]).abs() < 1e-9,
                    "quads {quads} {formulation:?}: {point:?}"
                );
            }
        }
    }
}

/// A patch test driven entirely by boundary edge loads (no prescribed boundary displacements)
/// reproduces the exact constant stress field, for Quad4 (full and enhanced) and Tri3.
#[test]
fn boundary_tractions_alone_reproduce_the_exact_constant_stress_field() {
    edge_load_patch(false, Quad4Formulation::Full);
    edge_load_patch(true, Quad4Formulation::Full);
    edge_load_patch(true, Quad4Formulation::Enhanced);
}
