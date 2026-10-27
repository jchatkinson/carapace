//! Plan `docs/csfm-plan.md` §4 I, "Native vs wire": the CSFM span-beam golden model (15 x 7 Q4
//! mesh, two embedded `CsfmRebar` paths, 20 load-control increments) converted to a
//! `CarapaceInputV1` and run through a `Session` must equal the same model built and stepped
//! natively with `carapace-core`, bit for bit, at every increment: node displacements, the bulk
//! Gauss strain/stress/`kc2` of the concrete, the bulk bar force/strain/stress/rupture, and the
//! controlling material-limit utilization with its location.
//!
//! The native model adds elements in the order the decoder does (triangles, quads, bars), so the
//! assembly sums in the same order and identical results are expected, not just close ones.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::continuum::embed_polyline;
use carapace_core::model::{
    CsfmConcreteLaw, Domain, Element, ElementId, LimitStateConfig, Material, Node, NodeId,
    PlaneMaterial, Quad4,
};
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, BulkQuantity, BulkTable, ConvergenceSpec, IntegratorSpec, RecorderSpec,
    SequenceSpec, StageSpec,
};
use carapace_wasm::input_v1::tables::ElementKind;
use carapace_wasm::input_v1::{decode, CarapaceInputV1, LimitKindDetail};
use serde_json::{json, Value};

use crate::common::empty_input;

const NAME: &str = "span-beam-5kip-gamma-0.001.json";

fn fixture() -> Value {
    let path = format!(
        "{}/../core/tests/csfm/fixtures/models/{NAME}",
        env!("CARGO_MANIFEST_DIR")
    );
    serde_json::from_str(&std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}")))
        .unwrap()
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

fn u(v: &Value) -> usize {
    v.as_u64().unwrap() as usize
}

struct Params {
    fck: f64,
    law: CsfmConcreteLaw,
    law_name: &'static str,
    softening: bool,
    gamma: Option<f64>,
    fctm: f64,
    ecm: f64,
    lambda: f64,
    tension_stiffening: bool,
    thickness: f64,
}

fn params(model: &Value) -> Params {
    let concrete = &model["concrete"];
    let (law, law_name) = match concrete["law"].as_str().unwrap() {
        "parabola-rectangle" => (CsfmConcreteLaw::ParabolaRectangle, "parabolaRectangle"),
        "elastic-ideally-plastic" => (
            CsfmConcreteLaw::ElasticIdeallyPlastic,
            "elasticIdeallyPlastic",
        ),
        other => panic!("{other}"),
    };
    let softening = concrete["softeningOn"].as_bool().unwrap();
    let gamma = concrete
        .get("tensileRegularization")
        .and_then(|r| r.get("gamma"))
        .map(f);
    let fck = f(&model["fck"]);
    let PlaneMaterial::CsfmConcrete(fields) =
        PlaneMaterial::csfm_concrete(fck, None, law, softening, gamma).unwrap()
    else {
        unreachable!()
    };
    Params {
        fck,
        law,
        law_name,
        softening,
        gamma,
        fctm: fields.fctm,
        ecm: fields.ecm,
        lambda: f(&model["lambda"]),
        tension_stiffening: model["tensionStiffening"].as_bool().unwrap(),
        thickness: f(&model["mesh"]["thickness"]),
    }
}

/// `(es, fy, ft, eps_u, diameter, bar area, rho_eff, polyline)` per rebar path.
type Path = (f64, f64, f64, f64, f64, f64, f64, Vec<[f64; 2]>);

fn paths(model: &Value) -> Vec<Path> {
    model["rebarPaths"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| {
            let area = f(&p["barAreaMm2"]) * f(&p["barCount"]);
            (
                f(&p["elasticModulusMpa"]),
                f(&p["fyMpa"]),
                f(&p["ftMpa"]),
                f(&p["ultimateStrain"]),
                f(&p["barDiameterMm"]),
                area,
                area / f(&p["effectiveAreaMm2"]),
                p["polyline"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| [f(&v["x"]), f(&v["y"])])
                    .collect(),
            )
        })
        .collect()
}

fn max_load(model: &Value) -> f64 {
    model["loads"]
        .as_array()
        .unwrap()
        .iter()
        .map(|l| f(&l["value"]).abs())
        .fold(0.0, f64::max)
}

/// The fixture model as a wire input: bulk recorders first (indices 0..7), then every node
/// displacement (7 + 2 * node + dof).
fn wire_input(fx: &Value) -> CarapaceInputV1 {
    let model = &fx["model"];
    let p = params(model);
    let mesh = &model["mesh"];
    let node_x = mesh["nodeX"].as_array().unwrap();
    let node_y = mesh["nodeY"].as_array().unwrap();
    let mut fixed = vec![0u8; node_x.len()];
    for bc in model["boundaryConditions"].as_array().unwrap() {
        fixed[u(&bc["node"])] |= 1 << u(&bc["dof"]);
    }
    let offsets: Vec<usize> = mesh["elementOffsets"]
        .as_array()
        .unwrap()
        .iter()
        .map(u)
        .collect();
    let nodes: Vec<usize> = mesh["elementNodes"]
        .as_array()
        .unwrap()
        .iter()
        .map(u)
        .collect();
    let mut quad_nodes = Vec::new();
    for (e, kind) in mesh["elementTypes"].as_array().unwrap().iter().enumerate() {
        assert_eq!(kind, "Q4", "the span beam is all Q4");
        quad_nodes.extend_from_slice(&nodes[offsets[e]..offsets[e] + 4]);
    }
    let quads = quad_nodes.len() / 4;

    let mut input = serde_json::to_value(empty_input(2)).unwrap();
    input["nodes"] = json!({
        "coords": node_x.iter().zip(node_y).flat_map(|(x, y)| [f(x), f(y)]).collect::<Vec<_>>(),
        "fixed": fixed, "massNodeIndex": [], "mass": [],
    });
    input["planeMaterials"] = json!([{
        "kind": "csfmConcrete", "fck": p.fck, "law": p.law_name, "softening": p.softening,
        "tensileRegularization": p.gamma,
    }]);
    input["quads"] = json!({
        "nodeIds": quad_nodes, "thickness": vec![p.thickness; quads], "material": vec![0; quads],
        "density": vec![0.0; quads], "formulation": vec!["full"; quads],
    });
    let rebar_paths = paths(model);
    input["materials"] = Value::Array(
        rebar_paths
            .iter()
            .map(|(es, fy, ft, eps_u, d, _, rho, _)| {
                json!({
                    "kind": "csfmRebar", "es": es, "fy": fy, "ft": ft, "epsU": eps_u,
                    "fctm": p.fctm, "ecm": p.ecm, "diameter": d, "rhoEff": rho,
                    "lambda": p.lambda, "tensionStiffening": p.tension_stiffening,
                })
            })
            .collect(),
    );
    input["embeddedBars"] = Value::Array(
        rebar_paths
            .iter()
            .enumerate()
            .map(|(m, (_, _, _, _, _, area, _, line))| {
                json!({
                    "polyline": line.iter().flat_map(|v| *v).collect::<Vec<_>>(),
                    "area": area, "material": m,
                })
            })
            .collect(),
    );
    let loads = model["loads"].as_array().unwrap();
    input["loadPatterns"] =
        json!({ "series": [{ "kind": "linear", "slope": 1.0 }], "scaleFactor": [1.0] });
    input["nodalLoads"] = json!({
        "pattern": vec![0; loads.len()],
        "node": loads.iter().map(|l| u(&l["node"])).collect::<Vec<_>>(),
        "dof": loads.iter().map(|l| u(&l["dof"])).collect::<Vec<_>>(),
        "value": loads.iter().map(|l| f(&l["value"])).collect::<Vec<_>>(),
        "stage": vec![0; loads.len()],
    });
    let mut input: CarapaceInputV1 = serde_json::from_value(input).unwrap();

    let schedule = fx["schedule"].as_array().unwrap();
    let tolerance = (1e-10 * max_load(model)).max(1e-10);
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "load".to_string(),
            steps: schedule.len() as u32,
            integrator: IntegratorSpec::LoadControl {
                increment: f(&schedule[0]),
            },
            algorithm: AlgorithmSpec::Newton {
                tangent: Default::default(),
                line_search: None,
            },
            convergence: Some(ConvergenceSpec::Combined {
                force_tol: tolerance,
                moment_tol: Some(tolerance),
                relative_tol: Some(1e-10),
                displacement_tol: None,
                max_iter: fx["solver"]["maxIterationsPerIncrement"].as_u64().unwrap() as u32,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: {
            let mut recorders = vec![
                RecorderSpec::Bulk {
                    table: BulkTable::Quads,
                    quantity: BulkQuantity::Strain,
                },
                RecorderSpec::Bulk {
                    table: BulkTable::Quads,
                    quantity: BulkQuantity::Stress,
                },
                RecorderSpec::Bulk {
                    table: BulkTable::Quads,
                    quantity: BulkQuantity::Kc2,
                },
                RecorderSpec::Bulk {
                    table: BulkTable::EmbeddedBars,
                    quantity: BulkQuantity::AxialForce,
                },
                RecorderSpec::Bulk {
                    table: BulkTable::EmbeddedBars,
                    quantity: BulkQuantity::Strain,
                },
                RecorderSpec::Bulk {
                    table: BulkTable::EmbeddedBars,
                    quantity: BulkQuantity::Stress,
                },
                RecorderSpec::Bulk {
                    table: BulkTable::EmbeddedBars,
                    quantity: BulkQuantity::Ruptured,
                },
            ];
            for node in 0..node_x.len() as u32 {
                for dof in 0..2 {
                    recorders.push(RecorderSpec::NodeDisp { node, dof });
                }
            }
            recorders
        },
    };
    input
}

#[derive(Default)]
struct Frame {
    displacements: Vec<f64>,
    strain: Vec<f64>,
    stress: Vec<f64>,
    kc2: Vec<f64>,
    bar_force: Vec<f64>,
    bar_strain: Vec<f64>,
    bar_stress: Vec<f64>,
    bar_ruptured: Vec<f64>,
    limit: Option<(f64, LimitKindDetail, ElementKind, u32)>,
}

fn run_native(fx: &Value) -> Vec<Frame> {
    let model = &fx["model"];
    let p = params(model);
    let mesh = &model["mesh"];
    let node_x = mesh["nodeX"].as_array().unwrap();
    let node_y = mesh["nodeY"].as_array().unwrap();
    let mut fixed: Vec<Vec<usize>> = vec![Vec::new(); node_x.len()];
    for bc in model["boundaryConditions"].as_array().unwrap() {
        fixed[u(&bc["node"])].push(u(&bc["dof"]));
    }
    let mut domain: Domain = Domain::new();
    let node_ids: Vec<NodeId> = (0..node_x.len())
        .map(|i| {
            let mut node = Node::new([f(&node_x[i]), f(&node_y[i])]);
            for &dof in &fixed[i] {
                node = node.fix(dof);
            }
            domain.add_node(node)
        })
        .collect();
    let concrete = PlaneMaterial::csfm_concrete(p.fck, None, p.law, p.softening, p.gamma).unwrap();
    let offsets: Vec<usize> = mesh["elementOffsets"]
        .as_array()
        .unwrap()
        .iter()
        .map(u)
        .collect();
    let nodes: Vec<usize> = mesh["elementNodes"]
        .as_array()
        .unwrap()
        .iter()
        .map(u)
        .collect();
    let quads: Vec<ElementId> = (0..mesh["elementTypes"].as_array().unwrap().len())
        .map(|e| {
            let n: [NodeId; 4] = std::array::from_fn(|i| node_ids[nodes[offsets[e] + i]]);
            domain.add_element(Element::Quad4(Quad4::new(n, p.thickness, concrete.clone())))
        })
        .collect();
    let mut bars: Vec<ElementId> = Vec::new();
    for (es, fy, ft, eps_u, d, area, rho, line) in paths(model) {
        let material = Material::csfm_rebar(
            es,
            fy,
            ft,
            eps_u,
            p.fctm,
            p.ecm,
            d,
            rho,
            p.lambda,
            p.tension_stiffening,
        )
        .unwrap();
        for segment in embed_polyline(&domain, &line, area, &material).unwrap() {
            bars.push(domain.add_element(Element::EmbeddedBar(segment)));
        }
    }
    let pattern = domain.add_load_pattern(carapace_core::model::LoadSeries::Linear { slope: 1.0 });
    for load in model["loads"].as_array().unwrap() {
        domain.add_nodal_load(
            pattern,
            node_ids[u(&load["node"])],
            u(&load["dof"]),
            f(&load["value"]),
        );
    }

    let schedule = fx["schedule"].as_array().unwrap();
    let tolerance = (1e-10 * max_load(model)).max(1e-10);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl {
            increment: f(&schedule[0]),
        })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::Combined {
            force_tol: tolerance,
            moment_tol: tolerance,
            relative_tol: 1e-10,
            displacement_tol: None,
            max_iter: fx["solver"]["maxIterationsPerIncrement"].as_u64().unwrap() as usize,
        })
        .build(domain);

    let mut frames = Vec::new();
    for _ in schedule {
        analysis.step().expect("native span beam converges");
        let d = analysis.domain();
        let mut frame = Frame::default();
        for &node in &node_ids {
            frame
                .displacements
                .extend_from_slice(&d.node(node).displacement[..2]);
        }
        for &q in &quads {
            for (g, response) in d.element_gauss_responses(q).unwrap().iter().enumerate() {
                frame.strain.extend_from_slice(&response.strain);
                frame.stress.extend_from_slice(&response.stress);
                frame.kc2.push(d.element_gauss_kc2(q, g).unwrap());
            }
        }
        for &b in &bars {
            let (force, strain, stress, ruptured) = d.element_embedded_bar_response(b).unwrap();
            frame.bar_force.push(force);
            frame.bar_strain.push(strain);
            frame.bar_stress.push(stress);
            frame.bar_ruptured.push(f64::from(ruptured));
        }
        frame.limit = d
            .controlling_limit_utilization(&LimitStateConfig::default())
            .map(|(id, utilization, kind)| {
                let (element_kind, index) = match quads.iter().position(|&q| q == id) {
                    Some(i) => (ElementKind::Quad4, i),
                    None => (
                        ElementKind::EmbeddedBar,
                        bars.iter().position(|&b| b == id).unwrap(),
                    ),
                };
                let kind = match kind {
                    carapace_core::model::LimitKind::ConcreteCrushing => {
                        LimitKindDetail::ConcreteCrushing
                    }
                    carapace_core::model::LimitKind::ConcreteTensileStrain => {
                        LimitKindDetail::ConcreteTensileStrain
                    }
                    carapace_core::model::LimitKind::RebarRupture => LimitKindDetail::RebarRupture,
                    carapace_core::model::LimitKind::RebarStrainLimit => {
                        LimitKindDetail::RebarStrainLimit
                    }
                };
                (utilization, kind, element_kind, index as u32)
            });
        frames.push(frame);
    }
    frames
}

/// One `advance` budget of the whole run, and 3-step budgets: the result must not depend on how
/// the host slices the run.
fn run_wire(fx: &Value, budget: u32) -> (Vec<Frame>, f64) {
    let mut session = decode(wire_input(fx)).expect("span beam decodes");
    let mut bulk: Vec<Vec<Vec<f64>>> = vec![Vec::new(); 7]; // recorder -> frames
    let mut displacements: Vec<Vec<f64>> = Vec::new();
    let mut limits = Vec::new();
    let mut disp_samples: Vec<Vec<f64>> = vec![Vec::new(); 256];
    let load_factor;
    loop {
        let outcome = session.advance(budget);
        assert!(outcome.error.is_none(), "{:?}", outcome.error);
        for batch in &outcome.bulk_batches {
            let width = (batch.elements * batch.points * batch.components) as usize;
            for frame in batch.values.chunks(width) {
                bulk[batch.recorder_index].push(frame.to_vec());
            }
        }
        for batch in &outcome.recorder_batches {
            disp_samples[batch.recorder_index - 7]
                .extend(batch.samples.iter().map(|&(_, value)| value));
        }
        limits.extend(outcome.limit_history.clone());
        if outcome.done {
            load_factor = outcome.load_factor;
            break;
        }
    }
    let steps = bulk[0].len();
    for k in 0..steps {
        displacements.push(disp_samples.iter().map(|series| series[k]).collect());
    }
    let frames = (0..steps)
        .map(|k| Frame {
            displacements: displacements[k].clone(),
            strain: bulk[0][k].clone(),
            stress: bulk[1][k].clone(),
            kc2: bulk[2][k].clone(),
            bar_force: bulk[3][k].clone(),
            bar_strain: bulk[4][k].clone(),
            bar_stress: bulk[5][k].clone(),
            bar_ruptured: bulk[6][k].clone(),
            limit: limits
                .get(k)
                .map(|s| (s.utilization, s.kind, s.element_kind, s.element_index)),
        })
        .collect();
    (frames, load_factor)
}

fn assert_identical(native: &[Frame], wire: &[Frame]) {
    assert_eq!(native.len(), wire.len(), "increment count");
    for (i, (n, w)) in native.iter().zip(wire).enumerate() {
        let fields = [
            ("displacements", &n.displacements, &w.displacements),
            ("strain", &n.strain, &w.strain),
            ("stress", &n.stress, &w.stress),
            ("kc2", &n.kc2, &w.kc2),
            ("bar force", &n.bar_force, &w.bar_force),
            ("bar strain", &n.bar_strain, &w.bar_strain),
            ("bar stress", &n.bar_stress, &w.bar_stress),
            ("bar ruptured", &n.bar_ruptured, &w.bar_ruptured),
        ];
        for (name, a, b) in fields {
            assert_eq!(a.len(), b.len(), "increment {i}: {name} length");
            for (k, (x, y)) in a.iter().zip(b).enumerate() {
                assert!(
                    x == y,
                    "increment {i}: {name}[{k}]: native {x:e} vs wire {y:e}"
                );
            }
        }
        assert_eq!(n.limit, w.limit, "increment {i}: controlling limit");
    }
}

#[test]
fn span_beam_through_the_wire_equals_the_native_run_including_recorders() {
    let fx = fixture();
    let native = run_native(&fx);
    assert_eq!(native.len(), 20);
    // The model really is nonlinear and loaded: the bars carry force and the concrete cracks.
    let last = native.last().unwrap();
    assert!(last.bar_force.iter().any(|&force| force > 1e3));
    assert!(last
        .limit
        .is_some_and(|(utilization, ..)| utilization > 0.0));

    let (wire, load_factor) = run_wire(&fx, 1000);
    assert!((load_factor - 1.0).abs() < 1e-12);
    assert_identical(&native, &wire);

    // Slicing the run into small `advance` budgets changes nothing.
    let (sliced, _) = run_wire(&fx, 3);
    assert_identical(&native, &sliced);
}
