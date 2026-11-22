//! The element interface is not tied to two-node elements or to a full
//! per-node DOF layout: an element reports its own nodes and pushes its own
//! contributions through the sinks. These tests drive a user-defined catalog
//! with a three-node element (compact layout: `ux` only) and a four-node
//! element (`ux, uy`) through the real `Domain`/`Analysis` stack and check
//! tangent, resistance, loads, lumped mass, reactions and recorded forces
//! against hand calculations.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    DofMask, DofRef, Domain, ElementForce, ElementLoadComponents, ElementOps, Node, NodeId,
    NodeList, NodeView, TangentSink, VectorSink,
};
use nalgebra::{SMatrix, SVector};
use slotmap::new_key_type;

new_key_type! { struct MockId; }

/// A constant nodal force applied to each node of a `Block4`, per unit factor.
#[derive(Clone, Copy, Debug, PartialEq)]
struct NodalPush {
    fx: f64,
    fy: f64,
}

impl std::ops::Add for NodalPush {
    type Output = NodalPush;
    fn add(self, o: NodalPush) -> NodalPush {
        NodalPush {
            fx: self.fx + o.fx,
            fy: self.fy + o.fy,
        }
    }
}

impl std::ops::Mul<f64> for NodalPush {
    type Output = NodalPush;
    fn mul(self, f: f64) -> NodalPush {
        NodalPush {
            fx: self.fx * f,
            fy: self.fy * f,
        }
    }
}

impl ElementLoadComponents for NodalPush {
    const COUNT: usize = 2;

    fn component(&self, index: usize) -> f64 {
        [self.fx, self.fy].get(index).copied().unwrap_or(0.0)
    }
}

#[derive(Clone)]
enum Mock {
    /// Three nodes in a line along x, two springs of stiffness `k` in series,
    /// acting on `ux` only: `K = k * [[1,-1,0],[-1,2,-1],[0,-1,1]]`.
    Chain3 { nodes: [NodeId; 3], k: f64 },
    /// Four nodes, `ux`/`uy` only, no stiffness: carries a lumped mass
    /// `mass / 4` per node per translation, and a `NodalPush` load.
    Block4 { nodes: [NodeId; 4], mass: f64 },
    /// Runtime-sized chain of springs, with compact `ux` layout.
    DynChain {
        nodes: Vec<NodeId>,
        k: f64,
        mass: f64,
    },
}

fn dofs<const N: usize>(nodes: &[NodeId], slots: &[u8]) -> [DofRef<NodeId>; N] {
    let mut out = Vec::new();
    for &n in nodes {
        for &s in slots {
            out.push((n, s));
        }
    }
    out.try_into().expect("N matches nodes x slots")
}

impl ElementOps<2, 3, NodeId> for Mock {
    type Load = NodalPush;
    type Id = MockId;

    fn nodes(&self) -> NodeList<'_, NodeId> {
        match self {
            Mock::Chain3 { nodes, .. } => nodes.iter().copied().collect(),
            Mock::Block4 { nodes, .. } => nodes.iter().copied().collect(),
            Mock::DynChain { nodes, .. } => NodeList::Borrowed(nodes),
        }
    }

    fn dof_mask(&self) -> DofMask {
        match self {
            Mock::DynChain { .. } => DofMask::none().with(0),
            _ => DofMask::all(3),
        }
    }

    fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        view: &NodeView<'_, 2, 3, NodeId>,
        _: Option<&NodalPush>,
        sink: &mut S,
    ) {
        if let Mock::Chain3 { nodes, k } = self {
            let u = SVector::<f64, 3>::from_fn(|a, _| view.get(nodes[a]).displacement[0]);
            let kmat =
                SMatrix::<f64, 3, 3>::new(1.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 1.0) * *k;
            sink.add(&dofs::<3>(nodes, &[0]), &kmat, &(kmat * u));
        }
        if let Mock::DynChain { nodes, k, .. } = self {
            let n = nodes.len();
            let mut matrix = vec![0.0; n * n];
            for i in 0..n - 1 {
                matrix[i * n + i] += k;
                matrix[i * n + i + 1] -= k;
                matrix[(i + 1) * n + i] -= k;
                matrix[(i + 1) * n + i + 1] += k;
            }
            let u: Vec<_> = nodes
                .iter()
                .map(|&id| view.get(id).displacement[0])
                .collect();
            let resistance: Vec<_> = (0..n)
                .map(|i| (0..n).map(|j| matrix[i * n + j] * u[j]).sum())
                .collect();
            let refs: Vec<_> = nodes.iter().map(|&id| (id, 0)).collect();
            sink.add_dyn(&refs, &matrix, &resistance);
        }
    }

    fn assemble_load<S: VectorSink<NodeId>>(
        &self,
        _: &NodeView<'_, 2, 3, NodeId>,
        load: Option<&NodalPush>,
        sink: &mut S,
    ) {
        if let (Mock::Block4 { nodes, .. }, Some(push)) = (self, load) {
            let v = SVector::<f64, 8>::from_fn(|a, _| if a % 2 == 0 { push.fx } else { push.fy });
            sink.add(&dofs::<8>(nodes, &[0, 1]), &v);
        }
        if let (Mock::DynChain { nodes, .. }, Some(push)) = (self, load) {
            let refs: Vec<_> = nodes.iter().map(|&id| (id, 0)).collect();
            let values = vec![push.fx; nodes.len()];
            sink.add_dyn(&refs, &values);
        }
    }

    fn assemble_mass<S: VectorSink<NodeId>>(&self, _: &NodeView<'_, 2, 3, NodeId>, sink: &mut S) {
        if let Mock::Block4 { nodes, mass } = self {
            sink.add(
                &dofs::<8>(nodes, &[0, 1]),
                &SVector::<f64, 8>::repeat(mass / 4.0),
            );
        }
        if let Mock::DynChain { nodes, mass, .. } = self {
            let refs: Vec<_> = nodes.iter().map(|&id| (id, 0)).collect();
            let values = vec![*mass; nodes.len()];
            sink.add_dyn(&refs, &values);
        }
    }

    fn commit(&mut self, _: &NodeView<'_, 2, 3, NodeId>, _: Option<&NodalPush>) {}

    fn local_force_width(&self) -> usize {
        match self {
            Mock::Chain3 { .. } => 3,
            Mock::Block4 { .. } => 8,
            Mock::DynChain { nodes, .. } => nodes.len(),
        }
    }

    fn local_force(&self, view: &NodeView<'_, 2, 3, NodeId>) -> ElementForce {
        match self {
            Mock::Chain3 { nodes, k } => {
                let u = SVector::<f64, 3>::from_fn(|a, _| view.get(nodes[a]).displacement[0]);
                let kmat =
                    SMatrix::<f64, 3, 3>::new(1.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 1.0) * *k;
                (kmat * u).into()
            }
            Mock::Block4 { .. } => SVector::<f64, 8>::zeros().into(),
            Mock::DynChain { nodes, k, .. } => {
                let mut force = vec![0.0; nodes.len()];
                for i in 0..nodes.len() - 1 {
                    let delta = k
                        * (view.get(nodes[i]).displacement[0]
                            - view.get(nodes[i + 1]).displacement[0]);
                    force[i] += delta;
                    force[i + 1] -= delta;
                }
                ElementForce::from_vec(force)
            }
        }
    }

    fn local_load_force(
        &self,
        _: &NodeView<'_, 2, 3, NodeId>,
        load: Option<&NodalPush>,
    ) -> ElementForce {
        match (self, load) {
            (Mock::Block4 { .. }, Some(p)) => {
                SVector::<f64, 8>::from_fn(|a, _| if a % 2 == 0 { p.fx } else { p.fy }).into()
            }
            (Mock::DynChain { nodes, .. }, Some(p)) => {
                ElementForce::from_vec(vec![p.fx; nodes.len()])
            }
            _ => SVector::<f64, 8>::zeros().into(),
        }
    }

    fn fiber_responses(&self, _: &NodeView<'_, 2, 3, NodeId>) -> Option<Vec<Vec<(f64, f64)>>> {
        None
    }
}

type MockDomain = Domain<2, 3, NodeId, Mock>;

fn linear_analysis(domain: MockDomain) -> carapace_core::analysis::Analysis<2, 3, NodeId, Mock> {
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

/// Two springs in series, fixed at one end and loaded at the other:
/// `u1 = P/k`, `u2 = 2P/k`. The element has three nodes and only a compact
/// `ux` layout, and its recorded local force is its own three-entry vector.
#[test]
fn three_node_element_with_compact_layout_assembles_and_solves() {
    let (k, p) = (50.0, 20.0);
    let mut domain = MockDomain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let mid = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
    let tip = domain.add_node(Node::new([2.0, 0.0]).fix(1).fix(2));
    let element = domain.add_element(Mock::Chain3 {
        nodes: [fixed, mid, tip],
        k,
    });
    domain.load_node(tip, 0, p);

    let mut analysis = linear_analysis(domain);
    analysis.step().expect("series springs solve");

    let d = analysis.domain();
    assert!((d.node(mid).displacement[0] - p / k).abs() < 1e-12);
    assert!((d.node(tip).displacement[0] - 2.0 * p / k).abs() < 1e-12);

    let f = d.element_local_force(element);
    assert_eq!(f.len(), 3);
    let expected = [-p, 0.0, p];
    for (i, e) in expected.iter().enumerate() {
        assert!((f[i] - e).abs() < 1e-9, "component {i}: {} vs {e}", f[i]);
    }
    // The support reaction balances the load.
    assert!((d.reaction(fixed, 0, 1.0) + p).abs() < 1e-9);
}

/// A four-node element's lumped mass lands `mass / 4` on `ux` and `uy` of
/// every node, and its equivalent load is scaled by the pattern factor and
/// seen by `reaction` only at the nodes it touches.
#[test]
fn four_node_element_mass_load_and_reaction() {
    let mut domain = MockDomain::new();
    let nodes: Vec<NodeId> = (0..4)
        .map(|i| domain.add_node(Node::new([i as f64, 0.0]).fix(2)))
        .collect();
    let untouched = domain.add_node(Node::new([9.0, 9.0]).fix(0).fix(1).fix(2));
    let block = domain.add_element(Mock::Block4 {
        nodes: [nodes[0], nodes[1], nodes[2], nodes[3]],
        mass: 8.0,
    });
    let pattern = domain.default_pattern();
    domain.add_element_load(pattern, block, NodalPush { fx: 3.0, fy: -1.0 });

    let mass = domain.assemble_mass_diagonal();
    assert_eq!(mass.len(), 8);
    assert!(
        mass.iter().all(|&m| (m - 2.0).abs() < 1e-15),
        "mass = {mass:?}"
    );

    // The default pattern is Linear(slope 1): at pseudo-time 0.5 the factor
    // is 0.5. A free node's reaction is `resistance - applied`.
    for &node in &nodes {
        assert!((domain.reaction(node, 0, 0.5) + 1.5).abs() < 1e-12);
        assert!((domain.reaction(node, 1, 0.5) - 0.5).abs() < 1e-12);
    }
    assert_eq!(
        domain.reaction(untouched, 0, 0.5),
        0.0,
        "untouched node sees nothing"
    );

    // End force of the loaded element: local_force (zero) minus factor * load.
    let end = domain.element_end_force(block, 0.5);
    assert_eq!(end.len(), 8);
    for c in 0..8 {
        let expected = if c % 2 == 0 { -1.5 } else { 0.5 };
        assert!((end[c] - expected).abs() < 1e-12, "component {c}");
    }
}

/// Elements of different sizes share nodes in one assembly.
#[test]
fn elements_of_different_sizes_assemble_into_one_system() {
    let (k, p) = (10.0, 5.0);
    let mut domain = MockDomain::new();
    let a = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let b = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
    let c = domain.add_node(Node::new([2.0, 0.0]).fix(1).fix(2));
    // The block adds mass but no stiffness, so its other nodes are fixed.
    let d = domain.add_node(Node::new([3.0, 0.0]).fix(0).fix(1).fix(2));
    let e = domain.add_node(Node::new([4.0, 0.0]).fix(0).fix(1).fix(2));
    domain.add_element(Mock::Chain3 {
        nodes: [a, b, c],
        k,
    });
    domain.add_element(Mock::Block4 {
        nodes: [b, c, d, e],
        mass: 4.0,
    });
    domain.load_node(c, 0, p);

    // Free DOFs are b.ux and c.ux; the block puts mass / 4 = 1.0 on each.
    let mass = domain.assemble_mass_diagonal();
    assert_eq!(mass.len(), 2);
    assert!((mass[0] - 1.0).abs() < 1e-15 && (mass[1] - 1.0).abs() < 1e-15);

    let mut analysis = linear_analysis(domain);
    analysis.step().expect("mixed-size system solves");
    let d = analysis.domain();
    // Springs in series, tip load on c: u_c = 2P/k, u_b = P/k.
    assert!((d.node(c).displacement[0] - 2.0 * p / k).abs() < 1e-12);
    assert!((d.node(b).displacement[0] - p / k).abs() < 1e-12);
}

fn variable_chain_model(extra_parallel_elements: bool) -> (MockDomain, Vec<NodeId>) {
    let mut domain = MockDomain::new();
    let nodes: Vec<_> = (0..7)
        .map(|i| {
            let node = Node::new([i as f64, 0.0]);
            domain.add_node(if i == 0 { node.fix(0) } else { node })
        })
        .collect();
    let element = domain.add_element(Mock::DynChain {
        nodes: nodes.clone(),
        k: 17.0,
        mass: 0.75,
    });
    let pattern = domain.default_pattern();
    domain.add_element_load(pattern, element, NodalPush { fx: 2.5, fy: 0.0 });
    if extra_parallel_elements {
        let filler = [0, 1, 2, 3]
            .map(|i| domain.add_node(Node::new([100.0 + i as f64, 0.0]).fix(0).fix(1).fix(2)));
        for _ in 0..255 {
            domain.add_element(Mock::Block4 {
                nodes: filler,
                mass: 0.0,
            });
        }
    }
    (domain, nodes)
}

fn equivalent_chain_model() -> (MockDomain, Vec<NodeId>) {
    let mut domain = MockDomain::new();
    let nodes: Vec<_> = (0..7)
        .map(|i| {
            let node = Node::new([i as f64, 0.0]).fix(1).fix(2).with_mass(0, 0.75);
            domain.add_node(if i == 0 { node.fix(0) } else { node })
        })
        .collect();
    for i in [0, 2, 4] {
        domain.add_element(Mock::Chain3 {
            nodes: [nodes[i], nodes[i + 1], nodes[i + 2]],
            k: 17.0,
        });
    }
    for &node in &nodes {
        domain.load_node(node, 0, 2.5);
    }
    (domain, nodes)
}

/// A seven-node runtime-sized element scatters tangent, resistance, load
/// and mass through slice sinks. Its solve and reactions match three
/// existing three-node chains plus nodal loads and masses.
#[test]
fn variable_connectivity_slice_sinks_match_fixed_size_reference() {
    let (mut variable, variable_nodes) = variable_chain_model(false);
    let (mut reference, reference_nodes) = equivalent_chain_model();
    let variable_mass = variable.assemble_mass_diagonal();
    let reference_mass = reference.assemble_mass_diagonal();
    assert_eq!(
        variable.num_free_dofs(),
        6,
        "nodes() and the ux-only dof_mask activate exactly six unknowns"
    );
    assert_eq!(variable_mass.as_slice(), reference_mass.as_slice());

    let mut variable = linear_analysis(variable);
    let mut reference = linear_analysis(reference);
    variable.step().expect("variable element solves");
    reference.step().expect("fixed-size reference solves");
    for (&a, &b) in variable_nodes.iter().zip(&reference_nodes) {
        assert!(
            (variable.domain().node(a).displacement[0]
                - reference.domain().node(b).displacement[0])
                .abs()
                < 1e-12
        );
        assert!(
            (variable.domain().reaction(a, 0, 1.0) - reference.domain().reaction(b, 0, 1.0)).abs()
                < 1e-10
        );
    }
}

#[cfg(feature = "parallel")]
#[test]
fn variable_connectivity_serial_and_parallel_scatter_are_bit_identical() {
    use carapace_core::analysis::{set_parallelism, Parallelism};

    let (serial, nodes) = variable_chain_model(true);
    let (parallel, parallel_nodes) = variable_chain_model(true);
    set_parallelism(Parallelism::Serial);
    let mut serial = linear_analysis(serial);
    serial.step().expect("serial solve");
    set_parallelism(Parallelism::Threads(2));
    let mut parallel = linear_analysis(parallel);
    parallel.step().expect("parallel solve");
    for (&a, &b) in nodes.iter().zip(&parallel_nodes) {
        assert_eq!(
            serial.domain().node(a).displacement[0].to_bits(),
            parallel.domain().node(b).displacement[0].to_bits()
        );
        assert_eq!(
            serial.domain().reaction(a, 0, 1.0).to_bits(),
            parallel.domain().reaction(b, 0, 1.0).to_bits()
        );
    }
    set_parallelism(Parallelism::Serial);
}
