//! `Quad4`: the 4-node bilinear isoparametric quadrilateral, 2x2 Gauss
//! integration, `[ux, uy]` per node (8 DOFs), plane stress or strain through
//! its `PlaneMaterial`. Each Gauss point owns an independent copy of the
//! material and commits it individually.

use nalgebra::{SMatrix, SVector};

use super::plane_common::{
    b_matrix, check_section, plane_coords, plane_displacement, plane_dofs, PlaneView,
};
use super::{DofMask, ElementForce, GaussResponse, TangentSink, VectorSink};
use crate::model::continuum::{
    physical_derivatives, quad4_dshape, quad4_shape, validate_quad, Jacobian, QUAD4_GAUSS_2X2,
};
use crate::model::{NodeId, PlaneMaterial, PlaneVector};

#[derive(Debug, Clone)]
pub struct Quad4 {
    pub nodes: [NodeId; 4],
    pub thickness: f64,
    pub density: f64,
    materials: [PlaneMaterial; 4],
    cache: Option<[PointGeometry; 4]>,
}

/// What depends only on the (fixed) coordinates at one Gauss point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct PointGeometry {
    pub b: SMatrix<f64, 3, 8>,
    /// `det J * w`.
    pub weight: f64,
    /// Shape-function values `N_i`.
    pub n: [f64; 4],
}

pub(super) fn point_geometry(coords: &[[f64; 2]; 4]) -> [PointGeometry; 4] {
    std::array::from_fn(|g| {
        let (xi, eta, w) = QUAD4_GAUSS_2X2[g];
        let dn = quad4_dshape(xi, eta);
        let jacobian = Jacobian::new(coords, &dn);
        PointGeometry {
            b: b_matrix(&physical_derivatives(&jacobian, &dn)),
            weight: jacobian.det * w,
            n: quad4_shape(xi, eta),
        }
    })
}

impl Quad4 {
    pub fn new(nodes: [NodeId; 4], thickness: f64, material: PlaneMaterial) -> Self {
        Quad4 {
            nodes,
            thickness,
            density: 0.0,
            materials: std::array::from_fn(|_| material.clone()),
            cache: None,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// The Gauss points' materials, in `QUAD4_GAUSS_2X2` order.
    pub fn materials(&self) -> &[PlaneMaterial; 4] {
        &self.materials
    }

    fn geometry(&self, view: &PlaneView<'_>) -> [PointGeometry; 4] {
        self.cache
            .unwrap_or_else(|| point_geometry(&plane_coords(&self.nodes, view)))
    }

    pub(super) fn dof_mask(&self) -> DofMask {
        DofMask::none().with(0).with(1)
    }

    pub(super) fn validate(&self, view: &PlaneView<'_>) -> Result<(), &'static str> {
        check_section(self.thickness, self.density, &self.materials[0])?;
        validate_quad(&plane_coords(&self.nodes, view))
    }

    pub(super) fn prepare(&mut self, view: &PlaneView<'_>) {
        self.cache = Some(point_geometry(&plane_coords(&self.nodes, view)));
    }

    pub(super) fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        view: &PlaneView<'_>,
        sink: &mut S,
    ) {
        let u = plane_displacement::<4, 8>(&self.nodes, view);
        let mut k = SMatrix::<f64, 8, 8>::zeros();
        let mut r = SVector::<f64, 8>::zeros();
        for (g, material) in self.geometry(view).iter().zip(&self.materials) {
            let (stress, d) = material.trial_stress_tangent(&(g.b * u));
            let scale = self.thickness * g.weight;
            k += scale * g.b.transpose() * d * g.b;
            r += scale * g.b.transpose() * stress;
        }
        sink.add(&plane_dofs::<4, 8>(&self.nodes), &k, &r);
    }

    /// Row-sum lumped mass `m_i = sum_g rho t N_i(g) det J(g) w(g)`, applied to
    /// each translational DOF of node `i`. Exact for the consistent mass here
    /// (bilinear `N_i` times linear `det J`), and equal to `rho t A / 4` only
    /// for parallelograms.
    pub(super) fn assemble_mass<S: VectorSink<NodeId>>(&self, view: &PlaneView<'_>, sink: &mut S) {
        let mut m = SVector::<f64, 8>::zeros();
        for g in self.geometry(view) {
            for i in 0..4 {
                let share = self.density * self.thickness * g.n[i] * g.weight;
                m[2 * i] += share;
                m[2 * i + 1] += share;
            }
        }
        sink.add(&plane_dofs::<4, 8>(&self.nodes), &m);
    }

    pub(super) fn commit(&mut self, view: &PlaneView<'_>) {
        let u = plane_displacement::<4, 8>(&self.nodes, view);
        let geometry = self.geometry(view);
        for (material, g) in self.materials.iter_mut().zip(&geometry) {
            *material = material.commit(&(g.b * u));
        }
    }

    pub(super) fn local_force(&self, view: &PlaneView<'_>) -> ElementForce {
        let u = plane_displacement::<4, 8>(&self.nodes, view);
        let mut r = SVector::<f64, 8>::zeros();
        for (g, material) in self.geometry(view).iter().zip(&self.materials) {
            let (stress, _) = material.trial_stress_tangent(&(g.b * u));
            r += self.thickness * g.weight * g.b.transpose() * stress;
        }
        r.into()
    }

    pub(super) fn gauss_responses(&self, view: &PlaneView<'_>) -> Vec<GaussResponse> {
        let u = plane_displacement::<4, 8>(&self.nodes, view);
        self.geometry(view)
            .iter()
            .zip(&self.materials)
            .map(|(g, material)| {
                let strain: PlaneVector = g.b * u;
                let (stress, _) = material.trial_stress_tangent(&strain);
                GaussResponse {
                    strain: strain.into(),
                    stress: stress.into(),
                }
            })
            .collect()
    }

    #[cfg(test)]
    pub(super) fn cache(&self) -> Option<[PointGeometry; 4]> {
        self.cache
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::SymmetricEigen;
    use slotmap::SlotMap;

    use super::*;
    use crate::model::{DofRef, Node};

    struct Collect {
        k: SMatrix<f64, 8, 8>,
        r: SVector<f64, 8>,
    }

    impl TangentSink<NodeId> for Collect {
        fn add<const N: usize>(
            &mut self,
            _dofs: &[DofRef<NodeId>; N],
            k: &SMatrix<f64, N, N>,
            r: &SVector<f64, N>,
        ) {
            self.k = SMatrix::from_fn(|i, j| k[(i, j)]);
            self.r = SVector::from_fn(|i, _| r[i]);
        }
    }

    struct Mass(SVector<f64, 8>);

    impl VectorSink<NodeId> for Mass {
        fn add<const N: usize>(&mut self, _dofs: &[DofRef<NodeId>; N], v: &SVector<f64, N>) {
            self.0 = SVector::from_fn(|i, _| v[i]);
        }
    }

    fn setup(coords: [[f64; 2]; 4]) -> (SlotMap<NodeId, Node>, [NodeId; 4]) {
        let mut nodes = SlotMap::with_key();
        let ids = coords.map(|c| nodes.insert(Node::new(c)));
        (nodes, ids)
    }

    fn assemble(el: &Quad4, nodes: &SlotMap<NodeId, Node>) -> Collect {
        let mut sink = Collect {
            k: SMatrix::zeros(),
            r: SVector::zeros(),
        };
        el.assemble_tangent(&PlaneView::new(nodes), &mut sink);
        sink
    }

    fn masses(el: &Quad4, nodes: &SlotMap<NodeId, Node>) -> SVector<f64, 8> {
        let mut sink = Mass(SVector::zeros());
        el.assemble_mass(&PlaneView::new(nodes), &mut sink);
        sink.0
    }

    const SKEW: [[f64; 2]; 4] = [[0.0, 0.0], [2.0, 0.2], [2.3, 1.7], [-0.1, 1.4]];

    fn material() -> PlaneMaterial {
        PlaneMaterial::plane_stress(30e3, 0.2).unwrap()
    }

    #[test]
    fn stiffness_matches_the_closed_form_for_the_unit_square() {
        // E = 1, nu = 0, t = 1, unit square (hand integration of the bilinear field):
        // K[ux0, ux0] = 1/3 + 1/6 = 1/2, K[ux0, uy0] = 1/8, K[ux0, ux1] = -1/3 + 1/12 = -1/4.
        let (nodes, ids) = setup([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let el = Quad4::new(ids, 1.0, PlaneMaterial::plane_stress(1.0, 0.0).unwrap());
        let k = assemble(&el, &nodes).k;
        assert!((k[(0, 0)] - 0.5).abs() < 1e-14);
        assert!((k[(0, 1)] - 0.125).abs() < 1e-14);
        assert!((k[(0, 2)] + 0.25).abs() < 1e-14);
        // Rigid translation: every row sums (per direction) to zero.
        for row in 0..8 {
            let sum_x: f64 = (0..4).map(|n| k[(row, 2 * n)]).sum();
            let sum_y: f64 = (0..4).map(|n| k[(row, 2 * n + 1)]).sum();
            assert!(sum_x.abs() < 1e-14 && sum_y.abs() < 1e-14);
        }
    }

    #[test]
    fn stiffness_is_symmetric_with_exactly_three_rigid_body_modes() {
        let (nodes, ids) = setup(SKEW);
        let el = Quad4::new(ids, 0.4, material());
        let k = assemble(&el, &nodes).k;
        assert!((k - k.transpose()).abs().max() < 1e-9 * k.abs().max());
        let eigen = SymmetricEigen::new(k).eigenvalues;
        let largest = eigen.max();
        assert_eq!(
            eigen.iter().filter(|&&l| l.abs() < 1e-9 * largest).count(),
            3,
            "{eigen}"
        );
        assert!(eigen.iter().all(|&l| l > -1e-9 * largest), "{eigen}");
    }

    #[test]
    fn a_rigid_motion_is_stress_free() {
        let (mut nodes, ids) = setup(SKEW);
        let el = Quad4::new(ids, 0.4, material());
        let (t, theta) = ([0.01, -0.02], 1e-5);
        for id in ids {
            let [x, y] = nodes[id].coords;
            nodes[id].displacement = [t[0] - theta * y, t[1] + theta * x, 0.0];
        }
        assert!(assemble(&el, &nodes).r.abs().max() < 1e-9);
    }

    #[test]
    fn cyclic_relabeling_of_the_nodes_gives_the_same_stiffness() {
        let (nodes, ids) = setup(SKEW);
        let base = assemble(&Quad4::new(ids, 0.4, material()), &nodes).k;
        for shift in 1..4 {
            let rotated: [NodeId; 4] = std::array::from_fn(|i| ids[(i + shift) % 4]);
            let k = assemble(&Quad4::new(rotated, 0.4, material()), &nodes).k;
            // Local node i of the rotated element is global node (i + shift) % 4.
            let perm = |a: usize| 2 * ((a / 2 + shift) % 4) + a % 2;
            for a in 0..8 {
                for b in 0..8 {
                    assert!(
                        (k[(a, b)] - base[(perm(a), perm(b))]).abs() < 1e-9 * base.abs().max(),
                        "shift {shift} ({a},{b})"
                    );
                }
            }
        }
    }

    #[test]
    fn lumped_masses_are_the_row_sums_of_the_consistent_mass() {
        // Trapezoid with rho t = 1: nodal masses 5/12, 5/12, 1/3, 1/3 (not A / 4 = 0.375).
        let (nodes, ids) = setup([[0.0, 0.0], [2.0, 0.0], [1.5, 1.0], [0.5, 1.0]]);
        let m = masses(&Quad4::new(ids, 1.0, material()).with_density(1.0), &nodes);
        for (node, expected) in [5.0 / 12.0, 5.0 / 12.0, 1.0 / 3.0, 1.0 / 3.0]
            .into_iter()
            .enumerate()
        {
            assert!(
                (m[2 * node] - expected).abs() < 1e-14
                    && (m[2 * node + 1] - expected).abs() < 1e-14
            );
        }
        // Parallelogram: A / 4 each.
        let (nodes, ids) = setup([[0.0, 0.0], [2.0, 0.0], [3.0, 1.5], [1.0, 1.5]]);
        let m = masses(&Quad4::new(ids, 0.5, material()).with_density(4.0), &nodes);
        assert!(m.iter().all(|&v| (v - 4.0 * 0.5 * 3.0 / 4.0).abs() < 1e-14));
    }

    #[test]
    fn mass_sums_to_rho_t_a_in_each_direction() {
        let (nodes, ids) = setup(SKEW);
        // Shoelace area of SKEW.
        let area = 0.5
            * (0..4)
                .map(|i| SKEW[i][0] * SKEW[(i + 1) % 4][1] - SKEW[(i + 1) % 4][0] * SKEW[i][1])
                .sum::<f64>();
        let (rho, t) = (2.5, 0.2);
        let m = masses(&Quad4::new(ids, t, material()).with_density(rho), &nodes);
        let sum_x: f64 = (0..4).map(|i| m[2 * i]).sum();
        let sum_y: f64 = (0..4).map(|i| m[2 * i + 1]).sum();
        assert!((sum_x - rho * t * area).abs() < 1e-12 && (sum_y - rho * t * area).abs() < 1e-12);
    }

    #[test]
    fn prepare_is_idempotent_and_the_cache_equals_recomputation() {
        let (nodes, ids) = setup(SKEW);
        let mut el = Quad4::new(ids, 0.4, material());
        let before = assemble(&el, &nodes).k;
        assert!(el.cache().is_none());
        let view = PlaneView::new(&nodes);
        el.prepare(&view);
        el.prepare(&view);
        assert_eq!(before, assemble(&el, &nodes).k);
        assert_eq!(
            el.cache(),
            Some(point_geometry(&plane_coords(&el.nodes, &view)))
        );
    }

    #[test]
    fn invalid_geometry_is_reported_with_the_offending_corner() {
        let (nodes, ids) = setup([[0.0, 0.0], [1.0, 0.0], [0.4, 0.4], [0.0, 1.0]]);
        let view = PlaneView::new(&nodes);
        let error = Quad4::new(ids, 1.0, material())
            .validate(&view)
            .unwrap_err();
        assert!(error.contains("corner 2"), "{error}");
        let (nodes, ids) = setup(SKEW);
        assert!(Quad4::new(ids, 1.0, material())
            .validate(&PlaneView::new(&nodes))
            .is_ok());
    }

    #[test]
    fn each_gauss_point_commits_its_own_material_history() {
        let (mut nodes, ids) = setup([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let mut el = Quad4::new(ids, 1.0, PlaneMaterial::Probe { last: [0.0; 3] });
        // Pure bending-like field: u_x = k x y gives gamma_xy = k x varying across Gauss points.
        for id in ids {
            let [x, y] = nodes[id].coords;
            nodes[id].displacement = [0.01 * x * y, 0.0, 0.0];
        }
        el.commit(&PlaneView::new(&nodes));
        let lasts: Vec<[f64; 3]> = el
            .materials()
            .iter()
            .map(|m| match m {
                PlaneMaterial::Probe { last } => *last,
                _ => unreachable!(),
            })
            .collect();
        assert!(lasts.windows(2).any(|w| w[0] != w[1]), "{lasts:?}");
        // The committed strain at each point equals what the element reports there.
        let reported = el.gauss_responses(&PlaneView::new(&nodes));
        for (last, point) in lasts.iter().zip(&reported) {
            assert!((PlaneVector::from(*last) - PlaneVector::from(point.strain)).norm() < 1e-15);
        }
    }
}
