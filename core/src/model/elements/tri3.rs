//! `Tri3`: the 3-node constant-strain triangle. One Gauss point, 6 DOFs
//! `[ux, uy]` per node, plane stress or strain through its `PlaneMaterial`.

use nalgebra::{SMatrix, SVector};

use super::plane_common::{
    b_matrix, check_section, load_vector, plane_coords, plane_displacement, plane_dofs, PlaneView,
};
use super::{DofMask, ElementForce, GaussResponse, TangentSink, VectorSink};
use crate::model::continuum::{physical_derivatives, tri3_dshape, validate_triangle, Jacobian};
use crate::model::{ElementLoad, NodeId, PlaneMaterial, PlaneVector};

#[derive(Debug, Clone)]
pub struct Tri3 {
    pub nodes: [NodeId; 3],
    pub thickness: f64,
    pub density: f64,
    /// The one Gauss point's material (an independent copy of the prototype).
    material: PlaneMaterial,
    cache: Option<Geometry>,
}

/// `B` and the area, which depend only on the (fixed) coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Geometry {
    b: SMatrix<f64, 3, 6>,
    area: f64,
}

impl Tri3 {
    pub fn new(
        n1: NodeId,
        n2: NodeId,
        n3: NodeId,
        thickness: f64,
        material: PlaneMaterial,
    ) -> Self {
        Tri3 {
            nodes: [n1, n2, n3],
            thickness,
            density: 0.0,
            material,
            cache: None,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    pub fn material(&self) -> &PlaneMaterial {
        &self.material
    }

    fn compute_geometry(&self, view: &PlaneView<'_>) -> Geometry {
        let coords = plane_coords(&self.nodes, view);
        let dn = tri3_dshape();
        let jacobian = Jacobian::new(&coords, &dn);
        Geometry {
            b: b_matrix(&physical_derivatives(&jacobian, &dn)),
            area: 0.5 * jacobian.det,
        }
    }

    fn geometry(&self, view: &PlaneView<'_>) -> Geometry {
        self.cache.unwrap_or_else(|| self.compute_geometry(view))
    }

    pub(super) fn dof_mask(&self) -> DofMask {
        DofMask::none().with(0).with(1)
    }

    pub(super) fn validate(&self, view: &PlaneView<'_>) -> Result<(), &'static str> {
        check_section(self.thickness, self.density, &self.material)?;
        validate_triangle(&plane_coords(&self.nodes, view))
    }

    pub(super) fn prepare(&mut self, view: &PlaneView<'_>) {
        self.cache = Some(self.compute_geometry(view));
    }

    fn strain(&self, g: &Geometry, view: &PlaneView<'_>) -> PlaneVector {
        g.b * plane_displacement::<3, 6>(&self.nodes, view)
    }

    pub(super) fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        view: &PlaneView<'_>,
        sink: &mut S,
    ) {
        let g = self.geometry(view);
        let (stress, d) = self.material.trial_stress_tangent(&self.strain(&g, view));
        let scale = self.thickness * g.area;
        let k = scale * g.b.transpose() * d * g.b;
        let r = scale * g.b.transpose() * stress;
        sink.add(&plane_dofs::<3, 6>(&self.nodes), &k, &r);
    }

    /// Exact lumped mass `rho t A / 3` on each translational DOF of each node.
    pub(super) fn assemble_mass<S: VectorSink<NodeId>>(&self, view: &PlaneView<'_>, sink: &mut S) {
        let m = self.density * self.thickness * self.geometry(view).area / 3.0;
        sink.add(
            &plane_dofs::<3, 6>(&self.nodes),
            &SVector::<f64, 6>::repeat(m),
        );
    }

    /// Consistent nodal forces of `load` (global DOF order): body force
    /// `b t A / 3` per node, plus edge tractions and pressures.
    pub(super) fn load_vector(&self, view: &PlaneView<'_>, load: &ElementLoad) -> SVector<f64, 6> {
        let share = self.thickness * self.geometry(view).area / 3.0;
        load_vector(
            &plane_coords(&self.nodes, view),
            &[share; 3],
            self.thickness,
            load,
        )
    }

    pub(super) fn assemble_load<S: VectorSink<NodeId>>(
        &self,
        view: &PlaneView<'_>,
        load: &ElementLoad,
        sink: &mut S,
    ) {
        sink.add(
            &plane_dofs::<3, 6>(&self.nodes),
            &self.load_vector(view, load),
        );
    }

    pub(super) fn commit(&mut self, view: &PlaneView<'_>) {
        let g = self.geometry(view);
        self.material = self.material.commit(&self.strain(&g, view));
    }

    /// The nodal resisting forces in global DOF order.
    pub(super) fn local_force(&self, view: &PlaneView<'_>) -> ElementForce {
        let g = self.geometry(view);
        let (stress, _) = self.material.trial_stress_tangent(&self.strain(&g, view));
        (self.thickness * g.area * g.b.transpose() * stress).into()
    }

    pub(super) fn gauss_responses(&self, view: &PlaneView<'_>) -> Vec<GaussResponse> {
        let g = self.geometry(view);
        let strain = self.strain(&g, view);
        let (stress, _) = self.material.trial_stress_tangent(&strain);
        vec![GaussResponse {
            strain: strain.into(),
            stress: stress.into(),
        }]
    }

    #[cfg(test)]
    pub(crate) fn cached(&self) -> bool {
        self.cache.is_some()
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::SymmetricEigen;
    use slotmap::SlotMap;

    use super::*;
    use crate::model::{DofRef, ElementLoad, Node};

    struct Collect {
        k: SMatrix<f64, 6, 6>,
        r: SVector<f64, 6>,
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

    struct Mass(SVector<f64, 6>);

    impl VectorSink<NodeId> for Mass {
        fn add<const N: usize>(&mut self, _dofs: &[DofRef<NodeId>; N], v: &SVector<f64, N>) {
            self.0 = SVector::from_fn(|i, _| v[i]);
        }
    }

    fn setup(coords: [[f64; 2]; 3]) -> (SlotMap<NodeId, Node>, [NodeId; 3]) {
        let mut nodes = SlotMap::with_key();
        let ids = coords.map(|c| nodes.insert(Node::new(c)));
        (nodes, ids)
    }

    fn assemble(el: &Tri3, nodes: &SlotMap<NodeId, Node>) -> Collect {
        let mut sink = Collect {
            k: SMatrix::zeros(),
            r: SVector::zeros(),
        };
        el.assemble_tangent(&PlaneView::new(nodes), &mut sink);
        sink
    }

    fn unit_triangle(e: f64, nu: f64, t: f64) -> (Tri3, SlotMap<NodeId, Node>) {
        let (nodes, [a, b, c]) = setup([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
        (
            Tri3::new(a, b, c, t, PlaneMaterial::plane_stress(e, nu).unwrap()),
            nodes,
        )
    }

    #[test]
    fn stiffness_matches_the_closed_form_for_the_unit_right_triangle() {
        // E = 1, nu = 0, t = 1: K = A B^T D B with D = diag(1, 1, 1/2) and A = 1/2; rows sum to zero
        // (rigid translation) and e.g. K[ux1, ux1] = A * 1 = 1/2, K[ux0, ux0] = A * (1 + 1/2) = 3/4.
        let (el, nodes) = unit_triangle(1.0, 0.0, 1.0);
        let k = assemble(&el, &nodes).k;
        #[rustfmt::skip]
        let expected = SMatrix::<f64, 6, 6>::from_row_slice(&[
             0.75,  0.25, -0.5,  -0.25, -0.25,  0.0,
             0.25,  0.75,  0.0,  -0.25, -0.25, -0.5,
            -0.5,   0.0,   0.5,   0.0,   0.0,   0.0,
            -0.25, -0.25,  0.0,   0.25,  0.25,  0.0,
            -0.25, -0.25,  0.0,   0.25,  0.25,  0.0,
             0.0,  -0.5,   0.0,   0.0,   0.0,   0.5,
        ]);
        assert!((k - expected).abs().max() < 1e-14, "{k}");
    }

    #[test]
    fn stiffness_is_symmetric_with_exactly_three_rigid_body_modes() {
        let (nodes, [a, b, c]) = setup([[0.3, 0.1], [2.1, 0.4], [0.9, 1.7]]);
        let el = Tri3::new(
            a,
            b,
            c,
            0.4,
            PlaneMaterial::plane_strain(30e3, 0.2).unwrap(),
        );
        let k = assemble(&el, &nodes).k;
        assert!((k - k.transpose()).abs().max() < 1e-9 * k.abs().max());
        let eigen = SymmetricEigen::new(k).eigenvalues;
        let largest = eigen.max();
        let zero = eigen.iter().filter(|&&l| l.abs() < 1e-9 * largest).count();
        assert_eq!(zero, 3, "{eigen}");
        assert!(
            eigen.iter().all(|&l| l > -1e-9 * largest),
            "not positive semidefinite: {eigen}"
        );
    }

    #[test]
    fn a_rigid_motion_is_stress_free() {
        let (mut nodes, ids) = setup([[0.3, 0.1], [2.1, 0.4], [0.9, 1.7]]);
        let el = Tri3::new(
            ids[0],
            ids[1],
            ids[2],
            0.4,
            PlaneMaterial::plane_stress(30e3, 0.3).unwrap(),
        );
        // Small rotation about the origin plus a translation: u = t + theta * [-y, x].
        let (t, theta) = ([0.01, -0.02], 1e-5);
        for id in ids {
            let [x, y] = nodes[id].coords;
            nodes[id].displacement = [t[0] - theta * y, t[1] + theta * x, 0.0];
        }
        let r = assemble(&el, &nodes).r;
        assert!(r.abs().max() < 1e-9, "{r}");
    }

    #[test]
    fn lumped_mass_sums_to_rho_t_a_in_each_direction() {
        let (nodes, [a, b, c]) = setup([[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]]);
        let (rho, t) = (2.5, 0.2);
        let el =
            Tri3::new(a, b, c, t, PlaneMaterial::plane_stress(1.0, 0.0).unwrap()).with_density(rho);
        let mut sink = Mass(SVector::zeros());
        el.assemble_mass(&PlaneView::new(&nodes), &mut sink);
        let area = 6.0;
        let (sum_x, sum_y): (f64, f64) = (
            (0..3).map(|i| sink.0[2 * i]).sum(),
            (0..3).map(|i| sink.0[2 * i + 1]).sum(),
        );
        assert!((sum_x - rho * t * area).abs() < 1e-12 && (sum_y - rho * t * area).abs() < 1e-12);
        assert!(sink
            .0
            .iter()
            .all(|&m| (m - rho * t * area / 3.0).abs() < 1e-12));
    }

    #[test]
    fn prepare_is_idempotent_and_the_cache_equals_recomputation() {
        let (nodes, [a, b, c]) = setup([[0.3, 0.1], [2.1, 0.4], [0.9, 1.7]]);
        let mut el = Tri3::new(
            a,
            b,
            c,
            0.4,
            PlaneMaterial::plane_stress(30e3, 0.3).unwrap(),
        );
        let before = assemble(&el, &nodes).k;
        assert!(!el.cached());
        let view = PlaneView::new(&nodes);
        el.prepare(&view);
        el.prepare(&view);
        assert!(el.cached());
        assert_eq!(before, assemble(&el, &nodes).k);
        assert_eq!(el.cache, Some(el.compute_geometry(&view)));
    }

    #[test]
    fn invalid_geometry_and_section_are_reported() {
        let (nodes, [a, b, c]) = setup([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
        let view = PlaneView::new(&nodes);
        let material = PlaneMaterial::plane_stress(1.0, 0.0).unwrap();
        assert!(Tri3::new(a, b, c, 1.0, material.clone())
            .validate(&view)
            .is_ok());
        assert!(Tri3::new(a, c, b, 1.0, material.clone())
            .validate(&view)
            .unwrap_err()
            .contains("clockwise"));
        assert!(Tri3::new(a, b, c, 0.0, material.clone())
            .validate(&view)
            .is_err());
        assert!(Tri3::new(a, b, c, 1.0, material.clone())
            .with_density(-1.0)
            .validate(&view)
            .is_err());
        let bad = PlaneMaterial::PlaneStress { e: -1.0, nu: 0.0 };
        assert!(Tri3::new(a, b, c, 1.0, bad).validate(&view).is_err());
    }

    #[test]
    fn gauss_response_reports_strain_and_the_stress_it_implies() {
        let (mut nodes, ids) = setup([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
        let material = PlaneMaterial::plane_stress(100.0, 0.25).unwrap();
        let el = Tri3::new(ids[0], ids[1], ids[2], 1.0, material.clone());
        // u = (a x + b y, c x + d y): constant strain [a, d, b + c].
        let (a, b, c, d) = (1e-3, 2e-3, -3e-4, 5e-4);
        for id in ids {
            let [x, y] = nodes[id].coords;
            nodes[id].displacement = [a * x + b * y, c * x + d * y, 0.0];
        }
        let response = el.gauss_responses(&PlaneView::new(&nodes));
        assert_eq!(response.len(), 1);
        let strain = PlaneVector::new(a, d, b + c);
        assert!((PlaneVector::from(response[0].strain) - strain).norm() < 1e-15);
        assert!(
            (PlaneVector::from(response[0].stress) - material.d_matrix() * strain).norm() < 1e-12
        );
    }

    #[test]
    fn body_force_and_edge_loads_total_their_resultants() {
        let (nodes, [a, b, c]) = setup([[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]]);
        let el = Tri3::new(a, b, c, 0.2, PlaneMaterial::plane_stress(1.0, 0.0).unwrap());
        let view = PlaneView::new(&nodes);
        let sums = |v: SVector<f64, 6>| {
            [
                (0..3).map(|i| v[2 * i]).sum::<f64>(),
                (0..3).map(|i| v[2 * i + 1]).sum::<f64>(),
            ]
        };
        let body = sums(el.load_vector(&view, &ElementLoad::body(1.5, -4.0)));
        assert!(
            (body[0] - 1.5 * 0.2 * 6.0).abs() < 1e-14 && (body[1] + 4.0 * 0.2 * 6.0).abs() < 1e-14
        );
        // Edge 0 is the base (length 4), edge 1 runs (4,0) -> (1,3) (length sqrt(18)).
        let base = sums(el.load_vector(&view, &ElementLoad::edge_traction(0, 2.0, 3.0)));
        assert!(
            (base[0] - 2.0 * 0.2 * 4.0).abs() < 1e-13 && (base[1] - 3.0 * 0.2 * 4.0).abs() < 1e-13
        );
        // Pressure on the base pushes up (+y, the inward normal of a counter-clockwise triangle).
        let pressure = sums(el.load_vector(&view, &ElementLoad::edge_pressure(0, 5.0)));
        assert!(pressure[0].abs() < 1e-13 && (pressure[1] - 5.0 * 0.2 * 4.0).abs() < 1e-13);
        // Edge 1's pressure resultant is its length times the unit inward normal, i.e. (-3, -3) * p t / ... .
        let slanted = sums(el.load_vector(&view, &ElementLoad::edge_pressure(1, 1.0)));
        // Tangent (-3, 3), inward (left) normal (-3, -3)/|.|; resultant = p t |edge| n = p t (-3, -3).
        assert!((slanted[0] + 0.6).abs() < 1e-13 && (slanted[1] + 0.6).abs() < 1e-13);
    }
}
