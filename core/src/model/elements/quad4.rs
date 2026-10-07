//! `Quad4`: the 4-node bilinear isoparametric quadrilateral, 2x2 Gauss
//! integration, `[ux, uy]` per node (8 DOFs), plane stress or strain through
//! its `PlaneMaterial`. Each Gauss point owns an independent copy of the
//! material and commits it individually.
//!
//! # Formulations
//!
//! [`Quad4Formulation::Full`] is the plain element. It shear-locks in bending.
//!
//! [`Quad4Formulation::Enhanced`] is the Wilson-Taylor incompatible-modes
//! element (QM6/Q6). For each displacement component it adds the two shapes
//! `1 - xi^2` and `1 - eta^2`, four internal parameters `alpha` per element,
//! discontinuous across edges and therefore condensed out locally. The
//! enhanced strain is `G~(xi, eta) alpha` with
//! `G~ = (detJ0 / detJ(xi, eta)) * [J0^{-1} dM/dxi]` expanded into strain
//! rows; `J0` is the Jacobian at the element center, and the `detJ0 / detJ`
//! factor (Taylor-Beresford-Wilson) is what preserves the constant-strain
//! patch test on distorted elements. The element stiffness is
//! `K_uu - K_ua K_aa^{-1} K_au`, the resistance is that stiffness times the
//! nodal displacements (the element is linear), and strains and stresses at
//! the Gauss points are recovered as `B u + G~ alpha` and `D` times that, with
//! `alpha = -K_aa^{-1} K_au u`.
//!
//! Scope: `Enhanced` accepts only linear `PlaneMaterial`s (for a nonlinear one
//! `alpha` would be element internal state needing its own local Newton
//! iteration). It is designed to cure shear (bending) locking. It is **not** a
//! mixed or B-bar formulation and carries no volumetric-locking guarantee for
//! nearly incompressible plane strain (`nu -> 0.5`); the plain element locks
//! completely there, and in the tests (Cook's membrane, a pressurized thick
//! cylinder at `nu = 0.4999`) the enhanced element was observed to do much
//! better, but that is an observation about those problems, not a property of
//! the formulation. Incompatible-mode elements also remain somewhat sensitive
//! to trapezoidal distortion.

use nalgebra::{SMatrix, SVector};

use super::plane_common::{
    b_matrix, check_section, load_vector, plane_coords, plane_displacement, plane_dofs, PlaneView,
};
use super::{DofMask, ElementForce, GaussResponse, TangentSink, VectorSink};
use crate::model::continuum::{
    physical_derivatives, quad4_dshape, quad4_shape, validate_quad, Jacobian, QUAD4_GAUSS_2X2,
};
use crate::model::{ElementLoad, NodeId, PlaneMaterial, PlaneMatrix, PlaneVector};

/// Which displacement field the element uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Quad4Formulation {
    /// Plain bilinear isoparametric.
    #[default]
    Full,
    /// Wilson-Taylor incompatible modes, condensed (linear materials only).
    Enhanced,
}

#[derive(Debug, Clone)]
pub struct Quad4 {
    pub nodes: [NodeId; 4],
    pub thickness: f64,
    pub density: f64,
    formulation: Quad4Formulation,
    materials: [PlaneMaterial; 4],
    cache: Option<Cache>,
}

/// Geometry-only quantities, computed once by `prepare`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Cache {
    points: [PointGeometry; 4],
    enhanced: Option<Enhanced>,
}

/// The condensed enhanced element: `k = K_uu - K_ua K_aa^{-1} K_au` and
/// `alpha = recovery * u`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Enhanced {
    k: SMatrix<f64, 8, 8>,
    recovery: SMatrix<f64, 4, 8>,
    /// `G~` at each Gauss point.
    g: [SMatrix<f64, 3, 4>; 4],
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

/// `G~` at each Gauss point (see the module doc comment).
fn enhanced_strain_matrices(coords: &[[f64; 2]; 4]) -> [SMatrix<f64, 3, 4>; 4] {
    let j0 = Jacobian::new(coords, &quad4_dshape(0.0, 0.0));
    std::array::from_fn(|g| {
        let (xi, eta, _) = QUAD4_GAUSS_2X2[g];
        let det = Jacobian::new(coords, &quad4_dshape(xi, eta)).det;
        // Reference derivatives of the modes 1 - xi^2 and 1 - eta^2, mapped by J0.
        let dm = [[-2.0 * xi, 0.0], [0.0, -2.0 * eta]];
        let gm = physical_derivatives(&j0, &dm);
        let scale = j0.det / det;
        let (g1, g2) = (gm[0], gm[1]);
        scale
            * SMatrix::<f64, 3, 4>::new(
                g1[0], g2[0], 0.0, 0.0, //
                0.0, 0.0, g1[1], g2[1], //
                g1[1], g2[1], g1[0], g2[0],
            )
    })
}

/// Condenses the incompatible modes out of the element.
fn condense(
    points: &[PointGeometry; 4],
    d: &[PlaneMatrix; 4],
    thickness: f64,
    coords: &[[f64; 2]; 4],
) -> Result<Enhanced, &'static str> {
    let g = enhanced_strain_matrices(coords);
    let mut k_uu = SMatrix::<f64, 8, 8>::zeros();
    let mut k_ua = SMatrix::<f64, 8, 4>::zeros();
    let mut k_aa = SMatrix::<f64, 4, 4>::zeros();
    for i in 0..4 {
        let scale = thickness * points[i].weight;
        let (b, gi) = (&points[i].b, &g[i]);
        k_uu += scale * b.transpose() * d[i] * b;
        k_ua += scale * b.transpose() * d[i] * gi;
        k_aa += scale * gi.transpose() * d[i] * gi;
    }
    let k_aa_inverse = k_aa
        .try_inverse()
        .filter(|inverse| inverse.iter().all(|v| v.is_finite()))
        .ok_or("enhanced Quad4 mode matrix is singular")?;
    let recovery = -k_aa_inverse * k_ua.transpose();
    Ok(Enhanced {
        k: k_uu + k_ua * recovery,
        recovery,
        g,
    })
}

impl Quad4 {
    pub fn new(nodes: [NodeId; 4], thickness: f64, material: PlaneMaterial) -> Self {
        Quad4 {
            nodes,
            thickness,
            density: 0.0,
            formulation: Quad4Formulation::Full,
            materials: std::array::from_fn(|_| material.clone()),
            cache: None,
        }
    }

    /// Selects the formulation; `Enhanced` requires a linear material
    /// (checked by `validate`).
    pub fn with_formulation(mut self, formulation: Quad4Formulation) -> Self {
        self.formulation = formulation;
        self
    }

    pub fn formulation(&self) -> Quad4Formulation {
        self.formulation
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// The Gauss points' materials, in `QUAD4_GAUSS_2X2` order.
    pub fn materials(&self) -> &[PlaneMaterial; 4] {
        &self.materials
    }

    fn compute_cache(&self, view: &PlaneView<'_>) -> Result<Cache, &'static str> {
        let coords = plane_coords(&self.nodes, view);
        let points = point_geometry(&coords);
        let enhanced = match self.formulation {
            Quad4Formulation::Full => None,
            Quad4Formulation::Enhanced => {
                let d: [PlaneMatrix; 4] = std::array::from_fn(|g| self.materials[g].d_matrix());
                Some(condense(&points, &d, self.thickness, &coords)?)
            }
        };
        Ok(Cache { points, enhanced })
    }

    /// The cache, or a fresh computation (an unvalidated enhanced element with a
    /// singular mode matrix gets zero stiffness rather than a panic here;
    /// `validate` is where that is reported).
    fn cache_or_compute(&self, view: &PlaneView<'_>) -> Cache {
        self.cache.unwrap_or_else(|| {
            self.compute_cache(view).unwrap_or_else(|_| Cache {
                points: point_geometry(&plane_coords(&self.nodes, view)),
                enhanced: Some(Enhanced {
                    k: SMatrix::zeros(),
                    recovery: SMatrix::zeros(),
                    g: [SMatrix::zeros(); 4],
                }),
            })
        })
    }

    pub(super) fn dof_mask(&self) -> DofMask {
        DofMask::none().with(0).with(1)
    }

    pub(super) fn validate(&self, view: &PlaneView<'_>) -> Result<(), &'static str> {
        check_section(self.thickness, self.density, &self.materials[0])?;
        validate_quad(&plane_coords(&self.nodes, view))?;
        if self.formulation == Quad4Formulation::Enhanced {
            if !self.materials.iter().all(PlaneMaterial::is_linear) {
                return Err("enhanced Quad4 requires a linear plane material");
            }
            self.compute_cache(view)?;
        }
        Ok(())
    }

    pub(super) fn prepare(&mut self, view: &PlaneView<'_>) {
        self.cache = self.compute_cache(view).ok();
    }

    /// Strain at Gauss point `g`: `B u`, plus `G~ alpha` for the enhanced element.
    fn strains(&self, cache: &Cache, u: &SVector<f64, 8>) -> [PlaneVector; 4] {
        let alpha = cache.enhanced.as_ref().map(|e| e.recovery * u);
        std::array::from_fn(|g| {
            let mut strain = cache.points[g].b * u;
            if let (Some(e), Some(alpha)) = (&cache.enhanced, &alpha) {
                strain += e.g[g] * alpha;
            }
            strain
        })
    }

    pub(super) fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        view: &PlaneView<'_>,
        sink: &mut S,
    ) {
        let u = plane_displacement::<4, 8>(&self.nodes, view);
        let cache = self.cache_or_compute(view);
        if let Some(enhanced) = &cache.enhanced {
            sink.add(
                &plane_dofs::<4, 8>(&self.nodes),
                &enhanced.k,
                &(enhanced.k * u),
            );
            return;
        }
        let mut k = SMatrix::<f64, 8, 8>::zeros();
        let mut r = SVector::<f64, 8>::zeros();
        for (g, material) in cache.points.iter().zip(&self.materials) {
            let (stress, d) = material.trial_stress_tangent(&(g.b * u));
            let scale = self.thickness * g.weight;
            k += scale * g.b.transpose() * d * g.b;
            r += scale * g.b.transpose() * stress;
        }
        sink.add(&plane_dofs::<4, 8>(&self.nodes), &k, &r);
    }

    /// `t * integral(N_i dA)` per node: the body-force share, and (times density) the
    /// row-sum lumped mass `m_i = sum_g rho t N_i(g) det J(g) w(g)`. Exact for the
    /// consistent mass here (bilinear `N_i` times linear `det J`), and equal to
    /// `t A / 4` only for parallelograms.
    fn nodal_weights(&self, view: &PlaneView<'_>) -> [f64; 4] {
        let mut weights = [0.0; 4];
        for g in self.cache_or_compute(view).points {
            for (i, w) in weights.iter_mut().enumerate() {
                *w += self.thickness * g.n[i] * g.weight;
            }
        }
        weights
    }

    /// Row-sum lumped mass, applied to each translational DOF of node `i`.
    pub(super) fn assemble_mass<S: VectorSink<NodeId>>(&self, view: &PlaneView<'_>, sink: &mut S) {
        let weights = self.nodal_weights(view);
        let m = SVector::<f64, 8>::from_fn(|a, _| self.density * weights[a / 2]);
        sink.add(&plane_dofs::<4, 8>(&self.nodes), &m);
    }

    /// Consistent nodal forces of `load` (global DOF order), sharing the mass weights
    /// for the body force.
    pub(super) fn load_vector(&self, view: &PlaneView<'_>, load: &ElementLoad) -> SVector<f64, 8> {
        load_vector(
            &plane_coords(&self.nodes, view),
            &self.nodal_weights(view),
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
            &plane_dofs::<4, 8>(&self.nodes),
            &self.load_vector(view, load),
        );
    }

    pub(super) fn commit(&mut self, view: &PlaneView<'_>) {
        let u = plane_displacement::<4, 8>(&self.nodes, view);
        let strains = self.strains(&self.cache_or_compute(view), &u);
        for (material, strain) in self.materials.iter_mut().zip(&strains) {
            *material = material.commit(strain);
        }
    }

    pub(super) fn local_force(&self, view: &PlaneView<'_>) -> ElementForce {
        let u = plane_displacement::<4, 8>(&self.nodes, view);
        let cache = self.cache_or_compute(view);
        if let Some(enhanced) = &cache.enhanced {
            return (enhanced.k * u).into();
        }
        let mut r = SVector::<f64, 8>::zeros();
        for (g, material) in cache.points.iter().zip(&self.materials) {
            let (stress, _) = material.trial_stress_tangent(&(g.b * u));
            r += self.thickness * g.weight * g.b.transpose() * stress;
        }
        r.into()
    }

    /// Committed strain and stress at the four Gauss points; for the enhanced
    /// element these include the internal-mode contribution.
    pub(super) fn gauss_responses(&self, view: &PlaneView<'_>) -> Vec<GaussResponse> {
        let u = plane_displacement::<4, 8>(&self.nodes, view);
        let strains = self.strains(&self.cache_or_compute(view), &u);
        strains
            .iter()
            .zip(&self.materials)
            .map(|(strain, material)| {
                let (stress, _) = material.trial_stress_tangent(strain);
                GaussResponse {
                    strain: (*strain).into(),
                    stress: stress.into(),
                }
            })
            .collect()
    }

    #[cfg(test)]
    pub(super) fn cache(&self) -> Option<Cache> {
        self.cache
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::SymmetricEigen;
    use slotmap::SlotMap;

    use super::*;
    use crate::model::{DofRef, ElementLoad, Node};

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
        assert_eq!(el.cache(), Some(el.compute_cache(&view).unwrap()));
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

    fn enhanced_on_skew() -> (Quad4, SlotMap<NodeId, Node>) {
        let (mut nodes, ids) = setup(SKEW);
        let el = Quad4::new(ids, 0.4, material()).with_formulation(Quad4Formulation::Enhanced);
        // A generic displacement state, bending-like so the internal modes engage.
        for (k, id) in ids.into_iter().enumerate() {
            let [x, y] = nodes[id].coords;
            nodes[id].displacement = [1e-3 * x * y + 1e-4 * k as f64, -5e-4 * x * x, 0.0];
        }
        (el, nodes)
    }

    #[test]
    fn recovered_stress_includes_the_internal_modes_and_differs_from_d_b_u() {
        let (el, nodes) = enhanced_on_skew();
        let view = PlaneView::new(&nodes);
        let cache = el.compute_cache(&view).unwrap();
        let enhanced = cache.enhanced.unwrap();
        let u = plane_displacement::<4, 8>(&el.nodes, &view);
        let alpha = enhanced.recovery * u;
        assert!(alpha.norm() > 1e-6, "internal modes should engage");
        let d = material().d_matrix();
        let responses = el.gauss_responses(&view);
        for g in 0..4 {
            let strain = cache.points[g].b * u + enhanced.g[g] * alpha;
            let stress = d * strain;
            assert!((PlaneVector::from(responses[g].stress) - stress).norm() < 1e-9);
            assert!((PlaneVector::from(responses[g].strain) - strain).norm() < 1e-15);
            let compatible_only = d * (cache.points[g].b * u);
            assert!((stress - compatible_only).norm() > 1e-3 * stress.norm());
        }
    }

    #[test]
    fn condensation_reproduces_the_resistance_from_the_recovered_stresses() {
        let (el, nodes) = enhanced_on_skew();
        let view = PlaneView::new(&nodes);
        let mut sink = Collect {
            k: SMatrix::zeros(),
            r: SVector::zeros(),
        };
        el.assemble_tangent(&view, &mut sink);
        let cache = el.compute_cache(&view).unwrap();
        let mut from_stress = SVector::<f64, 8>::zeros();
        for (g, response) in cache.points.iter().zip(el.gauss_responses(&view)) {
            from_stress +=
                el.thickness * g.weight * g.b.transpose() * PlaneVector::from(response.stress);
        }
        assert!((from_stress - sink.r).norm() < 1e-9 * sink.r.norm());
        assert_eq!(
            sink.r,
            sink.k * plane_displacement::<4, 8>(&el.nodes, &view)
        );
        // Condensation made K = K_uu - K_ua K_aa^{-1} K_au: symmetric, with three rigid-body modes.
        assert!((sink.k - sink.k.transpose()).abs().max() < 1e-9 * sink.k.abs().max());
        let eigen = SymmetricEigen::new(sink.k).eigenvalues;
        let largest = eigen.max();
        assert_eq!(
            eigen.iter().filter(|&&l| l.abs() < 1e-9 * largest).count(),
            3,
            "{eigen}"
        );
    }

    #[test]
    fn a_nonlinear_material_is_rejected_by_the_enhanced_element_only() {
        let (nodes, ids) = setup(SKEW);
        let view = PlaneView::new(&nodes);
        let probe = PlaneMaterial::Probe { last: [0.0; 3] };
        let full = Quad4::new(ids, 1.0, probe.clone());
        assert!(full.validate(&view).is_ok());
        let enhanced = full.with_formulation(Quad4Formulation::Enhanced);
        assert_eq!(
            enhanced.validate(&view),
            Err("enhanced Quad4 requires a linear plane material")
        );
    }

    #[test]
    fn the_enhanced_stiffness_is_softer_than_the_full_one_in_bending_only() {
        let (nodes, ids) = setup([[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let full = Quad4::new(ids, 1.0, material());
        let enhanced = full.clone().with_formulation(Quad4Formulation::Enhanced);
        let (kf, ke) = (assemble(&full, &nodes).k, assemble(&enhanced, &nodes).k);
        // A bending displacement pattern (u_x = y * x-independent shear-free bending shape).
        let bend = SVector::<f64, 8>::from_row_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            + SVector::from_fn(|i, _| {
                let node = i / 2;
                let [x, y] = [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]][node];
                if i % 2 == 0 {
                    x * (y - 0.5)
                } else {
                    -0.5 * x * x
                }
            });
        assert!((bend.transpose() * ke * bend)[0] < 0.5 * (bend.transpose() * kf * bend)[0]);
        // A constant-strain pattern costs the same energy.
        let stretch = SVector::<f64, 8>::from_fn(|i, _| {
            if i % 2 == 0 {
                1e-3 * [0.0, 2.0, 2.0, 0.0][i / 2]
            } else {
                0.0
            }
        });
        let (ef, ee) = (
            (stretch.transpose() * kf * stretch)[0],
            (stretch.transpose() * ke * stretch)[0],
        );
        assert!((ef - ee).abs() < 1e-12 * ef);
    }

    fn force_sums(v: &SVector<f64, 8>) -> [f64; 2] {
        [
            (0..4).map(|i| v[2 * i]).sum(),
            (0..4).map(|i| v[2 * i + 1]).sum(),
        ]
    }

    #[test]
    fn body_force_totals_b_t_a_and_edge_loads_total_traction_times_length() {
        let (nodes, ids) = setup([[0.0, 0.0], [2.0, 0.0], [1.5, 1.0], [0.5, 1.0]]);
        let el = Quad4::new(ids, 0.4, material());
        let view = PlaneView::new(&nodes);
        let area = 1.5; // trapezoid with parallel sides 2 and 1, height 1
        let body = force_sums(&el.load_vector(&view, &ElementLoad::body(3.0, -2.0)));
        assert!(
            (body[0] - 3.0 * 0.4 * area).abs() < 1e-14
                && (body[1] + 2.0 * 0.4 * area).abs() < 1e-14
        );
        // Edge 1 joins node 1 (2, 0) to node 2 (1.5, 1): length sqrt(1.25).
        let length = 1.25_f64.sqrt();
        let traction =
            force_sums(&el.load_vector(&view, &ElementLoad::edge_traction(1, 5.0, -7.0)));
        assert!((traction[0] - 5.0 * 0.4 * length).abs() < 1e-13);
        assert!((traction[1] + 7.0 * 0.4 * length).abs() < 1e-13);
        // Only the edge's own two nodes receive it.
        let v = el.load_vector(&view, &ElementLoad::edge_traction(1, 5.0, -7.0));
        assert!(v[0] == 0.0 && v[1] == 0.0 && v[6] == 0.0 && v[7] == 0.0);
    }

    #[test]
    fn pressure_acts_along_the_inward_normal() {
        let (nodes, ids) = setup([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let el = Quad4::new(ids, 1.0, material());
        let view = PlaneView::new(&nodes);
        // Counter-clockwise square: the inward normals of edges 0..4 are +y, -x, -y, +x.
        for (edge, expected) in [[0.0, 1.0], [-1.0, 0.0], [0.0, -1.0], [1.0, 0.0]]
            .into_iter()
            .enumerate()
        {
            let sum = force_sums(&el.load_vector(&view, &ElementLoad::edge_pressure(edge, 2.0)));
            assert!(
                (sum[0] - 2.0 * expected[0]).abs() < 1e-14,
                "edge {edge}: {sum:?}"
            );
            assert!(
                (sum[1] - 2.0 * expected[1]).abs() < 1e-14,
                "edge {edge}: {sum:?}"
            );
        }
    }

    #[test]
    fn gauss_responses_equal_d_b_u_for_the_plain_element() {
        let (mut nodes, ids) = setup(SKEW);
        let mat = material();
        let el = Quad4::new(ids, 1.0, mat.clone());
        // A bilinear displacement field: strains vary across the four Gauss points.
        for id in ids {
            let [x, y] = nodes[id].coords;
            nodes[id].displacement = [1e-3 * x * y, -2e-3 * x + 5e-4 * y, 0.0];
        }
        let view = PlaneView::new(&nodes);
        let u = plane_displacement::<4, 8>(&el.nodes, &view);
        let geometry = point_geometry(&plane_coords(&el.nodes, &view));
        for (g, response) in el.gauss_responses(&view).iter().enumerate() {
            let expected = mat.d_matrix() * (geometry[g].b * u);
            assert!((PlaneVector::from(response.stress) - expected).norm() < 1e-12);
        }
    }
}
