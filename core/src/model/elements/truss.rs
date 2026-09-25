use nalgebra::{SMatrix, SVector};

use super::super::{
    Material, Node, Node3, Node3Id, NodeId, SpatialDof, ELEMENT_DOF, SPATIAL_ELEMENT_DOF,
    SPATIAL_NDF, SPATIAL_NDIM,
};

/// Fixed-size spatial two-node element tangent/resistance pair types.
pub type SpatialElementMatrix = SMatrix<f64, SPATIAL_ELEMENT_DOF, SPATIAL_ELEMENT_DOF>;
pub type SpatialElementVector = SVector<f64, SPATIAL_ELEMENT_DOF>;

/// A 2-node axial truss: fixed-size element-local linear algebra (§2.4), no
/// heap allocation in the hot path. Only ever populates the translational
/// DOF entries of the (now 3-DOF-per-node) local/global matrices — the
/// rotation entries stay zero, so a node connected only to `Truss`/
/// `ZeroLength` elements needs its rotation DOF fixed by the model, or the
/// global system is singular in that DOF.
#[derive(Debug, Clone)]
pub struct Truss {
    pub node_i: NodeId,
    pub node_j: NodeId,
    pub area: f64,
    pub material: Material,
    /// Mass per unit volume. Zero (the default) means massless — existing
    /// models are unaffected unless they opt in via `with_density`.
    pub density: f64,
}

impl Truss {
    pub fn new(node_i: NodeId, node_j: NodeId, area: f64, material: Material) -> Self {
        Truss {
            node_i,
            node_j,
            area,
            material,
            density: 0.0,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    fn geometry(&self, node_i: &Node, node_j: &Node) -> (f64, f64, f64) {
        let dx = node_j.coords[0] - node_i.coords[0];
        let dy = node_j.coords[1] - node_i.coords[1];
        let length = (dx * dx + dy * dy).sqrt();
        (length, dx / length, dy / length)
    }

    /// Axial elongation from current nodal displacements, projected onto
    /// the (undeformed) bar axis — small-displacement theory. Shared by
    /// `form_tangent_and_resistance` (trial) and `commit`.
    fn strain(&self, node_i: &Node, node_j: &Node) -> f64 {
        let (length, cx, cy) = self.geometry(node_i, node_j);
        ((node_j.displacement[0] - node_i.displacement[0]) * cx
            + (node_j.displacement[1] - node_i.displacement[1]) * cy)
            / length
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        let (length, cx, cy) = self.geometry(node_i, node_j);
        let strain = self.strain(node_i, node_j);
        let (stress, tangent_modulus) = self.material.trial_stress_tangent(strain);

        // Local DOF order [ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]; rotation
        // entries (2, 5) stay zero — a truss carries no moment.
        let b = SVector::<f64, ELEMENT_DOF>::from_column_slice(&[-cx, -cy, 0.0, cx, cy, 0.0]);

        let k = (tangent_modulus * self.area / length) * (b * b.transpose());
        let resistance = (stress * self.area) * b;

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        let strain = self.strain(node_i, node_j);
        self.material = self.material.commit(strain);
    }

    /// A truss has no distinct local-axis frame to rotate into (its
    /// direction cosines are already baked into `form_tangent_and_resistance`'s
    /// `b`, unlike a beam-column's fixed `t`/current-chord transform) — its
    /// resistance vector already *is* the local nodal force, so this is
    /// that same value, recomputed fresh from the current committed strain
    /// (cheap: one material evaluation, no history to cache).
    pub(super) fn local_force(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        self.form_tangent_and_resistance(node_i, node_j).1
    }

    /// Lumped mass: half the element's total mass (`density * area * length`)
    /// at each node, split equally between that node's translational DOFs —
    /// a point mass has no directional preference. Zero rotational
    /// contribution (index 2, 5) — a truss carries no moment either.
    pub(super) fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        let (length, _cx, _cy) = self.geometry(node_i, node_j);
        let half = self.density * self.area * length / 2.0;
        SVector::<f64, ELEMENT_DOF>::from_column_slice(&[half, half, 0.0, half, half, 0.0])
    }
}

/// A small-displacement axial 3D truss.
///
/// The bar has three translational components at each end; its six rotational
/// rows/columns remain zero. A model made solely of `Truss3` therefore must
/// restrain unused rotations, exactly as a planar truss-only model restrains
/// its unused `rz` DOFs. The future `Domain3` will enforce this at solve time.
#[derive(Debug, Clone)]
pub struct Truss3 {
    pub node_i: Node3Id,
    pub node_j: Node3Id,
    pub area: f64,
    pub material: Material,
    /// Mass per unit volume; zero means massless.
    pub density: f64,
}

impl Truss3 {
    pub fn new(node_i: Node3Id, node_j: Node3Id, area: f64, material: Material) -> Self {
        Truss3 {
            node_i,
            node_j,
            area,
            material,
            density: 0.0,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// Initial length and global direction cosines. A zero-length bar is a
    /// model construction error, not a valid spatial element.
    fn geometry(&self, node_i: &Node3, node_j: &Node3) -> (f64, [f64; SPATIAL_NDIM]) {
        let dx = node_j.coords[0] - node_i.coords[0];
        let dy = node_j.coords[1] - node_i.coords[1];
        let dz = node_j.coords[2] - node_i.coords[2];
        let length = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(length > 0.0, "Truss3 endpoints must not coincide");
        (length, [dx / length, dy / length, dz / length])
    }

    fn axial_strain(&self, node_i: &Node3, node_j: &Node3) -> f64 {
        let (length, [cx, cy, cz]) = self.geometry(node_i, node_j);
        let relative = (node_j.displacement[SpatialDof::Ux as usize]
            - node_i.displacement[SpatialDof::Ux as usize])
            * cx
            + (node_j.displacement[SpatialDof::Uy as usize]
                - node_i.displacement[SpatialDof::Uy as usize])
                * cy
            + (node_j.displacement[SpatialDof::Uz as usize]
                - node_i.displacement[SpatialDof::Uz as usize])
                * cz;
        relative / length
    }

    /// Element tangent and internal resistance in global DOF order
    /// `[ux_i, uy_i, uz_i, rx_i, ry_i, rz_i, ux_j, ..., rz_j]`.
    pub fn form_tangent_and_resistance(
        &self,
        node_i: &Node3,
        node_j: &Node3,
    ) -> (SpatialElementMatrix, SpatialElementVector) {
        let (length, [cx, cy, cz]) = self.geometry(node_i, node_j);
        let strain = self.axial_strain(node_i, node_j);
        let (stress, tangent_modulus) = self.material.trial_stress_tangent(strain);

        let mut b = SpatialElementVector::zeros();
        b[SpatialDof::Ux as usize] = -cx;
        b[SpatialDof::Uy as usize] = -cy;
        b[SpatialDof::Uz as usize] = -cz;
        b[SPATIAL_NDF + SpatialDof::Ux as usize] = cx;
        b[SPATIAL_NDF + SpatialDof::Uy as usize] = cy;
        b[SPATIAL_NDF + SpatialDof::Uz as usize] = cz;

        let tangent = (tangent_modulus * self.area / length) * (b.clone() * b.transpose());
        let resistance = (stress * self.area) * b;
        (tangent, resistance)
    }

    pub fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        self.material = self.material.commit(self.axial_strain(node_i, node_j));
    }

    /// See `Truss::local_force`'s doc comment — same reasoning, no separate
    /// local frame to rotate into.
    pub fn local_force(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        self.form_tangent_and_resistance(node_i, node_j).1
    }

    /// Diagonal lumped mass: half of the member mass at each node, applied to
    /// each of its three translational DOFs. Rotational inertia remains zero
    /// until a spatial beam/mass formulation supplies it explicitly.
    pub fn form_mass(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        let (length, _) = self.geometry(node_i, node_j);
        let half = self.density * self.area * length / 2.0;
        let mut mass = SpatialElementVector::zeros();
        for dof in [SpatialDof::Ux, SpatialDof::Uy, SpatialDof::Uz] {
            mass[dof as usize] = half;
            mass[SPATIAL_NDF + dof as usize] = half;
        }
        mass
    }
}

#[cfg(test)]
mod spatial_tests {
    use super::*;

    #[test]
    fn skew_truss_matches_direction_cosine_stiffness_and_resistance() {
        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let mut node_j = Node3::new([2.0, 3.0, 6.0]);
        node_j.displacement[SpatialDof::Ux as usize] = 0.004;
        node_j.displacement[SpatialDof::Uy as usize] = 0.006;
        node_j.displacement[SpatialDof::Uz as usize] = 0.012;

        let truss = Truss3::new(
            Node3Id::default(),
            Node3Id::default(),
            2.0,
            Material::Elastic { e: 30_000.0 },
        );
        let (k, r) = truss.form_tangent_and_resistance(&node_i, &node_j);

        assert!((r[6] - 240.0 / 7.0).abs() < 1e-12);
        assert!((r[7] - 360.0 / 7.0).abs() < 1e-12);
        assert!((r[8] - 720.0 / 7.0).abs() < 1e-12);
        assert_eq!(r[3], 0.0);
        assert_eq!(r[9], 0.0);

        let mut u = SpatialElementVector::zeros();
        u[6] = 0.004;
        u[7] = 0.006;
        u[8] = 0.012;
        assert!((k * u - r).norm() < 1e-12);
    }

    #[test]
    fn mass_is_lumped_over_three_translations_only() {
        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let node_j = Node3::new([0.0, 0.0, 4.0]);
        let truss = Truss3::new(
            Node3Id::default(),
            Node3Id::default(),
            3.0,
            Material::Elastic { e: 1.0 },
        )
        .with_density(2.0);
        let mass = truss.form_mass(&node_i, &node_j);

        for i in [0, 1, 2, 6, 7, 8] {
            assert_eq!(mass[i], 12.0);
        }
        for i in [3, 4, 5, 9, 10, 11] {
            assert_eq!(mass[i], 0.0);
        }
    }
}
