use nalgebra::{SMatrix, SVector};

use super::super::{
    FiberSection, FiberSection3, Material, Node, Node3, Node3Id, NodeId, ELEMENT_DOF, NDF, SPATIAL_NDF,
};
use super::truss::{SpatialElementMatrix, SpatialElementVector};

/// A 2-node, zero-length connector: no geometry or integration, just direct
/// per-DOF material evaluation (§3.1) — each direction with a material
/// assigned independently relates that DOF's relative displacement between
/// the two nodes to a force along that same (global) direction. No
/// orientation vectors (unlike OpenSees' general `ZeroLength`, which can
/// evaluate materials along arbitrary local axes) — directions are the
/// global DOF axes (including rotation, since M3's DOF bump), which is all
/// the current scope needs.
#[derive(Debug, Clone)]
pub struct ZeroLength {
    pub node_i: NodeId,
    pub node_j: NodeId,
    materials: [Option<Material>; NDF],
}

impl ZeroLength {
    pub fn new(node_i: NodeId, node_j: NodeId) -> Self {
        ZeroLength {
            node_i,
            node_j,
            materials: std::array::from_fn(|_| None),
        }
    }

    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        let mut k = SMatrix::<f64, ELEMENT_DOF, ELEMENT_DOF>::zeros();
        let mut resistance = SVector::<f64, ELEMENT_DOF>::zeros();

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.trial_stress_tangent(relative);

            // Local DOF order [ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]; this
            // direction only couples node_i's and node_j's copy of the same
            // dof.
            let mut b = SVector::<f64, ELEMENT_DOF>::zeros();
            b[dof] = -1.0;
            b[NDF + dof] = 1.0;

            k += tangent_modulus * (b * b.transpose());
            resistance += force * b;
        }

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }
    }

    /// No orientation vectors, no separate local frame — see `Truss::
    /// local_force`'s doc comment for why this is just the resistance
    /// vector, recomputed fresh from each direction's current committed
    /// material state.
    pub(super) fn local_force(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        self.form_tangent_and_resistance(node_i, node_j).1
    }
}

/// `ZeroLength`'s spatial counterpart: independent per-DOF materials along
/// the six global directions `[ux, uy, uz, rx, ry, rz]` — no orientation
/// vectors, same simplification `ZeroLength` already makes relative to
/// OpenSees' general `ZeroLength` (which can evaluate materials along
/// arbitrary local axes via a user-supplied orientation). Arbitrarily
/// oriented local springs are a later, explicit transform feature (see
/// `docs/spatial-architecture.md`'s "Elements and transforms" section).
#[derive(Debug, Clone)]
pub struct ZeroLength3 {
    pub node_i: Node3Id,
    pub node_j: Node3Id,
    materials: [Option<Material>; SPATIAL_NDF],
}

impl ZeroLength3 {
    pub fn new(node_i: Node3Id, node_j: Node3Id) -> Self {
        ZeroLength3 {
            node_i,
            node_j,
            materials: std::array::from_fn(|_| None),
        }
    }

    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    pub(super) fn form_tangent_and_resistance(&self, node_i: &Node3, node_j: &Node3) -> (SpatialElementMatrix, SpatialElementVector) {
        let mut k = SpatialElementMatrix::zeros();
        let mut resistance = SpatialElementVector::zeros();

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.trial_stress_tangent(relative);

            let mut b = SpatialElementVector::zeros();
            b[dof] = -1.0;
            b[SPATIAL_NDF + dof] = 1.0;

            k += tangent_modulus * (b * b.transpose());
            resistance += force * b;
        }

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }
    }

    /// See `ZeroLength::local_force`'s doc comment.
    pub(super) fn local_force(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        self.form_tangent_and_resistance(node_i, node_j).1
    }
}

/// A 2-node, zero-length connector driven by a coupled `FiberSection`
/// instead of `ZeroLength`'s independent per-DOF materials — the axial
/// (`ux`) and flexural (`rz`) response come out of the same fiber
/// discretization, so axial force and bending moment interact through
/// shared fiber strain, the way a concentrated-plasticity "fiber hinge"
/// needs (e.g. at a beam-column end) and `ZeroLength` structurally cannot
/// produce. Modeled on OpenSees/Xara's `ZeroLengthSection`
/// (`SRC/element/Point/ZeroLengthSection.cpp`): the section's basic
/// deformation is the direct nodal DOF difference (no shape functions or
/// integration — there's no length to integrate over), and `K`/`P` are
/// `Bᵀ*k_section*B`/`Bᵀ*q` with `B` selecting `[ux, rz]`, the same pattern
/// `DispBeamColumn::form_tangent_and_resistance` uses per integration point
/// (here with a single "point" and unit scale).
///
/// `uy` — the one DOF `FiberSection` has no resultant for — can still carry
/// an independent uniaxial `Material` spring via `with_material`, exactly
/// like `ZeroLength`, additive with (and uncoupled from) the section's
/// axial-moment response.
#[derive(Debug, Clone)]
pub struct ZeroLengthSection {
    pub node_i: NodeId,
    pub node_j: NodeId,
    section: FiberSection,
    materials: [Option<Material>; NDF],
}

impl ZeroLengthSection {
    pub fn new(node_i: NodeId, node_j: NodeId, section: FiberSection) -> Self {
        ZeroLengthSection {
            node_i,
            node_j,
            section,
            materials: std::array::from_fn(|_| None),
        }
    }

    /// Independent spring for a DOF the section doesn't drive — intended
    /// for `dof = 1` (`uy`, shear), since the section already owns `ux`
    /// (dof 0) and `rz` (dof 2).
    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        let eps0 = node_j.displacement[0] - node_i.displacement[0];
        let kappa = node_j.displacement[2] - node_i.displacement[2];
        let (n, m, k_section) = self.section.trial(eps0, kappa);

        let b_eps0 = SVector::<f64, ELEMENT_DOF>::from_column_slice(&[-1.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
        let b_kappa = SVector::<f64, ELEMENT_DOF>::from_column_slice(&[0.0, 0.0, -1.0, 0.0, 0.0, 1.0]);

        let mut k = k_section[0][0] * (b_eps0 * b_eps0.transpose())
            + k_section[0][1] * (b_eps0 * b_kappa.transpose())
            + k_section[1][0] * (b_kappa * b_eps0.transpose())
            + k_section[1][1] * (b_kappa * b_kappa.transpose());
        let mut resistance = b_eps0 * n + b_kappa * m;

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.trial_stress_tangent(relative);

            let mut b = SVector::<f64, ELEMENT_DOF>::zeros();
            b[dof] = -1.0;
            b[NDF + dof] = 1.0;

            k += tangent_modulus * (b * b.transpose());
            resistance += force * b;
        }

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        let eps0 = node_j.displacement[0] - node_i.displacement[0];
        let kappa = node_j.displacement[2] - node_i.displacement[2];
        self.section.commit(eps0, kappa);

        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }
    }

    /// See `ZeroLength::local_force`'s doc comment.
    pub(super) fn local_force(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        self.form_tangent_and_resistance(node_i, node_j).1
    }
}

/// `ZeroLengthSection`'s spatial counterpart: coupled axial-biaxial-moment
/// response (`ux`, `ry`, `rz`) from a `FiberSection3`, plus optional
/// independent springs (via `with_material`) for `uy`, `uz` (shear) and
/// `rx` (torsion) — the three DOFs the section has no resultant for.
/// `FiberSection3` deliberately excludes torsion from its fiber loop (see
/// its doc comment), and there's no meaningful `G*J/length` term at zero
/// length the way `DispBeamColumn3` has, so an explicit `rx` material is
/// the only way to give this element torsional stiffness at all.
#[derive(Debug, Clone)]
pub struct ZeroLengthSection3 {
    pub node_i: Node3Id,
    pub node_j: Node3Id,
    section: FiberSection3,
    materials: [Option<Material>; SPATIAL_NDF],
}

impl ZeroLengthSection3 {
    pub fn new(node_i: Node3Id, node_j: Node3Id, section: FiberSection3) -> Self {
        ZeroLengthSection3 {
            node_i,
            node_j,
            section,
            materials: std::array::from_fn(|_| None),
        }
    }

    /// Independent spring for a DOF the section doesn't drive — `uy`/`uz`
    /// (shear, dofs 1/2) or `rx` (torsion, dof 3) — since the section
    /// already owns `ux` (dof 0), `ry` (dof 4) and `rz` (dof 5).
    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    pub(super) fn form_tangent_and_resistance(&self, node_i: &Node3, node_j: &Node3) -> (SpatialElementMatrix, SpatialElementVector) {
        let eps0 = node_j.displacement[0] - node_i.displacement[0];
        let kappa_z = node_j.displacement[5] - node_i.displacement[5];
        let kappa_y = node_j.displacement[4] - node_i.displacement[4];
        let (n, mz, my, k_section) = self.section.trial(eps0, kappa_z, kappa_y);

        let mut b_eps0 = SpatialElementVector::zeros();
        b_eps0[0] = -1.0;
        b_eps0[SPATIAL_NDF] = 1.0;

        let mut b_kappa_z = SpatialElementVector::zeros();
        b_kappa_z[5] = -1.0;
        b_kappa_z[SPATIAL_NDF + 5] = 1.0;

        let mut b_kappa_y = SpatialElementVector::zeros();
        b_kappa_y[4] = -1.0;
        b_kappa_y[SPATIAL_NDF + 4] = 1.0;

        let mut k = SpatialElementMatrix::zeros();
        let mut resistance = b_eps0 * n + b_kappa_z * mz + b_kappa_y * my;

        let b = [&b_eps0, &b_kappa_z, &b_kappa_y];
        for (bi, row_i) in b.iter().enumerate() {
            for (bj, row_j) in b.iter().enumerate() {
                k += k_section[bi][bj] * (*row_i * row_j.transpose());
            }
        }

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.trial_stress_tangent(relative);

            let mut bm = SpatialElementVector::zeros();
            bm[dof] = -1.0;
            bm[SPATIAL_NDF + dof] = 1.0;

            k += tangent_modulus * (bm * bm.transpose());
            resistance += force * bm;
        }

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        let eps0 = node_j.displacement[0] - node_i.displacement[0];
        let kappa_z = node_j.displacement[5] - node_i.displacement[5];
        let kappa_y = node_j.displacement[4] - node_i.displacement[4];
        self.section.commit(eps0, kappa_z, kappa_y);

        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }
    }

    /// See `ZeroLength::local_force`'s doc comment.
    pub(super) fn local_force(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        self.form_tangent_and_resistance(node_i, node_j).1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Fiber, Fiber3, Material, Node};

    /// Proves the `Material` enum's dispatch generalizes to `ZeroLength`
    /// across every M2 variant and regime (elastic, plastic plateau, open
    /// gap, closed gap, tension, compression) by driving the element
    /// directly off manually-set node displacements — independent of
    /// `Analysis`/`Algorithm`, which for M2 is still linear-only (no
    /// Newton iteration until M4) and so can't itself resolve equilibrium
    /// across a material's nonlinear regimes in one step.
    #[test]
    fn zero_length_dispatches_across_material_regimes() {
        let node_i = Node::new([0.0, 0.0]);

        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = -0.02; // compression, past both EPP yield and Gap closure

        // `resistance = force * b` with `b[dof_i] = -1`, so node_i's
        // component of the resistance vector is `-force`.
        let epp = ZeroLength::new(NodeId::default(), NodeId::default())
            .with_material(0, Material::elastic_pp(100.0, 0.01));
        let (_, r) = epp.form_tangent_and_resistance(&node_i, &node_j);
        assert_eq!(r[0], 100.0 * 0.01, "should be clamped to yield force");

        let gap = ZeroLength::new(NodeId::default(), NodeId::default())
            .with_material(0, Material::Gap { e: 100.0, gap: 0.01 });
        let (_, r) = gap.form_tangent_and_resistance(&node_i, &node_j);
        assert_eq!(r[0], 100.0 * (0.02 - 0.01), "gap engaged past closure");

        let ent = ZeroLength::new(NodeId::default(), NodeId::default())
            .with_material(0, Material::Ent { e: 100.0 });
        let (_, r) = ent.form_tangent_and_resistance(&node_i, &node_j);
        assert_eq!(r[0], 100.0 * 0.02, "ENT carries full compression");

        let mut node_tension = Node::new([0.0, 0.0]);
        node_tension.displacement[0] = 0.02; // tension
        let (k, r) = ent.form_tangent_and_resistance(&node_i, &node_tension);
        assert_eq!(r[0], 0.0, "ENT carries no tension");
        assert_eq!(k[(0, 0)], 0.0, "no tangent stiffness while open");
    }

    /// `ZeroLength3`'s spatial counterpart, exercising a rotational DOF
    /// (index 3, `SpatialDof::Rx`) in addition to a translational one — the
    /// planar `ZeroLength` test above only ever touches dof 0.
    #[test]
    fn zero_length3_dispatches_independently_per_global_direction() {
        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let mut node_j = Node3::new([0.0, 0.0, 0.0]);
        node_j.displacement[0] = 0.01; // Ux
        node_j.displacement[3] = -0.02; // Rx

        let zl = ZeroLength3::new(Node3Id::default(), Node3Id::default())
            .with_material(0, Material::Elastic { e: 100.0 })
            .with_material(3, Material::elastic_pp(50.0, 0.01));

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);
        // `resistance = force * b` with `b[dof_i] = -1`, so node_i's
        // component of the resistance vector is `-force` (see the planar
        // `ZeroLength` test's doc comment).
        assert_eq!(r[0], -(100.0 * 0.01), "Ux: linear elastic");
        assert_eq!(r[SPATIAL_NDF], 100.0 * 0.01);
        assert_eq!(r[3], 50.0 * 0.01, "Rx: clamped to yield force (sign follows relative disp)");
        assert_eq!(r[SPATIAL_NDF + 3], -(50.0 * 0.01));
        assert_eq!(k[(0, 0)], 100.0, "Ux tangent unaffected by Rx material");
        assert_eq!(k[(3, 3)], 0.0, "Rx tangent is zero past yield");
        for dof in [1, 2, 4, 5] {
            assert_eq!(r[dof], 0.0, "unassigned directions carry no force");
        }
    }

    /// `ZeroLengthSection`: two symmetric elastic fibers reproduce the
    /// closed-form `EA`/`EI` (already validated directly against
    /// `FiberSection::trial` in `fiber_section.rs`'s
    /// `two_symmetric_elastic_fibers_reproduce_ea_ei_exactly`), now routed
    /// through the element's `B`-matrix mapping onto `[ux, rz]`.
    #[test]
    fn zero_length_section_reproduces_symmetric_ea_ei() {
        let (e, area, iz): (f64, f64, f64) = (30000.0, 2.0, 1000.0);
        let h = (iz / area).sqrt();
        let section = FiberSection::new(vec![
            Fiber::new(h, area / 2.0, Material::Elastic { e }),
            Fiber::new(-h, area / 2.0, Material::Elastic { e }),
        ]);

        let node_i = Node::new([0.0, 0.0]);
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 0.001; // ux -> eps0
        node_j.displacement[2] = 0.0002; // rz -> kappa

        let zl = ZeroLengthSection::new(NodeId::default(), NodeId::default(), section);
        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((k[(0, 0)] - e * area).abs() < 1e-6, "EA mismatch: {}", k[(0, 0)]);
        assert!((k[(2, 2)] - e * iz).abs() < 1e-6, "EI mismatch: {}", k[(2, 2)]);
        assert!(k[(0, 2)].abs() < 1e-9, "no axial-moment coupling for a symmetric section");

        let n = e * area * 0.001;
        let m = e * iz * 0.0002;
        assert!((r[3] - n).abs() < 1e-6, "node_j axial force");
        assert!((r[0] + n).abs() < 1e-6, "node_i axial force");
        assert!((r[5] - m).abs() < 1e-6, "node_j moment");
        assert!((r[2] + m).abs() < 1e-6, "node_i moment");
        assert_eq!(r[1], 0.0, "no material assigned to uy");
    }

    /// An off-centroid single fiber couples axial and moment response
    /// (`k_section[0][1] != 0`) — something `ZeroLength`'s independent
    /// per-DOF materials can never produce, since each of its directions is
    /// evaluated against its own material in isolation.
    #[test]
    fn zero_length_section_couples_axial_and_moment_for_an_asymmetric_section() {
        let (e, area, y): (f64, f64, f64) = (1000.0, 2.0, 3.0);
        let section = FiberSection::new(vec![Fiber::new(y, area, Material::Elastic { e })]);

        let (_n, _m, k_section) = section.trial(0.001, 0.0002);
        assert!(k_section[0][1].abs() > 1e-9, "test setup: section must actually couple");

        let node_i = Node::new([0.0, 0.0]);
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 0.001;
        node_j.displacement[2] = 0.0002;

        let zl = ZeroLengthSection::new(NodeId::default(), NodeId::default(), section);
        let (k, _r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        // b_eps0[0] = -1, b_kappa[2] = -1, so k(0,2) = k_section[0][1].
        assert!(
            (k[(0, 2)] - k_section[0][1]).abs() < 1e-9,
            "axial-moment coupling should carry through to the element tangent"
        );
    }

    /// `ZeroLengthSection3`: four corner elastic fibers reproduce `EA`,
    /// `EIy`, `EIz` (see `fiber_section.rs`'s
    /// `four_corner_elastic_fibers_reproduce_ea_eiy_eiz_with_zero_cross_coupling`)
    /// routed through the element's `B`-matrix mapping onto `[ux, ry, rz]`;
    /// `uy`/`uz`/`rx` — DOFs the section has no resultant for — stay zero.
    #[test]
    fn zero_length_section3_reproduces_biaxial_ea_eiy_eiz() {
        let (e, area, iy, iz): (f64, f64, f64, f64) = (30000.0, 4.0, 500.0, 2000.0);
        let hz = (iy / area).sqrt();
        let hy = (iz / area).sqrt();
        let a4 = area / 4.0;
        let section = FiberSection3::new(vec![
            Fiber3::new(hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(hy, -hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, -hz, a4, Material::Elastic { e }),
        ]);

        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let mut node_j = Node3::new([0.0, 0.0, 0.0]);
        node_j.displacement[0] = 0.001; // ux -> eps0
        node_j.displacement[5] = 0.0002; // rz -> kappa_z
        node_j.displacement[4] = 0.0003; // ry -> kappa_y

        let zl = ZeroLengthSection3::new(Node3Id::default(), Node3Id::default(), section);
        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((k[(0, 0)] - e * area).abs() < 1e-6, "EA mismatch");
        assert!((k[(5, 5)] - e * iz).abs() < 1e-6, "EIz mismatch");
        assert!((k[(4, 4)] - e * iy).abs() < 1e-6, "EIy mismatch");

        for dof in [1usize, 2, 3] {
            for local in [dof, SPATIAL_NDF + dof] {
                assert_eq!(r[local], 0.0, "unmodeled DOF carries no force");
                for col in 0..2 * SPATIAL_NDF {
                    assert_eq!(k[(local, col)], 0.0, "unmodeled DOF has no stiffness");
                }
            }
        }
    }

    /// A section commit propagates into the element's subsequent trial
    /// response, the same commit/trial split `Material` and `ZeroLength`
    /// already rely on — proves `ZeroLengthSection::commit` actually drives
    /// `FiberSection::commit` (not just recomputing `eps0`/`kappa` and
    /// discarding them).
    #[test]
    fn zero_length_section_commit_propagates_fiber_history() {
        let mat = || Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);

        // Pure axial fiber (y = 0) isolates the material's own history from
        // any bending coupling.
        let section = FiberSection::new(vec![Fiber::new(0.0, 1.0, mat())]);
        let mut zl = ZeroLengthSection::new(NodeId::default(), NodeId::default(), section);

        // Independent "by hand" copy, committed in lockstep, to compute the
        // expected post-yield residual stress from scratch (same pattern as
        // `fiber_section.rs`'s RC test).
        let mut hand = mat();

        let node_i = Node::new([0.0, 0.0]);
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 0.01; // well past yield (epsy = 60/29000 ~= 0.00207)

        zl.commit(&node_i, &node_j);
        hand = hand.commit(0.01);

        // Fully unload back to zero relative displacement — with committed
        // plastic strain, elastic unload from 0.01 leaves residual stress.
        node_j.displacement[0] = 0.0;
        let (_k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        let (expected_n, _tangent) = hand.trial_stress_tangent(0.0);
        assert!(expected_n.abs() > 1e-6, "test setup: material should have residual stress after yield");
        assert!(
            (r[3] - expected_n).abs() < 1e-9,
            "committed plastic history should carry into subsequent trial"
        );
    }

    /// A `uy` spring alongside the section is additive and uncoupled: it
    /// contributes independently to `k`/`r` on top of (not instead of) the
    /// section's coupled axial-moment response.
    #[test]
    fn zero_length_section_combines_with_an_independent_shear_spring() {
        let (e, area, iz): (f64, f64, f64) = (30000.0, 2.0, 1000.0);
        let h = (iz / area).sqrt();
        let section = FiberSection::new(vec![
            Fiber::new(h, area / 2.0, Material::Elastic { e }),
            Fiber::new(-h, area / 2.0, Material::Elastic { e }),
        ]);
        let shear_e = 500.0;

        let zl = ZeroLengthSection::new(NodeId::default(), NodeId::default(), section)
            .with_material(1, Material::Elastic { e: shear_e });

        let node_i = Node::new([0.0, 0.0]);
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 0.001;
        node_j.displacement[1] = 0.002; // uy -> shear spring
        node_j.displacement[2] = 0.0002;

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((k[(1, 1)] - shear_e).abs() < 1e-9, "shear spring stiffness");
        assert_eq!(k[(1, 0)], 0.0, "shear spring uncoupled from axial");
        assert_eq!(k[(1, 2)], 0.0, "shear spring uncoupled from moment");
        assert!((r[4] - shear_e * 0.002).abs() < 1e-9, "node_j shear force");

        // Section response unchanged by the added spring.
        assert!((k[(0, 0)] - e * area).abs() < 1e-6, "EA unaffected by shear spring");
        assert!((k[(2, 2)] - e * iz).abs() < 1e-6, "EI unaffected by shear spring");
    }

    /// `ZeroLengthSection3` analogue: an `rx` torsion spring alongside the
    /// section's biaxial axial-moment response — the only way this element
    /// gets torsional stiffness at all (see the type's doc comment).
    #[test]
    fn zero_length_section3_combines_with_an_independent_torsion_spring() {
        let (e, area, iy, iz): (f64, f64, f64, f64) = (30000.0, 4.0, 500.0, 2000.0);
        let hz = (iy / area).sqrt();
        let hy = (iz / area).sqrt();
        let a4 = area / 4.0;
        let section = FiberSection3::new(vec![
            Fiber3::new(hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(hy, -hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, -hz, a4, Material::Elastic { e }),
        ]);
        let torsion_e = 700.0;

        let zl = ZeroLengthSection3::new(Node3Id::default(), Node3Id::default(), section)
            .with_material(3, Material::Elastic { e: torsion_e });

        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let mut node_j = Node3::new([0.0, 0.0, 0.0]);
        node_j.displacement[0] = 0.001;
        node_j.displacement[3] = 0.004; // rx -> torsion spring
        node_j.displacement[5] = 0.0002;
        node_j.displacement[4] = 0.0003;

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((k[(3, 3)] - torsion_e).abs() < 1e-9, "torsion spring stiffness");
        assert!((r[9] - torsion_e * 0.004).abs() < 1e-9, "node_j torque");
        assert!((k[(0, 0)] - e * area).abs() < 1e-6, "EA unaffected by torsion spring");
        assert!((k[(5, 5)] - e * iz).abs() < 1e-6, "EIz unaffected by torsion spring");
        assert!((k[(4, 4)] - e * iy).abs() < 1e-6, "EIy unaffected by torsion spring");
    }
}
