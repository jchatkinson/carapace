use nalgebra::{Matrix3, SMatrix, SVector, Vector3};

use super::{Node, Node3, SpatialElementMatrix};

/// Coordinate-transformation strategy for planar frame elements: relates
/// an element's local (basic) stiffness/forces to the global system. The
/// co-rotational option follows the deformed chord and removes rigid-body
/// translation/rotation from the element's basic deformations.
///
/// `Linear` and `PDelta` select the existing small-displacement paths.
/// `Corotational2d` centralizes the current-chord kinematics used by the
/// planar elastic, displacement-based, and force-based frame elements.
/// `GeomTransf3` (below) remains a separate spatial type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomTransf {
    /// Small-displacement: local/global relationship fixed at the element's
    /// original (undeformed) orientation, no correction for the effect of
    /// axial force on bending stiffness.
    Linear,
    /// Same small-displacement geometry as `Linear`, plus a linearized
    /// geometric-stiffness correction from the current axial force (a
    /// first-order P-Delta effect) — not full corotational tracking.
    PDelta,
    /// Finite-rotation planar transformation. The current chord defines the
    /// element axes; axial extension and end rotations relative to that
    /// chord are the three objective basic deformations.
    Corotational,
}

/// Current co-rotational kinematics for one planar two-node frame member.
/// Basic displacement order is `[axial extension, theta_i - chord_rotation,
/// theta_j - chord_rotation]`, matching Xara/OpenSees `CorotCrdTransf2d`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Corotational2d {
    initial_length: f64,
    deformed_length: f64,
    basic: SVector<f64, 3>,
    /// Derivative of basic deformations with respect to
    /// `[ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]`.
    jacobian: SMatrix<f64, 3, 6>,
    /// Block-diagonal global-to-current-local rotation, used for local loads.
    global_to_local: SMatrix<f64, 6, 6>,
    /// Current chord direction in global coordinates.
    direction: SVector<f64, 2>,
}

impl Corotational2d {
    pub(crate) fn new(node_i: &Node, node_j: &Node) -> Self {
        let reference_dx = node_j.coords[0] - node_i.coords[0];
        let reference_dy = node_j.coords[1] - node_i.coords[1];
        let initial_length = reference_dx.hypot(reference_dy);
        assert!(
            initial_length > 0.0,
            "Corotational2d: element endpoints must not coincide"
        );
        let reference_angle = reference_dy.atan2(reference_dx);

        let dx = reference_dx + node_j.displacement[0] - node_i.displacement[0];
        let dy = reference_dy + node_j.displacement[1] - node_i.displacement[1];
        let deformed_length = dx.hypot(dy);
        assert!(
            deformed_length > 1e-12,
            "Corotational2d: deformed element length must be nonzero"
        );

        let direction = SVector::<f64, 2>::new(dx / deformed_length, dy / deformed_length);
        let relative_angle = dy.atan2(dx) - reference_angle;
        let chord_rotation = relative_angle.sin().atan2(relative_angle.cos());
        let basic = SVector::<f64, 3>::new(
            deformed_length - initial_length,
            node_i.displacement[2] - chord_rotation,
            node_j.displacement[2] - chord_rotation,
        );

        // d(theta_chord)/d(relative_position) = [-dy, dx] / length^2.
        let inv_l2 = 1.0 / (deformed_length * deformed_length);
        let alpha_dx = -dy * inv_l2;
        let alpha_dy = dx * inv_l2;
        let mut jacobian = SMatrix::<f64, 3, 6>::zeros();
        jacobian[(0, 0)] = -direction[0];
        jacobian[(0, 1)] = -direction[1];
        jacobian[(0, 3)] = direction[0];
        jacobian[(0, 4)] = direction[1];
        for row in [1, 2] {
            jacobian[(row, 0)] = alpha_dx;
            jacobian[(row, 1)] = alpha_dy;
            jacobian[(row, 3)] = -alpha_dx;
            jacobian[(row, 4)] = -alpha_dy;
        }
        jacobian[(1, 2)] = 1.0;
        jacobian[(2, 5)] = 1.0;

        let cx = direction[0];
        let cy = direction[1];
        #[rustfmt::skip]
        let block = SMatrix::<f64, 3, 3>::new(
             cx,  cy, 0.0,
            -cy,  cx, 0.0,
            0.0, 0.0, 1.0,
        );
        let mut global_to_local = SMatrix::<f64, 6, 6>::zeros();
        global_to_local.fixed_view_mut::<3, 3>(0, 0).copy_from(&block);
        global_to_local.fixed_view_mut::<3, 3>(3, 3).copy_from(&block);

        Self {
            initial_length,
            deformed_length,
            basic,
            jacobian,
            global_to_local,
            direction,
        }
    }

    pub(crate) fn initial_length(&self) -> f64 {
        self.initial_length
    }

    pub(crate) fn basic_deformation(&self) -> SVector<f64, 3> {
        self.basic
    }

    /// Map work-conjugate basic forces into the global nodal residual.
    pub(crate) fn global_resistance(&self, basic_force: &SVector<f64, 3>) -> SVector<f64, 6> {
        self.jacobian.transpose() * basic_force
    }

    /// Same nodal force as `global_resistance`, expressed in the *current*
    /// (deformed-chord) local frame instead of global coordinates — the
    /// corotational counterpart of `Linear`/`PDelta`'s fixed-orientation
    /// `t * r_global` recovery, used for element-force recording so a
    /// corotational member reports axial/shear/moment against its own
    /// current axis, not raw global components.
    pub(crate) fn local_resistance(&self, basic_force: &SVector<f64, 3>) -> SVector<f64, 6> {
        self.global_to_local * self.global_resistance(basic_force)
    }

    /// Consistent tangent: `Bᵀ k_basic B` plus the geometric Hessian of the
    /// basic deformation map (axial-force and end-moment contributions).
    pub(crate) fn global_tangent(
        &self,
        basic_tangent: &SMatrix<f64, 3, 3>,
        basic_force: &SVector<f64, 3>,
    ) -> SMatrix<f64, 6, 6> {
        let mut tangent = self.jacobian.transpose() * basic_tangent * self.jacobian;
        let dx = self.direction[0] * self.deformed_length;
        let dy = self.direction[1] * self.deformed_length;
        let l = self.deformed_length;
        let inv_l2 = 1.0 / (l * l);
        let inv_l4 = inv_l2 * inv_l2;

        let h_length =
            (SMatrix::<f64, 2, 2>::identity() - self.direction * self.direction.transpose()) / l;
        let h_angle = SMatrix::<f64, 2, 2>::new(
            2.0 * dx * dy * inv_l4,
            (dy * dy - dx * dx) * inv_l4,
            (dy * dy - dx * dx) * inv_l4,
            -2.0 * dx * dy * inv_l4,
        );
        let h_relative = basic_force[0] * h_length - (basic_force[1] + basic_force[2]) * h_angle;
        for i in 0..2 {
            for j in 0..2 {
                tangent[(i, j)] += h_relative[(i, j)];
                tangent[(i, j + 3)] -= h_relative[(i, j)];
                tangent[(i + 3, j)] -= h_relative[(i, j)];
                tangent[(i + 3, j + 3)] += h_relative[(i, j)];
            }
        }
        tangent
    }

    /// Follower uniform transverse load in the current local frame, integrated
    /// over the undeformed member length like the reference formulation.
    pub(crate) fn global_uniform_transverse_load(&self, load: f64) -> SVector<f64, 6> {
        let local = SVector::<f64, 6>::from_row_slice(&[
            0.0,
            load * self.initial_length / 2.0,
            load * self.initial_length.powi(2) / 12.0,
            0.0,
            load * self.initial_length / 2.0,
            -load * self.initial_length.powi(2) / 12.0,
        ]);
        self.global_to_local.transpose() * local
    }

    /// Local nodal modes for the three basic deformations, used to integrate
    /// displacement-based fiber section response in the objective basic system.
    pub(crate) fn local_basic_modes() -> SMatrix<f64, 6, 3> {
        #[rustfmt::skip]
        let modes = SMatrix::<f64, 6, 3>::from_row_slice(&[
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
        ]);
        modes
    }
}

/// `GeomTransf`'s spatial counterpart: builds a member's local basis and
/// local↔global rotation from its two end-node coordinates and an
/// orientation vector (`vec_xz`, per `docs/spatial-architecture.md`'s
/// "Coordinate frames and orientation" — `x = normalize(pJ - pI)`, `y =
/// normalize(vec_xz × x)`, `z = x × y`). Unlike planar `GeomTransf`, this
/// isn't just a behavior flag: the local-axis construction is real,
/// nontrivial, and shared by every spatial frame element (`ElasticBeamColumn3`
/// today; `DispBeamColumn3`/`ForceBeamColumn3` will reuse it once they
/// exist), so it's a real type rather than something each element inlines a
/// copy of. Verified against Xara's `LinearCrdTransf3d`.
///
/// `Linear3` and `PDelta3` are both "cheap" per the same reasoning as
/// planar `GeomTransf`'s doc comment — full spatial corotational geometry
/// is its own later milestone (M18), not a third variant here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeomTransf3 {
    /// Small-displacement: local/global relationship fixed at the member's
    /// original (undeformed) orientation, no correction for the effect of
    /// axial force on bending stiffness.
    Linear3 {
        /// A vector not parallel to the member axis, fixing the local y/z
        /// orientation (roll) — not a physical direction. E.g. `[0.0, 0.0,
        /// 1.0]` (global z) for a horizontal member picks local z as the
        /// in-plane bending axis for gravity loads.
        vec_xz: [f64; 3],
    },
    /// Same small-displacement geometry as `Linear3`, plus a linearized
    /// geometric-stiffness correction from the current axial force in
    /// *both* bending planes (a first-order P-Delta effect) — not full
    /// corotational tracking. See `ElasticBeamColumn3::geometric_stiffness`
    /// for the actual correction; this variant only selects it.
    PDelta3 { vec_xz: [f64; 3] },
}

impl GeomTransf3 {
    pub fn linear(vec_xz: [f64; 3]) -> Self {
        GeomTransf3::Linear3 { vec_xz }
    }

    pub fn p_delta(vec_xz: [f64; 3]) -> Self {
        GeomTransf3::PDelta3 { vec_xz }
    }

    /// Member length and the local-axis rotation `R` (rows = local x/y/z
    /// unit vectors, each expressed in global coordinates), so `R *
    /// global_vec` is that same vector in local coordinates.
    pub(crate) fn local_axes(&self, node_i: &Node3, node_j: &Node3) -> (f64, Matrix3<f64>) {
        let (GeomTransf3::Linear3 { vec_xz } | GeomTransf3::PDelta3 { vec_xz }) = self;

        let p_i = Vector3::from(node_i.coords);
        let p_j = Vector3::from(node_j.coords);
        let d = p_j - p_i;
        let length = d.norm();
        assert!(length > 0.0, "GeomTransf3: element endpoints must not coincide");
        let x_axis = d / length;

        let v = Vector3::from(*vec_xz);
        let y_unnormalized = v.cross(&x_axis);
        let y_norm = y_unnormalized.norm();
        assert!(y_norm > 1e-12, "GeomTransf3: vec_xz must not be parallel to the member axis");
        let y_axis = y_unnormalized / y_norm;
        let z_axis = x_axis.cross(&y_axis);

        #[rustfmt::skip]
        let r = Matrix3::new(
            x_axis[0], x_axis[1], x_axis[2],
            y_axis[0], y_axis[1], y_axis[2],
            z_axis[0], z_axis[1], z_axis[2],
        );
        (length, r)
    }

    /// Block-diagonal rotation from global to local DOFs — four copies of
    /// `r` (translation@i, rotation@i, translation@j, rotation@j).
    pub(crate) fn rotation_matrix(r: &Matrix3<f64>) -> SpatialElementMatrix {
        let mut t = SpatialElementMatrix::zeros();
        for offset in [0, 3, 6, 9] {
            t.fixed_view_mut::<3, 3>(offset, offset).copy_from(r);
        }
        t
    }
}

#[cfg(test)]
mod geom_transf3_tests {
    use super::*;

    #[test]
    fn axis_aligned_member_with_z_up_gives_identity_rotation() {
        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let node_j = Node3::new([5.0, 0.0, 0.0]);
        let transform = GeomTransf3::linear([0.0, 0.0, 1.0]);
        let (length, r) = transform.local_axes(&node_i, &node_j);
        assert_eq!(length, 5.0);
        assert!((r - Matrix3::identity()).norm() < 1e-12);
    }

    #[test]
    #[should_panic(expected = "must not coincide")]
    fn coincident_endpoints_panics() {
        let node = Node3::new([1.0, 2.0, 3.0]);
        GeomTransf3::linear([0.0, 0.0, 1.0]).local_axes(&node, &node);
    }

    #[test]
    #[should_panic(expected = "must not be parallel")]
    fn vec_xz_parallel_to_axis_panics() {
        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let node_j = Node3::new([5.0, 0.0, 0.0]);
        GeomTransf3::linear([1.0, 0.0, 0.0]).local_axes(&node_i, &node_j);
    }
}
