/// Coordinate-transformation strategy for frame elements: relates an
/// element's local (basic) stiffness/forces to the global system. Closed
/// enum (§2.1) — `Linear` and `PDelta` are both "cheap" per the
/// implementation plan §3.4; `Corotational` (large-displacement) is its own
/// milestone (M9) if/when actually needed, not a third variant here.
///
/// Purely a behavior-selecting flag: the 2×2 direction-cosine rotation
/// itself is cheap enough that each planar frame element inlines its own
/// copy (`ElasticBeamColumn`/`DispBeamColumn`/`ForceBeamColumn` each have
/// their own `geometry`/`transformation` methods) rather than sharing one
/// through this type. `GeomTransf3` (below) can't get away with that — see
/// its doc comment.
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
}

use nalgebra::{Matrix3, Vector3};

use super::{Node3, SpatialElementMatrix};

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
