//! `EmbeddedBar`: one segment of an embedded rebar polyline (plan
//! `docs/csfm-plan.md` §2.3), 2D only. Its nodes are its **host**
//! element's own 3 or 4 nodes (never extra DOFs), so it fits
//! `MAX_ELEMENT_NODES` and the existing `TangentSink`/`VectorSink` path.
//!
//! The segment measures the *average* axial strain between its two end
//! points, which are generally not host nodes: `eps = d . (sum N_b,i u_i -
//! sum N_a,i u_i) / L`, where `N_a`/`N_b` are the host's own shape
//! functions evaluated at the segment's start/end natural coordinates, `d`
//! is the segment's (undeformed) unit direction and `L` its length. Both
//! are baked into one strain-displacement coefficient per host node at
//! construction (`coeffs`), exactly as `Truss::geometry` bakes its own
//! direction cosines from the undeformed configuration — small-displacement
//! theory, so the coefficients never change after that.
//!
//! Discretizing a polyline into these segments and finding each one's host
//! is `embed_polyline` (`model/continuum/embed.rs`), not this module: this
//! type only knows how to assemble *given* a host and its own shape values.

use nalgebra::{SMatrix, SVector};

use super::plane_common::{plane_dofs, PlaneView};
use super::{DofMask, ElementForce, TangentSink};
use crate::model::{LimitKind, LimitStateConfig, Material, NodeId};

/// `EmbeddedBar`'s host, carrying the host's own node list and this
/// segment's strain-displacement coefficient per host node (`coeffs[i] =
/// (N_b[i] - N_a[i]) * direction / length`, so `strain = sum_i coeffs[i] .
/// host_displacement[i]`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EmbeddedBarHost {
    Tri3 {
        nodes: [NodeId; 3],
        coeffs: [[f64; 2]; 3],
    },
    Quad4 {
        nodes: [NodeId; 4],
        coeffs: [[f64; 2]; 4],
    },
}

#[derive(Debug, Clone)]
pub struct EmbeddedBar {
    host: EmbeddedBarHost,
    length: f64,
    area: f64,
    material: Material,
}

impl EmbeddedBar {
    /// Builds a segment bound to a `Tri3` host. `n_a`/`n_b` are the host's
    /// own shape-function values at the segment's start/end natural
    /// coordinates; `direction` is the segment's unit direction and
    /// `length` its length, both in the undeformed configuration.
    pub fn in_tri3(
        nodes: [NodeId; 3],
        n_a: [f64; 3],
        n_b: [f64; 3],
        direction: [f64; 2],
        length: f64,
        area: f64,
        material: Material,
    ) -> Self {
        let coeffs = std::array::from_fn(|i| {
            let a = n_b[i] - n_a[i];
            [a * direction[0] / length, a * direction[1] / length]
        });
        EmbeddedBar {
            host: EmbeddedBarHost::Tri3 { nodes, coeffs },
            length,
            area,
            material,
        }
    }

    /// Builds a segment bound to a `Quad4` host. See [`Self::in_tri3`].
    pub fn in_quad4(
        nodes: [NodeId; 4],
        n_a: [f64; 4],
        n_b: [f64; 4],
        direction: [f64; 2],
        length: f64,
        area: f64,
        material: Material,
    ) -> Self {
        let coeffs = std::array::from_fn(|i| {
            let a = n_b[i] - n_a[i];
            [a * direction[0] / length, a * direction[1] / length]
        });
        EmbeddedBar {
            host: EmbeddedBarHost::Quad4 { nodes, coeffs },
            length,
            area,
            material,
        }
    }

    pub fn host(&self) -> &EmbeddedBarHost {
        &self.host
    }

    pub fn material(&self) -> &Material {
        &self.material
    }

    pub fn area(&self) -> f64 {
        self.area
    }

    pub(super) fn dof_mask(&self) -> DofMask {
        DofMask::none().with(0).with(1)
    }

    fn strain<const NN: usize>(
        nodes: &[NodeId; NN],
        coeffs: &[[f64; 2]; NN],
        view: &PlaneView<'_>,
    ) -> f64 {
        (0..NN)
            .map(|i| {
                let displacement = view.get(nodes[i]).displacement;
                coeffs[i][0] * displacement[0] + coeffs[i][1] * displacement[1]
            })
            .sum()
    }

    fn committed_strain(&self, view: &PlaneView<'_>) -> f64 {
        match &self.host {
            EmbeddedBarHost::Tri3 { nodes, coeffs } => Self::strain(nodes, coeffs, view),
            EmbeddedBarHost::Quad4 { nodes, coeffs } => Self::strain(nodes, coeffs, view),
        }
    }

    /// Tangent and resisting force for an `NN`-node host (`N = 2 * NN`).
    /// Virtual work over the segment is `sigma * area * L * delta_eps`
    /// with `delta_eps = c . delta_u` (`c` flattening `coeffs`, which
    /// already carries its own `1/L` factor — see the module doc
    /// comment), so `R = sigma * area * L * c` and, since `d(sigma)/d(u)
    /// = Et * c`, `K = Et * area * L * (c ⊗ c)`.
    fn tangent_and_resistance<const NN: usize, const N: usize>(
        nodes: &[NodeId; NN],
        coeffs: &[[f64; 2]; NN],
        length: f64,
        area: f64,
        material: &Material,
        view: &PlaneView<'_>,
    ) -> (SMatrix<f64, N, N>, SVector<f64, N>) {
        debug_assert_eq!(N, 2 * NN);
        let strain = Self::strain(nodes, coeffs, view);
        let (stress, tangent_modulus) = material.trial_stress_tangent(strain);
        let mut c = SVector::<f64, N>::zeros();
        for i in 0..NN {
            c[2 * i] = coeffs[i][0];
            c[2 * i + 1] = coeffs[i][1];
        }
        let k = (tangent_modulus * area * length) * (c * c.transpose());
        let r = (stress * area * length) * c;
        (k, r)
    }

    pub(super) fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        view: &PlaneView<'_>,
        sink: &mut S,
    ) {
        match &self.host {
            EmbeddedBarHost::Tri3 { nodes, coeffs } => {
                let (k, r) = Self::tangent_and_resistance::<3, 6>(
                    nodes,
                    coeffs,
                    self.length,
                    self.area,
                    &self.material,
                    view,
                );
                sink.add(&plane_dofs::<3, 6>(nodes), &k, &r);
            }
            EmbeddedBarHost::Quad4 { nodes, coeffs } => {
                let (k, r) = Self::tangent_and_resistance::<4, 8>(
                    nodes,
                    coeffs,
                    self.length,
                    self.area,
                    &self.material,
                    view,
                );
                sink.add(&plane_dofs::<4, 8>(nodes), &k, &r);
            }
        }
    }

    pub(super) fn commit(&mut self, view: &PlaneView<'_>) {
        let strain = self.committed_strain(view);
        self.material = self.material.commit(strain);
    }

    pub(super) fn local_force_width(&self) -> usize {
        match &self.host {
            EmbeddedBarHost::Tri3 { .. } => 6,
            EmbeddedBarHost::Quad4 { .. } => 8,
        }
    }

    pub(super) fn local_force(&self, view: &PlaneView<'_>) -> ElementForce {
        match &self.host {
            EmbeddedBarHost::Tri3 { nodes, coeffs } => {
                let (_, r) = Self::tangent_and_resistance::<3, 6>(
                    nodes,
                    coeffs,
                    self.length,
                    self.area,
                    &self.material,
                    view,
                );
                r.into()
            }
            EmbeddedBarHost::Quad4 { nodes, coeffs } => {
                let (_, r) = Self::tangent_and_resistance::<4, 8>(
                    nodes,
                    coeffs,
                    self.length,
                    self.area,
                    &self.material,
                    view,
                );
                r.into()
            }
        }
    }

    /// Axial force, strain, stress at the current committed state, and
    /// whether this segment's material has ruptured — recorded as one
    /// `GaussResponse`-shaped point (`strain`/`stress` populate only the
    /// first component; the rest are `0.0`), matching the pattern
    /// `gauss_responses` already uses for continuum points.
    pub fn local_response(&self, view: &PlaneView<'_>) -> (f64, f64, f64, bool) {
        let strain = self.committed_strain(view);
        let (stress, _) = self.material.trial_stress_tangent(strain);
        let axial_force = stress * self.area;
        let ruptured = matches!(&self.material, Material::CsfmRebar(_))
            && self.material.csfm_rebar_ruptured(strain);
        (axial_force, strain, stress, ruptured)
    }

    /// The `CsfmRebar` rupture/strain-limit utilization at the current
    /// (trial) strain (`None` if the material isn't `CsfmRebar`).
    pub(super) fn limit_utilization(
        &self,
        view: &PlaneView<'_>,
        limits: &LimitStateConfig,
    ) -> Option<(f64, LimitKind)> {
        if !matches!(&self.material, Material::CsfmRebar(_)) {
            return None;
        }
        let strain = self.committed_strain(view);
        Some(
            self.material
                .csfm_rebar_limit_utilization(strain, limits.rebar_strain_limit),
        )
    }
}
