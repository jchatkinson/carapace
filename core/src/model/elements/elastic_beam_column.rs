use nalgebra::{SMatrix, SVector};

use super::super::transform::Corotational2d;
use super::super::{GeomTransf, GeomTransf3, Node, Node3, Node3Id, NodeId};
use super::truss::{SpatialElementMatrix, SpatialElementVector};

/// A 2-node, prismatic, linear-elastic 2D beam-column (Euler-Bernoulli, no
/// shear deformation): closed-form stiffness, no iteration (§3.1) — unlike
/// `Truss`/`ZeroLength`, there's no `Material` dispatch here, since the
/// element's response is fully determined by `e`/`a`/`iz` with no nonlinear
/// stress-strain law to evaluate.
#[derive(Debug, Clone)]
pub struct ElasticBeamColumn {
    pub node_i: NodeId,
    pub node_j: NodeId,
    pub e: f64,
    pub a: f64,
    pub iz: f64,
    pub transform: GeomTransf,
    /// Mass per unit volume. Zero (the default) means massless — existing
    /// models are unaffected unless they opt in via `with_density`.
    pub density: f64,
}

impl ElasticBeamColumn {
    pub fn new(node_i: NodeId, node_j: NodeId, e: f64, a: f64, iz: f64, transform: GeomTransf) -> Self {
        ElasticBeamColumn {
            node_i,
            node_j,
            e,
            a,
            iz,
            transform,
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

    /// Local-DOF-order `[u1, v1, θ1, u2, v2, θ2]` elastic stiffness —
    /// standard Euler-Bernoulli prismatic beam (Przemieniecki / Cook §2).
    fn local_elastic_stiffness(&self, length: f64) -> SMatrix<f64, 6, 6> {
        let (e, a, i) = (self.e, self.a, self.iz);
        let l = length;
        let ea_l = e * a / l;
        let ei = e * i;

        #[rustfmt::skip]
        let k = SMatrix::<f64, 6, 6>::from_row_slice(&[
             ea_l,             0.0,             0.0,  -ea_l,             0.0,             0.0,
              0.0,  12.0*ei/l.powi(3),  6.0*ei/l.powi(2),   0.0, -12.0*ei/l.powi(3),  6.0*ei/l.powi(2),
              0.0,   6.0*ei/l.powi(2),      4.0*ei/l,       0.0,  -6.0*ei/l.powi(2),      2.0*ei/l,
            -ea_l,             0.0,             0.0,   ea_l,             0.0,             0.0,
              0.0, -12.0*ei/l.powi(3), -6.0*ei/l.powi(2),   0.0,  12.0*ei/l.powi(3), -6.0*ei/l.powi(2),
              0.0,   6.0*ei/l.powi(2),      2.0*ei/l,       0.0,  -6.0*ei/l.powi(2),      4.0*ei/l,
        ]);
        k
    }

    /// First-order (linearized) geometric-stiffness correction from the
    /// current axial force `p` (tension positive), added to the elastic
    /// bending block for `GeomTransf::PDelta` — the standard consistent
    /// geometric stiffness matrix (e.g. Przemieniecki §2), *not* a full
    /// corotational update. Only the tangent used to solve for the next
    /// displacement increment is corrected this way; the resisting-force
    /// recovery (`form_tangent_and_resistance`'s second return value)
    /// stays purely elastic. Consistent path-following P-Delta (where the
    /// resisting force itself must reflect the same correction) needs
    /// Newton iteration to track correctly and is deferred alongside M4.
    fn geometric_stiffness(&self, p: f64, length: f64) -> SMatrix<f64, 6, 6> {
        let l = length;
        let mut kg = SMatrix::<f64, 6, 6>::zeros();
        #[rustfmt::skip]
        let block = [
            [ 6.0/5.0,       l/10.0,      -6.0/5.0,       l/10.0],
            [   l/10.0,  2.0*l*l/15.0,      -l/10.0,   -l*l/30.0],
            [-6.0/5.0,      -l/10.0,       6.0/5.0,      -l/10.0],
            [   l/10.0,    -l*l/30.0,      -l/10.0,  2.0*l*l/15.0],
        ];
        let indices = [1, 2, 4, 5];
        for (bi, &gi) in indices.iter().enumerate() {
            for (bj, &gj) in indices.iter().enumerate() {
                kg[(gi, gj)] = (p / l) * block[bi][bj];
            }
        }
        kg
    }

    /// Block-diagonal rotation from global to local DOFs: `d_local = T * d_global`.
    fn transformation(&self, cx: f64, cy: f64) -> SMatrix<f64, 6, 6> {
        #[rustfmt::skip]
        let block = SMatrix::<f64, 3, 3>::new(
             cx,  cy, 0.0,
            -cy,  cx, 0.0,
            0.0, 0.0, 1.0,
        );
        let mut t = SMatrix::<f64, 6, 6>::zeros();
        t.fixed_view_mut::<3, 3>(0, 0).copy_from(&block);
        t.fixed_view_mut::<3, 3>(3, 3).copy_from(&block);
        t
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, 6, 6>, SVector<f64, 6>) {
        if self.transform == GeomTransf::Corotational {
            let transform = Corotational2d::new(node_i, node_j);
            let l = transform.initial_length();
            let ea_l = self.e * self.a / l;
            let eiz_l = self.e * self.iz / l;
            let k_basic = SMatrix::<f64, 3, 3>::new(
                ea_l,
                0.0,
                0.0,
                0.0,
                4.0 * eiz_l,
                2.0 * eiz_l,
                0.0,
                2.0 * eiz_l,
                4.0 * eiz_l,
            );
            let q = k_basic * transform.basic_deformation();
            return (
                transform.global_tangent(&k_basic, &q),
                transform.global_resistance(&q),
            );
        }

        let (length, cx, cy) = self.geometry(node_i, node_j);
        let t = self.transformation(cx, cy);
        let k_local = self.local_elastic_stiffness(length);

        let d_global = SVector::<f64, 6>::from_row_slice(&[
            node_i.displacement[0],
            node_i.displacement[1],
            node_i.displacement[2],
            node_j.displacement[0],
            node_j.displacement[1],
            node_j.displacement[2],
        ]);
        let d_local = t * d_global;

        let resistance_local = k_local * d_local;
        let resistance = t.transpose() * resistance_local;

        let k_total_local = match self.transform {
            GeomTransf::Linear => k_local,
            GeomTransf::PDelta => {
                let axial_force = resistance_local[3]; // tension-positive axial force at node j
                k_local + self.geometric_stiffness(axial_force, length)
            }
            GeomTransf::Corotational => unreachable!("handled above"),
        };
        let k = t.transpose() * k_total_local * t;

        (k, resistance)
    }

    /// Local nodal force — `resistance_local`/`Corotational2d::
    /// local_resistance`, the pre-transform value `form_tangent_and_
    /// resistance` already computes and discards, recomputed fresh here
    /// (stateless/elastic, so this is cheap and always exactly consistent
    /// with the current committed displacement — no cache to keep in
    /// sync).
    pub(super) fn local_force(&self, node_i: &Node, node_j: &Node) -> SVector<f64, 6> {
        if self.transform == GeomTransf::Corotational {
            let transform = Corotational2d::new(node_i, node_j);
            let l = transform.initial_length();
            let ea_l = self.e * self.a / l;
            let eiz_l = self.e * self.iz / l;
            let k_basic = SMatrix::<f64, 3, 3>::new(
                ea_l, 0.0, 0.0, 0.0, 4.0 * eiz_l, 2.0 * eiz_l, 0.0, 2.0 * eiz_l, 4.0 * eiz_l,
            );
            let q = k_basic * transform.basic_deformation();
            return transform.local_resistance(&q);
        }

        let (length, cx, cy) = self.geometry(node_i, node_j);
        let t = self.transformation(cx, cy);
        let k_local = self.local_elastic_stiffness(length);
        let d_global = SVector::<f64, 6>::from_row_slice(&[
            node_i.displacement[0],
            node_i.displacement[1],
            node_i.displacement[2],
            node_j.displacement[0],
            node_j.displacement[1],
            node_j.displacement[2],
        ]);
        k_local * (t * d_global)
    }

    /// Equivalent nodal load (global coordinates) from a uniform transverse
    /// load `w` (force/length, local +y direction) applied to this element
    /// by whichever `LoadPattern` is currently being assembled (§3.4;
    /// `Domain::assemble_reference_load` passes `w` in — it's no longer a
    /// field on the element itself, since a pattern-scoped load can't live
    /// on `Element`, same reasoning as `Node`'s load leaving `Node`), via
    /// consistent (virtual-work) Hermite cubic shape-function integration.
    pub(super) fn form_load_vector(&self, node_i: &Node, node_j: &Node, w: f64) -> SVector<f64, 6> {
        if w == 0.0 {
            return SVector::<f64, 6>::zeros();
        }
        if self.transform == GeomTransf::Corotational {
            return Corotational2d::new(node_i, node_j).global_uniform_transverse_load(w);
        }
        let (length, cx, cy) = self.geometry(node_i, node_j);
        let t = self.transformation(cx, cy);
        let l = length;
        let local = SVector::<f64, 6>::from_row_slice(&[
            0.0,
            w * l / 2.0,
            w * l * l / 12.0,
            0.0,
            w * l / 2.0,
            -w * l * l / 12.0,
        ]);
        t.transpose() * local
    }

    /// Lumped mass: half the element's total mass (`density * a * length`)
    /// at each node's translational DOFs, zero rotational contribution —
    /// the simplest standard lumped-mass model (§3.4: "lumped, to start");
    /// a consistent (non-diagonal) mass matrix or a nonzero rotational
    /// lumped inertia is a further refinement, not built until needed.
    pub(super) fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, 6> {
        let (length, _cx, _cy) = self.geometry(node_i, node_j);
        let half = self.density * self.a * length / 2.0;
        SVector::<f64, 6>::from_column_slice(&[half, half, 0.0, half, half, 0.0])
    }
}

/// `ElasticBeamColumn`'s spatial counterpart: a 2-node, prismatic,
/// linear-elastic 3D Euler-Bernoulli beam-column (`E, G, A, J, Iy, Iz`),
/// closed-form, `GeomTransf3::Linear3` only — see this struct's
/// module-level context (spatial-architecture plan, "Elements and
/// transforms") for why spatial P-Delta isn't implemented yet: it needs
/// force/tangent recovery already consistent in the linear case first.
///
/// Local DOF order per node `[u, v, w, rx, ry, rz]` (axial, then bending
/// about local z causing local-y displacement `v`/`rz`, then bending about
/// local y causing local-z displacement `w`/`ry`, torsion `rx`) — matches
/// `SpatialDof`'s `[Ux, Uy, Uz, Rx, Ry, Rz]` order exactly, so no per-node
/// DOF reordering is needed between local and global, only `transform`'s
/// rotation. Verified against Xara/OpenSees's `ElasticBeam3d` (basic-
/// stiffness terms, local DOF order).
#[derive(Debug, Clone)]
pub struct ElasticBeamColumn3 {
    pub node_i: Node3Id,
    pub node_j: Node3Id,
    pub e: f64,
    pub g: f64,
    pub a: f64,
    /// Torsion constant (St. Venant), not a polar moment of inertia for
    /// non-circular sections.
    pub j: f64,
    /// Second moment of area about local y — governs bending that
    /// displaces local z (`w`/`ry`).
    pub iy: f64,
    /// Second moment of area about local z — governs bending that
    /// displaces local y (`v`/`rz`), the plane a 2D `ElasticBeamColumn`
    /// implicitly bends in.
    pub iz: f64,
    pub transform: GeomTransf3,
    /// Mass per unit volume; zero means massless.
    pub density: f64,
}

impl ElasticBeamColumn3 {
    #[allow(clippy::too_many_arguments)]
    pub fn new(node_i: Node3Id, node_j: Node3Id, e: f64, g: f64, a: f64, j: f64, iy: f64, iz: f64, transform: GeomTransf3) -> Self {
        ElasticBeamColumn3 {
            node_i,
            node_j,
            e,
            g,
            a,
            j,
            iy,
            iz,
            transform,
            density: 0.0,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// Local-DOF-order 3D Euler-Bernoulli stiffness (no shear deformation).
    /// Axial and torsion are simple 2×2 blocks; the two bending planes are
    /// each the same 2D beam 4×4 block `ElasticBeamColumn::local_elastic_
    /// stiffness` uses, keyed to `[v, rz]` (local z bending, `Iz`) or `[w,
    /// ry]` (local y bending, `Iy`) — the `Iy` block's off-diagonal
    /// (shear/moment coupling) terms carry the opposite sign from `Iz`'s,
    /// the standard result of `y`/`z` forming a right-handed triad with
    /// `x`: positive rotation about `+y` produces a *negative* `w` slope,
    /// where positive rotation about `+z` produces a *positive* `v` slope.
    fn local_elastic_stiffness(&self, length: f64) -> SpatialElementMatrix {
        let (e, g, a, j, iy, iz) = (self.e, self.g, self.a, self.j, self.iy, self.iz);
        let l = length;
        let mut k = SpatialElementMatrix::zeros();

        let ea_l = e * a / l;
        k[(0, 0)] = ea_l;
        k[(6, 6)] = ea_l;
        k[(0, 6)] = -ea_l;
        k[(6, 0)] = -ea_l;

        let gj_l = g * j / l;
        k[(3, 3)] = gj_l;
        k[(9, 9)] = gj_l;
        k[(3, 9)] = -gj_l;
        k[(9, 3)] = -gj_l;

        // Bending about local z: [v1, rz1, v2, rz2] = indices [1, 5, 7, 11].
        let eiz = e * iz;
        let z_dof = [1, 5, 7, 11];
        #[rustfmt::skip]
        let z_block = [
            [ 12.0*eiz/l.powi(3),   6.0*eiz/l.powi(2), -12.0*eiz/l.powi(3),   6.0*eiz/l.powi(2)],
            [  6.0*eiz/l.powi(2),       4.0*eiz/l,       -6.0*eiz/l.powi(2),      2.0*eiz/l],
            [-12.0*eiz/l.powi(3),  -6.0*eiz/l.powi(2),  12.0*eiz/l.powi(3),  -6.0*eiz/l.powi(2)],
            [  6.0*eiz/l.powi(2),       2.0*eiz/l,       -6.0*eiz/l.powi(2),      4.0*eiz/l],
        ];
        for (bi, &gi) in z_dof.iter().enumerate() {
            for (bj, &gj) in z_dof.iter().enumerate() {
                k[(gi, gj)] = z_block[bi][bj];
            }
        }

        // Bending about local y: [w1, ry1, w2, ry2] = indices [2, 4, 8, 10];
        // off-diagonal signs flipped relative to the z block (see doc comment).
        let eiy = e * iy;
        let y_dof = [2, 4, 8, 10];
        #[rustfmt::skip]
        let y_block = [
            [ 12.0*eiy/l.powi(3),  -6.0*eiy/l.powi(2), -12.0*eiy/l.powi(3),  -6.0*eiy/l.powi(2)],
            [ -6.0*eiy/l.powi(2),       4.0*eiy/l,        6.0*eiy/l.powi(2),      2.0*eiy/l],
            [-12.0*eiy/l.powi(3),   6.0*eiy/l.powi(2),  12.0*eiy/l.powi(3),   6.0*eiy/l.powi(2)],
            [ -6.0*eiy/l.powi(2),       2.0*eiy/l,        6.0*eiy/l.powi(2),      4.0*eiy/l],
        ];
        for (bi, &gi) in y_dof.iter().enumerate() {
            for (bj, &gj) in y_dof.iter().enumerate() {
                k[(gi, gj)] = y_block[bi][bj];
            }
        }

        k
    }

    /// First-order (linearized) geometric-stiffness correction from the
    /// current axial force `p` (tension positive), added to the elastic
    /// bending blocks for `GeomTransf3::PDelta3` — the spatial analogue of
    /// planar `ElasticBeamColumn::geometric_stiffness` (see its doc comment
    /// for the "linearized, not full corotational" caveat, which applies
    /// identically here). Corrects *both* bending planes: the standard
    /// consistent geometric-stiffness block depends only on the transverse
    /// shape functions and axial force, not on `EI`, so it's the exact same
    /// numeric block in each plane — but expressed in `[w, ry]` rather than
    /// `[v, rz]`, its off-diagonal (shear/rotation coupling) entries flip
    /// sign for the same reason `local_elastic_stiffness`'s `y_block` does
    /// relative to `z_block` (see that method's doc comment): `ry = -dw/dx`
    /// while `rz = +dv/dx`, and this block's off-diagonal terms are each
    /// linear in exactly one rotation, so substituting the sign-flipped
    /// rotation flips exactly those terms, leaving the pure-translation and
    /// pure-rotation entries unchanged.
    fn geometric_stiffness(&self, p: f64, length: f64) -> SpatialElementMatrix {
        let l = length;
        let mut kg = SpatialElementMatrix::zeros();

        // Bending about local z: [v1, rz1, v2, rz2] = indices [1, 5, 7, 11].
        let z_dof = [1, 5, 7, 11];
        #[rustfmt::skip]
        let z_block = [
            [ 6.0/5.0,       l/10.0,      -6.0/5.0,       l/10.0],
            [   l/10.0,  2.0*l*l/15.0,      -l/10.0,   -l*l/30.0],
            [-6.0/5.0,      -l/10.0,       6.0/5.0,      -l/10.0],
            [   l/10.0,    -l*l/30.0,      -l/10.0,  2.0*l*l/15.0],
        ];
        for (bi, &gi) in z_dof.iter().enumerate() {
            for (bj, &gj) in z_dof.iter().enumerate() {
                kg[(gi, gj)] = (p / l) * z_block[bi][bj];
            }
        }

        // Bending about local y: [w1, ry1, w2, ry2] = indices [2, 4, 8, 10];
        // off-diagonal (w-ry coupling) signs flipped relative to the z
        // block — see this method's doc comment.
        let y_dof = [2, 4, 8, 10];
        #[rustfmt::skip]
        let y_block = [
            [ 6.0/5.0,      -l/10.0,      -6.0/5.0,      -l/10.0],
            [  -l/10.0,  2.0*l*l/15.0,       l/10.0,   -l*l/30.0],
            [-6.0/5.0,       l/10.0,       6.0/5.0,       l/10.0],
            [  -l/10.0,    -l*l/30.0,       l/10.0,  2.0*l*l/15.0],
        ];
        for (bi, &gi) in y_dof.iter().enumerate() {
            for (bj, &gj) in y_dof.iter().enumerate() {
                kg[(gi, gj)] = (p / l) * y_block[bi][bj];
            }
        }

        kg
    }

    pub(super) fn form_tangent_and_resistance(&self, node_i: &Node3, node_j: &Node3) -> (SpatialElementMatrix, SpatialElementVector) {
        let (length, r) = self.transform.local_axes(node_i, node_j);
        let t = GeomTransf3::rotation_matrix(&r);
        let k_local = self.local_elastic_stiffness(length);

        let d_global = SpatialElementVector::from_iterator(node_i.displacement.iter().chain(node_j.displacement.iter()).copied());
        let d_local = t * d_global;

        let resistance_local = k_local * d_local;
        let resistance = t.transpose() * resistance_local;

        let k_total_local = match self.transform {
            GeomTransf3::Linear3 { .. } => k_local,
            GeomTransf3::PDelta3 { .. } => {
                let axial_force = resistance_local[6]; // tension-positive axial force at node j
                k_local + self.geometric_stiffness(axial_force, length)
            }
        };
        let k = t.transpose() * k_total_local * t;

        (k, resistance)
    }

    /// Local nodal force, spatial counterpart of the planar `local_force`
    /// above — same `resistance_local` value `form_tangent_and_resistance`
    /// already computes and discards, recomputed fresh (cheap, stateless).
    pub(super) fn local_force(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        let (length, r) = self.transform.local_axes(node_i, node_j);
        let t = GeomTransf3::rotation_matrix(&r);
        let k_local = self.local_elastic_stiffness(length);
        let d_global = SpatialElementVector::from_iterator(node_i.displacement.iter().chain(node_j.displacement.iter()).copied());
        k_local * (t * d_global)
    }

    /// Equivalent nodal load (global coordinates) from local transverse
    /// loads `wy`/`wz` (force/length, local `+y`/`+z` directions) — the
    /// direct spatial generalization of planar `ElasticBeamColumn::
    /// form_load_vector`'s consistent (virtual-work) Hermite-cubic
    /// equivalent load, applied independently in each bending plane. The
    /// `wy`-induced terms (`[v, rz]`, indices `[1,5,7,11]`) are numerically
    /// identical to the planar formula; the `wz`-induced terms (`[w, ry]`,
    /// indices `[2,4,8,10]`) carry the same rotation-coefficient sign flip
    /// as `local_elastic_stiffness`'s `y_block` relative to `z_block` (`ry
    /// = -dw/dx` vs `rz = +dv/dx`) — verified against Xara/OpenSees's
    /// `ElasticBeam3d::addLoad`'s `Beam3dUniformLoad` fixed-end-force
    /// derivation (`q0`: `MI=-Mz, MJ=Mz` for the `wy`-driven `z`-bending
    /// pair and `MI=My, MJ=-My` for the `wz`-driven `y`-bending pair, each
    /// negated here since a fixed-end restraint force is the negative of
    /// the equivalent nodal load it corresponds to), not just asserted by
    /// analogy to the planar case.
    pub(super) fn form_load_vector(&self, node_i: &Node3, node_j: &Node3, wy: f64, wz: f64) -> SpatialElementVector {
        if wy == 0.0 && wz == 0.0 {
            return SpatialElementVector::zeros();
        }
        let (length, r) = self.transform.local_axes(node_i, node_j);
        let t = GeomTransf3::rotation_matrix(&r);
        let l = length;

        let mut local = SpatialElementVector::zeros();
        local[1] = wy * l / 2.0;
        local[5] = wy * l * l / 12.0;
        local[7] = wy * l / 2.0;
        local[11] = -wy * l * l / 12.0;

        local[2] = wz * l / 2.0;
        local[4] = -wz * l * l / 12.0;
        local[8] = wz * l / 2.0;
        local[10] = wz * l * l / 12.0;

        t.transpose() * local
    }

    /// Lumped mass: half the element's total mass at each node's three
    /// translational DOFs, zero rotational contribution — same
    /// simplification as planar `ElasticBeamColumn::form_mass`.
    pub(super) fn form_mass(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        let (length, _r) = self.transform.local_axes(node_i, node_j);
        let half = self.density * self.a * length / 2.0;
        let mut mass = SpatialElementVector::zeros();
        for dof in [0, 1, 2, 6, 7, 8] {
            mass[dof] = half;
        }
        mass
    }
}

#[cfg(test)]
mod spatial_tests {
    use nalgebra::Matrix3;

    use super::*;

    /// A cantilever along global x with `vec_xz = [0, 0, 1]` (global z),
    /// which makes local `[x, y, z] = global [x, y, z]` exactly (`R` is the
    /// identity) — isolates the local stiffness formula itself from the
    /// `GeomTransf3` rotation, which `transform.rs`'s and `truss.rs`'s own
    /// tests already cover for the transform-only case. With `node_i` fully
    /// fixed, the free-DOF system reduces to `local_elastic_stiffness`'s
    /// `[6..12, 6..12]` block applied directly (no rotation needed), so tip
    /// displacements can be checked against textbook cantilever closed
    /// forms for axial, biaxial bending, and torsion simultaneously (all
    /// four modes are decoupled in a prismatic doubly-symmetric section).
    #[test]
    fn axis_aligned_cantilever_matches_closed_form_axial_biaxial_bending_and_torsion() {
        let (e, g, a, j, iy, iz) = (30_000.0, 12_000.0, 20.0, 5.0, 400.0, 800.0);
        let length = 100.0;
        let transform = GeomTransf3::linear([0.0, 0.0, 1.0]);
        let beam = ElasticBeamColumn3::new(Node3Id::default(), Node3Id::default(), e, g, a, j, iy, iz, transform);

        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let node_j = Node3::new([length, 0.0, 0.0]);
        let (geom_length, r) = beam.transform.local_axes(&node_i, &node_j);
        assert_eq!(geom_length, length);
        assert!((r - Matrix3::identity()).norm() < 1e-12, "vec_xz=[0,0,1] with a member along x must give R = identity");

        let k_local = beam.local_elastic_stiffness(length);

        // node_i fixed at zero, so the free (node_j) system is exactly
        // k_local's bottom-right 6x6 block.
        let k_free = k_local.fixed_view::<6, 6>(6, 6).into_owned();
        let (axial_force, fy, fz, torque) = (100.0, 50.0, 30.0, 20.0);
        // [ux, uy, uz, rx, ry, rz] at node_j; no applied end moments (my, mz).
        let f = SVector::<f64, 6>::from_row_slice(&[axial_force, fy, fz, torque, 0.0, 0.0]);
        let d = k_free.lu().solve(&f).expect("k_free is nonsingular for a real prismatic section");

        let expected_ux = axial_force * length / (e * a);
        let expected_uy = fy * length.powi(3) / (3.0 * e * iz);
        let expected_rz = fy * length.powi(2) / (2.0 * e * iz);
        let expected_uz = fz * length.powi(3) / (3.0 * e * iy);
        let expected_ry = -fz * length.powi(2) / (2.0 * e * iy);
        let expected_rx = torque * length / (g * j);

        let tol = 1e-9;
        assert!((d[0] - expected_ux).abs() < tol, "axial: {} vs {}", d[0], expected_ux);
        assert!((d[1] - expected_uy).abs() < tol, "z-bending deflection: {} vs {}", d[1], expected_uy);
        assert!((d[5] - expected_rz).abs() < tol, "z-bending rotation: {} vs {}", d[5], expected_rz);
        assert!((d[2] - expected_uz).abs() < tol, "y-bending deflection: {} vs {}", d[2], expected_uz);
        assert!((d[4] - expected_ry).abs() < tol, "y-bending rotation: {} vs {}", d[4], expected_ry);
        assert!((d[3] - expected_rx).abs() < tol, "torsion: {} vs {}", d[3], expected_rx);
    }

    /// The same cantilever driven through the full `form_tangent_and_
    /// resistance` (transform included, not just the local block) with a
    /// skew member and a non-trivial orientation, checked by round-tripping:
    /// applying the closed-form local displacement (rotated into global
    /// coordinates) must reproduce the applied global load exactly, proving
    /// `GeomTransf3`/the local stiffness compose correctly rather than
    /// being separately-plausible but inconsistent.
    #[test]
    fn skew_cantilever_transform_round_trips_a_known_local_solution() {
        let (e, g, a, j, iy, iz) = (30_000.0, 12_000.0, 20.0, 5.0, 400.0, 800.0);
        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let node_j = Node3::new([3.0, 4.0, 12.0]); // length 13, skew in all three axes
        let transform = GeomTransf3::linear([1.0, 0.0, 0.0]);
        let beam = ElasticBeamColumn3::new(Node3Id::default(), Node3Id::default(), e, g, a, j, iy, iz, transform);

        let (length, r) = beam.transform.local_axes(&node_i, &node_j);
        assert!((length - 13.0).abs() < 1e-12);

        // A purely local axial elongation at node_j: local d = [delta,0,0,0,0,0].
        let delta = 0.005;
        let mut local_d = SVector::<f64, 6>::zeros();
        local_d[0] = delta;
        let global_d_j = r.transpose() * local_d.fixed_rows::<3>(0).into_owned();
        let global_rot_j = r.transpose() * local_d.fixed_rows::<3>(3).into_owned();

        let mut node_j_displaced = node_j.clone();
        node_j_displaced.displacement[0] = global_d_j[0];
        node_j_displaced.displacement[1] = global_d_j[1];
        node_j_displaced.displacement[2] = global_d_j[2];
        node_j_displaced.displacement[3] = global_rot_j[0];
        node_j_displaced.displacement[4] = global_rot_j[1];
        node_j_displaced.displacement[5] = global_rot_j[2];

        let (_, resistance) = beam.form_tangent_and_resistance(&node_i, &node_j_displaced);

        // Pure axial elongation should produce a pure axial force pair,
        // recoverable via the same projection `Truss3`'s skew test uses.
        let expected_axial_force = e * a * delta / length;
        let axial_disp_of_force = resistance.fixed_rows::<3>(6).into_owned().dot(&r.row(0).transpose());
        assert!((axial_disp_of_force - expected_axial_force).abs() < 1e-9);
    }
}
