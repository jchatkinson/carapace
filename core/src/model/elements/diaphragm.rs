//! Horizontal 3D diaphragm made from independent shear fields and axial grid members.

use arrayvec::ArrayVec;
use nalgebra::{SMatrix, SVector};

use super::stud_wall::shear_law_response;
use super::{DofMask, NodeView, TangentSink, VectorSink};
use crate::model::{ElementForce, Material, Node3Id, NodeKey, WallShearLaw};

/// Plan grid and properties for a horizontal diaphragm.
#[derive(Debug, Clone)]
pub struct DiaphragmLayout {
    /// Global X coordinates, strictly increasing.
    pub xs: Vec<f64>,
    /// Global Z coordinates, strictly increasing.
    pub zs: Vec<f64>,
    /// Global elevation (Y).
    pub y: f64,
    /// Shear law per bay, indexed `i + j * (nx - 1)`; `None` is an opening.
    pub bays: Vec<Option<WallShearLaw>>,
    /// Axial stiffnesses on constant-Z lines, indexed by Z line.
    pub ea_x: Vec<f64>,
    /// Axial stiffnesses on constant-X lines, indexed by X line.
    pub ea_z: Vec<f64>,
}

/// Recovered shear and axial response (tension positive).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DiaphragmResponse {
    /// Per-bay engineering shear strain and shear flow.
    pub bays: Vec<DiaphragmBayResponse>,
    /// Axial forces on X-directed segments, by Z line then X segment.
    pub axial_x: Vec<Vec<f64>>,
    /// Axial forces on Z-directed segments, by X line then Z segment.
    pub axial_z: Vec<Vec<f64>>,
    /// Mid-span chord force on each edge parallel to X (constant-Z boundary lines).
    pub chord_forces_x: [f64; 2],
    /// Maximum absolute boundary shear flow on edges parallel to X and Z.
    pub max_boundary_shear_flow: [f64; 4],
}

/// One bay's shear response.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct DiaphragmBayResponse {
    /// Engineering shear strain.
    pub gamma: f64,
    /// Shear flow.
    pub q: f64,
}

/// Horizontal diaphragm element. Nodes are ordered `i + j * nx` and carry only `ux, uz`.
#[derive(Debug, Clone)]
pub struct Diaphragm<K: NodeKey = Node3Id> {
    nodes: Vec<K>,
    layout: DiaphragmLayout,
    materials: Vec<Option<Material>>,
}

impl Diaphragm<Node3Id> {
    /// Builds a diaphragm; axis-aligned global X-Z grid with global Y up. Geometry is checked by `validate`.
    pub fn new(nodes: Vec<Node3Id>, layout: DiaphragmLayout) -> Result<Self, String> {
        validate_layout(&layout)?;
        let count = layout.bays.len();
        let materials = layout
            .bays
            .iter()
            .map(|law| match law {
                Some(WallShearLaw::Material(m)) => Some(m.clone()),
                _ => None,
            })
            .collect::<Vec<_>>();
        if nodes.len() != layout.xs.len() * layout.zs.len() {
            return Err("node count must equal nx * nz".into());
        }
        if nodes
            .iter()
            .copied()
            .collect::<std::collections::HashSet<_>>()
            .len()
            != nodes.len()
        {
            return Err("diaphragm nodes must be distinct".into());
        }
        debug_assert_eq!(materials.len(), count);
        Ok(Self {
            nodes,
            layout,
            materials,
        })
    }
}

impl<K: NodeKey> Diaphragm<K> {
    /// Connected grid nodes, ordered by `i + j * nx`.
    pub fn nodes(&self) -> &[K] {
        &self.nodes
    }
    /// Element layout.
    pub fn layout(&self) -> &DiaphragmLayout {
        &self.layout
    }
    /// Exact active global plan translations: `ux` and `uz` only.
    pub fn dof_mask(&self) -> DofMask {
        DofMask::none().with(0).with(2)
    }

    /// Checks axis-aligned grid coordinates against `(xs[i], y, zs[j])` at relative tolerance `1e-9`.
    pub fn validate<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> Result<(), &'static str> {
        let tol = self
            .layout
            .xs
            .iter()
            .chain(&self.layout.zs)
            .chain(std::iter::once(&self.layout.y))
            .fold(1.0_f64, |a, b| a.max(b.abs()))
            * 1e-9;
        for (j, z) in self.layout.zs.iter().enumerate() {
            for (i, x) in self.layout.xs.iter().enumerate() {
                let p = nodes.get(self.nodes[i + j * self.layout.xs.len()]).coords;
                if (p[0] - x)
                    .abs()
                    .max((p[1] - self.layout.y).abs())
                    .max((p[2] - z).abs())
                    > tol
                {
                    return Err(
                        "diaphragm nodes must lie on its axis-aligned X-Z grid at elevation Y",
                    );
                }
            }
        }
        Ok(())
    }

    pub fn prepare<const NDIM: usize, const NDOF: usize>(
        &mut self,
        _nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) {
    }
    pub fn accepts_load(&self, _load: &crate::model::ElementLoad3) -> bool {
        false
    }
    pub fn assemble_load<const NDIM: usize, const NDOF: usize, S: VectorSink<K>>(
        &self,
        _nodes: &NodeView<'_, NDIM, NDOF, K>,
        _load: Option<&crate::model::ElementLoad3>,
        _sink: &mut S,
    ) {
    }

    pub fn assemble_tangent<const NDIM: usize, const NDOF: usize, S: TangentSink<K>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        sink: &mut S,
    ) {
        let nx = self.layout.xs.len();
        for j in 0..self.layout.zs.len() - 1 {
            for i in 0..nx - 1 {
                let index = i + j * (nx - 1);
                let ids = [
                    self.nodes[i + j * nx],
                    self.nodes[i + 1 + j * nx],
                    self.nodes[i + (j + 1) * nx],
                    self.nodes[i + 1 + (j + 1) * nx],
                ];
                let dx = self.layout.xs[i + 1] - self.layout.xs[i];
                let dz = self.layout.zs[j + 1] - self.layout.zs[j];
                if self.layout.bays[index].is_some() {
                    let (_, q, tangent) = self.shear(index, ids, dx, dz, nodes);
                    let b = shear_b(dx, dz);
                    let scale = dx * dz * tangent;
                    let k = b * b.transpose() * scale;
                    let r = b * (dx * dz * q);
                    let mut dofs = ArrayVec::<(K, u8), 8>::new();
                    dofs.extend(ids.iter().flat_map(|id| [(*id, 0), (*id, 2)]));
                    let mut flat = ArrayVec::<f64, 64>::new();
                    for row in 0..8 {
                        for col in 0..8 {
                            flat.push(k[(row, col)]);
                        }
                    }
                    sink.add_dyn(&dofs, &flat, r.as_slice());
                }
            }
        }
        for j in 0..self.layout.zs.len() {
            for i in 0..nx - 1 {
                let a = self.nodes[i + j * nx];
                let b = self.nodes[i + 1 + j * nx];
                add_bar(
                    sink,
                    a,
                    b,
                    0,
                    self.layout.ea_x[j] / (self.layout.xs[i + 1] - self.layout.xs[i]),
                    nodes.get(a).displacement[0],
                    nodes.get(b).displacement[0],
                );
            }
        }
        for i in 0..nx {
            for j in 0..self.layout.zs.len() - 1 {
                let a = self.nodes[i + j * nx];
                let b = self.nodes[i + (j + 1) * nx];
                add_bar(
                    sink,
                    a,
                    b,
                    2,
                    self.layout.ea_z[i] / (self.layout.zs[j + 1] - self.layout.zs[j]),
                    nodes.get(a).displacement[2],
                    nodes.get(b).displacement[2],
                );
            }
        }
    }

    pub fn commit<const NDIM: usize, const NDOF: usize>(
        &mut self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        _load: Option<&crate::model::ElementLoad3>,
    ) {
        let nx = self.layout.xs.len();
        for j in 0..self.layout.zs.len() - 1 {
            for i in 0..nx - 1 {
                let index = i + j * (nx - 1);
                if self.materials[index].is_none() {
                    continue;
                }
                let ids = [
                    self.nodes[i + j * nx],
                    self.nodes[i + 1 + j * nx],
                    self.nodes[i + (j + 1) * nx],
                    self.nodes[i + 1 + (j + 1) * nx],
                ];
                let gamma = shear_b(
                    self.layout.xs[i + 1] - self.layout.xs[i],
                    self.layout.zs[j + 1] - self.layout.zs[j],
                )
                .dot(&displacements(ids, nodes));
                if let Some(m) = &mut self.materials[index] {
                    *m = m.commit(gamma);
                }
            }
        }
    }
    pub fn local_force_width(&self) -> usize {
        self.nodes.len() * 2
    }
    pub fn local_force<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> ElementForce {
        let mut out = vec![0.; self.local_force_width()];
        let nx = self.layout.xs.len();
        for j in 0..self.layout.zs.len() - 1 {
            for i in 0..nx - 1 {
                let ix = i + j * (nx - 1);
                if self.layout.bays[ix].is_none() {
                    continue;
                }
                let ids = [
                    self.nodes[i + j * nx],
                    self.nodes[i + 1 + j * nx],
                    self.nodes[i + (j + 1) * nx],
                    self.nodes[i + 1 + (j + 1) * nx],
                ];
                let dx = self.layout.xs[i + 1] - self.layout.xs[i];
                let dz = self.layout.zs[j + 1] - self.layout.zs[j];
                let (_, q, _) = self.shear(ix, ids, dx, dz, nodes);
                let f = shear_b(dx, dz) * (dx * dz * q);
                for (a, id) in ids.iter().enumerate() {
                    let n = self.nodes.iter().position(|x| x == id).unwrap();
                    out[2 * n] += f[2 * a];
                    out[2 * n + 1] += f[2 * a + 1];
                }
            }
        }
        for j in 0..self.layout.zs.len() {
            for i in 0..nx - 1 {
                let (a, b) = (i + j * nx, i + 1 + j * nx);
                let k = self.layout.ea_x[j] / (self.layout.xs[i + 1] - self.layout.xs[i]);
                let n = k
                    * (nodes.get(self.nodes[b]).displacement[0]
                        - nodes.get(self.nodes[a]).displacement[0]);
                out[2 * a] -= n;
                out[2 * b] += n;
            }
        }
        for i in 0..nx {
            for j in 0..self.layout.zs.len() - 1 {
                let (a, b) = (i + j * nx, i + (j + 1) * nx);
                let k = self.layout.ea_z[i] / (self.layout.zs[j + 1] - self.layout.zs[j]);
                let n = k
                    * (nodes.get(self.nodes[b]).displacement[2]
                        - nodes.get(self.nodes[a]).displacement[2]);
                out[2 * a + 1] -= n;
                out[2 * b + 1] += n;
            }
        }
        ElementForce::from_vec(out)
    }
    pub fn local_load_force<const NDIM: usize, const NDOF: usize>(
        &self,
        _nodes: &NodeView<'_, NDIM, NDOF, K>,
        _load: Option<&crate::model::ElementLoad3>,
    ) -> ElementForce {
        ElementForce::from_vec(vec![0.; self.local_force_width()])
    }
    pub fn response<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> DiaphragmResponse {
        let nx = self.layout.xs.len();
        let nz = self.layout.zs.len();
        let bays = (0..self.layout.bays.len())
            .map(|ix| {
                let i = ix % (nx - 1);
                let j = ix / (nx - 1);
                if self.layout.bays[ix].is_none() {
                    return DiaphragmBayResponse::default();
                }
                let ids = [
                    self.nodes[i + j * nx],
                    self.nodes[i + 1 + j * nx],
                    self.nodes[i + (j + 1) * nx],
                    self.nodes[i + 1 + (j + 1) * nx],
                ];
                let (g, q, _) = self.shear(
                    ix,
                    ids,
                    self.layout.xs[i + 1] - self.layout.xs[i],
                    self.layout.zs[j + 1] - self.layout.zs[j],
                    nodes,
                );
                DiaphragmBayResponse { gamma: g, q }
            })
            .collect();
        let axial_x = (0..nz)
            .map(|j| {
                (0..nx - 1)
                    .map(|i| {
                        let (a, b) = (self.nodes[i + j * nx], self.nodes[i + 1 + j * nx]);
                        self.layout.ea_x[j] / (self.layout.xs[i + 1] - self.layout.xs[i])
                            * (nodes.get(b).displacement[0] - nodes.get(a).displacement[0])
                    })
                    .collect()
            })
            .collect();
        let axial_z = (0..nx)
            .map(|i| {
                (0..nz - 1)
                    .map(|j| {
                        let (a, b) = (self.nodes[i + j * nx], self.nodes[i + (j + 1) * nx]);
                        self.layout.ea_z[i] / (self.layout.zs[j + 1] - self.layout.zs[j])
                            * (nodes.get(b).displacement[2] - nodes.get(a).displacement[2])
                    })
                    .collect()
            })
            .collect();
        let mut out = DiaphragmResponse {
            bays,
            axial_x,
            axial_z,
            ..Default::default()
        };
        for i in 0..nx - 1 {
            let q0 = out.bays[i].q.abs();
            let q1 = out.bays[i + (nz - 2) * (nx - 1)].q.abs();
            out.max_boundary_shear_flow[0] = out.max_boundary_shear_flow[0].max(q0);
            out.max_boundary_shear_flow[1] = out.max_boundary_shear_flow[1].max(q1);
        }
        for j in 0..nz - 1 {
            out.max_boundary_shear_flow[2] =
                out.max_boundary_shear_flow[2].max(out.bays[j * (nx - 1)].q.abs());
            out.max_boundary_shear_flow[3] =
                out.max_boundary_shear_flow[3].max(out.bays[j * (nx - 1) + nx - 2].q.abs());
        }
        let mid = (0..nx - 1)
            .min_by(|&a, &b| {
                let ca = (self.layout.xs[a] + self.layout.xs[a + 1]) * 0.5;
                let cb = (self.layout.xs[b] + self.layout.xs[b + 1]) * 0.5;
                (ca - (self.layout.xs[0] + self.layout.xs[nx - 1]) * 0.5)
                    .abs()
                    .total_cmp(&(cb - (self.layout.xs[0] + self.layout.xs[nx - 1]) * 0.5).abs())
            })
            .unwrap();
        out.chord_forces_x = [out.axial_x[0][mid], out.axial_x[nz - 1][mid]];
        out
    }
    fn shear<const NDIM: usize, const NDOF: usize>(
        &self,
        index: usize,
        ids: [K; 4],
        dx: f64,
        dz: f64,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> (f64, f64, f64) {
        let gamma = shear_b(dx, dz).dot(&displacements(ids, nodes));
        let law = self.layout.bays[index].as_ref().unwrap();
        let (q, tangent) = shear_law_response(law, self.materials[index].as_ref(), gamma);
        (gamma, q, tangent)
    }
}

fn displacements<const NDIM: usize, const NDOF: usize, K: NodeKey>(
    ids: [K; 4],
    nodes: &NodeView<'_, NDIM, NDOF, K>,
) -> SVector<f64, 8> {
    SVector::from_fn(|i, _| nodes.get(ids[i / 2]).displacement[if i % 2 == 0 { 0 } else { 2 }])
}
fn shear_b(dx: f64, dz: f64) -> SVector<f64, 8> {
    SVector::from([
        -0.5 / dz,
        -0.5 / dx,
        -0.5 / dz,
        0.5 / dx,
        0.5 / dz,
        -0.5 / dx,
        0.5 / dz,
        0.5 / dx,
    ])
}
fn add_bar<K: NodeKey, S: TangentSink<K>>(
    sink: &mut S,
    a: K,
    b: K,
    slot: u8,
    k: f64,
    ua: f64,
    ub: f64,
) {
    let mut dofs = ArrayVec::<(K, u8), 2>::new();
    dofs.extend([(a, slot), (b, slot)]);
    let m = SMatrix::<f64, 2, 2>::new(k, -k, -k, k);
    let r = SVector::<f64, 2>::new(k * (ua - ub), k * (ub - ua));
    sink.add_dyn(
        &dofs,
        &[m[(0, 0)], m[(0, 1)], m[(1, 0)], m[(1, 1)]],
        r.as_slice(),
    );
}
fn validate_layout(l: &DiaphragmLayout) -> Result<(), String> {
    let (nx, nz) = (l.xs.len(), l.zs.len());
    if nx < 2 || nz < 2 {
        return Err("xs and zs must each contain at least two coordinates".into());
    }
    if l.xs
        .iter()
        .chain(&l.zs)
        .chain(std::iter::once(&l.y))
        .any(|x| !x.is_finite())
        || l.xs.windows(2).chain(l.zs.windows(2)).any(|w| w[1] <= w[0])
    {
        return Err("grid coordinates must be finite and strictly increasing".into());
    }
    if l.bays.len() != (nx - 1) * (nz - 1) {
        return Err("bays length must equal (nx - 1) * (nz - 1)".into());
    }
    if l.ea_x.len() != nz
        || l.ea_z.len() != nx
        || l.ea_x
            .iter()
            .chain(&l.ea_z)
            .any(|x| !x.is_finite() || *x < 0.)
    {
        return Err("EA arrays must match grid lines and be finite nonnegative".into());
    }
    if l.ea_x
        .iter()
        .enumerate()
        .any(|(j, x)| j > 0 && j + 1 < nz && *x == 0.)
        || l.ea_z
            .iter()
            .enumerate()
            .any(|(i, x)| i > 0 && i + 1 < nx && *x == 0.)
    {
        return Err("interior grid lines require positive EA".into());
    }
    for law in l.bays.iter().flatten() {
        match law {
            WallShearLaw::Linear { ga } if !ga.is_finite() || *ga <= 0. => {
                return Err("linear ga must be finite and positive".into())
            }
            WallShearLaw::O86 {
                bv,
                nail_spacing,
                nail_diameter,
            } if [bv, nail_spacing, nail_diameter]
                .iter()
                .any(|x| !x.is_finite() || **x <= 0.) =>
            {
                return Err("O86 parameters must be finite and positive".into())
            }
            WallShearLaw::Material(m)
                if !m.initial_tangent().is_finite() || m.initial_tangent() <= 0. =>
            {
                return Err("material initial tangent must be finite and positive".into())
            }
            _ => {}
        }
    }
    Ok(())
}
