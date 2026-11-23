//! Horizontal 3D diaphragm made from independent shear fields and axial grid members.

use arrayvec::ArrayVec;
use nalgebra::{SMatrix, SVector};

use super::stud_wall::shear_law_response;
use super::{DofMask, NodeView, TangentSink, VectorSink};
use crate::model::{ElementForce, Material, Node3Id, NodeKey, WallShearLaw};

/// Grid coordinates follow the cyclic in-plane axes after its normal: `(n+1)%3`, `(n+2)%3`.
#[derive(Debug, Clone)]
pub struct DiaphragmLayout {
    /// First in-plane coordinate, strictly increasing.
    pub u: Vec<f64>,
    /// Second in-plane coordinate, strictly increasing.
    pub v: Vec<f64>,
    /// Shared coordinate along the normal axis.
    pub level: f64,
    /// Shear law per bay, indexed `i + j * (nx - 1)`; `None` is an opening.
    pub bays: Vec<Option<WallShearLaw>>,
    /// Axial stiffnesses on constant-v lines, indexed by v line.
    pub ea_u: Vec<f64>,
    /// Axial stiffnesses on constant-u lines, indexed by u line.
    pub ea_v: Vec<f64>,
}

/// Recovered shear and axial response (tension positive).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DiaphragmResponse {
    /// Per-bay engineering shear strain and shear flow.
    pub bays: Vec<DiaphragmBayResponse>,
    /// Axial forces on X-directed segments, by v line then u segment.
    pub axial_u: Vec<Vec<f64>>,
    /// Axial forces on Z-directed segments, by u line then v segment.
    pub axial_v: Vec<Vec<f64>>,
    /// Mid-span chord force on each edge parallel to u (constant-v boundary lines).
    pub chord_forces_u: [f64; 2],
    /// Maximum absolute boundary shear flow on edges parallel to u and v.
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

/// Axis-aligned diaphragm element. Nodes are ordered `i + j * nu` and carry the two in-plane translations.
#[derive(Debug, Clone)]
pub struct Diaphragm<K: NodeKey = Node3Id> {
    nodes: Vec<K>,
    layout: DiaphragmLayout,
    materials: Vec<Option<Material>>,
}

impl Diaphragm<Node3Id> {
    /// Builds a diaphragm whose axis-aligned plane is derived from its nodes. Geometry is checked by `validate`.
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
        if nodes.len() != layout.u.len() * layout.v.len() {
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
    pub fn dof_mask<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> DofMask {
        let (a, b) = self.axes(nodes).map(|(_, a, b)| (a, b)).unwrap_or((0, 2));
        DofMask::none().with(a).with(b)
    }

    fn axes<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> Result<(usize, usize, usize), &'static str> {
        let first = nodes.get(self.nodes[0]).coords;
        let scale = self.nodes.iter().fold(1.0_f64, |m, id| {
            nodes.get(*id).coords.iter().fold(m, |s, x| s.max(x.abs()))
        });
        let tol = scale * 1e-9;
        let normal = (0..3)
            .find(|&axis| {
                self.nodes
                    .iter()
                    .all(|id| (nodes.get(*id).coords[axis] - first[axis]).abs() <= tol)
            })
            .ok_or("diaphragm grid nodes do not share a global X, Y or Z coordinate")?;
        Ok((normal, (normal + 1) % 3, (normal + 2) % 3))
    }

    /// Checks the grid against `(u[i], v[j], level)` on its geometry-derived axis-aligned plane.
    pub fn validate<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> Result<(), &'static str> {
        let (n, a, b) = self.axes(nodes)?;
        let tol = self
            .layout
            .u
            .iter()
            .chain(&self.layout.v)
            .chain(std::iter::once(&self.layout.level))
            .fold(1.0_f64, |m, x| m.max(x.abs()))
            * 1e-9;
        for (j, v) in self.layout.v.iter().enumerate() {
            for (i, u) in self.layout.u.iter().enumerate() {
                let p = nodes.get(self.nodes[i + j * self.layout.u.len()]).coords;
                if (p[a] - u)
                    .abs()
                    .max((p[b] - v).abs())
                    .max((p[n] - self.layout.level).abs())
                    > tol
                {
                    return Err("diaphragm nodes must lie on the axis-aligned grid described by u, v and level");
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
        let (_, axis_u, axis_v) = self.axes(nodes).expect("validated diaphragm geometry");
        let nx = self.layout.u.len();
        for j in 0..self.layout.v.len() - 1 {
            for i in 0..nx - 1 {
                let index = i + j * (nx - 1);
                let ids = [
                    self.nodes[i + j * nx],
                    self.nodes[i + 1 + j * nx],
                    self.nodes[i + (j + 1) * nx],
                    self.nodes[i + 1 + (j + 1) * nx],
                ];
                let dx = self.layout.u[i + 1] - self.layout.u[i];
                let dz = self.layout.v[j + 1] - self.layout.v[j];
                if self.layout.bays[index].is_some() {
                    let (_, q, tangent) = self.shear(index, ids, [dx, dz], [axis_u, axis_v], nodes);
                    let b = shear_b(dx, dz);
                    let scale = dx * dz * tangent;
                    let k = b * b.transpose() * scale;
                    let r = b * (dx * dz * q);
                    let mut dofs = ArrayVec::<(K, u8), 8>::new();
                    dofs.extend(
                        ids.iter()
                            .flat_map(|id| [(*id, axis_u as u8), (*id, axis_v as u8)]),
                    );
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
        for j in 0..self.layout.v.len() {
            for i in 0..nx - 1 {
                let a = self.nodes[i + j * nx];
                let b = self.nodes[i + 1 + j * nx];
                add_bar(
                    sink,
                    a,
                    b,
                    axis_u as u8,
                    self.layout.ea_u[j] / (self.layout.u[i + 1] - self.layout.u[i]),
                    nodes.get(a).displacement[axis_u],
                    nodes.get(b).displacement[axis_u],
                );
            }
        }
        for i in 0..nx {
            for j in 0..self.layout.v.len() - 1 {
                let a = self.nodes[i + j * nx];
                let b = self.nodes[i + (j + 1) * nx];
                add_bar(
                    sink,
                    a,
                    b,
                    axis_v as u8,
                    self.layout.ea_v[i] / (self.layout.v[j + 1] - self.layout.v[j]),
                    nodes.get(a).displacement[axis_v],
                    nodes.get(b).displacement[axis_v],
                );
            }
        }
    }

    pub fn commit<const NDIM: usize, const NDOF: usize>(
        &mut self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        _load: Option<&crate::model::ElementLoad3>,
    ) {
        let (_, axis_u, axis_v) = self.axes(nodes).expect("validated diaphragm geometry");
        let nx = self.layout.u.len();
        for j in 0..self.layout.v.len() - 1 {
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
                    self.layout.u[i + 1] - self.layout.u[i],
                    self.layout.v[j + 1] - self.layout.v[j],
                )
                .dot(&displacements(ids, axis_u, axis_v, nodes));
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
        let (_, axis_u, axis_v) = self.axes(nodes).expect("validated diaphragm geometry");
        let mut out = vec![0.; self.local_force_width()];
        let nx = self.layout.u.len();
        for j in 0..self.layout.v.len() - 1 {
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
                let dx = self.layout.u[i + 1] - self.layout.u[i];
                let dz = self.layout.v[j + 1] - self.layout.v[j];
                let (_, q, _) = self.shear(ix, ids, [dx, dz], [axis_u, axis_v], nodes);
                let f = shear_b(dx, dz) * (dx * dz * q);
                for (a, id) in ids.iter().enumerate() {
                    let n = self.nodes.iter().position(|x| x == id).unwrap();
                    out[2 * n] += f[2 * a];
                    out[2 * n + 1] += f[2 * a + 1];
                }
            }
        }
        for j in 0..self.layout.v.len() {
            for i in 0..nx - 1 {
                let (a, b) = (i + j * nx, i + 1 + j * nx);
                let k = self.layout.ea_u[j] / (self.layout.u[i + 1] - self.layout.u[i]);
                let n = k
                    * (nodes.get(self.nodes[b]).displacement[axis_u]
                        - nodes.get(self.nodes[a]).displacement[axis_u]);
                out[2 * a] -= n;
                out[2 * b] += n;
            }
        }
        for i in 0..nx {
            for j in 0..self.layout.v.len() - 1 {
                let (a, b) = (i + j * nx, i + (j + 1) * nx);
                let k = self.layout.ea_v[i] / (self.layout.v[j + 1] - self.layout.v[j]);
                let n = k
                    * (nodes.get(self.nodes[b]).displacement[axis_v]
                        - nodes.get(self.nodes[a]).displacement[axis_v]);
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
        let (_, axis_u, axis_v) = self.axes(nodes).expect("validated diaphragm geometry");
        let nx = self.layout.u.len();
        let nz = self.layout.v.len();
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
                    [
                        self.layout.u[i + 1] - self.layout.u[i],
                        self.layout.v[j + 1] - self.layout.v[j],
                    ],
                    [axis_u, axis_v],
                    nodes,
                );
                DiaphragmBayResponse { gamma: g, q }
            })
            .collect();
        let axial_u = (0..nz)
            .map(|j| {
                (0..nx - 1)
                    .map(|i| {
                        let (a, b) = (self.nodes[i + j * nx], self.nodes[i + 1 + j * nx]);
                        self.layout.ea_u[j] / (self.layout.u[i + 1] - self.layout.u[i])
                            * (nodes.get(b).displacement[axis_u]
                                - nodes.get(a).displacement[axis_u])
                    })
                    .collect()
            })
            .collect();
        let axial_v = (0..nx)
            .map(|i| {
                (0..nz - 1)
                    .map(|j| {
                        let (a, b) = (self.nodes[i + j * nx], self.nodes[i + (j + 1) * nx]);
                        self.layout.ea_v[i] / (self.layout.v[j + 1] - self.layout.v[j])
                            * (nodes.get(b).displacement[axis_v]
                                - nodes.get(a).displacement[axis_v])
                    })
                    .collect()
            })
            .collect();
        let mut out = DiaphragmResponse {
            bays,
            axial_u,
            axial_v,
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
                let ca = (self.layout.u[a] + self.layout.u[a + 1]) * 0.5;
                let cb = (self.layout.u[b] + self.layout.u[b + 1]) * 0.5;
                (ca - (self.layout.u[0] + self.layout.u[nx - 1]) * 0.5)
                    .abs()
                    .total_cmp(&(cb - (self.layout.u[0] + self.layout.u[nx - 1]) * 0.5).abs())
            })
            .unwrap();
        out.chord_forces_u = [out.axial_u[0][mid], out.axial_u[nz - 1][mid]];
        out
    }
    fn shear<const NDIM: usize, const NDOF: usize>(
        &self,
        index: usize,
        ids: [K; 4],
        span: [f64; 2],
        axes: [usize; 2],
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> (f64, f64, f64) {
        let gamma = shear_b(span[0], span[1]).dot(&displacements(ids, axes[0], axes[1], nodes));
        let law = self.layout.bays[index].as_ref().unwrap();
        let (q, tangent) = shear_law_response(law, self.materials[index].as_ref(), gamma);
        (gamma, q, tangent)
    }
}

fn displacements<const NDIM: usize, const NDOF: usize, K: NodeKey>(
    ids: [K; 4],
    a: usize,
    b: usize,
    nodes: &NodeView<'_, NDIM, NDOF, K>,
) -> SVector<f64, 8> {
    SVector::from_fn(|i, _| nodes.get(ids[i / 2]).displacement[if i % 2 == 0 { a } else { b }])
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
    let (nx, nz) = (l.u.len(), l.v.len());
    if nx < 2 || nz < 2 {
        return Err("u and v must each contain at least two coordinates".into());
    }
    if l.u
        .iter()
        .chain(&l.v)
        .chain(std::iter::once(&l.level))
        .any(|x| !x.is_finite())
        || l.u.windows(2).chain(l.v.windows(2)).any(|w| w[1] <= w[0])
    {
        return Err("grid coordinates must be finite and strictly increasing".into());
    }
    if l.bays.len() != (nx - 1) * (nz - 1) {
        return Err("bays length must equal (nx - 1) * (nz - 1)".into());
    }
    if l.ea_u.len() != nz
        || l.ea_v.len() != nx
        || l.ea_u
            .iter()
            .chain(&l.ea_v)
            .any(|x| !x.is_finite() || *x < 0.)
    {
        return Err("EA arrays must match grid lines and be finite nonnegative".into());
    }
    if l.ea_u
        .iter()
        .enumerate()
        .any(|(j, x)| j > 0 && j + 1 < nz && *x == 0.)
        || l.ea_v
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
