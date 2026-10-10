//! `Shell3`: the 3-node flat shell (OpenSees `ShellDKGT`): a DKT plate-bending triangle (Kirchhoff, so no transverse
//! shear) combined with an Allman-type membrane triangle that carries the drilling rotation, six DOFs per node
//! (18 total), the 4-point rule `(1/3,1/3,1/3)` with weight `-9/16` and `(1/5,3/5,1/5)` permutations with `25/48`.
//!
//! Like [`Shell4`](super::Shell4) it reproduces OpenSees term by term: the local frame is `e1 = x1 - x0`, `e2` its
//! Gram-Schmidt completion with `x2 - x0`, `e3 = e1 x e2`; the section ordering, sign convention and bending-row
//! negation are the same, and the mass is a row-sum lumped translational mass (OpenSees' is consistent). The DKGT has
//! no drilling penalty and no shear strains: its `gamma_xz`, `gamma_yz` rows are zero.

use nalgebra::{SMatrix, SVector, Vector3};

use super::{DofMask, DofRef, ElementForce, NodeView, ShellResponse, TangentSink, VectorSink};
use crate::model::continuum::validate_triangle;
use crate::model::{ElementLoad3, Node3Id, ShellSection, ShellVector, SPATIAL_NDF, SPATIAL_NDIM};

type View<'a> = NodeView<'a, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>;
type Vec18 = SVector<f64, 18>;
type Mat18 = SMatrix<f64, 18, 18>;
type Bmat = SMatrix<f64, 8, 18>;

/// Area coordinates `(L1, L2, L3)` of the four Gauss points and their weights.
const GAUSS: [([f64; 3], f64); 4] = [
    ([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], -9.0 / 16.0),
    ([0.2, 0.6, 0.2], 25.0 / 48.0),
    ([0.6, 0.2, 0.2], 25.0 / 48.0),
    ([0.2, 0.2, 0.6], 25.0 / 48.0),
];

#[derive(Debug, Clone)]
pub struct Shell3 {
    pub nodes: [Node3Id; 3],
    sections: [ShellSection; 4],
    cache: Option<Cache>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Cache {
    /// Local axes `e1`, `e2`, `e3` in global components.
    axes: [Vector3<f64>; 3],
    points: [PointGeometry; 4],
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct PointGeometry {
    /// Strain-displacement matrix of the generalized strains, in global DOFs.
    b: Bmat,
    /// The same with the bending-bending block negated, for the internal force.
    b_force: Bmat,
    /// `0.5 w xsj` (`w` summing to 1, so the weights sum to the area).
    dvol: f64,
    /// Area coordinates.
    n: [f64; 3],
}

type Frame = ([Vector3<f64>; 3], [[f64; 3]; 2]);

fn frame(coords: &[[f64; 3]; 3]) -> Result<Frame, &'static str> {
    let x: [Vector3<f64>; 3] = std::array::from_fn(|i| Vector3::from(coords[i]));
    let mut v1 = x[1] - x[0];
    let mut v2 = x[2] - x[0];
    let n1 = v1.norm();
    if !n1.is_finite() || n1 == 0.0 {
        return Err("degenerate shell geometry");
    }
    v1 /= n1;
    v2 -= v2.dot(&v1) * v1;
    let n2 = v2.norm();
    if !n2.is_finite() || n2 == 0.0 {
        return Err("degenerate shell geometry");
    }
    v2 /= n2;
    let v3 = v1.cross(&v2);
    Ok((
        [v1, v2, v3],
        [
            std::array::from_fn(|i| x[i].dot(&v1)),
            std::array::from_fn(|i| x[i].dot(&v2)),
        ],
    ))
}

/// The DKT bending shape-function derivatives, OpenSees' `shapeBend`: rows `Hx,x`, `Hx,y`, `Hy,x`, `Hy,y`
/// over the nine DOFs `(w, theta1, theta2)` per node.
fn shape_bend(l: [f64; 3], x: &[[f64; 3]; 2], area: f64) -> [[f64; 9]; 4] {
    let (x12, x23, x31) = (x[0][0] - x[0][1], x[0][1] - x[0][2], x[0][2] - x[0][0]);
    let (y12, y23, y31) = (x[1][0] - x[1][1], x[1][1] - x[1][2], x[1][2] - x[1][0]);
    let (l12, l23, l31) = (
        x12 * x12 + y12 * y12,
        x23 * x23 + y23 * y23,
        x31 * x31 + y31 * y31,
    );
    // Per mid-side 4 (23), 5 (31), 6 (12): a, b, c, d, e.
    let side = |xs: f64, ys: f64, len: f64| {
        (
            -xs / len,
            0.25 * 3.0 * xs * ys / len,
            0.25 * (xs * xs - 2.0 * ys * ys) / len,
            -ys / len,
            0.25 * (ys * ys - 2.0 * xs * xs) / len,
        )
    };
    let (a4, b4, _c4, d4, _e4) = side(x23, y23, l23);
    let (a5, b5, _c5, d5, _e5) = side(x31, y31, l31);
    let (a6, b6, _c6, d6, _e6) = side(x12, y12, l12);
    let (p4, p5, p6) = (6.0 * a4, 6.0 * a5, 6.0 * a6);
    let (t4, t5, t6) = (6.0 * d4, 6.0 * d5, 6.0 * d6);
    let (q4, q5, q6) = (4.0 * b4, 4.0 * b5, 4.0 * b6);
    let (r4, r5, r6) = (
        3.0 * y23 * y23 / l23,
        3.0 * y31 * y31 / l31,
        3.0 * y12 * y12 / l12,
    );

    let (l1, l2) = (l[1], l[2]);
    let hx_xi = [
        p6 * (1.0 - 2.0 * l1) + (p5 - p6) * l2,
        q6 * (1.0 - 2.0 * l1) - (q5 + q6) * l2,
        -4.0 + 6.0 * (l1 + l2) + r6 * (1.0 - 2.0 * l1) - l2 * (r5 + r6),
        -p6 * (1.0 - 2.0 * l1) + l2 * (p4 + p6),
        q6 * (1.0 - 2.0 * l1) - l2 * (q6 - q4),
        -2.0 + 6.0 * l1 + r6 * (1.0 - 2.0 * l1) + l2 * (r4 - r6),
        -l2 * (p5 + p4),
        l2 * (q4 - q5),
        -l2 * (r5 - r4),
    ];
    let hx_eta = [
        -p5 * (1.0 - 2.0 * l2) - l1 * (p6 - p5),
        q5 * (1.0 - 2.0 * l2) - l1 * (q5 + q6),
        -4.0 + 6.0 * (l1 + l2) + r5 * (1.0 - 2.0 * l2) - l1 * (r5 + r6),
        l1 * (p4 + p6),
        l1 * (q4 - q6),
        -l1 * (r6 - r4),
        p5 * (1.0 - 2.0 * l2) - l1 * (p4 + p5),
        q5 * (1.0 - 2.0 * l2) + l1 * (q4 - q5),
        -2.0 + 6.0 * l2 + r5 * (1.0 - 2.0 * l2) + l1 * (r4 - r5),
    ];
    let hy_xi = [
        t6 * (1.0 - 2.0 * l1) + l2 * (t5 - t6),
        1.0 + r6 * (1.0 - 2.0 * l1) - l2 * (r5 + r6),
        -q6 * (1.0 - 2.0 * l1) + l2 * (q5 + q6),
        -t6 * (1.0 - 2.0 * l1) + l2 * (t4 + t6),
        -1.0 + r6 * (1.0 - 2.0 * l1) + l2 * (r4 - r6),
        -q6 * (1.0 - 2.0 * l1) - l2 * (q4 - q6),
        -l2 * (t4 + t5),
        l2 * (r4 - r5),
        -l2 * (q4 - q5),
    ];
    let hy_eta = [
        -t5 * (1.0 - 2.0 * l2) - l1 * (t6 - t5),
        1.0 + r5 * (1.0 - 2.0 * l2) - l1 * (r5 + r6),
        -q5 * (1.0 - 2.0 * l2) + l1 * (q5 + q6),
        l1 * (t4 + t6),
        l1 * (r4 - r6),
        -l1 * (q4 - q6),
        t5 * (1.0 - 2.0 * l2) - l1 * (t4 + t5),
        -1.0 + r5 * (1.0 - 2.0 * l2) + l1 * (r4 - r5),
        -q5 * (1.0 - 2.0 * l2) - l1 * (q4 - q5),
    ];
    let s = 1.0 / (2.0 * area);
    let mut out = [[0.0; 9]; 4];
    for i in 0..9 {
        out[0][i] = (hx_xi[i] * y31 + hx_eta[i] * y12) * s;
        out[1][i] = (hx_xi[i] * -x31 + hx_eta[i] * -x12) * s;
        out[2][i] = (hy_xi[i] * y31 + hy_eta[i] * y12) * s;
        out[3][i] = (hy_xi[i] * -x31 + hy_eta[i] * -x12) * s;
    }
    out
}

fn point_geometry(axes: &[Vector3<f64>; 3], x: &[[f64; 3]; 2]) -> [PointGeometry; 4] {
    let area = 0.5
        * (x[0][0] * x[1][1] + x[0][1] * x[1][2] + x[0][2] * x[1][0]
            - x[0][0] * x[1][2]
            - x[0][1] * x[1][0]
            - x[0][2] * x[1][1]);
    let b: [f64; 3] = [x[1][1] - x[1][2], x[1][2] - x[1][0], x[1][0] - x[1][1]];
    let c: [f64; 3] = [-x[0][1] + x[0][2], -x[0][2] + x[0][0], -x[0][0] + x[0][1]];
    let xsj = (x[0][1] - x[0][0]) * (x[1][2] - x[1][0]) - (x[1][1] - x[1][0]) * (x[0][2] - x[0][0]);

    std::array::from_fn(|g| {
        let (l, w) = GAUSS[g];
        let bend = shape_bend(l, x, area);
        let mut local = Bmat::zeros();
        for j in 0..3 {
            let (nx, ny) = (b[j] / (2.0 * area), c[j] / (2.0 * area));
            let d = drill_row(j, l, &b, &c, area);
            let col = 6 * j;
            // Membrane rows: [N,x 0 Nu,x; 0 N,y Nv,y; N,y N,x Nv,x + Nu,y] on (u, v, theta_z).
            local[(0, col)] = nx;
            local[(0, col + 5)] = d[0];
            local[(1, col + 1)] = ny;
            local[(1, col + 5)] = d[3];
            local[(2, col)] = ny;
            local[(2, col + 1)] = nx;
            local[(2, col + 5)] = d[1] + d[2];
            // Bending rows on (w, theta_x, theta_y) = local columns 2, 3, 4; OpenSees negates Bbend.
            for k in 0..3 {
                local[(3, col + 2 + k)] = -bend[0][3 * j + k];
                local[(4, col + 2 + k)] = -bend[3][3 * j + k];
                local[(5, col + 2 + k)] = -(bend[1][3 * j + k] + bend[2][3 * j + k]);
            }
        }
        // Local -> global DOFs: each node's translations and rotations rotate with the axes.
        let mut t = SMatrix::<f64, 18, 18>::zeros();
        for j in 0..3 {
            for block in 0..2 {
                for (r, axis) in axes.iter().enumerate() {
                    for k in 0..3 {
                        t[(6 * j + 3 * block + r, 6 * j + 3 * block + k)] = axis[k];
                    }
                }
            }
        }
        let b_global = local * t;
        let mut b_force = b_global;
        // Negate the bending-bending block of the *local* B before rotating: rows 3..6, local columns 2..5.
        let mut local_force = local;
        for j in 0..3 {
            for r in 3..6 {
                for k in 2..5 {
                    local_force[(r, 6 * j + k)] *= -1.0;
                }
            }
        }
        b_force.copy_from(&(local_force * t));
        PointGeometry {
            b: b_global,
            b_force,
            dvol: 0.5 * w * xsj,
            n: l,
        }
    })
}

/// `[shpDrill[0..4][i]]` of OpenSees at area coordinates `l`.
fn drill_row(i: usize, l: [f64; 3], b: &[f64; 3], c: &[f64; 3], area: f64) -> [f64; 4] {
    let s = 1.0 / (4.0 * area);
    match i {
        0 => {
            let pb = b[2] * l[1] - b[1] * l[2];
            let pc = c[2] * l[1] - c[1] * l[2];
            [b[0] * pb * s, c[0] * pb * s, b[0] * pc * s, c[0] * pc * s]
        }
        1 => {
            let pb = b[0] * l[2] - b[2] * l[0];
            let pc = c[0] * l[2] - c[2] * l[0];
            [b[1] * pb * s, c[1] * pb * s, b[1] * pc * s, c[1] * pc * s]
        }
        _ => {
            let pb = b[1] * l[0] - b[0] * l[1];
            let pc = c[1] * l[0] - c[0] * l[1];
            [b[2] * pb * s, c[2] * pb * s, b[2] * pc * s, c[2] * pc * s]
        }
    }
}

fn dofs(nodes: &[Node3Id; 3]) -> [DofRef<Node3Id>; 18] {
    std::array::from_fn(|a| (nodes[a / 6], (a % 6) as u8))
}

impl Shell3 {
    pub fn new(nodes: [Node3Id; 3], section: ShellSection) -> Self {
        Shell3 {
            nodes,
            sections: std::array::from_fn(|_| section.clone()),
            cache: None,
        }
    }

    /// The Gauss points' sections, in OpenSees' Gauss-point order.
    pub fn sections(&self) -> &[ShellSection; 4] {
        &self.sections
    }

    fn coords(&self, view: &View<'_>) -> [[f64; 3]; 3] {
        std::array::from_fn(|i| view.get(self.nodes[i]).coords)
    }

    fn compute_cache(&self, view: &View<'_>) -> Result<Cache, &'static str> {
        let (axes, xl) = frame(&self.coords(view))?;
        Ok(Cache {
            axes,
            points: point_geometry(&axes, &xl),
        })
    }

    fn cache_or_compute(&self, view: &View<'_>) -> Cache {
        self.cache.unwrap_or_else(|| {
            self.compute_cache(view)
                .expect("unvalidated Shell3 geometry")
        })
    }

    pub(super) fn dof_mask(&self) -> DofMask {
        DofMask::all(SPATIAL_NDF)
    }

    pub(super) fn validate(&self, view: &View<'_>) -> Result<(), &'static str> {
        self.sections[0].validate().map_err(|e| e.message())?;
        let (_, x) = frame(&self.coords(view))?;
        validate_triangle(&[[x[0][0], x[1][0]], [x[0][1], x[1][1]], [x[0][2], x[1][2]]])
    }

    pub(super) fn prepare(&mut self, view: &View<'_>) {
        self.cache = self.compute_cache(view).ok();
    }

    fn displacement(&self, view: &View<'_>) -> Vec18 {
        Vec18::from_fn(|a, _| view.get(self.nodes[a / 6]).displacement[a % 6])
    }

    pub(super) fn assemble_tangent<S: TangentSink<Node3Id>>(&self, view: &View<'_>, sink: &mut S) {
        let cache = self.cache_or_compute(view);
        let u = self.displacement(view);
        let mut k = Mat18::zeros();
        let mut r = Vec18::zeros();
        for (g, section) in cache.points.iter().zip(&self.sections) {
            let (stress, d) = section.trial_resultant(&(g.b * u));
            k += g.dvol * g.b_force.transpose() * d * g.b;
            r += g.dvol * g.b_force.transpose() * stress;
        }
        sink.add(&dofs(&self.nodes), &k, &r);
    }

    /// `integral(N_i dA)` per node (`area / 3`, exactly): the pressure and body-force share and the lumped mass.
    fn nodal_weights(&self, cache: &Cache) -> [f64; 3] {
        let mut w = [0.0; 3];
        for g in &cache.points {
            for (i, wi) in w.iter_mut().enumerate() {
                *wi += g.n[i] * g.dvol;
            }
        }
        w
    }

    pub(super) fn assemble_mass<S: VectorSink<Node3Id>>(&self, view: &View<'_>, sink: &mut S) {
        let w = self.nodal_weights(&self.cache_or_compute(view));
        let rho_h = self.sections[0].rho_h();
        let m = Vec18::from_fn(|a, _| if a % 6 < 3 { rho_h * w[a / 6] } else { 0.0 });
        sink.add(&dofs(&self.nodes), &m);
    }

    /// Equivalent nodal forces of `load` in global axes (body acceleration times `rho h`, pressure along `e3`).
    pub(super) fn load_vector(&self, view: &View<'_>, load: &ElementLoad3) -> Vec18 {
        let cache = self.cache_or_compute(view);
        let w = self.nodal_weights(&cache);
        let rho_h = self.sections[0].rho_h();
        let normal = cache.axes[2];
        Vec18::from_fn(|a, _| {
            let (node, slot) = (a / 6, a % 6);
            if slot < 3 {
                w[node] * (rho_h * load.body[slot] + load.pressure * normal[slot])
            } else {
                0.0
            }
        })
    }

    pub(super) fn assemble_load<S: VectorSink<Node3Id>>(
        &self,
        view: &View<'_>,
        load: &ElementLoad3,
        sink: &mut S,
    ) {
        sink.add(&dofs(&self.nodes), &self.load_vector(view, load));
    }

    pub(super) fn commit(&mut self, view: &View<'_>) {
        let cache = self.cache_or_compute(view);
        let u = self.displacement(view);
        for (section, g) in self.sections.iter_mut().zip(&cache.points) {
            let strain: ShellVector = g.b * u;
            *section = section.commit(&strain);
        }
    }

    /// The element's resisting force in global axes (18).
    pub(super) fn local_force(&self, view: &View<'_>) -> ElementForce {
        let cache = self.cache_or_compute(view);
        let u = self.displacement(view);
        let mut r = Vec18::zeros();
        for (g, section) in cache.points.iter().zip(&self.sections) {
            let (stress, _) = section.trial_resultant(&(g.b * u));
            r += g.dvol * g.b_force.transpose() * stress;
        }
        r.into()
    }

    /// Generalized strain and resultants at each Gauss point (local axes, section ordering).
    pub(super) fn shell_responses(&self, view: &View<'_>) -> Vec<ShellResponse> {
        let cache = self.cache_or_compute(view);
        let u = self.displacement(view);
        cache
            .points
            .iter()
            .zip(&self.sections)
            .map(|(g, section)| {
                let strain = g.b * u;
                ShellResponse {
                    strain: strain.into(),
                    resultant: section.trial_resultant(&strain).0.into(),
                }
            })
            .collect()
    }
}
