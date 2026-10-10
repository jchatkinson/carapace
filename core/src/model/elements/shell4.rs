//! `Shell4`: the 4-node MITC4 flat shell (OpenSees `ShellMITC4`), small-displacement, 2x2 Gauss
//! integration, six DOFs per node (24 total).
//!
//! The formulation follows OpenSees so results agree to solver precision:
//!
//! * Local frame: `e1` from the edge midpoints (`0.5 (x2 + x1 - x3 - x0)`), `e2` its
//!   Gram-Schmidt completion with `0.5 (x3 + x2 - x1 - x0)`, `e3 = e1 x e2`. The shell is
//!   treated as flat in that frame (nodes are projected on `e1`, `e2`); warped quads are
//!   tolerated the same way.
//! * Membrane: bilinear displacement field. Bending: Reissner-Mindlin rotations, with the
//!   transverse shear strains interpolated from the four edge midpoints (MITC4) to avoid
//!   shear locking. Drilling: a penalty on the difference between the in-plane rotation
//!   and the skew part of the displacement gradient, stiffness the smallest eigenvalue of the
//!   section's membrane tangent.
//! * Mass: OpenSees uses a consistent translational mass with no rotational inertia; here
//!   `assemble_mass` is a row-sum lumped translational mass (same total, same centroid).
//!
//! Generalized strains and resultants use the section's ordering and sign convention
//! (`ShellSection`): `[eps_x, eps_y, gamma_xy, kappa_x, kappa_y, 2 kappa_xy, gamma_xz,
//! gamma_yz]`. The bending rows of the `B` matrix used for the internal force are negated
//! because the section's bending tangent is negative (the OpenSees convention).

use nalgebra::{Matrix2, SMatrix, SVector, Vector3};

use super::{DofMask, DofRef, ElementForce, NodeView, ShellResponse, TangentSink, VectorSink};
use crate::model::continuum::validate_quad;
use crate::model::{ElementLoad3, Node3Id, ShellSection, ShellVector, SPATIAL_NDF, SPATIAL_NDIM};

type View<'a> = NodeView<'a, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>;
type Vec24 = SVector<f64, 24>;
type Mat24 = SMatrix<f64, 24, 24>;
type Bmat = SMatrix<f64, 8, 24>;

const ONE_OVER_ROOT3: f64 = 0.577_350_269_189_625_7;
/// Gauss points in OpenSees' order: `(-,-), (+,-), (+,+), (-,+)`.
const SG: [f64; 4] = [
    -ONE_OVER_ROOT3,
    ONE_OVER_ROOT3,
    ONE_OVER_ROOT3,
    -ONE_OVER_ROOT3,
];
const TG: [f64; 4] = [
    -ONE_OVER_ROOT3,
    -ONE_OVER_ROOT3,
    ONE_OVER_ROOT3,
    ONE_OVER_ROOT3,
];

#[derive(Debug, Clone)]
pub struct Shell4 {
    pub nodes: [Node3Id; 4],
    sections: [ShellSection; 4],
    cache: Option<Cache>,
}

/// Geometry-only quantities, computed once by `prepare`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Cache {
    /// Local axes `e1`, `e2`, `e3` in global components.
    axes: [Vector3<f64>; 3],
    points: [PointGeometry; 4],
    /// Drilling stiffness penalty.
    ktt: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct PointGeometry {
    /// Strain-displacement matrix of the generalized strains.
    b: Bmat,
    /// The same with the bending-bending block negated, for the internal force.
    b_force: Bmat,
    /// Drilling strain-displacement row.
    drill: Vec24,
    /// `det J * w`.
    dvol: f64,
    /// Shape-function values.
    n: [f64; 4],
}

type Frame = ([Vector3<f64>; 3], [[f64; 2]; 4]);

/// The local frame and the nodes' coordinates in it.
fn frame(coords: &[[f64; 3]; 4]) -> Result<Frame, &'static str> {
    let x: [Vector3<f64>; 4] = std::array::from_fn(|i| Vector3::from(coords[i]));
    let mut v1 = 0.5 * (x[2] + x[1] - x[3] - x[0]);
    let mut v2 = 0.5 * (x[3] + x[2] - x[1] - x[0]);
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
        std::array::from_fn(|i| [x[i].dot(&v1), x[i].dot(&v2)]),
    ))
}

/// OpenSees' `shape2d`: `(dN/dx, dN/dy, N)` per node at `(s, t)` and the Jacobian determinant.
fn shape2d(s: f64, t: f64, xl: &[[f64; 2]; 4]) -> ([[f64; 4]; 3], f64) {
    const SS: [f64; 4] = [-0.5, 0.5, 0.5, -0.5];
    const TT: [f64; 4] = [-0.5, -0.5, 0.5, 0.5];
    let mut shp = [[0.0; 4]; 3];
    for (i, (&si, &ti)) in SS.iter().zip(&TT).enumerate() {
        shp[2][i] = (0.5 + si * s) * (0.5 + ti * t);
        shp[0][i] = si * (0.5 + ti * t);
        shp[1][i] = ti * (0.5 + si * s);
    }
    let mut xs = [[0.0; 2]; 2];
    for (a, row) in xs.iter_mut().enumerate() {
        for (b, v) in row.iter_mut().enumerate() {
            *v = (0..4).map(|k| xl[k][a] * shp[b][k]).sum();
        }
    }
    let det = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];
    let inv = 1.0 / det;
    let sx = [
        [xs[1][1] * inv, -xs[0][1] * inv],
        [-xs[1][0] * inv, xs[0][0] * inv],
    ];
    #[allow(clippy::needless_range_loop)] // reads and writes two rows of `shp` at column `i`
    for i in 0..4 {
        let dx = shp[0][i] * sx[0][0] + shp[1][i] * sx[1][0];
        shp[1][i] = shp[0][i] * sx[0][1] + shp[1][i] * sx[1][1];
        shp[0][i] = dx;
    }
    (shp, det)
}

fn point_geometry(axes: &[Vector3<f64>; 3], xl: &[[f64; 2]; 4]) -> [PointGeometry; 4] {
    let [g1, g2, g3] = axes;
    let (dx34, dy34) = (xl[2][0] - xl[3][0], xl[2][1] - xl[3][1]);
    let (dx21, dy21) = (xl[1][0] - xl[0][0], xl[1][1] - xl[0][1]);
    let (dx32, dy32) = (xl[2][0] - xl[1][0], xl[2][1] - xl[1][1]);
    let (dx41, dy41) = (xl[3][0] - xl[0][0], xl[3][1] - xl[0][1]);
    let q = 0.25;
    let mut g = SMatrix::<f64, 4, 12>::zeros();
    for (row, entries) in [
        [
            (0, -0.5),
            (1, -dy41 * q),
            (2, dx41 * q),
            (9, 0.5),
            (10, -dy41 * q),
            (11, dx41 * q),
        ],
        [
            (0, -0.5),
            (1, -dy21 * q),
            (2, dx21 * q),
            (3, 0.5),
            (4, -dy21 * q),
            (5, dx21 * q),
        ],
        [
            (3, -0.5),
            (4, -dy32 * q),
            (5, dx32 * q),
            (6, 0.5),
            (7, -dy32 * q),
            (8, dx32 * q),
        ],
        [
            (6, 0.5),
            (7, -dy34 * q),
            (8, dx34 * q),
            (9, -0.5),
            (10, -dy34 * q),
            (11, dx34 * q),
        ],
    ]
    .into_iter()
    .enumerate()
    {
        for (col, v) in entries {
            g[(row, col)] = v;
        }
    }

    let ax = -xl[0][0] + xl[1][0] + xl[2][0] - xl[3][0];
    let bx = xl[0][0] - xl[1][0] + xl[2][0] - xl[3][0];
    let cx = -xl[0][0] - xl[1][0] + xl[2][0] + xl[3][0];
    let ay = -xl[0][1] + xl[1][1] + xl[2][1] - xl[3][1];
    let by = xl[0][1] - xl[1][1] + xl[2][1] - xl[3][1];
    let cy = -xl[0][1] - xl[1][1] + xl[2][1] + xl[3][1];
    let alpha = (ay / ax).atan();
    let beta = std::f64::consts::FRAC_PI_2 - (cx / cy).atan();
    let rot = Matrix2::new(beta.sin(), -alpha.sin(), -beta.cos(), alpha.cos());

    std::array::from_fn(|i| {
        let (s, t) = (SG[i], TG[i]);
        let r1 = (cx + s * bx).hypot(cy + s * by);
        let r2 = (ax + t * bx).hypot(ay + t * by);
        let (shp, det) = shape2d(s, t, xl);

        let mut ms = SMatrix::<f64, 2, 4>::zeros();
        ms[(1, 0)] = 1.0 - s;
        ms[(0, 1)] = 1.0 - t;
        ms[(1, 2)] = 1.0 + s;
        ms[(0, 3)] = 1.0 + t;
        let mut bsv = ms * g;
        for j in 0..12 {
            bsv[(0, j)] *= r1 / (8.0 * det);
            bsv[(1, j)] *= r2 / (8.0 * det);
        }
        let bs = rot * bsv;

        let mut b = Bmat::zeros();
        let mut drill = Vec24::zeros();
        for j in 0..4 {
            let (nx, ny, n) = (shp[0][j], shp[1][j], shp[2][j]);
            let c = 6 * j;
            // Membrane: [Nx 0; 0 Ny; Ny Nx] times [g1; g2] on the translations.
            // Bending: [0 -Nx; Ny 0; Nx -Ny] times [g1; g2] on the rotations.
            for k in 0..3 {
                b[(0, c + k)] = nx * g1[k];
                b[(1, c + k)] = ny * g2[k];
                b[(2, c + k)] = ny * g1[k] + nx * g2[k];
                b[(3, c + 3 + k)] = -nx * g2[k];
                b[(4, c + 3 + k)] = ny * g1[k];
                b[(5, c + 3 + k)] = nx * g1[k] - ny * g2[k];
                // Shear: [w, theta1, theta2] -> [g3.u, g1.theta, g2.theta].
                for r in 0..2 {
                    b[(6 + r, c + k)] = bs[(r, 3 * j)] * g3[k];
                    b[(6 + r, c + 3 + k)] = bs[(r, 3 * j + 1)] * g1[k] + bs[(r, 3 * j + 2)] * g2[k];
                }
                let (b1, b2, b6) = (-0.5 * ny, 0.5 * nx, -n);
                drill[c + k] = b1 * g1[k] + b2 * g2[k];
                drill[c + 3 + k] = b6 * g3[k];
            }
        }
        let mut b_force = b;
        for j in 0..4 {
            for r in 3..6 {
                for k in 0..3 {
                    b_force[(r, 6 * j + 3 + k)] *= -1.0;
                }
            }
        }
        PointGeometry {
            b,
            b_force,
            drill,
            dvol: det,
            n: shp[2],
        }
    })
}

fn dofs(nodes: &[Node3Id; 4]) -> [DofRef<Node3Id>; 24] {
    std::array::from_fn(|a| (nodes[a / 6], (a % 6) as u8))
}

impl Shell4 {
    pub fn new(nodes: [Node3Id; 4], section: ShellSection) -> Self {
        Shell4 {
            nodes,
            sections: std::array::from_fn(|_| section.clone()),
            cache: None,
        }
    }

    /// The Gauss points' sections, in OpenSees' Gauss-point order.
    pub fn sections(&self) -> &[ShellSection; 4] {
        &self.sections
    }

    fn coords(&self, view: &View<'_>) -> [[f64; 3]; 4] {
        std::array::from_fn(|i| view.get(self.nodes[i]).coords)
    }

    fn compute_cache(&self, view: &View<'_>) -> Result<Cache, &'static str> {
        let (axes, xl) = frame(&self.coords(view))?;
        Ok(Cache {
            axes,
            points: point_geometry(&axes, &xl),
            ktt: self.sections[0].drilling_stiffness(),
        })
    }

    fn cache_or_compute(&self, view: &View<'_>) -> Cache {
        self.cache.unwrap_or_else(|| {
            self.compute_cache(view)
                .expect("unvalidated Shell4 geometry")
        })
    }

    pub(super) fn dof_mask(&self) -> DofMask {
        DofMask::all(SPATIAL_NDF)
    }

    pub(super) fn validate(&self, view: &View<'_>) -> Result<(), &'static str> {
        self.sections[0].validate().map_err(|e| e.message())?;
        let (_, xl) = frame(&self.coords(view))?;
        validate_quad(&xl)
    }

    pub(super) fn prepare(&mut self, view: &View<'_>) {
        self.cache = self.compute_cache(view).ok();
    }

    fn displacement(&self, view: &View<'_>) -> Vec24 {
        Vec24::from_fn(|a, _| view.get(self.nodes[a / 6]).displacement[a % 6])
    }

    /// Generalized strain at each Gauss point.
    fn strains(&self, cache: &Cache, u: &Vec24) -> [ShellVector; 4] {
        std::array::from_fn(|g| cache.points[g].b * u)
    }

    pub(super) fn assemble_tangent<S: TangentSink<Node3Id>>(&self, view: &View<'_>, sink: &mut S) {
        let cache = self.cache_or_compute(view);
        let u = self.displacement(view);
        let mut k = Mat24::zeros();
        let mut r = Vec24::zeros();
        for (g, section) in cache.points.iter().zip(&self.sections) {
            let (stress, d) = section.trial_resultant(&(g.b * u));
            k += g.dvol * g.b_force.transpose() * d * g.b;
            k += cache.ktt * g.dvol * g.drill * g.drill.transpose();
            r += g.dvol * g.b_force.transpose() * stress;
            r += cache.ktt * g.dvol * g.drill.dot(&u) * g.drill;
        }
        sink.add(&dofs(&self.nodes), &k, &r);
    }

    /// `integral(N_i dA)` per node: the pressure and body-force share, and (times `rho h`)
    /// the row-sum lumped mass.
    fn nodal_weights(&self, cache: &Cache) -> [f64; 4] {
        let mut w = [0.0; 4];
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
        let m = Vec24::from_fn(|a, _| if a % 6 < 3 { rho_h * w[a / 6] } else { 0.0 });
        sink.add(&dofs(&self.nodes), &m);
    }

    /// Equivalent nodal forces of `load` in global axes: the body force (a body
    /// acceleration scaled by `rho h`) and the pressure along `e3`.
    pub(super) fn load_vector(&self, view: &View<'_>, load: &ElementLoad3) -> Vec24 {
        let cache = self.cache_or_compute(view);
        let w = self.nodal_weights(&cache);
        let rho_h = self.sections[0].rho_h();
        let normal = cache.axes[2];
        Vec24::from_fn(|a, _| {
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
        let u = self.displacement(view);
        let strains = self.strains(&self.cache_or_compute(view), &u);
        for (section, strain) in self.sections.iter_mut().zip(&strains) {
            *section = section.commit(strain);
        }
    }

    /// The element's resisting force in global axes (24).
    pub(super) fn local_force(&self, view: &View<'_>) -> ElementForce {
        let cache = self.cache_or_compute(view);
        let u = self.displacement(view);
        let mut r = Vec24::zeros();
        for (g, section) in cache.points.iter().zip(&self.sections) {
            let (stress, _) = section.trial_resultant(&(g.b * u));
            r += g.dvol * g.b_force.transpose() * stress;
            r += cache.ktt * g.dvol * g.drill.dot(&u) * g.drill;
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
