use arrayvec::ArrayVec;
use nalgebra::DMatrix;

/// Where a wall-type element's local in-plane axes sit in the model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WallPlane {
    /// The model plane; its x axis follows the first two base nodes.
    Planar,
    /// A vertical plane with a unit plan direction and global Y up.
    Vertical { ex: f64, ez: f64 },
}

impl WallPlane {
    /// Creates a vertical plane after checking that its plan direction is unit length. A component below 1e-12
    /// (e.g. `cos 90°`) is snapped to zero, so an axis-aligned wall does not activate a near-zero-stiffness slot.
    pub fn vertical(ex: f64, ez: f64) -> Result<Self, String> {
        if !ex.is_finite() || !ez.is_finite() || ((ex * ex + ez * ez).sqrt() - 1.0).abs() > 1e-9 {
            return Err("wall plane direction must be finite and unit length".into());
        }
        let snap = |v: f64| if v.abs() < 1e-12 { 0. } else { v };
        let (ex, ez) = (snap(ex), snap(ez));
        let n = ex.hypot(ez);
        Ok(Self::Vertical {
            ex: ex / n,
            ez: ez / n,
        })
    }

    pub(crate) fn slots(&self) -> ArrayVec<u8, 3> {
        let mut out = ArrayVec::new();
        match *self {
            Self::Planar => out.extend([0, 1]),
            Self::Vertical { ex, ez } => {
                if ex != 0. {
                    out.push(0);
                }
                out.push(1);
                if ez != 0. {
                    out.push(2);
                }
            }
        }
        out
    }

    pub(crate) fn map<const NDIM: usize>(&self, p0: [f64; NDIM], p1: [f64; NDIM]) -> PlaneMap {
        let mut origin = [0.; 3];
        origin[..NDIM.min(3)].copy_from_slice(&p0[..NDIM.min(3)]);
        match *self {
            Self::Planar => {
                let dx = p1[0] - p0[0];
                let dy = p1[1] - p0[1];
                let l = dx.hypot(dy);
                let (c, s) = (dx / l, dy / l);
                let mut slots = ArrayVec::new();
                slots.extend([0, 1]);
                let mut t = ArrayVec::new();
                t.extend([[c, -s], [s, c]]);
                PlaneMap {
                    slots,
                    t,
                    origin,
                    ex: [c, s, 0.],
                    ey: [-s, c, 0.],
                }
            }
            Self::Vertical { ex, ez } => {
                let mut slots = ArrayVec::new();
                let mut t = ArrayVec::new();
                if ex != 0. {
                    slots.push(0);
                    t.push([ex, 0.]);
                }
                slots.push(1);
                t.push([0., 1.]);
                if ez != 0. {
                    slots.push(2);
                    t.push([ez, 0.]);
                }
                PlaneMap {
                    slots,
                    t,
                    origin,
                    ex: [ex, 0., ez],
                    ey: [0., 1., 0.],
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct PlaneMap {
    pub slots: ArrayVec<u8, 3>,
    pub t: ArrayVec<[f64; 2], 3>,
    pub origin: [f64; 3],
    pub ex: [f64; 3],
    pub ey: [f64; 3],
}

impl PlaneMap {
    pub fn local(&self, d: &[f64]) -> [f64; 2] {
        let mut out = [0.; 2];
        for (k, &slot) in self.slots.iter().enumerate() {
            out[0] += self.t[k][0] * d[slot as usize];
            out[1] += self.t[k][1] * d[slot as usize];
        }
        out
    }
    pub fn global(&self, f: [f64; 2]) -> ArrayVec<f64, 3> {
        self.slots
            .iter()
            .enumerate()
            .map(|(k, _)| self.t[k][0] * f[0] + self.t[k][1] * f[1])
            .collect()
    }
    pub fn point(&self, p: &[f64]) -> (f64, f64, f64) {
        let d = [
            p[0] - self.origin[0],
            p.get(1).copied().unwrap_or(0.) - self.origin[1],
            p.get(2).copied().unwrap_or(0.) - self.origin[2],
        ];
        let x = d.iter().zip(self.ex).map(|(a, b)| a * b).sum();
        let y = d.iter().zip(self.ey).map(|(a, b)| a * b).sum();
        let off = if p.len() == 2 {
            0.
        } else {
            (d[0] * (self.ex[1] * self.ey[2] - self.ex[2] * self.ey[1])
                + d[1] * (self.ex[2] * self.ey[0] - self.ex[0] * self.ey[2])
                + d[2] * (self.ex[0] * self.ey[1] - self.ex[1] * self.ey[0]))
                .abs()
        };
        (x, y, off)
    }
    pub fn to_global_matrix(&self, kl: &DMatrix<f64>, n_nodes: usize) -> Vec<f64> {
        let m = self.slots.len();
        let n = m * n_nodes;
        let mut out = vec![0.; n * n];
        for i in 0..n_nodes {
            for j in 0..n_nodes {
                for a in 0..m {
                    for b in 0..m {
                        let mut v = 0.;
                        for l in 0..2 {
                            for q in 0..2 {
                                v += self.t[a][l] * kl[(2 * i + l, 2 * j + q)] * self.t[b][q];
                            }
                        }
                        out[(i * m + a) * n + j * m + b] = v;
                    }
                }
            }
        }
        out
    }
    pub fn to_global_vector(&self, local: &[f64], n_nodes: usize) -> Vec<f64> {
        let m = self.slots.len();
        let mut out = vec![0.; m * n_nodes];
        for i in 0..n_nodes {
            for a in 0..m {
                out[i * m + a] = self.t[a][0] * local[2 * i] + self.t[a][1] * local[2 * i + 1];
            }
        }
        out
    }
    pub fn to_global_cols(&self, local: &DMatrix<f64>, n_nodes: usize) -> Vec<f64> {
        let m = self.slots.len();
        let p = local.ncols();
        let mut out = vec![0.; m * n_nodes * p];
        for i in 0..n_nodes {
            for a in 0..m {
                for j in 0..p {
                    out[(i * m + a) * p + j] =
                        self.t[a][0] * local[(2 * i, j)] + self.t[a][1] * local[(2 * i + 1, j)];
                }
            }
        }
        out
    }
    // This maps local strain rows back to the global translation columns, so the map is required.
    #[allow(clippy::wrong_self_convention)]
    pub fn from_global_rows(&self, local: &DMatrix<f64>, n_nodes: usize) -> Vec<f64> {
        let m = self.slots.len();
        let p = local.nrows();
        let mut out = vec![0.; p * m * n_nodes];
        for i in 0..n_nodes {
            for a in 0..p {
                for b in 0..m {
                    out[a * m * n_nodes + i * m + b] =
                        local[(a, 2 * i)] * self.t[b][0] + local[(a, 2 * i + 1)] * self.t[b][1];
                }
            }
        }
        out
    }
}
