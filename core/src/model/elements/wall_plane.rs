use arrayvec::ArrayVec;
use nalgebra::DMatrix;

/// Geometry-derived map between a wall's two local translations and global translations.
#[derive(Debug, Clone)]
pub(crate) struct PlaneMap {
    pub slots: ArrayVec<u8, 3>,
    pub t: ArrayVec<[f64; 2], 3>,
    pub origin: [f64; 3],
    pub ex: [f64; 3],
    pub ey: [f64; 3],
}

impl PlaneMap {
    /// Derives local x from `c0 -> c1` and local y from `c0 -> c3`.
    pub fn from_corners<const NDIM: usize>(
        p0: [f64; NDIM],
        p1: [f64; NDIM],
        p3: [f64; NDIM],
    ) -> Result<Self, &'static str> {
        let mut origin = [0.; 3];
        let mut dx = [0.; 3];
        let mut dy = [0.; 3];
        for i in 0..NDIM.min(3) {
            origin[i] = p0[i];
            dx[i] = p1[i] - p0[i];
            dy[i] = p3[i] - p0[i];
        }
        let lx = dx.iter().map(|v| v * v).sum::<f64>().sqrt();
        let ly = dy.iter().map(|v| v * v).sum::<f64>().sqrt();
        if !lx.is_finite() || !ly.is_finite() || lx <= f64::MIN_POSITIVE || ly <= f64::MIN_POSITIVE
        {
            return Err("wall axes are degenerate");
        }
        let ex = dx.map(|v| v / lx);
        let ey = dy.map(|v| v / ly);
        let dot = ex.iter().zip(ey).map(|(a, b)| a * b).sum::<f64>();
        if !dot.is_finite() || dot.abs() > 1e-9 {
            return Err("wall axes are not orthogonal");
        }
        let mut slots = ArrayVec::new();
        let mut t = ArrayVec::new();
        for i in 0..NDIM.min(3) {
            if ex[i].abs() > 1e-12 || ey[i].abs() > 1e-12 {
                slots.push(i as u8);
                t.push([ex[i], ey[i]]);
            }
        }
        Ok(Self {
            slots,
            t,
            origin,
            ex,
            ey,
        })
    }

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
