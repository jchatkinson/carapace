//! Top-plate loads discretized onto the internal frame.

use super::*;

pub(super) fn discretize_load(
    load: &WallLoad,
    slot: usize,
    l: &StudWallLayout,
    nodes: &[(f64, [usize; 3])],
    segments: &[ChordSegment],
    ext_top: &[(f64, usize)],
    p: &mut DVector<f64>,
) -> Result<(), StudWallError> {
    // Lever-rule split of a horizontal force onto the exterior top row (floor zone only).
    let lever = |p: &mut DVector<f64>, x: f64, f: f64| {
        let i = ext_top
            .windows(2)
            .position(|w| x <= w[1].0)
            .unwrap_or(ext_top.len() - 2);
        let (a, b) = (ext_top[i], ext_top[i + 1]);
        let t = ((x - a.0) / (b.0 - a.0)).clamp(0., 1.);
        p[a.1] += (1. - t) * f;
        p[b.1] += t * f;
    };
    let mut point = |x: f64, fx: f64, fy: f64| -> Result<(), StudWallError> {
        if ![x, fx, fy].iter().all(|v| v.is_finite()) || x < 0. || x > l.length {
            return Err(StudWallError(format!(
                "load case {slot} contains invalid coordinates or intensities"
            )));
        }
        let fx = if ext_top.is_empty() {
            fx
        } else {
            lever(p, x, fx);
            0.
        };
        if let Some((_, d)) = nodes
            .iter()
            .find(|(xx, _)| (*xx - x).abs() < 1e-12 * l.length)
        {
            p[d[0]] += fx;
            p[d[1]] += fy;
            return Ok(());
        }
        let pair = segments
            .iter()
            .find(|w| x > w[0].0 && x < w[1].0)
            .ok_or_else(|| {
                StudWallError(format!("load case {slot} point is outside the top chord"))
            })?;
        let le = pair[1].0 - pair[0].0;
        let t = (x - pair[0].0) / le;
        if l.top_distribution == ChordDistribution::Direct {
            p[pair[0].1[0]] += (1. - t) * fx;
            p[pair[1].1[0]] += t * fx;
            p[pair[0].1[1]] += (1. - t) * fy;
            p[pair[1].1[1]] += t * fy;
        } else {
            let h1 = 1. - 3. * t * t + 2. * t * t * t;
            let h2 = le * (t - 2. * t * t + t * t * t);
            let h3 = 3. * t * t - 2. * t * t * t;
            let h4 = le * (-t * t + t * t * t);
            for (d, v) in [
                (pair[0].1[0], (1. - t) * fx),
                (pair[1].1[0], t * fx),
                (pair[0].1[1], h1 * fy),
                (pair[0].1[2], h2 * fy),
                (pair[1].1[1], h3 * fy),
                (pair[1].1[2], h4 * fy),
            ] {
                p[d] += v;
            }
        }
        Ok(())
    };
    match *load {
        WallLoad::TopPoint { x, fx, fy } => point(x, fx, fy)?,
        WallLoad::TopLine {
            x0,
            x1,
            wx0,
            wx1,
            wy0,
            wy1,
        } => {
            if ![x0, x1, wx0, wx1, wy0, wy1].iter().all(|v| v.is_finite())
                || x0 < 0.
                || x1 > l.length
                || x1 <= x0
            {
                return Err(StudWallError(format!(
                    "load case {slot} has invalid line bounds or intensities"
                )));
            }
            for w in segments {
                let a = x0.max(w[0].0);
                let b = x1.min(w[1].0);
                if b <= a {
                    continue;
                }
                let gauss = [
                    -0.906179845938664,
                    -0.538469310105683,
                    0.,
                    0.538469310105683,
                    0.906179845938664,
                ];
                let weights = [
                    0.236926885056189,
                    0.478628670499366,
                    0.568888888888889,
                    0.478628670499366,
                    0.236926885056189,
                ];
                let dx = (b - a) / 2.;
                let mid = (a + b) / 2.;
                for (z, gw) in gauss.iter().zip(weights) {
                    let x = mid + dx * z;
                    let t = (x - x0) / (x1 - x0);
                    let wx = wx0 + (wx1 - wx0) * t;
                    let wy = wy0 + (wy1 - wy0) * t;
                    let s = (x - w[0].0) / (w[1].0 - w[0].0);
                    let fac = dx * gw;
                    let wx = if ext_top.is_empty() {
                        wx
                    } else {
                        lever(p, x, fac * wx);
                        0.
                    };
                    if l.top_distribution == ChordDistribution::Direct {
                        for (d, v) in [
                            (w[0].1[0], (1. - s) * wx),
                            (w[1].1[0], s * wx),
                            (w[0].1[1], (1. - s) * wy),
                            (w[1].1[1], s * wy),
                        ] {
                            p[d] += fac * v;
                        }
                    } else {
                        let le = w[1].0 - w[0].0;
                        for (d, v) in [
                            (w[0].1[0], (1. - s) * wx),
                            (w[1].1[0], s * wx),
                            (w[0].1[1], (1. - 3. * s * s + 2. * s.powi(3)) * wy),
                            (w[0].1[2], le * (s - 2. * s * s + s.powi(3)) * wy),
                            (w[1].1[1], (3. * s * s - 2. * s.powi(3)) * wy),
                            (w[1].1[2], le * (-s * s + s.powi(3)) * wy),
                        ] {
                            p[d] += fac * v;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

pub(super) fn segment_load(
    loads: &[WallLoad],
    x0: f64,
    x1: f64,
    length: f64,
    beam: bool,
) -> [f64; 6] {
    let mut out = [0.; 6];
    let le = x1 - x0;
    let mut add = |x: f64, fx: f64, fy: f64, scale: f64| {
        if x <= x0 + 1e-12 * length.max(1.) || x >= x1 - 1e-12 * length.max(1.) {
            return;
        }
        let t = (x - x0) / le;
        out[0] += scale * (1. - t) * fx;
        out[3] += scale * t * fx;
        if beam {
            let n = [
                1. - 3. * t * t + 2. * t.powi(3),
                le * (t - 2. * t * t + t.powi(3)),
                3. * t * t - 2. * t.powi(3),
                le * (-t * t + t.powi(3)),
            ];
            for (j, d) in [1, 2, 4, 5].into_iter().enumerate() {
                out[d] += scale * n[j] * fy;
            }
        } else {
            out[1] += scale * (1. - t) * fy;
            out[4] += scale * t * fy;
        }
    };
    for load in loads {
        match *load {
            WallLoad::TopPoint { x, fx, fy } => add(x, fx, fy, 1.),
            WallLoad::TopLine {
                x0: a,
                x1: b,
                wx0,
                wx1,
                wy0,
                wy1,
            } => {
                let a = a.max(x0);
                let b = b.min(x1);
                if b <= a {
                    continue;
                }
                let z = [
                    -0.906179845938664,
                    -0.538469310105683,
                    0.,
                    0.538469310105683,
                    0.906179845938664,
                ];
                let w = [
                    0.236926885056189,
                    0.478628670499366,
                    0.568888888888889,
                    0.478628670499366,
                    0.236926885056189,
                ];
                let mid = (a + b) * 0.5;
                let half = (b - a) * 0.5;
                for i in 0..5 {
                    let x = mid + half * z[i];
                    let t =
                        (x - (match *load {
                            WallLoad::TopLine { x0, .. } => x0,
                            _ => 0.,
                        })) / (match *load {
                            WallLoad::TopLine { x0, x1, .. } => x1 - x0,
                            _ => 1.,
                        });
                    add(x, wx0 + (wx1 - wx0) * t, wy0 + (wy1 - wy0) * t, half * w[i]);
                }
            }
        }
    }
    out
}
