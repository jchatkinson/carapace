//! Recovered internal response (stud, jack, header, plate, pier, base and hold-down forces) in local axes.

use super::*;

/// Full-height stud axial force (tension positive). Inside a sheathed pier the two ends differ by the shear
/// lag `(q_left − q_right)·H`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallAxialResponse {
    pub x: f64,
    pub bottom: f64,
    pub top: f64,
}
/// Jack or cripple axial force (tension positive).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallMemberResponse {
    pub x: f64,
    pub axial: f64,
}
/// Header moment (positive sagging) and shear at a header node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallHeaderResponse {
    pub x: f64,
    pub moment: f64,
    pub shear: f64,
}
/// Plate axial force (tension positive), moment (positive sagging) and shear at a plate node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallChordResponse {
    pub x: f64,
    pub axial: f64,
    pub moment: f64,
    pub shear: f64,
}
/// One full-height pier: its shear field and overturning.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallPierResponse {
    /// First full-height stud of the pier.
    pub xa: f64,
    /// Last full-height stud of the pier.
    pub xb: f64,
    /// Shear flow (force per length); positive for a +x load at the top.
    pub q: f64,
    /// Compatible pier shear strain (sheathing strain plus shear lag).
    pub gamma: f64,
    /// Moment of the pier's stud axial forces at their bottom ends about the pier mid-width, CCW positive: the
    /// overturning moment of everything above the base (negative for a +x load).
    pub moment_base: f64,
    /// As `moment_base`, at the stud top ends; `moment_top = moment_base + q·b·H`.
    pub moment_top: f64,
    /// End-post centroid spacing, `xb - xa`.
    pub jd: f64,
    /// Overturning share only, `[−moment_base/jd, +moment_base/jd]`, tension positive. This is the
    /// CSA O86 / SDPWS M/jd design force; gravity share comes from each post's stud force or gravity analysis.
    pub chords: [f64; 2],
}
/// Force on the wall at one supported bottom position (local axes; gravity gives positive `fy`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallBaseResponse {
    pub x: f64,
    pub fx: f64,
    pub fy: f64,
}
/// A pier-end hold-down of a bearing wall: the post's rise relative to the surface and the law's force
/// (tension positive).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallHoldDownResponse {
    pub x: f64,
    pub deformation: f64,
    pub force: f64,
}
/// Recovered internal response, in local axes, tension positive. See `docs/stud-wall.md` §7 for design use.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StudWallResponse {
    /// Full-height studs, in `x` order.
    pub studs: Vec<WallAxialResponse>,
    /// Jacks, two per opening.
    pub jacks: Vec<WallMemberResponse>,
    /// Cripples above openings.
    pub cripples: Vec<WallMemberResponse>,
    /// Header end reactions `[left, right]` per opening (the jack forces).
    pub header_reactions: Vec<[f64; 2]>,
    /// Header moment and shear at each header node.
    pub header: Vec<WallHeaderResponse>,
    /// Top plate results at each top-plate node.
    pub top_chord: Vec<WallChordResponse>,
    /// Bottom plate results at each bottom-plate node.
    pub bottom_chord: Vec<WallChordResponse>,
    /// Full-height piers, in `x` order (empty without sheathing).
    pub piers: Vec<WallPierResponse>,
    /// Bottom exterior nodes (`c0`, `c1`, bottom stations) and, with `foundation`, every internal support, by `x`.
    /// With `bearing`: every bearing position (spring and hold-down forces) instead.
    pub base_nodes: Vec<WallBaseResponse>,
    /// Pier-end hold-downs of a bearing wall, by `x` (empty otherwise).
    pub hold_downs: Vec<WallHoldDownResponse>,
    /// `[Fx, Fy, M]`: sums over `base_nodes`, `M = Σ (x − L/2)·fy`.
    pub base: [f64; 3],
    /// Floor-zone bays between consecutive top positions (empty without a zone).
    pub floor_bays: Vec<WallFloorBay>,
    /// Rim axial force per bay (tension positive; empty without a zone).
    pub floor_rim: Vec<f64>,
}
/// One floor-zone bay's shear field. `q` is zero for a bay over an opening.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallFloorBay {
    pub x0: f64,
    pub x1: f64,
    /// Shear strain.
    pub gamma: f64,
    /// Shear flow `ga·gamma`.
    pub q: f64,
}

impl<K: NodeKey> StudWall<K> {
    pub fn response<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        cases: &[WallCaseLoad],
    ) -> StudWallResponse {
        let ue = self.local_displacement(nodes);
        let mut ui = -(&self.recovery * &ue);
        for c in cases {
            ui += &self.load_recovery[c.case as usize] * c.factor;
        }
        let nonlinear_state = self
            .nonlinear
            .as_ref()
            .and_then(|_| self.nonlinear_trial(ue.as_slice(), cases));
        let nonlinear_q = nonlinear_state.as_ref().map(|(_, q, _)| q.clone());
        let nonlinear_gamma = nonlinear_state.as_ref().map(|(g, _, _)| g.clone());
        if let (Some(cache), Some((_, q, _))) = (&self.nonlinear, &nonlinear_state) {
            for i in 0..ui.len() {
                for j in 0..q.len() {
                    ui[i] -= cache.internal_recovery[(i, j)] * q[j];
                }
            }
        }
        let mut ufull = ue.as_slice().to_vec();
        ufull.extend(ui.iter().copied());
        let q: Vec<f64> = if self.pier_comps == 0 {
            self.piers
                .iter()
                .map(|p| {
                    let gamma = p.beta.dot(&DVector::from_column_slice(&ufull));
                    let WallShearLaw::Linear { ga } = &self.layout.sheathing.as_ref().unwrap().law
                    else {
                        unreachable!("nonlinear laws use the local pier solve")
                    };
                    ga * p.area * gamma / (p.area + ga * p.c)
                })
                .collect()
        } else if let Some(q) = nonlinear_q.as_ref() {
            q[..self.pier_comps].to_vec()
        } else {
            vec![f64::NAN; self.piers.len()]
        };
        let disp = |d: usize| {
            if d < ue.len() {
                ue[d]
            } else {
                ui[d - ue.len()]
            }
        };
        let mut r = StudWallResponse::default();
        for b in &self.bars {
            let n = b.k * (disp(b.b) - disp(b.a));
            if b.stud {
                let q_left = self.field_side(b.x, true, &q);
                let q_right = self.field_side(b.x, false, &q);
                let delta = q_left - q_right;
                r.studs.push(WallAxialResponse {
                    x: b.x,
                    bottom: n - delta * b.length / 2.,
                    top: n + delta * b.length / 2.,
                });
            } else if b.jack {
                r.jacks.push(WallMemberResponse { x: b.x, axial: n });
            } else {
                r.cripples.push(WallMemberResponse { x: b.x, axial: n });
            }
        }
        r.header_reactions = r.jacks.iter().map(|j| [j.x, j.axial]).collect();
        for h in &self.headers {
            let l = h.x1 - h.x0;
            let [v0, t0, v1, t1] = h.dofs;
            let d = [disp(v0), disp(t0), disp(v1), disp(t1)];
            let k = beam4(h.ei, l);
            let f: [f64; 4] = std::array::from_fn(|a| (0..4).map(|b| k[a][b] * d[b]).sum());
            for (x, m, shear) in [(h.x0, f[1], f[0]), (h.x1, f[3], -f[2])] {
                r.header.push(WallHeaderResponse {
                    x,
                    moment: m,
                    shear,
                });
            }
        }
        r.top_chord = self.chord_response(&self.top_nodes, &ue, &ui, WallChord::Top, &q, cases);
        r.bottom_chord =
            self.chord_response(&self.bottom_nodes, &ue, &ui, WallChord::Bottom, &q, cases);
        for (i, p) in self.piers.iter().enumerate() {
            let gamma = p.beta.dot(&DVector::from_column_slice(&ufull));
            let mid = (p.xa + p.xb) * 0.5;
            let mut moment_base = 0.;
            let mut moment_top = 0.;
            for stud in r.studs.iter().filter(|s| s.x >= p.xa && s.x <= p.xb) {
                moment_base += (stud.x - mid) * stud.bottom;
                moment_top += (stud.x - mid) * stud.top;
            }
            let jd = p.xb - p.xa;
            r.piers.push(WallPierResponse {
                xa: p.xa,
                xb: p.xb,
                q: q[i],
                gamma,
                moment_base,
                moment_top,
                jd,
                chords: [-moment_base / jd, moment_base / jd],
            });
        }
        let mut resisting = &self.stiffness * &ue;
        if let (Some(cache), Some(q)) = (&self.nonlinear, &nonlinear_q) {
            for i in 0..resisting.len() {
                for j in 0..q.len() {
                    resisting[i] += cache.coupling[(i, j)] * q[j];
                }
            }
        }
        for c in cases {
            resisting -= &self.load_slots[c.case as usize] * c.factor;
        }
        for (i, id) in self.nodes.iter().enumerate() {
            if self.layout.bearing.is_some() {
                break;
            }
            if i < 2
                || self
                    .stations
                    .iter()
                    .any(|s| s.node == *id && s.chord == WallChord::Bottom)
            {
                let fx = resisting[2 * i];
                let fy = resisting[2 * i + 1];
                let x = if i == 0 {
                    0.
                } else if i == 1 {
                    self.layout.length
                } else {
                    self.stations
                        .iter()
                        .find(|s| s.node == *id)
                        .map(|s| s.x)
                        .unwrap_or(0.)
                };
                r.base[0] += fx;
                r.base[1] += fy;
                r.base[2] += (x - self.layout.length / 2.) * fy;
                r.base_nodes.push(WallBaseResponse { x, fx, fy });
            }
        }
        if let Some(bl) = &self.layout.bearing {
            // Spring forces on the wall (gravity positive): `-k (plate - surface)`; a hold-down's law force.
            let rel = |d: usize, slot: usize, base: &[(usize, f64)]| {
                disp(d) - base.iter().map(|&(e, w)| w * ue[2 * e + slot]).sum::<f64>()
            };
            for b in &self.bearing {
                let fx = -bl.kh * rel(b.u, 0, &b.base);
                let fy = match b.hold {
                    Some(c) => {
                        let (deformation, force) = match (&nonlinear_gamma, &nonlinear_q) {
                            (Some(g), Some(q)) => (g[c], q[c]),
                            _ => (f64::NAN, f64::NAN),
                        };
                        r.hold_downs.push(WallHoldDownResponse {
                            x: b.x,
                            deformation,
                            force,
                        });
                        -force
                    }
                    None => -bl.kv * rel(b.v, 1, &b.base),
                };
                r.base[0] += fx;
                r.base[1] += fy;
                r.base[2] += (b.x - self.layout.length / 2.) * fy;
                r.base_nodes.push(WallBaseResponse { x: b.x, fx, fy });
            }
        }
        if let Some(z) = self.layout.floor {
            for b in &self.zone_bays {
                let u: [f64; 8] = std::array::from_fn(|i| disp(b.dofs[i]));
                let (gamma, q, rim) = bay_response(b.x1 - b.x0, z.depth, b.ga, z.rim_ea, &u);
                r.floor_bays.push(WallFloorBay {
                    x0: b.x0,
                    x1: b.x1,
                    gamma,
                    q,
                });
                r.floor_rim.push(rim);
            }
        }
        if !self.foundation_supports.is_empty() {
            let mut residual =
                self.full_stiffness.as_ref().unwrap() * DVector::from_column_slice(&ufull);
            if let Some(q) = &nonlinear_q {
                for (i, pier) in self.comps.iter().enumerate() {
                    residual += &pier.beta * (pier.area * q[i]);
                }
            }
            for c in cases {
                residual -= &self.full_load_slots[c.case as usize] * c.factor;
            }
            for &(x, u, v) in &self.foundation_supports {
                let fx = residual[u];
                let fy = residual[v];
                r.base[0] += fx;
                r.base[1] += fy;
                r.base[2] += (x - self.layout.length / 2.) * fy;
                r.base_nodes.push(WallBaseResponse { x, fx, fy });
            }
            r.base_nodes.sort_by(|a, b| a.x.total_cmp(&b.x));
        }
        r
    }
    pub(super) fn chord_response(
        &self,
        _nodes: &[(f64, [usize; 3])],
        ue: &DVector<f64>,
        ui: &DVector<f64>,
        chord: WallChord,
        q: &[f64],
        cases: &[WallCaseLoad],
    ) -> Vec<WallChordResponse> {
        let get = |d: usize| {
            if d < ue.len() {
                ue[d]
            } else {
                ui[d - ue.len()]
            }
        };
        let mut out = Vec::with_capacity(2 * self.axial_chords.len());
        for (ai, b) in self
            .axial_chords
            .iter()
            .enumerate()
            .filter(|(_, b)| b.chord == chord)
        {
            let l = b.x1 - b.x0;
            let n = b.ea * (get(b.b) - get(b.a)) / l;
            let p0 = cases
                .iter()
                .map(|c| c.factor * self.axial_load_slots[c.case as usize][ai][0])
                .sum::<f64>();
            let p1 = cases
                .iter()
                .map(|c| c.factor * self.axial_load_slots[c.case as usize][ai][1])
                .sum::<f64>();
            let qseg = self
                .piers
                .iter()
                .enumerate()
                .find(|(_, p)| (b.x0 + b.x1) * 0.5 >= p.xa && (b.x0 + b.x1) * 0.5 <= p.xb)
                .map(|(i, _)| q[i])
                .unwrap_or(0.);
            let beam_i = self.beams.iter().position(|m| {
                m.chord == chord && (m.x0 - b.x0).abs() < 1e-10 && (m.x1 - b.x1).abs() < 1e-10
            });
            let (m0, v0, m1, v1) = if let Some(bi) = beam_i {
                let m = &self.beams[bi];
                let w = m.dofs;
                let d = [get(w[1]), get(w[2]), get(w[4]), get(w[5])];
                let k = beam4(m.ei, l);
                let f: [f64; 4] = std::array::from_fn(|a| (0..4).map(|j| k[a][j] * d[j]).sum());
                let p: [f64; 6] = std::array::from_fn(|j| {
                    cases
                        .iter()
                        .map(|c| c.factor * self.beam_load_slots[c.case as usize][bi][j])
                        .sum()
                });
                (f[1] - p[2], f[0] - p[1], f[3] - p[5], -(f[2] - p[4]))
            } else {
                (0., 0., 0., 0.)
            };
            for (x, force, moment, shear) in [(b.x0, n + p0, m0, v0), (b.x1, n - p1, m1, v1)] {
                let field = if chord == WallChord::Top {
                    qseg * (x - (b.x0 + b.x1) / 2.)
                } else {
                    -qseg * (x - (b.x0 + b.x1) / 2.)
                };
                out.push(WallChordResponse {
                    x,
                    axial: force + field,
                    moment,
                    shear,
                });
            }
        }
        out
    }
    pub(super) fn field_side(&self, x: f64, left: bool, q: &[f64]) -> f64 {
        self.piers
            .iter()
            .enumerate()
            .filter(|(_, p)| {
                if left {
                    (p.xb - x).abs() < 1e-10
                } else {
                    (p.xa - x).abs() < 1e-10
                }
            })
            .map(|(i, _)| q[i])
            .sum()
    }
}
