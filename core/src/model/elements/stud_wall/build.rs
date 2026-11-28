//! Construction of a stud wall: the internal frame (plates, studs, headers, jacks, sheathing, floor zone and
//! bearing springs), its static condensation onto the exterior nodes, and the prepared load cases.

use super::*;

impl<K: NodeKey> StudWall<K> {
    pub(super) fn build(
        corners: [K; 4],
        stations: Vec<WallStation<K>>,
        layout: StudWallLayout,
        cases: Vec<Vec<WallLoad>>,
    ) -> Result<Self, StudWallError> {
        validate_layout(&layout, &stations)?;
        let is_nonlinear = layout
            .sheathing
            .as_ref()
            .is_some_and(|s| !matches!(&s.law, WallShearLaw::Linear { .. }));
        let potential_piers = if is_nonlinear {
            pier_bounds(&layout).len()
        } else {
            0
        };
        if potential_piers > MAX_WALL_PIERS {
            return Err(StudWallError(format!(
                "pier count {potential_piers} exceeds nonlinear stud-wall limit of {MAX_WALL_PIERS}"
            )));
        }
        let mut nodes = corners.to_vec();
        nodes.extend(stations.iter().map(|s| s.node));
        if nodes
            .iter()
            .copied()
            .collect::<std::collections::HashSet<_>>()
            .len()
            != nodes.len()
        {
            return Err(StudWallError(
                "stations and corners must use distinct exterior nodes".into(),
            ));
        }
        if nodes.len() > MAX_EXTERIOR_NODES {
            return Err(StudWallError(format!(
                "exterior node count {} exceeds stud-wall limit of {MAX_EXTERIOR_NODES}",
                nodes.len()
            )));
        }
        let ne = nodes.len() * 2;
        let mut model = FrameBuild::new(ne);
        let mut top_x: Vec<f64> = layout.studs.iter().map(|s| s.x).collect();
        for o in &layout.openings {
            for s in &layout.studs {
                if s.x > o.x0 && s.x < o.x1 {
                    top_x.push(s.x);
                }
            }
        }
        top_x.extend(
            stations
                .iter()
                .filter(|s| s.chord == WallChord::Top)
                .map(|s| s.x),
        );
        if layout.top_distribution == ChordDistribution::Beam {
            for o in &layout.openings {
                top_x.extend([o.x0, o.x1]);
            }
        }
        top_x.sort_by(f64::total_cmp);
        dedup_x(&mut top_x, layout.length);
        let bottom_studs: Vec<f64> = layout
            .studs
            .iter()
            .filter(|s| !layout.openings.iter().any(|o| s.x > o.x0 && s.x < o.x1))
            .map(|s| s.x)
            .collect();
        let mut bottom_x = bottom_studs.clone();
        for o in &layout.openings {
            bottom_x.extend([o.x0, o.x1]);
        }
        bottom_x.extend(
            stations
                .iter()
                .filter(|s| s.chord == WallChord::Bottom)
                .map(|s| s.x),
        );
        bottom_x.sort_by(f64::total_cmp);
        dedup_x(&mut bottom_x, layout.length);
        let top_dofs = if layout.top_distribution == ChordDistribution::Beam {
            3
        } else {
            2
        };
        let top_nodes = make_chord_nodes(
            &mut model,
            &top_x,
            layout.height,
            WallChord::Top,
            top_dofs,
            layout.floor.is_some(),
            &corners,
            &stations,
            &nodes,
        );
        let bottom_nodes = make_chord_nodes(
            &mut model,
            &bottom_x,
            0.,
            WallChord::Bottom,
            3,
            false,
            &corners,
            &stations,
            &nodes,
        );
        // Floor zone: bays between consecutive top positions (corners and top stations), rim included.
        let mut zone_bays = Vec::new();
        if let Some(z) = layout.floor {
            let mut pos: Vec<(f64, usize)> = vec![(0., 3), (layout.length, 2)];
            for (j, s) in stations.iter().enumerate() {
                if s.chord == WallChord::Top {
                    pos.push((s.x, 4 + j));
                }
            }
            pos.sort_by(|a, b| a.0.total_cmp(&b.0));
            let tol = 1e-9 * layout.length;
            for w in pos.windows(2) {
                let (a, b) = (w[0], w[1]);
                let fa = find_chord(&top_nodes, a.0)
                    .ok_or_else(|| StudWallError("missing top position node".into()))?;
                let fb = find_chord(&top_nodes, b.0)
                    .ok_or_else(|| StudWallError("missing top position node".into()))?;
                let ga = if layout
                    .openings
                    .iter()
                    .any(|o| a.0 >= o.x0 - tol && b.0 <= o.x1 + tol)
                {
                    0.
                } else {
                    z.ga
                };
                let dofs = [
                    fa[0],
                    fa[1],
                    fb[0],
                    fb[1],
                    2 * a.1,
                    2 * a.1 + 1,
                    2 * b.1,
                    2 * b.1 + 1,
                ];
                let k = bay_stiffness(b.0 - a.0, z.depth, ga, z.rim_ea);
                for i in 0..8 {
                    for j in 0..8 {
                        model.k[(dofs[i], dofs[j])] += k[(i, j)];
                    }
                }
                zone_bays.push(ZoneBayRec {
                    x0: a.0,
                    x1: b.0,
                    dofs,
                    ga,
                });
            }
        }
        let mut header_nodes = Vec::new();
        for o in &layout.openings {
            let mut xs = vec![o.x0, o.x1];
            xs.extend(
                layout
                    .studs
                    .iter()
                    .filter(|s| s.x > o.x0 && s.x < o.x1)
                    .map(|s| s.x),
            );
            xs.sort_by(f64::total_cmp);
            dedup_x(&mut xs, layout.length);
            let group: Vec<_> = xs
                .iter()
                .map(|&x| {
                    let d = model.add_internal_node(x, o.head_height, 2);
                    let r = [d[0].unwrap(), d[1].unwrap()];
                    header_nodes.push((x, r));
                    r
                })
                .collect();
            for i in 0..xs.len() - 1 {
                let d = xs[i + 1] - xs[i];
                model.add_beam4(
                    [group[i][0], group[i][1], group[i + 1][0], group[i + 1][1]],
                    o.header_ei,
                    d,
                );
            }
        }
        let mut bars = Vec::new();
        let mut beams = Vec::new();
        let mut headers = Vec::new();
        let mut axial_chords = Vec::new();
        // Chord elements, with rotations retained internally and exterior translations mapped directly.
        let top_segments = add_chord_elements(
            &mut model,
            &top_nodes,
            layout.top_ea,
            layout.top_ei.unwrap_or(0.),
            layout.top_distribution,
            WallChord::Top,
            &layout.openings,
            &mut beams,
            &mut axial_chords,
        );
        for pair in bottom_nodes.windows(2) {
            if layout
                .openings
                .iter()
                .any(|o| o.door && pair[0].0 >= o.x0 - 1e-12 && pair[1].0 <= o.x1 + 1e-12)
            {
                continue;
            }
            let [a, b] = released_segment(&mut model, pair, &layout.openings);
            add_beam3(
                &mut model,
                a.0,
                b.0,
                a.1,
                b.1,
                layout.bottom_ea,
                layout.bottom_ei,
                WallChord::Bottom,
                &mut beams,
            );
            axial_chords.push(ChordAxialRec {
                x0: pair[0].0,
                x1: pair[1].0,
                a: pair[0].1[0],
                b: pair[1].1[0],
                ea: layout.bottom_ea,
                chord: WallChord::Bottom,
            });
        }
        for (x, stud) in layout.studs.iter().map(|s| (s.x, s)) {
            if layout.openings.iter().any(|o| x > o.x0 && x < o.x1) {
                continue;
            }
            let a = find_chord(&bottom_nodes, x)
                .ok_or_else(|| StudWallError("missing bottom stud node".into()))?;
            let b = find_chord(&top_nodes, x)
                .ok_or_else(|| StudWallError("missing top stud node".into()))?;
            let k = stud.ea / layout.height;
            model.add_bar(a[1], b[1], k);
            bars.push(BarRec {
                x,
                length: layout.height,
                k,
                a: a[1],
                b: b[1],
                stud: true,
                jack: false,
            });
        }
        for o in &layout.openings {
            let hnode = header_nodes
                .iter()
                .find(|(x, _)| (*x - o.x0).abs() < 1e-12)
                .unwrap()
                .1;
            let hright = header_nodes
                .iter()
                .find(|(x, _)| (*x - o.x1).abs() < 1e-12)
                .unwrap()
                .1;
            for (x, h) in [(o.x0, hnode[0]), (o.x1, hright[0])] {
                let b = find_chord(&bottom_nodes, x).unwrap()[1];
                let k = o.jack_ea / o.head_height;
                model.add_bar(b, h, k);
                bars.push(BarRec {
                    x,
                    length: o.head_height,
                    k,
                    a: b,
                    b: h,
                    stud: false,
                    jack: true,
                });
            }
            for s in layout.studs.iter().filter(|s| s.x > o.x0 && s.x < o.x1) {
                let a = header_nodes
                    .iter()
                    .find(|(x, _)| (*x - s.x).abs() < 1e-12)
                    .unwrap()
                    .1[0];
                let b = find_chord(&top_nodes, s.x).unwrap()[1];
                let len = layout.height - o.head_height;
                let k = s.ea / len;
                model.add_bar(a, b, k);
                bars.push(BarRec {
                    x: s.x,
                    length: len,
                    k,
                    a,
                    b,
                    stud: false,
                    jack: false,
                });
            }
        }
        for o in &layout.openings {
            let mut xs: Vec<_> = header_nodes
                .iter()
                .filter(|(x, _)| *x >= o.x0 && *x <= o.x1)
                .copied()
                .collect();
            xs.sort_by(|a, b| a.0.total_cmp(&b.0));
            for w in xs.windows(2) {
                headers.push(HeaderRec {
                    x0: w[0].0,
                    x1: w[1].0,
                    dofs: [w[0].1[0], w[0].1[1], w[1].1[0], w[1].1[1]],
                    ei: o.header_ei,
                });
            }
        }
        let piers = add_sheathing(&mut model, &layout, &top_nodes, &bottom_nodes)?;
        let (bearing, holds) = add_bearing(&mut model, &layout, &stations, &bottom_nodes)?;
        let pier_comps = if is_nonlinear { piers.len() } else { 0 };
        let mut comps: Vec<PierRec> = piers[..pier_comps].to_vec();
        comps.extend(holds);
        if comps.len() > MAX_WALL_PIERS {
            return Err(StudWallError(format!(
                "{} nonlinear components (piers and hold-downs) exceed the stud-wall limit of {MAX_WALL_PIERS}",
                comps.len()
            )));
        }
        let bearing: Vec<BearingRec> = bearing
            .into_iter()
            .map(|b| BearingRec {
                hold: b.hold.map(|h| pier_comps + h),
                ..b
            })
            .collect();
        let mut plan_ties = Vec::new();
        if layout.bearing.is_some() {
            let tol = 1e-9 * layout.length.max(1.);
            let bottoms = [(0, 0.), (1, layout.length)].into_iter().chain(
                stations
                    .iter()
                    .enumerate()
                    .filter(|(_, s)| s.chord == WallChord::Bottom)
                    .map(|(j, s)| (4 + j, s.x)),
            );
            for (e, x) in bottoms {
                if let Some(b) = bearing.iter().find(|b| (b.x - x).abs() <= tol) {
                    plan_ties.push((e, b.base.clone()));
                }
            }
        }
        let mut plan_tied = vec![0u8; 4 + stations.len()];
        for (e, w) in &plan_ties {
            plan_tied[*e] = 1;
            for &(b, _) in w {
                if plan_tied[b] == 0 {
                    plan_tied[b] = 2;
                }
            }
        }
        let mut foundation_supports = Vec::new();
        if layout.foundation {
            // A bottom station replaces the internal support at its position (for a hold-down there).
            for &(x, d) in &bottom_nodes {
                if x > 0. && x < layout.length && d[1] >= ne {
                    model.fixed_internal.extend([d[0], d[1]]);
                    foundation_supports.push((x, d[0], d[1]));
                }
            }
        }
        let full_stiffness = layout.foundation.then(|| model.k.clone());
        let mut nonlinear = if !comps.is_empty() {
            Some(make_nonlinear_cache(&model, ne, &comps)?)
        } else {
            None
        };
        let (kc, recovery, _ki) = model.condense(ne)?;
        // Horizontal top loads act on the exterior top row when there is a floor zone.
        let ext_top: Vec<(f64, usize)> = if layout.floor.is_some() {
            let mut v = vec![(0., 6), (layout.length, 4)];
            v.extend(
                stations
                    .iter()
                    .enumerate()
                    .filter(|(_, s)| s.chord == WallChord::Top)
                    .map(|(j, s)| (s.x, 2 * (4 + j))),
            );
            v.sort_by(|a, b| a.0.total_cmp(&b.0));
            v
        } else {
            Vec::new()
        };
        let mut load_slots = Vec::with_capacity(cases.len());
        let mut load_recovery = Vec::with_capacity(cases.len());
        let mut full_load_slots = Vec::with_capacity(cases.len());
        for (slot, slot_loads) in cases.iter().enumerate() {
            let mut p = DVector::zeros(model.ndof);
            for l in slot_loads {
                discretize_load(
                    l,
                    slot,
                    &layout,
                    &top_nodes,
                    &top_segments,
                    &ext_top,
                    &mut p,
                )?;
            }
            if let Some(cache) = &mut nonlinear {
                cache
                    .load_strain
                    .push(pier_load_strain(&cache.load_projection, ne, &p));
            }
            let (f, w) = model.condense_load(ne, &p)?;
            if layout.foundation {
                full_load_slots.push(p);
            }
            load_slots.push(f);
            load_recovery.push(w);
        }
        let beam_load_slots = (0..cases.len())
            .map(|slot| {
                beams
                    .iter()
                    .map(|b| {
                        if b.chord == WallChord::Top {
                            segment_load(&cases[slot], b.x0, b.x1, layout.length, true)
                        } else {
                            [0.; 6]
                        }
                    })
                    .collect()
            })
            .collect();
        let axial_load_slots = (0..cases.len())
            .map(|slot| {
                axial_chords
                    .iter()
                    .map(|b| {
                        let v = if b.chord == WallChord::Top {
                            segment_load(
                                &cases[slot],
                                b.x0,
                                b.x1,
                                layout.length,
                                layout.top_distribution == ChordDistribution::Beam,
                            )
                        } else {
                            [0.; 6]
                        };
                        // A floor zone takes the horizontal load on the exterior row, not on the plate.
                        if layout.floor.is_some() {
                            [0., 0.]
                        } else {
                            [v[0], v[3]]
                        }
                    })
                    .collect()
            })
            .collect();
        let material_states = (0..comps.len())
            .map(|i| {
                if i >= pier_comps {
                    layout.bearing.as_ref().and_then(|b| b.hold_down.clone())
                } else if let Some(WallShearLaw::Material(m)) =
                    layout.sheathing.as_ref().map(|s| &s.law)
                {
                    Some(m.clone())
                } else {
                    None
                }
            })
            .collect();
        Ok(Self {
            corners,
            stations,
            nodes,
            layout,
            stiffness: kc,
            load_slots,
            beam_load_slots,
            axial_load_slots,
            recovery,
            load_recovery,
            bars,
            axial_chords,
            beams,
            headers,
            piers,
            top_nodes,
            bottom_nodes,
            zone_bays,
            prepared: None,
            prepared_nonlinear: None,
            full_stiffness,
            full_load_slots,
            foundation_supports,
            bearing,
            plan_ties,
            plan_tied,
            comps,
            pier_comps,
            nonlinear,
            material_states,
            committed_gamma: ArrayVec::new(),
            committed_q: ArrayVec::new(),
            cases,
            pending: None,
            offset: None,
        })
    }
}

// Internal frame builder uses exterior DOFs first and internal DOFs after them.
pub(super) struct FrameBuild {
    pub(super) ndof: usize,
    pub(super) k: DMatrix<f64>,
    pub(super) fixed_internal: Vec<usize>,
}
impl FrameBuild {
    pub(super) fn new(ne: usize) -> Self {
        Self {
            ndof: ne,
            k: DMatrix::zeros(ne, ne),
            fixed_internal: Vec::new(),
        }
    }
    pub(super) fn grow(&mut self) {
        let n = self.ndof + 1;
        let mut k = DMatrix::zeros(n, n);
        k.view_mut((0, 0), (self.ndof, self.ndof))
            .copy_from(&self.k);
        self.k = k;
        self.ndof = n;
    }
    pub(super) fn add_internal_node(&mut self, _x: f64, _y: f64, n: usize) -> [Option<usize>; 3] {
        let mut d = [None; 3];
        for slot in d.iter_mut().take(n) {
            self.grow();
            *slot = Some(self.ndof - 1);
        }
        d
    }
    pub(super) fn add_bar(&mut self, a: usize, b: usize, k: f64) {
        self.k[(a, a)] += k;
        self.k[(b, b)] += k;
        self.k[(a, b)] -= k;
        self.k[(b, a)] -= k;
    }
    pub(super) fn add_beam3(&mut self, a: [usize; 3], b: [usize; 3], ea: f64, ei: f64, l: f64) {
        let k = beam6(ea, ei, l);
        let d = [a[0], a[1], a[2], b[0], b[1], b[2]];
        for i in 0..6 {
            for j in 0..6 {
                self.k[(d[i], d[j])] += k[i][j];
            }
        }
    }
    pub(super) fn add_beam4(&mut self, d: [usize; 4], ei: f64, l: f64) {
        let k = beam4(ei, l);
        for i in 0..4 {
            for j in 0..4 {
                self.k[(d[i], d[j])] += k[i][j];
            }
        }
    }
    pub(super) fn condense(&self, ne: usize) -> Result<CondensedSystem, StudWallError> {
        let ni = self.ndof - ne;
        let kee = self.k.view((0, 0), (ne, ne)).into_owned();
        if ni == 0 {
            return Ok((kee, DMatrix::zeros(0, ne), DMatrix::zeros(0, 0)));
        }
        let free: Vec<_> = (ne..self.ndof)
            .filter(|d| !self.fixed_internal.contains(d))
            .collect();
        let mut recovery = DMatrix::zeros(ni, ne);
        if free.is_empty() {
            return Ok((kee, recovery, DMatrix::zeros(0, 0)));
        }
        let kei = DMatrix::from_fn(ne, free.len(), |i, j| self.k[(i, free[j])]);
        let kie = DMatrix::from_fn(free.len(), ne, |i, j| self.k[(free[i], j)]);
        let kii = DMatrix::from_fn(free.len(), free.len(), |i, j| self.k[(free[i], free[j])]);
        let lu = kii.clone().lu();
        let Some(r) = lu.solve(&kie) else {
            return Err(StudWallError("internal stiffness K_ii is singular".into()));
        };
        let kc = kee - &kei * &r;
        if !kc.iter().all(|x| x.is_finite()) {
            return Err(StudWallError("internal stiffness K_ii is singular".into()));
        }
        for (j, &d) in free.iter().enumerate() {
            recovery.row_mut(d - ne).copy_from(&r.row(j));
        }
        Ok((kc, recovery, kii))
    }
    pub(super) fn condense_load(
        &self,
        ne: usize,
        p: &DVector<f64>,
    ) -> Result<(DVector<f64>, DVector<f64>), StudWallError> {
        let pe = p.rows(0, ne).into_owned();
        let free: Vec<_> = (ne..self.ndof)
            .filter(|d| !self.fixed_internal.contains(d))
            .collect();
        if free.is_empty() {
            return Ok((pe, DVector::zeros(0)));
        }
        let kii = DMatrix::from_fn(free.len(), free.len(), |i, j| self.k[(free[i], free[j])]);
        let pi = DVector::from_iterator(free.len(), free.iter().map(|&d| p[d]));
        let w_free = kii
            .lu()
            .solve(&pi)
            .ok_or_else(|| StudWallError("internal stiffness K_ii is singular".into()))?;
        let kei = DMatrix::from_fn(ne, free.len(), |i, j| self.k[(i, free[j])]);
        let mut w = DVector::zeros(self.ndof - ne);
        for (j, &d) in free.iter().enumerate() {
            w[d - ne] = w_free[j];
        }
        Ok((pe - kei * &w_free, w))
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn make_chord_nodes<K: NodeKey>(
    model: &mut FrameBuild,
    xs: &[f64],
    y: f64,
    chord: WallChord,
    nslot: usize,
    zone: bool,
    corners: &[K; 4],
    stations: &[WallStation<K>],
    nodes: &[K],
) -> Vec<(f64, [usize; 3])> {
    xs.iter()
        .map(|&x| {
            // Each exterior node contributes translation slots; beam rotations remain internal.
            let external = if x.abs() < 1e-12 {
                Some(if chord == WallChord::Top {
                    (corners[3], 3)
                } else {
                    (corners[0], 0)
                })
            } else if (x - xs.last().copied().unwrap_or(x)).abs() < 1e-12 {
                Some(if chord == WallChord::Top {
                    (corners[2], 2)
                } else {
                    (corners[1], 1)
                })
            } else {
                stations
                    .iter()
                    .find(|s| s.chord == chord && (s.x - x).abs() < 1e-12)
                    .and_then(|s| nodes.iter().position(|n| *n == s.node).map(|i| (s.node, i)))
            };
            let mut d = [None; 3];
            if let (Some((_, i)), true) = (external, zone && chord == WallChord::Top) {
                // Under a floor zone the framing node keeps its own `u` and shares the exterior node's `v`
                // (the rigid post); the exterior `u` is reached through the zone bays.
                model.grow();
                d[0] = Some(model.ndof - 1);
                d[1] = Some(2 * i + 1);
                if nslot == 3 {
                    model.grow();
                    d[2] = Some(model.ndof - 1);
                }
            } else if let Some((_, i)) = external {
                d[0] = Some(2 * i);
                d[1] = Some(2 * i + 1);
                if nslot == 3 {
                    model.grow();
                    d[2] = Some(model.ndof - 1);
                }
            } else {
                d = model.add_internal_node(x, y, nslot);
            }
            (
                x,
                [d[0].unwrap(), d[1].unwrap(), d[2].unwrap_or(usize::MAX)],
            )
        })
        .collect()
}
pub(super) type ChordSegment = [(f64, [usize; 3]); 2];

/// A plate segment over an opening is hinged at the opening edges: its end there gets its own rotation DOF,
/// so the plate ties the piers on either side axially only (no frame action through plate bending).
pub(super) fn released_segment(
    model: &mut FrameBuild,
    pair: &[(f64, [usize; 3])],
    openings: &[WallOpening],
) -> ChordSegment {
    let mut seg = [pair[0], pair[1]];
    let tol = 1e-12 * pair[1].0.abs().max(1.);
    for o in openings {
        if pair[0].0 < o.x0 - tol || pair[1].0 > o.x1 + tol || seg[0].1[2] == usize::MAX {
            continue;
        }
        for (end, edge) in [(0, o.x0), (1, o.x1)] {
            if (seg[end].0 - edge).abs() <= tol {
                model.grow();
                seg[end].1[2] = model.ndof - 1;
            }
        }
    }
    seg
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add_chord_elements(
    model: &mut FrameBuild,
    nodes: &[(f64, [usize; 3])],
    ea: f64,
    ei: f64,
    dist: ChordDistribution,
    chord: WallChord,
    openings: &[WallOpening],
    beams: &mut Vec<BeamRec>,
    axial: &mut Vec<ChordAxialRec>,
) -> Vec<ChordSegment> {
    let mut segments = Vec::new();
    for pair in nodes.windows(2) {
        let l = pair[1].0 - pair[0].0;
        if chord == WallChord::Top && dist == ChordDistribution::Direct {
            model.add_bar(pair[0].1[0], pair[1].1[0], ea / l);
            segments.push([pair[0], pair[1]]);
        } else {
            let [a, b] = released_segment(model, pair, openings);
            add_beam3(model, a.0, b.0, a.1, b.1, ea, ei, chord, beams);
            segments.push([a, b]);
        }
        axial.push(ChordAxialRec {
            x0: pair[0].0,
            x1: pair[1].0,
            a: pair[0].1[0],
            b: pair[1].1[0],
            ea,
            chord,
        });
    }
    segments
}
#[allow(clippy::too_many_arguments)]
pub(super) fn add_beam3(
    model: &mut FrameBuild,
    x0: f64,
    x1: f64,
    a: [usize; 3],
    b: [usize; 3],
    ea: f64,
    ei: f64,
    chord: WallChord,
    beams: &mut Vec<BeamRec>,
) {
    let dofs = [a[0], a[1], a[2], b[0], b[1], b[2]];
    model.add_beam3(a, b, ea, ei, x1 - x0);
    beams.push(BeamRec {
        x0,
        x1,
        dofs,
        ei,
        chord,
    });
}
pub(super) fn beam6(ea: f64, ei: f64, l: f64) -> [[f64; 6]; 6] {
    let mut k = [[0.; 6]; 6];
    let a = ea / l;
    for (i, j, v) in [(0, 0, a), (0, 3, -a), (3, 0, -a), (3, 3, a)] {
        k[i][j] += v;
    }
    let b = 12. * ei / l.powi(3);
    let c = 6. * ei / l.powi(2);
    let d = 4. * ei / l;
    let e = 2. * ei / l;
    let m = [[b, c, -b, c], [c, d, -c, e], [-b, -c, b, -c], [c, e, -c, d]];
    for i in 0..4 {
        for j in 0..4 {
            k[[1, 2, 4, 5][i]][[1, 2, 4, 5][j]] += m[i][j];
        }
    }
    k
}
pub(super) fn beam4(ei: f64, l: f64) -> [[f64; 4]; 4] {
    let b = 12. * ei / l.powi(3);
    let c = 6. * ei / l.powi(2);
    let d = 4. * ei / l;
    let e = 2. * ei / l;
    [[b, c, -b, c], [c, d, -c, e], [-b, -c, b, -c], [c, e, -c, d]]
}
pub(super) fn find_chord(nodes: &[(f64, [usize; 3])], x: f64) -> Option<[usize; 3]> {
    nodes
        .iter()
        .find(|(xx, _)| (*xx - x).abs() < 1e-10)
        .map(|(_, d)| *d)
}

pub(super) fn add_sheathing(
    model: &mut FrameBuild,
    l: &StudWallLayout,
    top: &[(f64, [usize; 3])],
    bottom: &[(f64, [usize; 3])],
) -> Result<Vec<PierRec>, StudWallError> {
    let Some(sh) = &l.sheathing else {
        return Ok(Vec::new());
    };
    let ga = match &sh.law {
        WallShearLaw::Linear { ga } => Some(*ga),
        WallShearLaw::O86 { .. } | WallShearLaw::Material(_) => None,
    };
    let ranges = pier_bounds(l);
    let mut out = Vec::new();
    for (xa, xb) in ranges {
        let width = xb - xa;
        if width <= 0. {
            continue;
        }
        let mut beta = DVector::zeros(model.ndof);
        for (nodes, sign) in [(top, 1.), (bottom, -1.)] {
            for pair in nodes.windows(2) {
                let a = pair[0].0.max(xa);
                let b = pair[1].0.min(xb);
                if b <= a {
                    continue;
                }
                let len = pair[1].0 - pair[0].0;
                let t0 = (a - pair[0].0) / len;
                let t1 = (b - pair[0].0) / len;
                let w0 = (b - a) * (2. - t0 - t1) / 2.;
                let w1 = (b - a) * (t0 + t1) / 2.;
                beta[pair[0].1[0]] += sign * w0 / (width * l.height);
                beta[pair[1].1[0]] += sign * w1 / (width * l.height);
            }
        }
        for (x, sgn) in [(xa, -1.), (xb, 1.)] {
            let t = find_chord(top, x)
                .ok_or_else(|| StudWallError("pier boundary has no top stud node".into()))?;
            let b = find_chord(bottom, x)
                .ok_or_else(|| StudWallError("pier boundary has no bottom stud node".into()))?;
            beta[t[1]] += sgn * 0.5 / width;
            beta[b[1]] += sgn * 0.5 / width;
        }
        let mut c = l.height.powi(3) / 12. * (1. / stud_at(l, xa).ea + 1. / stud_at(l, xb).ea);
        for nodes in [top, bottom] {
            let ea = if std::ptr::eq(nodes, top) {
                l.top_ea
            } else {
                l.bottom_ea
            };
            for p in nodes.windows(2) {
                let d = (p[1].0.min(xb) - p[0].0.max(xa)).max(0.);
                c += d.powi(3) / (12. * ea);
            }
        }
        let area = width * l.height;
        if let Some(ga) = ga {
            let fac = area * area / (area / ga + c);
            for i in 0..model.ndof {
                for j in 0..model.ndof {
                    model.k[(i, j)] += fac * beta[i] * beta[j];
                }
            }
        }
        out.push(PierRec {
            xa,
            xb,
            beta,
            c,
            area,
        });
    }
    Ok(out)
}
/// Springs from every bottom-plate bearing position to the interpolated base surface ([`WallBearing`]): `kh` and
/// `kv` into the frame stiffness, and at each pier end (with a hold-down law) a nonlinear component in place of
/// `kv`. Returns the bearing positions (`hold` indexing the returned components) and the components.
/// The global translation slots along a wall's unit normal (`ex × ey`), with the normal's component on each.
pub(super) fn normal_slots(map: &PlaneMap) -> ArrayVec<(u8, f64), 3> {
    let (a, b) = (map.ex, map.ey);
    let n = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    (0..3)
        .filter(|&i| n[i].abs() > 1e-12)
        .map(|i| (i as u8, n[i]))
        .collect()
}
pub(super) fn add_bearing<K: NodeKey>(
    model: &mut FrameBuild,
    l: &StudWallLayout,
    stations: &[WallStation<K>],
    bottom: &[(f64, [usize; 3])],
) -> Result<(Vec<BearingRec>, Vec<PierRec>), StudWallError> {
    let Some(b) = &l.bearing else {
        return Ok((Vec::new(), Vec::new()));
    };
    let tol = 1e-9 * l.length.max(1.);
    // Base positions: sorted, nodes at one `x` grouped (exterior node index 4 + station index).
    let mut base: Vec<(f64, Vec<usize>)> = Vec::new();
    let mut pts: Vec<(f64, usize)> = stations
        .iter()
        .enumerate()
        .filter(|(_, s)| s.chord == WallChord::Base)
        .map(|(j, s)| (s.x, 4 + j))
        .collect();
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    for (x, e) in pts {
        match base.last_mut() {
            Some((bx, g)) if (x - *bx).abs() <= tol => g.push(e),
            _ => base.push((x, vec![e])),
        }
    }
    let weights = |x: f64| -> Option<Vec<(usize, f64)>> {
        let share = |g: &[usize], w: f64| {
            g.iter()
                .map(|&e| (e, w / g.len() as f64))
                .collect::<Vec<_>>()
        };
        if let Some((_, g)) = base.iter().find(|(bx, _)| (x - bx).abs() <= tol) {
            return Some(share(g, 1.));
        }
        let i = base.windows(2).position(|w| w[0].0 < x && x < w[1].0)?;
        let ((x0, g0), (x1, g1)) = (&base[i], &base[i + 1]);
        let t = (x - x0) / (x1 - x0);
        let mut v = share(g0, 1. - t);
        v.extend(share(g1, t));
        Some(v)
    };
    let ends = l.pier_ends();
    let mut recs = Vec::new();
    let mut holds = Vec::new();
    for &(x, d) in bottom {
        let w = weights(x).ok_or_else(|| {
            StudWallError(format!(
                "bearing position x = {x} lies outside the base nodes' range"
            ))
        })?;
        let hold = b.hold_down.is_some() && ends.iter().any(|e| (e - x).abs() <= tol);
        for (slot, dof, k) in [(0, d[0], b.kh), (1, d[1], b.kv)] {
            let mut idx = vec![(dof, 1.)];
            idx.extend(w.iter().map(|&(e, wt)| (2 * e + slot, -wt)));
            if slot == 1 && hold {
                let mut beta = DVector::zeros(model.ndof);
                for &(i, v) in &idx {
                    beta[i] += v;
                }
                holds.push(PierRec {
                    xa: x,
                    xb: x,
                    beta,
                    c: 0.,
                    area: 1.,
                });
            } else {
                for &(i, vi) in &idx {
                    for &(j, vj) in &idx {
                        model.k[(i, j)] += k * vi * vj;
                    }
                }
            }
        }
        recs.push(BearingRec {
            x,
            u: d[0],
            v: d[1],
            base: w,
            hold: hold.then(|| holds.len() - 1),
        });
    }
    Ok((recs, holds))
}
pub(super) fn stud_at(l: &StudWallLayout, x: f64) -> &WallStud {
    l.studs.iter().find(|s| (s.x - x).abs() < 1e-10).unwrap()
}
pub(super) fn pier_bounds(l: &StudWallLayout) -> Vec<(f64, f64)> {
    let mut runs = Vec::new();
    let mut a = 0.;
    for o in &l.openings {
        runs.push((a, o.x0));
        a = o.x1;
    }
    runs.push((a, l.length));
    runs.into_iter()
        .filter_map(|(a, b)| {
            let xs: Vec<_> = l
                .studs
                .iter()
                .filter(|s| s.x >= a && s.x <= b)
                .map(|s| s.x)
                .collect();
            (xs.len() > 1).then(|| (*xs.first().unwrap(), *xs.last().unwrap()))
        })
        .collect()
}
