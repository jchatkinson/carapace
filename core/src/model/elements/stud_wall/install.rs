//! Staged sheathing installation (`docs/stud-wall-panel.md` §6.1, §6.4): a wall whose shear fields (pier sheathing and floor-zone
//! bays) carry nothing until installed, so gravity goes down through studs, plates, headers, jacks and posts alone (the stick
//! takedown), and the sheathing takes only what is applied after installation (lateral load).
//!
//! A staged wall is built twice: as used before installation, with every shear field at [`STAGED_SCALE`] of its stiffness (the
//! sheathing law's initial tangent as a linear law; the floor zone's `ga` scaled), which keeps the racking mode from being singular
//! without carrying gravity; and as installed. [`StudWall::install_sheathing`] swaps the installed build in at the committed state
//! with offsets that keep the exterior forces, the internal displacements and the hold-downs' strains continuous, and records
//! each shear field's strain `γ0` there: from then on a field responds to `γ - γ0`, stress-free at installation.

use super::*;
use nalgebra::SVector;

/// The stiffness factor of a staged wall's shear fields before installation.
pub const STAGED_SCALE: f64 = 1e-6;

/// What installation adds to the installed build so its state continues the uninstalled one's.
#[derive(Debug, Clone)]
pub(super) struct InstallOffset {
    /// Added to the local resisting force (`K u + C q`), loads included: `r_before(u0) - r_after(u0)`.
    pub force: DVector<f64>,
    /// Added to the recovered internal displacements.
    pub internal: DVector<f64>,
    /// Added to the components' imposed strains `y`.
    pub strain: Vec<f64>,
    /// Per component: the strain its law is measured from (a nonlinear pier's `γ0`; 0 for a hold-down).
    pub gamma0: Vec<f64>,
    /// Per pier: its shear strain at installation.
    pub pier_gamma0: Vec<f64>,
    /// Per floor-zone bay: its shear strain at installation.
    pub bay_gamma0: Vec<f64>,
    /// Added to the foundation's full-system residual (internal support reactions).
    pub full: Option<DVector<f64>>,
}

impl<K: NodeKey> StudWall<K> {
    /// This wall with its shear fields uninstalled until [`Self::install_sheathing`]. A wall with no sheathing and no floor
    /// zone is returned as it is.
    pub fn stage_sheathing(self) -> Result<Self, StudWallError> {
        if self.pending.is_some()
            || self.offset.is_some()
            || (self.layout.sheathing.is_none() && self.layout.floor.is_none())
        {
            return Ok(self);
        }
        let mut l = self.layout.clone();
        if let Some(sh) = &mut l.sheathing {
            let material = match &sh.law {
                WallShearLaw::Material(m) => Some(m),
                _ => None,
            };
            let k0 = shear_law_response(&sh.law, material, 0.).1;
            sh.law = WallShearLaw::Linear {
                ga: STAGED_SCALE * k0,
            };
        }
        if let Some(z) = &mut l.floor {
            z.ga *= STAGED_SCALE;
        }
        let mut before = Self::build(self.corners, self.stations.clone(), l, self.cases.clone())?;
        before.pending = Some(Box::new(self));
        Ok(before)
    }

    /// Whether this wall's shear fields wait for installation.
    pub fn is_staged(&self) -> bool {
        self.pending.is_some()
    }

    /// Installs a staged wall's shear fields at the nodes' committed state with the load cases `cases` in effect (see the
    /// module doc). Returns whether there was anything to install.
    pub fn install_sheathing<const NDIM: usize, const NDOF: usize>(
        &mut self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        cases: &[WallCaseLoad],
    ) -> bool {
        let Some(mut after) = self.pending.take() else {
            return false;
        };
        let u0 = self.local_displacement(nodes);
        // The uninstalled state: its hold-downs (its only components) and internal displacements.
        let (g_before, q_before) = match self
            .nonlinear
            .as_ref()
            .and_then(|_| self.nonlinear_trial(u0.as_slice(), cases))
        {
            Some((g, q, _)) => (g.to_vec(), q.to_vec()),
            None => (Vec::new(), Vec::new()),
        };
        let qb = self.nonlinear.is_some().then_some(q_before.as_slice());
        let ui0 = self.internal_displacement(&u0, cases, qb);
        let r_before = self.resisting_local(&u0, qb) - self.load_local(cases);
        let d0 = DVector::from_iterator(u0.len() + ui0.len(), u0.iter().chain(ui0.iter()).copied());

        // The installed build's components: its nonlinear piers (from their strain now, stress-free), then the hold-downs as they are.
        let pier_gamma0: Vec<f64> = after.piers.iter().map(|p| p.beta.dot(&d0)).collect();
        let pc = after.pier_comps;
        let mut g_after: Vec<f64> = pier_gamma0[..pc].to_vec();
        let mut q_after = vec![0.; pc];
        g_after.extend(g_before.iter().copied());
        q_after.extend(q_before.iter().copied());
        let count = after.comps.len();
        let strain = match &after.nonlinear {
            Some(cache) => {
                let y = after.imposed_strain(u0.as_slice(), cases);
                (0..count)
                    .map(|k| {
                        g_after[k]
                            + (0..count)
                                .map(|l| cache.compliance[(k, l)] * q_after[l])
                                .sum::<f64>()
                            - y[k]
                    })
                    .collect()
            }
            None => Vec::new(),
        };
        let qa = after.nonlinear.is_some().then_some(q_after.as_slice());
        let internal = &ui0 - after.internal_displacement(&u0, cases, qa);
        let force = r_before - (after.resisting_local(&u0, qa) - after.load_local(cases));
        let bay_gamma0 = after
            .zone_bays
            .iter()
            .map(|b| {
                let z = after.layout.floor.expect("bays come from a floor zone");
                let u = SVector::<f64, 8>::from_fn(|i, _| d0[b.dofs[i]]);
                floor_zone::bay_operator(b.x1 - b.x0, z.depth).dot(&u)
            })
            .collect();
        let full = after.full_stiffness.as_ref().map(|k_after| {
            let residual = |w: &Self, k: &DMatrix<f64>, q: &[f64]| {
                let mut r = k * &d0;
                for (i, c) in w.comps.iter().enumerate() {
                    r += &c.beta * (c.area * q[i]);
                }
                for c in cases {
                    r -= &w.full_load_slots[c.case as usize] * c.factor;
                }
                r
            };
            residual(
                self,
                self.full_stiffness.as_ref().expect("same foundation"),
                &q_before,
            ) - residual(&after, k_after, &q_after)
        });
        let hold_states: Vec<_> = self.material_states.clone();
        for (j, m) in hold_states.into_iter().enumerate() {
            after.material_states[pc + j] = m;
        }
        after.committed_gamma = g_after.iter().copied().collect();
        after.committed_q = q_after.iter().copied().collect();
        let mut gamma0 = pier_gamma0[..pc].to_vec();
        gamma0.resize(count, 0.);
        after.offset = Some(InstallOffset {
            force,
            internal,
            strain,
            gamma0,
            pier_gamma0,
            bay_gamma0,
            full,
        });
        after.prepare(nodes);
        *self = *after;
        true
    }

    /// The components' imposed strains `y` (exterior displacements and loads, plus an installation offset).
    pub(super) fn imposed_strain(
        &self,
        ue: &[f64],
        cases: &[WallCaseLoad],
    ) -> ArrayVec<f64, MAX_WALL_PIERS> {
        let mut y = ArrayVec::new();
        let Some(cache) = &self.nonlinear else {
            return y;
        };
        for i in 0..self.comps.len() {
            let mut yi = (0..ue.len())
                .map(|j| cache.strain_map[(i, j)] * ue[j])
                .sum::<f64>();
            for c in cases {
                yi += cache.load_strain[c.case as usize][i] * c.factor;
            }
            y.push(yi + self.offset.as_ref().map_or(0., |o| o.strain[i]));
        }
        y
    }

    /// The strain a component's law is measured from: a nonlinear pier's strain at installation, else 0.
    pub(super) fn law_origin(&self, i: usize) -> f64 {
        self.offset.as_ref().map_or(0., |o| o.gamma0[i])
    }

    /// Internal displacements for exterior displacements `ue`, load cases and (nonlinear walls) component forces `q`.
    pub(super) fn internal_displacement(
        &self,
        ue: &DVector<f64>,
        cases: &[WallCaseLoad],
        q: Option<&[f64]>,
    ) -> DVector<f64> {
        let mut ui = -(&self.recovery * ue);
        for c in cases {
            ui += &self.load_recovery[c.case as usize] * c.factor;
        }
        if let (Some(cache), Some(q)) = (&self.nonlinear, q) {
            for i in 0..ui.len() {
                for (j, qj) in q.iter().enumerate() {
                    ui[i] -= cache.internal_recovery[(i, j)] * qj;
                }
            }
        }
        if let Some(o) = &self.offset {
            ui += &o.internal;
        }
        ui
    }

    /// The local resisting force `K u + C q` (loads not subtracted), plus an installation offset.
    pub(super) fn resisting_local(&self, ue: &DVector<f64>, q: Option<&[f64]>) -> DVector<f64> {
        let mut r = &self.stiffness * ue;
        if let (Some(cache), Some(q)) = (&self.nonlinear, q) {
            for i in 0..r.len() {
                for (j, qj) in q.iter().enumerate() {
                    r[i] += cache.coupling[(i, j)] * qj;
                }
            }
        }
        if let Some(o) = &self.offset {
            r += &o.force;
        }
        r
    }

    /// The local equivalent load of `cases`.
    pub(super) fn load_local(&self, cases: &[WallCaseLoad]) -> DVector<f64> {
        let mut f = DVector::zeros(self.nodes.len() * 2);
        for c in cases {
            f += &self.load_slots[c.case as usize] * c.factor;
        }
        f
    }
}
