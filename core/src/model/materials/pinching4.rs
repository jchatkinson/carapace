use super::Material;

/// `Pinching4Material`'s `Cstate`/`Tstate` — an `int` 0..4 in the OpenSees
/// source, modelled here as a real enum.
///
/// Unlike every other material in this catalog, this is *not* a
/// positive/negative/interior split:
///
/// - [`NearOrigin`](Pinching4State::NearOrigin) (`0`) — before the first
///   real excursion in either direction, i.e. still inside the tiny
///   synthesized near-origin bracket `[-u, +u]` that `SetEnvelope` builds.
/// - [`Positive`](Pinching4State::Positive) (`1`) — riding the positive
///   (damaged) envelope.
/// - [`Negative`](Pinching4State::Negative) (`2`) — riding the negative
///   (damaged) envelope.
/// - [`Trilinear3`](Pinching4State::Trilinear3) (`3`) — the pinched
///   trilinear path connecting a negative excursion back *toward* (not
///   necessarily reaching) the positive side.
/// - [`Trilinear4`](Pinching4State::Trilinear4) (`4`) — the mirror image:
///   a positive excursion heading back toward the negative side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pinching4State {
    NearOrigin,
    Positive,
    Negative,
    Trilinear3,
    Trilinear4,
}

/// `Pinching4Material`'s `DmgCyc` flag (a raw `int` in OpenSees): whether
/// the second damage term is driven by the dissipated-energy ratio
/// (`DmgCyc == 0`) or by an accumulated half-cycle count (`DmgCyc == 1`).
///
/// Note the two are not symmetric in the source: the energy-based term is
/// additionally gated on `Tenergy > elasticStrainEnergy`, the cycle-based
/// one is not (see `update_dmg`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pinching4DmgCyc {
    EnergyBased,
    CycleBased,
}

/// The full field set behind `Material::Pinching4` (see that variant's doc
/// comment for why it's boxed). `Copy` for the same reason
/// `HystereticFields` is — every field is an `f64`, a small `f64` array, or
/// a fieldless enum, so one cheap copy out of the `Box` per `evaluate` call
/// beats borrow-juggling a destructure.
///
/// # What is stored vs. what is derived
///
/// OpenSees carries several *derived* arrays/scalars as members, refreshed
/// in `commitState`. None of them are stored here, because each is an exact
/// function of what is:
///
/// - `envlpPosDamgdStress`/`envlpNegDamgdStress` — always
///   `envlp*Stress[i] * (1 - gammaFUsed)` for whatever `gammaFUsed` was
///   current when they were last written, so one scalar per side recovers
///   them exactly. They *are* rebuilt as trial-local `[f64; 6]` arrays
///   inside `evaluate` (see the note there about the per-side lag).
/// - `kElasticPosDamgd`/`kElasticNegDamgd` — `kElastic* * (1 - gammaKUsed)`,
///   same argument.
/// - `uMaxDamgd`/`uMinDamgd` — `C*StrainDmnd * (1 + CgammaD)`; `commitState`
///   writes them from exactly those two committed fields, and nothing
///   mutates them mid-step, so they are recomputed at the top of `evaluate`.
///
/// The backbone is kept as the 6-point envelope arrays `SetEnvelope`
/// produces rather than the 16 raw `stress1p..strain4n` inputs: points 1-4
/// of each array *are* those inputs, and points 0/5 (the synthesized
/// near-origin point and the far extrapolated point) are the only things
/// the raw values are otherwise needed for.
///
/// History fields mirror the `C*` half of OpenSees' `T*`/`C*` pairs:
/// `state`/`strain`/`stress`/`strain_rate`, the current state's strain
/// bracket and the stresses at its ends
/// (`low_state_*`/`hgh_state_*`), the peak strain demand on each side
/// (`min_strain_dmnd`/`max_strain_dmnd`), cumulative dissipated `energy`,
/// the three damage factors `gamma_k`/`gamma_d`/`gamma_f`, the two
/// *last-transition* damage factors `gamma_k_used`/`gamma_f_used`, and the
/// accumulated half-cycle count `n_cycle`.
#[derive(Debug, Clone, Copy)]
pub struct Pinching4Fields {
    // --- envelope (fixed at construction) ---
    pub envlp_pos_strain: [f64; 6],
    pub envlp_pos_stress: [f64; 6],
    pub envlp_neg_strain: [f64; 6],
    pub envlp_neg_stress: [f64; 6],
    pub k_elastic_pos: f64,
    pub k_elastic_neg: f64,
    pub energy_capacity: f64,
    // --- pinching parameters (fixed) ---
    pub r_disp_p: f64,
    pub r_force_p: f64,
    pub u_force_p: f64,
    pub r_disp_n: f64,
    pub r_force_n: f64,
    pub u_force_n: f64,
    // --- damage parameters (fixed) ---
    pub gamma_k_params: [f64; 4],
    pub gamma_k_limit: f64,
    pub gamma_d_params: [f64; 4],
    pub gamma_d_limit: f64,
    pub gamma_f_params: [f64; 4],
    pub gamma_f_limit: f64,
    pub dmg_cyc: Pinching4DmgCyc,
    // --- committed history ---
    pub state: Pinching4State,
    pub strain: f64,
    pub stress: f64,
    pub strain_rate: f64,
    pub low_state_strain: f64,
    pub low_state_stress: f64,
    pub hgh_state_strain: f64,
    pub hgh_state_stress: f64,
    pub min_strain_dmnd: f64,
    pub max_strain_dmnd: f64,
    pub energy: f64,
    pub gamma_k: f64,
    pub gamma_d: f64,
    pub gamma_f: f64,
    pub gamma_k_used: f64,
    pub gamma_f_used: f64,
    pub n_cycle: f64,
}

/// Everything `Pinching4Material`'s `setTrialStrain` treats as a mutable
/// `T*` member for the duration of one call. Bundled into a struct purely
/// so `getstate`/`update_dmg` can take one `&mut` instead of a dozen.
///
/// `pos_damgd_stress`/`neg_damgd_stress` are trial-local copies of
/// OpenSees' `envlpPosDamgdStress`/`envlpNegDamgdStress`, seeded from the
/// committed `gamma_f_used` and refreshed *one side at a time* by
/// `getstate` — faithfully reproducing a real asymmetry in the source: a
/// transition out of state 1 rewrites only the negative side's damaged
/// stresses (and only `kElasticPosDamgd`), leaving the other side one
/// commit behind until `commitState` resynchronizes both. The same applies
/// to `k_elastic_pos_damgd`/`k_elastic_neg_damgd`.
struct Trial {
    state: Pinching4State,
    low_strain: f64,
    low_stress: f64,
    hgh_strain: f64,
    hgh_stress: f64,
    min_dmnd: f64,
    max_dmnd: f64,
    energy: f64,
    gamma_k: f64,
    gamma_d: f64,
    gamma_f: f64,
    gamma_k_used: f64,
    gamma_f_used: f64,
    k_elastic_pos_damgd: f64,
    k_elastic_neg_damgd: f64,
    pos_damgd_stress: [f64; 6],
    neg_damgd_stress: [f64; 6],
    n_cycle: f64,
}

impl Trial {
    /// `Pinching4Material::posEnvlpStress` — a piecewise search over the
    /// 6-point positive envelope, using the *damaged* stresses but the
    /// undamaged strains (strength degradation scales force, not
    /// deformation). Ported literally, `k == 0.0` sentinel and all.
    fn pos_envlp_stress(&self, p: &Pinching4Fields, u: f64) -> f64 {
        let mut k = 0.0;
        let mut f = 0.0;
        let mut i = 0;
        while k == 0.0 && i <= 4 {
            if u <= p.envlp_pos_strain[i + 1] {
                k = (self.pos_damgd_stress[i + 1] - self.pos_damgd_stress[i])
                    / (p.envlp_pos_strain[i + 1] - p.envlp_pos_strain[i]);
                f = self.pos_damgd_stress[i] + (u - p.envlp_pos_strain[i]) * k;
            }
            i += 1;
        }
        if k == 0.0 {
            k = (self.pos_damgd_stress[5] - self.pos_damgd_stress[4]) / (p.envlp_pos_strain[5] - p.envlp_pos_strain[4]);
            f = self.pos_damgd_stress[5] + k * (u - p.envlp_pos_strain[5]);
        }
        f
    }

    /// `Pinching4Material::posEnvlpTangent`.
    fn pos_envlp_tangent(&self, p: &Pinching4Fields, u: f64) -> f64 {
        let mut k = 0.0;
        let mut i = 0;
        while k == 0.0 && i <= 4 {
            if u <= p.envlp_pos_strain[i + 1] {
                k = (self.pos_damgd_stress[i + 1] - self.pos_damgd_stress[i])
                    / (p.envlp_pos_strain[i + 1] - p.envlp_pos_strain[i]);
            }
            i += 1;
        }
        if k == 0.0 {
            k = (self.pos_damgd_stress[5] - self.pos_damgd_stress[4]) / (p.envlp_pos_strain[5] - p.envlp_pos_strain[4]);
        }
        k
    }

    /// `Pinching4Material::negEnvlpStress`.
    fn neg_envlp_stress(&self, p: &Pinching4Fields, u: f64) -> f64 {
        let mut k = 0.0;
        let mut f = 0.0;
        let mut i = 0;
        while k == 0.0 && i <= 4 {
            if u >= p.envlp_neg_strain[i + 1] {
                k = (self.neg_damgd_stress[i] - self.neg_damgd_stress[i + 1])
                    / (p.envlp_neg_strain[i] - p.envlp_neg_strain[i + 1]);
                f = self.neg_damgd_stress[i + 1] + (u - p.envlp_neg_strain[i + 1]) * k;
            }
            i += 1;
        }
        if k == 0.0 {
            k = (self.neg_damgd_stress[4] - self.neg_damgd_stress[5]) / (p.envlp_neg_strain[4] - p.envlp_neg_strain[5]);
            f = self.neg_damgd_stress[5] + k * (u - p.envlp_neg_strain[5]);
        }
        f
    }

    /// `Pinching4Material::negEnvlpTangent`.
    fn neg_envlp_tangent(&self, p: &Pinching4Fields, u: f64) -> f64 {
        let mut k = 0.0;
        let mut i = 0;
        while k == 0.0 && i <= 4 {
            if u >= p.envlp_neg_strain[i + 1] {
                k = (self.neg_damgd_stress[i] - self.neg_damgd_stress[i + 1])
                    / (p.envlp_neg_strain[i] - p.envlp_neg_strain[i + 1]);
            }
            i += 1;
        }
        if k == 0.0 {
            k = (self.neg_damgd_stress[4] - self.neg_damgd_stress[5]) / (p.envlp_neg_strain[4] - p.envlp_neg_strain[5]);
        }
        k
    }
}

impl Material {
    /// Ported from `Pinching4Material`'s full (asymmetric-backbone)
    /// constructor. OpenSees' second overload — the one that mirrors the
    /// positive backbone onto the negative side and copies
    /// `rDispP`/`rForceP`/`uForceP` across — is not ported, for the same
    /// reason `Material::hysteretic` skipped its 2-point overload: a caller
    /// can pass the negated values directly.
    ///
    /// The four-term damage parameter groups are taken as `[f64; 4]` arrays
    /// (`gammaK1..gammaK4` etc.) rather than 12 separate scalars — a purely
    /// mechanical regrouping of OpenSees' flat signature that keeps this
    /// under 30 arguments instead of 39.
    ///
    /// Debug builds assert the backbone is one-to-one (all positive strains
    /// `> 0`, all negative strains `< 0`); OpenSees only prints a warning
    /// and carries on with a nonsensical envelope.
    #[allow(clippy::too_many_arguments)]
    pub fn pinching4(
        stress1p: f64,
        strain1p: f64,
        stress2p: f64,
        strain2p: f64,
        stress3p: f64,
        strain3p: f64,
        stress4p: f64,
        strain4p: f64,
        stress1n: f64,
        strain1n: f64,
        stress2n: f64,
        strain2n: f64,
        stress3n: f64,
        strain3n: f64,
        stress4n: f64,
        strain4n: f64,
        r_disp_p: f64,
        r_force_p: f64,
        u_force_p: f64,
        r_disp_n: f64,
        r_force_n: f64,
        u_force_n: f64,
        gamma_k_params: [f64; 4],
        gamma_k_limit: f64,
        gamma_d_params: [f64; 4],
        gamma_d_limit: f64,
        gamma_f_params: [f64; 4],
        gamma_f_limit: f64,
        gamma_e: f64,
        dmg_cyc: Pinching4DmgCyc,
    ) -> Self {
        debug_assert!(
            strain1p > 0.0 && strain2p > 0.0 && strain3p > 0.0 && strain4p > 0.0,
            "positive backbone strains must all be positive"
        );
        debug_assert!(
            strain1n < 0.0 && strain2n < 0.0 && strain3n < 0.0 && strain4n < 0.0,
            "negative backbone strains must all be negative"
        );

        // `Pinching4Material::SetEnvelope`.
        let k_pos = stress1p / strain1p;
        let k_neg = stress1n / strain1n;
        let k = if k_pos > k_neg { k_pos } else { k_neg };
        let u = if strain1p > -strain1n { 1e-4 * strain1p } else { -1e-4 * strain1n };

        let mut envlp_pos_strain = [u, strain1p, strain2p, strain3p, strain4p, 0.0];
        let mut envlp_pos_stress = [u * k, stress1p, stress2p, stress3p, stress4p, 0.0];
        let mut envlp_neg_strain = [-u, strain1n, strain2n, strain3n, strain4n, 0.0];
        let mut envlp_neg_stress = [-u * k, stress1n, stress2n, stress3n, stress4n, 0.0];

        let k1 = (stress4p - stress3p) / (strain4p - strain3p);
        let k2 = (stress4n - stress3n) / (strain4n - strain3n);

        envlp_pos_strain[5] = 1e6 * strain4p;
        envlp_pos_stress[5] = if k1 > 0.0 {
            stress4p + k1 * (envlp_pos_strain[5] - strain4p)
        } else {
            stress4p * 1.1
        };
        envlp_neg_strain[5] = 1e6 * strain4n;
        envlp_neg_stress[5] = if k2 > 0.0 {
            stress4n + k2 * (envlp_neg_strain[5] - strain4n)
        } else {
            stress4n * 1.1
        };

        let k_elastic_pos = envlp_pos_stress[1] / envlp_pos_strain[1];
        let k_elastic_neg = envlp_neg_stress[1] / envlp_neg_strain[1];

        let mut energy_pos = 0.5 * envlp_pos_strain[0] * envlp_pos_stress[0];
        for jt in 0..4 {
            energy_pos +=
                0.5 * (envlp_pos_stress[jt] + envlp_pos_stress[jt + 1]) * (envlp_pos_strain[jt + 1] - envlp_pos_strain[jt]);
        }
        let mut energy_neg = 0.5 * envlp_neg_strain[0] * envlp_neg_stress[0];
        for jy in 0..4 {
            energy_neg +=
                0.5 * (envlp_neg_stress[jy] + envlp_neg_stress[jy + 1]) * (envlp_neg_strain[jy + 1] - envlp_neg_strain[jy]);
        }
        let max_energy = if energy_pos > energy_neg { energy_pos } else { energy_neg };
        let energy_capacity = gamma_e * max_energy;

        // `Pinching4Material::revertToStart`.
        Material::Pinching4(Box::new(Pinching4Fields {
            envlp_pos_strain,
            envlp_pos_stress,
            envlp_neg_strain,
            envlp_neg_stress,
            k_elastic_pos,
            k_elastic_neg,
            energy_capacity,
            r_disp_p,
            r_force_p,
            u_force_p,
            r_disp_n,
            r_force_n,
            u_force_n,
            gamma_k_params,
            gamma_k_limit,
            gamma_d_params,
            gamma_d_limit,
            gamma_f_params,
            gamma_f_limit,
            dmg_cyc,
            state: Pinching4State::NearOrigin,
            strain: 0.0,
            stress: 0.0,
            strain_rate: 0.0,
            low_state_strain: envlp_neg_strain[0],
            low_state_stress: envlp_neg_stress[0],
            hgh_state_strain: envlp_pos_strain[0],
            hgh_state_stress: envlp_pos_stress[0],
            min_strain_dmnd: envlp_neg_strain[1],
            max_strain_dmnd: envlp_pos_strain[1],
            energy: 0.0,
            gamma_k: 0.0,
            gamma_d: 0.0,
            gamma_f: 0.0,
            gamma_k_used: 0.0,
            gamma_f_used: 0.0,
            n_cycle: 0.0,
        }))
    }
}

/// Ported from `Pinching4Material::setTrialStrain` (plus the parts of
/// `commitState` that only recompute derived quantities).
///
/// Structure follows the C++ exactly: seed the trial state from the
/// committed state, compute `dstrain` (with OpenSees' `1e-12` deadband),
/// run the transition function `getstate`, evaluate stress/tangent for
/// whichever state we ended up in, accumulate energy, then run `update_dmg`
/// unconditionally.
pub(super) fn evaluate(m: &Material, strain: f64) -> (f64, f64, Material) {
    let Material::Pinching4(fields) = m else {
        unreachable!()
    };
    let p = **fields;

    let mut dstrain = strain - p.strain;
    if dstrain < 1e-12 && dstrain > -1e-12 {
        dstrain = 0.0;
    }

    // `commitState` derives these three from committed fields only, and
    // nothing mutates them during a trial step — so they're recomputed here
    // rather than stored (see `Pinching4Fields`' doc comment).
    let u_max_damgd = p.max_strain_dmnd * (1.0 + p.gamma_d);
    let u_min_damgd = p.min_strain_dmnd * (1.0 + p.gamma_d);

    let mut t = Trial {
        state: p.state,
        low_strain: p.low_state_strain,
        low_stress: p.low_state_stress,
        hgh_strain: p.hgh_state_strain,
        hgh_stress: p.hgh_state_stress,
        min_dmnd: p.min_strain_dmnd,
        max_dmnd: p.max_strain_dmnd,
        energy: p.energy,
        gamma_k: p.gamma_k,
        gamma_d: p.gamma_d,
        gamma_f: p.gamma_f,
        gamma_k_used: p.gamma_k_used,
        gamma_f_used: p.gamma_f_used,
        k_elastic_pos_damgd: p.k_elastic_pos * (1.0 - p.gamma_k_used),
        k_elastic_neg_damgd: p.k_elastic_neg * (1.0 - p.gamma_k_used),
        pos_damgd_stress: scaled(&p.envlp_pos_stress, 1.0 - p.gamma_f_used),
        neg_damgd_stress: scaled(&p.envlp_neg_stress, 1.0 - p.gamma_f_used),
        n_cycle: p.n_cycle,
    };

    getstate(&mut t, &p, strain, dstrain, u_max_damgd, u_min_damgd);

    let (stress, tangent) = match t.state {
        Pinching4State::NearOrigin => {
            let tangent = p.envlp_pos_stress[0] / p.envlp_pos_strain[0];
            (tangent * strain, tangent)
        }
        Pinching4State::Positive => (t.pos_envlp_stress(&p, strain), t.pos_envlp_tangent(&p, strain)),
        Pinching4State::Negative => (t.neg_envlp_stress(&p, strain), t.neg_envlp_tangent(&p, strain)),
        Pinching4State::Trilinear3 => {
            let kunload = if t.hgh_strain < 0.0 {
                t.k_elastic_neg_damgd
            } else {
                t.k_elastic_pos_damgd
            };
            let mut s3_strain = [t.low_strain, 0.0, 0.0, t.hgh_strain];
            let mut s3_stress = [t.low_stress, 0.0, 0.0, t.hgh_stress];
            get_state3(&mut s3_strain, &mut s3_stress, kunload, &t, &p);
            (
                envlp34_stress(&s3_strain, &s3_stress, strain),
                envlp34_tangent(&s3_strain, &s3_stress, strain),
            )
        }
        Pinching4State::Trilinear4 => {
            let kunload = if t.low_strain < 0.0 {
                t.k_elastic_neg_damgd
            } else {
                t.k_elastic_pos_damgd
            };
            let mut s4_strain = [t.low_strain, 0.0, 0.0, t.hgh_strain];
            let mut s4_stress = [t.low_stress, 0.0, 0.0, t.hgh_stress];
            get_state4(&mut s4_strain, &mut s4_stress, kunload, &t, &p);
            (
                envlp34_stress(&s4_strain, &s4_stress, strain),
                envlp34_tangent(&s4_strain, &s4_stress, strain),
            )
        }
    };

    let denergy = 0.5 * (stress + p.stress) * dstrain;
    let elastic_strain_energy = if strain > 0.0 {
        0.5 * stress / t.k_elastic_pos_damgd * stress
    } else {
        0.5 * stress / t.k_elastic_neg_damgd * stress
    };
    t.energy = p.energy + denergy;

    update_dmg(&mut t, &p, strain, dstrain, elastic_strain_energy, p.n_cycle);

    // `commitState`: `CstrainRate` keeps its previous value across a
    // zero-increment step rather than latching zero.
    let strain_rate = if dstrain > 1e-12 || dstrain < -1e-12 { dstrain } else { p.strain_rate };

    (
        stress,
        tangent,
        Material::Pinching4(Box::new(Pinching4Fields {
            state: t.state,
            strain,
            stress,
            strain_rate,
            low_state_strain: t.low_strain,
            low_state_stress: t.low_stress,
            hgh_state_strain: t.hgh_strain,
            hgh_state_stress: t.hgh_stress,
            min_strain_dmnd: t.min_dmnd,
            max_strain_dmnd: t.max_dmnd,
            energy: t.energy,
            gamma_k: t.gamma_k,
            gamma_d: t.gamma_d,
            gamma_f: t.gamma_f,
            gamma_k_used: t.gamma_k_used,
            gamma_f_used: t.gamma_f_used,
            n_cycle: t.n_cycle,
            ..p
        })),
    )
}

fn scaled(v: &[f64; 6], factor: f64) -> [f64; 6] {
    let mut out = [0.0; 6];
    for i in 0..6 {
        out[i] = v[i] * factor;
    }
    out
}

/// Ported from `Pinching4Material::getstate` — the transition function.
///
/// A transition is only *considered* when the trial strain leaves the
/// current state's `[low, hgh]` bracket or the strain rate changes sign
/// (`cid`); which transition is taken then depends on the current state and
/// the sign of `du`. Note that `du == 0.0` always sets `cid`, which is
/// harmless: every inner branch that could fire requires a strictly signed
/// `du`, except state 0's, which tests the bracket itself.
fn getstate(t: &mut Trial, p: &Pinching4Fields, u: f64, du: f64, u_max_damgd: f64, u_min_damgd: f64) {
    // `CstrainRate`/`CgammaF`/`CgammaK`/`Cstrain`/`Cstress` are read from
    // the committed state; `t.*` are the trial values being built.
    let c_strain = p.strain;
    let c_stress = p.stress;

    let cid = du * p.strain_rate <= 0.0;
    let mut cis = false;
    let mut new_state = Pinching4State::NearOrigin;

    if u < t.low_strain || u > t.hgh_strain || cid {
        match t.state {
            Pinching4State::NearOrigin => {
                if u > t.hgh_strain {
                    cis = true;
                    new_state = Pinching4State::Positive;
                    t.low_strain = p.envlp_pos_strain[0];
                    t.low_stress = p.envlp_pos_stress[0];
                    t.hgh_strain = p.envlp_pos_strain[5];
                    t.hgh_stress = p.envlp_pos_stress[5];
                } else if u < t.low_strain {
                    cis = true;
                    new_state = Pinching4State::Negative;
                    t.low_strain = p.envlp_neg_strain[5];
                    t.low_stress = p.envlp_neg_stress[5];
                    t.hgh_strain = p.envlp_neg_strain[0];
                    t.hgh_stress = p.envlp_neg_stress[0];
                }
            }
            Pinching4State::Positive if du < 0.0 => {
                cis = true;
                if c_strain > t.max_dmnd {
                    t.max_dmnd = u - du;
                }
                if t.max_dmnd < u_max_damgd {
                    t.max_dmnd = u_max_damgd;
                }
                if u < u_min_damgd {
                    new_state = Pinching4State::Negative;
                    t.gamma_f_used = p.gamma_f;
                    t.neg_damgd_stress = scaled(&p.envlp_neg_stress, 1.0 - t.gamma_f_used);
                    t.low_strain = p.envlp_neg_strain[5];
                    t.low_stress = p.envlp_neg_stress[5];
                    t.hgh_strain = p.envlp_neg_strain[0];
                    t.hgh_stress = p.envlp_neg_stress[0];
                } else {
                    new_state = Pinching4State::Trilinear3;
                    t.low_strain = u_min_damgd;
                    t.gamma_f_used = p.gamma_f;
                    t.neg_damgd_stress = scaled(&p.envlp_neg_stress, 1.0 - t.gamma_f_used);
                    t.low_stress = t.neg_envlp_stress(p, u_min_damgd);
                    t.hgh_strain = c_strain;
                    t.hgh_stress = c_stress;
                }
                t.gamma_k_used = p.gamma_k;
                t.k_elastic_pos_damgd = p.k_elastic_pos * (1.0 - t.gamma_k_used);
            }
            Pinching4State::Negative if du > 0.0 => {
                cis = true;
                if c_strain < t.min_dmnd {
                    t.min_dmnd = c_strain;
                }
                if t.min_dmnd > u_min_damgd {
                    t.min_dmnd = u_min_damgd;
                }
                if u > u_max_damgd {
                    new_state = Pinching4State::Positive;
                    t.gamma_f_used = p.gamma_f;
                    t.pos_damgd_stress = scaled(&p.envlp_pos_stress, 1.0 - t.gamma_f_used);
                    t.low_strain = p.envlp_pos_strain[0];
                    t.low_stress = p.envlp_pos_stress[0];
                    t.hgh_strain = p.envlp_pos_strain[5];
                    t.hgh_stress = p.envlp_pos_stress[5];
                } else {
                    new_state = Pinching4State::Trilinear4;
                    t.low_strain = c_strain;
                    t.low_stress = c_stress;
                    t.hgh_strain = u_max_damgd;
                    t.gamma_f_used = p.gamma_f;
                    t.pos_damgd_stress = scaled(&p.envlp_pos_stress, 1.0 - t.gamma_f_used);
                    t.hgh_stress = t.pos_envlp_stress(p, u_max_damgd);
                }
                t.gamma_k_used = p.gamma_k;
                t.k_elastic_neg_damgd = p.k_elastic_neg * (1.0 - t.gamma_k_used);
            }
            Pinching4State::Trilinear3 => {
                if u < t.low_strain {
                    cis = true;
                    new_state = Pinching4State::Negative;
                    t.low_strain = p.envlp_neg_strain[5];
                    t.hgh_strain = p.envlp_neg_strain[0];
                    t.low_stress = t.neg_damgd_stress[5];
                    t.hgh_stress = t.neg_damgd_stress[0];
                } else if u > u_max_damgd && du > 0.0 {
                    cis = true;
                    new_state = Pinching4State::Positive;
                    t.low_strain = p.envlp_pos_strain[0];
                    t.low_stress = p.envlp_pos_stress[0];
                    t.hgh_strain = p.envlp_pos_strain[5];
                    t.hgh_stress = p.envlp_pos_stress[5];
                } else if du > 0.0 {
                    cis = true;
                    new_state = Pinching4State::Trilinear4;
                    t.low_strain = c_strain;
                    t.low_stress = c_stress;
                    t.hgh_strain = u_max_damgd;
                    t.gamma_f_used = p.gamma_f;
                    t.pos_damgd_stress = scaled(&p.envlp_pos_stress, 1.0 - t.gamma_f_used);
                    t.hgh_stress = t.pos_envlp_stress(p, u_max_damgd);
                    t.gamma_k_used = p.gamma_k;
                    t.k_elastic_neg_damgd = p.k_elastic_neg * (1.0 - t.gamma_k_used);
                }
            }
            Pinching4State::Trilinear4 => {
                if u > t.hgh_strain {
                    cis = true;
                    new_state = Pinching4State::Positive;
                    t.low_strain = p.envlp_pos_strain[0];
                    t.low_stress = t.pos_damgd_stress[0];
                    t.hgh_strain = p.envlp_pos_strain[5];
                    t.hgh_stress = t.pos_damgd_stress[5];
                } else if u < u_min_damgd && du < 0.0 {
                    cis = true;
                    new_state = Pinching4State::Negative;
                    t.low_strain = p.envlp_neg_strain[5];
                    t.low_stress = t.neg_damgd_stress[5];
                    t.hgh_strain = p.envlp_neg_strain[0];
                    t.hgh_stress = t.neg_damgd_stress[0];
                } else if du < 0.0 {
                    cis = true;
                    new_state = Pinching4State::Trilinear3;
                    t.low_strain = u_min_damgd;
                    t.gamma_f_used = p.gamma_f;
                    t.neg_damgd_stress = scaled(&p.envlp_neg_stress, 1.0 - t.gamma_f_used);
                    t.low_stress = t.neg_envlp_stress(p, u_min_damgd);
                    t.hgh_strain = c_strain;
                    t.hgh_stress = c_stress;
                    t.gamma_k_used = p.gamma_k;
                    t.k_elastic_pos_damgd = p.k_elastic_pos * (1.0 - t.gamma_k_used);
                }
            }
            _ => {}
        }
    }

    if cis {
        t.state = new_state;
    }
}

/// Ported from `Pinching4Material::getState3` — builds the 4-point
/// trilinear path for state 3. Points 0 and 3 are the state's bracket
/// (already filled in by the caller); points 1 and 2 are constructed from
/// the pinching parameters and then run through a cascade of
/// degenerate-geometry repairs.
///
/// Translated as a literal line-by-line trace of the C++, including the
/// several places where the "same" check is spelled differently between
/// `getState3` and `getState4` — those asymmetries are in the source and
/// are deliberately preserved rather than reconciled.
fn get_state3(strain: &mut [f64; 4], stress: &mut [f64; 4], kunload: f64, t: &Trial, p: &Pinching4Fields) {
    let kmax = if kunload > t.k_elastic_neg_damgd {
        kunload
    } else {
        t.k_elastic_neg_damgd
    };

    if strain[0] * strain[3] < 0.0 {
        // Trilinear unload/reload path expected; first define the point
        // for reloading.
        strain[1] = t.low_strain * p.r_disp_n;
        if p.r_force_n - p.u_force_n > 1e-8 {
            stress[1] = t.low_stress * p.r_force_n;
        } else if t.min_dmnd < p.envlp_neg_strain[3] {
            let st1 = t.low_stress * p.u_force_n * (1.0 + 1e-6);
            let st2 = t.neg_damgd_stress[4] * (1.0 + 1e-6);
            stress[1] = if st1 < st2 { st1 } else { st2 };
        } else {
            let st1 = t.neg_damgd_stress[3] * p.u_force_n * (1.0 + 1e-6);
            let st2 = t.neg_damgd_stress[4] * (1.0 + 1e-6);
            stress[1] = if st1 < st2 { st1 } else { st2 };
        }

        // If the reload stiffness exceeds the unload stiffness, reduce it
        // to match.
        if (stress[1] - stress[0]) / (strain[1] - strain[0]) > t.k_elastic_neg_damgd {
            strain[1] = t.low_strain + (stress[1] - stress[0]) / t.k_elastic_neg_damgd;
        }

        // Check that the reloading point is not behind point 4.
        if strain[1] > strain[3] {
            // Path taken to be a straight line between points 1 and 4.
            let du = strain[3] - strain[0];
            let df = stress[3] - stress[0];
            strain[1] = strain[0] + 0.33 * du;
            strain[2] = strain[0] + 0.67 * du;
            stress[1] = stress[0] + 0.33 * df;
            stress[2] = stress[0] + 0.67 * df;
        } else {
            if t.min_dmnd < p.envlp_neg_strain[3] {
                stress[2] = p.u_force_n * t.neg_damgd_stress[4];
            } else {
                stress[2] = p.u_force_n * t.neg_damgd_stress[3];
            }
            strain[2] = t.hgh_strain - (t.hgh_stress - stress[2]) / kunload;

            if strain[2] > strain[3] {
                // Point 3 should be along a line between 2 and 4.
                let du = strain[3] - strain[1];
                let df = stress[3] - stress[1];
                strain[2] = strain[1] + 0.5 * du;
                stress[2] = stress[1] + 0.5 * df;
            } else if (stress[2] - stress[1]) / (strain[2] - strain[1]) > kmax {
                // Linear unload-reload path expected.
                let du = strain[3] - strain[0];
                let df = stress[3] - stress[0];
                strain[1] = strain[0] + 0.33 * du;
                strain[2] = strain[0] + 0.67 * du;
                stress[1] = stress[0] + 0.33 * df;
                stress[2] = stress[0] + 0.67 * df;
            } else if (strain[2] < strain[1]) || ((stress[2] - stress[1]) / (strain[2] - strain[1]) < 0.0) {
                if strain[2] < 0.0 {
                    // Point 3 should be along a line between 2 and 4.
                    let du = strain[3] - strain[1];
                    let df = stress[3] - stress[1];
                    strain[2] = strain[1] + 0.5 * du;
                    stress[2] = stress[1] + 0.5 * df;
                } else if strain[1] > 0.0 {
                    // Point 2 should be along a line between 1 and 3.
                    let du = strain[2] - strain[0];
                    let df = stress[2] - stress[0];
                    strain[1] = strain[0] + 0.5 * du;
                    stress[1] = stress[0] + 0.5 * df;
                } else {
                    let avgforce = 0.5 * (stress[2] + stress[1]);
                    let dfr = if avgforce < 0.0 { -avgforce / 100.0 } else { avgforce / 100.0 };
                    let slope12 = (stress[1] - stress[0]) / (strain[1] - strain[0]);
                    let slope34 = (stress[3] - stress[2]) / (strain[3] - strain[2]);
                    stress[1] = avgforce - dfr;
                    stress[2] = avgforce + dfr;
                    strain[1] = strain[0] + (stress[1] - stress[0]) / slope12;
                    strain[2] = strain[3] - (stress[3] - stress[2]) / slope34;
                }
            }
        }
    } else {
        // Linear unload-reload path is expected.
        let du = strain[3] - strain[0];
        let df = stress[3] - stress[0];
        strain[1] = strain[0] + 0.33 * du;
        strain[2] = strain[0] + 0.67 * du;
        stress[1] = stress[0] + 0.33 * df;
        stress[2] = stress[0] + 0.67 * df;
    }

    final_check(strain, stress);
}

/// Ported from `Pinching4Material::getState4` — the mirror image of
/// [`get_state3`]; see its doc comment. Note the extra `uForceP == 0.0`
/// special case at the top, which `getState3` has no counterpart for: that
/// asymmetry is in the OpenSees source.
///
/// `clippy::if_same_then_else` fires on that very special case (`u_force_p
/// == 0.0` and `r_force_p - u_force_p > 1e-8` do the same thing), and is
/// allowed rather than collapsed for the same reason `hysteretic.rs`
/// allows it: the two conditions are genuinely different in the source,
/// and merging them would obscure the line-by-line correspondence this
/// whole cascade depends on for reviewability.
#[allow(clippy::if_same_then_else)]
fn get_state4(strain: &mut [f64; 4], stress: &mut [f64; 4], kunload: f64, t: &Trial, p: &Pinching4Fields) {
    let kmax = if kunload > t.k_elastic_pos_damgd {
        kunload
    } else {
        t.k_elastic_pos_damgd
    };

    if strain[0] * strain[3] < 0.0 {
        // Trilinear unload/reload path expected.
        strain[2] = t.hgh_strain * p.r_disp_p;
        if p.u_force_p == 0.0 {
            stress[2] = t.hgh_stress * p.r_force_p;
        } else if p.r_force_p - p.u_force_p > 1e-8 {
            stress[2] = t.hgh_stress * p.r_force_p;
        } else if t.max_dmnd > p.envlp_pos_strain[3] {
            let st1 = t.hgh_stress * p.u_force_p * (1.0 + 1e-6);
            let st2 = t.pos_damgd_stress[4] * (1.0 + 1e-6);
            stress[2] = if st1 > st2 { st1 } else { st2 };
        } else {
            let st1 = t.pos_damgd_stress[3] * p.u_force_p * (1.0 + 1e-6);
            let st2 = t.pos_damgd_stress[4] * (1.0 + 1e-6);
            stress[2] = if st1 > st2 { st1 } else { st2 };
        }

        // If the reload stiffness exceeds the unload stiffness, reduce it
        // to match.
        if (stress[3] - stress[2]) / (strain[3] - strain[2]) > t.k_elastic_pos_damgd {
            strain[2] = t.hgh_strain - (stress[3] - stress[2]) / t.k_elastic_pos_damgd;
        }

        // Check that the reloading point is not behind point 1.
        if strain[2] < strain[0] {
            // Path taken to be a straight line between points 1 and 4.
            let du = strain[3] - strain[0];
            let df = stress[3] - stress[0];
            strain[1] = strain[0] + 0.33 * du;
            strain[2] = strain[0] + 0.67 * du;
            stress[1] = stress[0] + 0.33 * df;
            stress[2] = stress[0] + 0.67 * df;
        } else {
            if t.max_dmnd > p.envlp_pos_strain[3] {
                stress[1] = p.u_force_p * t.pos_damgd_stress[4];
            } else {
                stress[1] = p.u_force_p * t.pos_damgd_stress[3];
            }
            strain[1] = t.low_strain + (-t.low_stress + stress[1]) / kunload;

            if strain[1] < strain[0] {
                // Point 2 should be along a line between 1 and 3.
                let du = strain[2] - strain[0];
                let df = stress[2] - stress[0];
                strain[1] = strain[0] + 0.5 * du;
                stress[1] = stress[0] + 0.5 * df;
            } else if (stress[2] - stress[1]) / (strain[2] - strain[1]) > kmax {
                // Linear unload-reload path expected.
                let du = strain[3] - strain[0];
                let df = stress[3] - stress[0];
                strain[1] = strain[0] + 0.33 * du;
                strain[2] = strain[0] + 0.67 * du;
                stress[1] = stress[0] + 0.33 * df;
                stress[2] = stress[0] + 0.67 * df;
            } else if (strain[2] < strain[1]) || ((stress[2] - stress[1]) / (strain[2] - strain[1]) < 0.0) {
                if strain[1] > 0.0 {
                    // Point 2 should be along a line between 1 and 3.
                    let du = strain[2] - strain[0];
                    let df = stress[2] - stress[0];
                    strain[1] = strain[0] + 0.5 * du;
                    stress[1] = stress[0] + 0.5 * df;
                } else if strain[2] < 0.0 {
                    // Point 2 should be along a line between 2 and 4.
                    let du = strain[3] - strain[1];
                    let df = stress[3] - stress[1];
                    strain[2] = strain[1] + 0.5 * du;
                    stress[2] = stress[1] + 0.5 * df;
                } else {
                    let avgforce = 0.5 * (stress[2] + stress[1]);
                    let dfr = if avgforce < 0.0 { -avgforce / 100.0 } else { avgforce / 100.0 };
                    let slope12 = (stress[1] - stress[0]) / (strain[1] - strain[0]);
                    let slope34 = (stress[3] - stress[2]) / (strain[3] - strain[2]);
                    stress[1] = avgforce - dfr;
                    stress[2] = avgforce + dfr;
                    strain[1] = strain[0] + (stress[1] - stress[0]) / slope12;
                    strain[2] = strain[3] - (stress[3] - stress[2]) / slope34;
                }
            }
        }
    } else {
        // Linear unload-reload path is expected.
        let du = strain[3] - strain[0];
        let df = stress[3] - stress[0];
        strain[1] = strain[0] + 0.33 * du;
        strain[2] = strain[0] + 0.67 * du;
        stress[1] = stress[0] + 0.33 * df;
        stress[2] = stress[0] + 0.67 * df;
    }

    final_check(strain, stress);
}

/// The identical "final check" tail shared verbatim by `getState3` and
/// `getState4`: if any segment of the path runs backwards in strain or
/// stress, collapse the whole path to the straight 33%/67% split; and —
/// only reachable once that collapse has happened, since `slope` is
/// otherwise left at its initial `0.0` — if the collapsed slope is
/// positive but flatter than the secant to point 1, pin point 2 to the
/// origin and put point 3 at the midpoint.
///
/// The `slope == 0.0` gating is an artifact of the C++ control flow rather
/// than an obviously intended design, so it's preserved exactly: both the
/// loop-local shadowing of `du`/`df` and the fact that the second `if` is
/// evaluated on every iteration are load-bearing for which branch fires.
fn final_check(strain: &mut [f64; 4], stress: &mut [f64; 4]) {
    let check_slope = stress[0] / strain[0];
    let mut slope = 0.0;

    let mut i = 0;
    while i < 3 {
        let du = strain[i + 1] - strain[i];
        let df = stress[i + 1] - stress[i];
        if du < 0.0 || df < 0.0 {
            let du = strain[3] - strain[0];
            let df = stress[3] - stress[0];
            strain[1] = strain[0] + 0.33 * du;
            strain[2] = strain[0] + 0.67 * du;
            stress[1] = stress[0] + 0.33 * df;
            stress[2] = stress[0] + 0.67 * df;
            slope = df / du;
            i = 3;
        }
        if slope > 1e-8 && slope < check_slope {
            strain[1] = 0.0;
            stress[1] = 0.0;
            strain[2] = strain[3] / 2.0;
            stress[2] = stress[3] / 2.0;
        }
        i += 1;
    }
}

/// `Pinching4Material::Envlp3Stress`/`Envlp4Stress` — identical bodies in
/// the source, so shared here. Piecewise-interpolates the 4-point state-3/4
/// path, extrapolating off the first or last segment when the strain falls
/// outside it.
fn envlp34_stress(strain: &[f64; 4], stress: &[f64; 4], u: f64) -> f64 {
    let mut k = 0.0;
    let mut f = 0.0;
    for i in 0..3 {
        if u >= strain[i] {
            k = (stress[i + 1] - stress[i]) / (strain[i + 1] - strain[i]);
            f = stress[i] + (u - strain[i]) * k;
        }
    }
    if k == 0.0 {
        let i = if u < strain[0] { 0 } else { 2 };
        k = (stress[i + 1] - stress[i]) / (strain[i + 1] - strain[i]);
        f = stress[i] + (u - strain[i]) * k;
    }
    f
}

/// `Pinching4Material::Envlp3Tangent`/`Envlp4Tangent`; see
/// [`envlp34_stress`].
fn envlp34_tangent(strain: &[f64; 4], stress: &[f64; 4], u: f64) -> f64 {
    let mut k = 0.0;
    for i in 0..3 {
        if u >= strain[i] {
            k = (stress[i + 1] - stress[i]) / (strain[i + 1] - strain[i]);
        }
    }
    if k == 0.0 {
        let i = if u < strain[0] { 0 } else { 2 };
        k = (stress[i + 1] - stress[i]) / (strain[i + 1] - strain[i]);
    }
    k
}

/// Ported from `Pinching4Material::updateDmg`. Runs unconditionally at the
/// end of every trial evaluation — *not* gated behind a state transition,
/// unlike `gamma_k_used`/`gamma_f_used`.
///
/// This produces the first of three damage "clocks": the values written
/// here become the *committed* `gamma_k`/`gamma_d`/`gamma_f` of the next
/// step, which in turn become `gamma_k_used`/`gamma_f_used` only if some
/// later step happens to transition state. Each damage factor combines a
/// ductility term (`umax_abs/uult_abs`) with either a dissipated-energy
/// ratio or a half-cycle count, then is clamped by its own limit and — for
/// the stiffness factor — by `gamma_k_lim_env`, an envelope-derived bound
/// that prevents damage from ever *increasing* the tangent.
fn update_dmg(
    t: &mut Trial,
    p: &Pinching4Fields,
    strain: f64,
    dstrain: f64,
    elastic_strain_energy: f64,
    c_n_cycle: f64,
) {
    let umax_abs = if t.max_dmnd > -t.min_dmnd { t.max_dmnd } else { -t.min_dmnd };
    let uult_abs = if p.envlp_pos_strain[4] > -p.envlp_neg_strain[4] {
        p.envlp_pos_strain[4]
    } else {
        -p.envlp_neg_strain[4]
    };
    t.n_cycle = c_n_cycle + dstrain.abs() / (4.0 * umax_abs);

    let within = strain < uult_abs && strain > -uult_abs;

    if within && t.energy < p.energy_capacity {
        t.gamma_k = p.gamma_k_params[0] * (umax_abs / uult_abs).powf(p.gamma_k_params[2]);
        t.gamma_d = p.gamma_d_params[0] * (umax_abs / uult_abs).powf(p.gamma_d_params[2]);
        t.gamma_f = p.gamma_f_params[0] * (umax_abs / uult_abs).powf(p.gamma_f_params[2]);

        if t.energy > elastic_strain_energy && p.dmg_cyc == Pinching4DmgCyc::EnergyBased {
            let tes = (t.energy - elastic_strain_energy) / p.energy_capacity;
            t.gamma_k += p.gamma_k_params[1] * tes.powf(p.gamma_k_params[3]);
            t.gamma_d += p.gamma_d_params[1] * tes.powf(p.gamma_d_params[3]);
            t.gamma_f += p.gamma_f_params[1] * tes.powf(p.gamma_f_params[3]);
        } else if p.dmg_cyc == Pinching4DmgCyc::CycleBased {
            t.gamma_k += p.gamma_k_params[1] * t.n_cycle.powf(p.gamma_k_params[3]);
            t.gamma_d += p.gamma_d_params[1] * t.n_cycle.powf(p.gamma_d_params[3]);
            t.gamma_f += p.gamma_f_params[1] * t.n_cycle.powf(p.gamma_f_params[3]);
        }

        let gamma_k_lim_env = gamma_k_lim_env(t, p);
        let k1 = if t.gamma_k < p.gamma_k_limit { t.gamma_k } else { p.gamma_k_limit };
        t.gamma_k = if k1 < gamma_k_lim_env { k1 } else { gamma_k_lim_env };
        if t.gamma_d >= p.gamma_d_limit {
            t.gamma_d = p.gamma_d_limit;
        }
        if t.gamma_f >= p.gamma_f_limit {
            t.gamma_f = p.gamma_f_limit;
        }
    } else if within {
        let gamma_k_lim_env = gamma_k_lim_env(t, p);
        t.gamma_k = if p.gamma_k_limit < gamma_k_lim_env {
            p.gamma_k_limit
        } else {
            gamma_k_lim_env
        };
        t.gamma_d = p.gamma_d_limit;
        t.gamma_f = p.gamma_f_limit;
    }
}

/// `updateDmg`'s `gammaKLimEnv`: how much the current peak-ductility-demand
/// point has already softened relative to the undamaged elastic stiffness,
/// which bounds how much *further* stiffness damage may be applied.
///
/// The two copies in the C++ differ only in `>` vs. `>=` when picking the
/// governing side — irrelevant (they pick numerically equal values on a
/// tie), so they're shared here.
fn gamma_k_lim_env(t: &Trial, p: &Pinching4Fields) -> f64 {
    let kmin_p = t.pos_envlp_stress(p, t.max_dmnd) / t.max_dmnd;
    let kmin_n = t.neg_envlp_stress(p, t.min_dmnd) / t.min_dmnd;
    let kmin = if (kmin_p / p.k_elastic_pos) > (kmin_n / p.k_elastic_neg) {
        kmin_p / p.k_elastic_pos
    } else {
        kmin_n / p.k_elastic_neg
    };
    if 0.0 > (1.0 - kmin) { 0.0 } else { 1.0 - kmin }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A symmetric 4-point-per-side backbone shared by the tests below:
    /// positive points `(0.01, 10)`, `(0.02, 15)`, `(0.03, 17)`,
    /// `(0.04, 10)` — an elastic slope of 1000, then 500, then 200, then a
    /// softening -700 — with the negative side its exact mirror.
    ///
    /// All twelve damage parameters and all three limits are zero, so
    /// `gamma_k`/`gamma_d`/`gamma_f` stay identically zero for the whole
    /// test: the damaged envelope equals the undamaged one and
    /// `u_max_damgd`/`u_min_damgd` equal the raw strain demands. That's the
    /// same tactic `Hysteretic`'s `beta = 0` test material used — it keeps
    /// the state machine and the geometric cascade (the actual subjects of
    /// these tests) hand-computable.
    fn pinching4_test_material(r_disp: f64, r_force: f64, u_force: f64) -> Material {
        Material::pinching4(
            10.0, 0.01, 15.0, 0.02, 17.0, 0.03, 10.0, 0.04, -10.0, -0.01, -15.0, -0.02, -17.0, -0.03, -10.0, -0.04,
            r_disp, r_force, u_force, r_disp, r_force, u_force, [0.0; 4], 0.0, [0.0; 4], 0.0, [0.0; 4], 0.0, 10.0,
            Pinching4DmgCyc::EnergyBased,
        )
    }

    /// Test 1: state 0 (`NearOrigin`) gives the synthesized near-origin
    /// point's slope exactly. `SetEnvelope` puts that point at
    /// `u = 1e-4*strain1p` with stress `u*k`, `k = max(stress1p/strain1p,
    /// stress1n/strain1n) = 1000`, so the state-0 tangent is exactly 1000 —
    /// the same slope as the first real backbone segment, but arrived at
    /// through an entirely different code path.
    ///
    /// "Exactly" up to one ulp: the tangent is recovered as
    /// `(u*k)/u` rather than as `k` itself, which for `u = 1e-6` lands on
    /// `1000.0000000000001`. Reconstructing `k` directly would be a
    /// deviation from the source, so the test carries the tolerance
    /// instead.
    #[test]
    fn pinching4_gives_exact_initial_tangent_at_zero_strain() {
        let m = pinching4_test_material(0.2, 0.4, 0.05);
        let (stress, tangent) = m.trial_stress_tangent(0.0);
        assert_eq!(stress, 0.0);
        assert!((tangent - 1000.0).abs() < 1e-9, "got {tangent}");
        assert_eq!(m.initial_tangent(), tangent);
    }

    /// Test 2: state 1 rides the 6-point positive envelope. Checks every
    /// real breakpoint (1-4), a midpoint of every real segment, and the
    /// extrapolated point-5 segment beyond `strain4p`.
    ///
    /// Note the point-0/point-1 segment: it runs from `(1e-6, 1e-3)` to
    /// `(0.01, 10)`, i.e. slope exactly 1000 and therefore
    /// indistinguishable from a straight elastic line through the origin —
    /// which is the whole point of the synthesized point.
    #[test]
    fn pinching4_matches_the_six_point_envelope_on_first_loading() {
        // One commit gets us out of state 0 and onto the positive
        // envelope; every trial below then stays in state 1 (all are
        // inside the state's `[1e-6, 4e4]` bracket, and all are increasing
        // so the strain rate never flips sign).
        let m = pinching4_test_material(0.2, 0.4, 0.05).commit(0.005);

        let close = |got: (f64, f64), want: (f64, f64)| {
            assert!((got.0 - want.0).abs() < 1e-9, "stress: got {}, want {}", got.0, want.0);
            assert!((got.1 - want.1).abs() < 1e-9, "tangent: got {}, want {}", got.1, want.1);
        };

        close(m.trial_stress_tangent(0.005), (5.0, 1000.0));
        close(m.trial_stress_tangent(0.01), (10.0, 1000.0)); // breakpoint 1
        close(m.trial_stress_tangent(0.015), (12.5, 500.0));
        close(m.trial_stress_tangent(0.02), (15.0, 500.0)); // breakpoint 2
        close(m.trial_stress_tangent(0.025), (16.0, 200.0));
        close(m.trial_stress_tangent(0.03), (17.0, 200.0)); // breakpoint 3
        close(m.trial_stress_tangent(0.035), (13.5, -700.0)); // softening
        close(m.trial_stress_tangent(0.04), (10.0, -700.0)); // breakpoint 4

        // Past breakpoint 4, the extrapolated point 5 takes over: since the
        // 3-4 segment is softening (k1 < 0), `SetEnvelope` puts point 5 at
        // `(1e6*strain4p, 1.1*stress4p) = (4e4, 11)` — a nearly flat
        // segment, not a continuation of the -700 softening slope.
        let k5 = (11.0 - 10.0) / (40000.0 - 0.04);
        close(m.trial_stress_tangent(0.05), (10.0 + 0.01 * k5, k5));
    }

    /// Test 3: a full reversal off the positive envelope into state 3,
    /// checked against a hand-walked `get_state3` trace that takes the
    /// *clean* trilinear path — no degenerate-geometry fallback fires.
    ///
    /// `r_disp <= r_force` (0.2 vs. 0.4) is what keeps the reload-stiffness
    /// cap from engaging: the secant from point 1 to point 2 comes out at
    /// `1000*(1-r_force)/(1-r_disp) = 750`, under the elastic 1000.
    ///
    /// Hand trace, after committing 0.005 then 0.02 (so `Cstrain = 0.02`,
    /// `Cstress = 15`, `max_strain_dmnd` still at its initial `0.01`):
    /// reversing to -0.005 sets `max_dmnd = 0.02`, and since -0.005 is not
    /// past `u_min_damgd = -0.01`, the new state is 3 with bracket
    /// `low = (-0.01, negEnvlpStress(-0.01) = -10)` and
    /// `hgh = (0.02, 15)`. Then:
    ///   point 1 = `(-0.01*0.2, -10*0.4) = (-0.002, -4)`;
    ///   point 2 = stress `u_force*envlp_neg[3] = 0.05*-17 = -0.85`,
    ///     strain `0.02 - (15+0.85)/1000 = 0.00415`.
    /// Both survive every check, so the path is
    /// `(-0.01,-10) -> (-0.002,-4) -> (0.00415,-0.85) -> (0.02,15)`.
    #[test]
    fn pinching4_reloads_along_a_clean_trilinear_path_after_a_reversal() {
        let m = pinching4_test_material(0.2, 0.4, 0.05).commit(0.005).commit(0.02);

        // On the first segment of the path: slope (-4+10)/(-0.002+0.01) =
        // 750, so stress = -10 + 0.005*750 = -6.25.
        let (stress, tangent) = m.trial_stress_tangent(-0.005);
        assert!((stress - -6.25).abs() < 1e-9, "got {stress}");
        assert!((tangent - 750.0).abs() < 1e-9, "got {tangent}");

        // On the second (pinched) segment, at zero strain: slope
        // (-0.85+4)/(0.00415+0.002) = 3.15/0.00615, so stress =
        // -4 + 0.002*that.
        let k = 3.15 / 0.00615;
        let (stress, tangent) = m.trial_stress_tangent(0.0);
        assert!((stress - (-4.0 + 0.002 * k)).abs() < 1e-9, "got {stress}");
        assert!((tangent - k).abs() < 1e-9, "got {tangent}");
        assert!(
            tangent < 750.0,
            "the pinched segment should be softer than the reload segment, got {tangent}"
        );
    }

    /// Test 4: the same reversal, but with pinching parameters chosen to
    /// drive one of `get_state3`'s degenerate-geometry fallbacks — the
    /// "linear unload-reload path expected" branch, taken when the
    /// constructed segment from point 1 to point 2 is stiffer than `kmax`
    /// (the greater of the unload stiffness and the damaged negative
    /// elastic stiffness, both 1000 here).
    ///
    /// `r_disp = 0.05` pulls point 1 in close to the origin while
    /// `r_force = 0.9` keeps its stress near the full envelope value, so
    /// the 1-to-2 secant comes out at `(−0.85+9)/(0.00415+0.0005) ≈ 1753`,
    /// over `kmax`. The whole path is then discarded and redrawn as the
    /// straight 33%/67% split between points 0 and 3:
    /// `(-0.01,-10) -> (-0.0001,-1.75) -> (0.0101, 6.75) -> (0.02, 15)`.
    #[test]
    fn pinching4_falls_back_to_a_linear_path_when_the_reload_segment_is_too_stiff() {
        let m = pinching4_test_material(0.05, 0.9, 0.05).commit(0.005).commit(0.02);

        // First segment of the redrawn path: from (-0.01,-10) to
        // (-0.0001,-1.75), slope 8.25/0.0099 = 2500/3.
        let k = 8.25 / 0.0099;
        let (stress, tangent) = m.trial_stress_tangent(-0.005);
        assert!((stress - (-10.0 + 0.005 * k)).abs() < 1e-9, "got {stress}");
        assert!((tangent - k).abs() < 1e-9, "got {tangent}");
        assert!((stress - -5.833333333333333).abs() < 1e-9, "got {stress}");

        // The fallback is genuinely a different path, not a relabelling of
        // the trilinear one: the un-fallen-back construction would have put
        // point 1 at (-0.0005, -9), a first-segment slope of 1/0.0095 ≈ 105
        // rather than ~833.
        assert!(
            (tangent - 1.0 / 0.0095).abs() > 1.0,
            "expected the fallback path, not the raw trilinear construction"
        );
    }

    /// The `Trilinear3 -> Trilinear4` transition, i.e. reversing *again*
    /// mid-reload. Exercises the one `getstate` branch that both refreshes
    /// the positive damaged envelope and rebuilds the bracket from the
    /// committed point, and confirms the resulting state-4 path is anchored
    /// at that reversal point.
    #[test]
    fn pinching4_reverses_from_state_three_into_state_four() {
        let m = pinching4_test_material(0.2, 0.4, 0.05)
            .commit(0.005)
            .commit(0.02)
            .commit(-0.005)
            .commit(0.0);
        let Material::Pinching4(f) = &m else { panic!() };
        assert_eq!(f.state, Pinching4State::Trilinear4);
        // The state-4 bracket runs from the reversal point (the previous
        // committed state, -0.005) up to `u_max_damgd` (0.02, the peak
        // strain demand, undamaged since gamma_d is zero).
        assert!((f.low_state_strain - -0.005).abs() < 1e-12);
        assert!((f.hgh_state_strain - 0.02).abs() < 1e-12);
        // ...anchored at the stress the state-3 path had reached there
        // (-6.25, the value test 3 hand-derived).
        assert!((f.low_state_stress - -6.25).abs() < 1e-9, "got {}", f.low_state_stress);
        // ...and at the *damaged* positive envelope stress at that peak
        // demand, 15 (undamaged here). Evaluating at the top of the
        // bracket walks the whole state-4 path and lands exactly there.
        assert!((f.hgh_state_stress - 15.0).abs() < 1e-9, "got {}", f.hgh_state_stress);
        let (stress, _) = m.trial_stress_tangent(0.02);
        assert!((stress - 15.0).abs() < 1e-9, "got {stress}");
    }
}
