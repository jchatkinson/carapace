use super::Material;

/// `HystereticMaterial`'s `CloadIndicator`/`TloadIndicator` — which
/// direction the material last loaded in (`None` only before the first
/// nonzero strain increment ever).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HystereticLoad {
    None,
    Positive,
    Negative,
}

/// The full field set behind `Material::Hysteretic` (see that variant's
/// doc comment for why it's boxed rather than inlined). `Copy` despite the
/// field count — every field is itself `f64`/`HystereticLoad`, so this is
/// cheap to copy out of its `Box` for the duration of one `evaluate` call,
/// which is simpler than borrow-juggling a destructure through `&Box<_>`.
///
/// The backbone (`mom1p..mom3n`/`rot1p..rot3n`) and pinching/damage
/// parameters (`pinch_x`, `pinch_y`, `damfc1`, `damfc2`, `beta`) are fixed
/// at construction, along with the envelope slopes/limit (`e1p..e3n`,
/// `eup`, `eun`, `energy_a`) `HystereticMaterial::setEnvelope` derives from
/// them once — recomputing those every `evaluate` call would be wasted
/// work, since they never change. History: `rot_max`/`rot_min` are the
/// peak strain reached on each side (`CrotMax`/`CrotMin`); `rot_pu`/
/// `rot_nu` are the pinching reference points captured at the last
/// reversal into negative/positive loading respectively (`CrotPu`/
/// `CrotNu`); `energy_d` is cumulative dissipated energy, used only to
/// compute the damage factor on a reversal; `load` mirrors
/// `CloadIndicator`.
#[derive(Debug, Clone, Copy)]
pub struct HystereticFields {
    pub mom1p: f64,
    pub rot1p: f64,
    pub mom2p: f64,
    pub rot2p: f64,
    pub mom3p: f64,
    pub rot3p: f64,
    pub mom1n: f64,
    pub rot1n: f64,
    pub mom2n: f64,
    pub rot2n: f64,
    pub mom3n: f64,
    pub rot3n: f64,
    pub pinch_x: f64,
    pub pinch_y: f64,
    pub damfc1: f64,
    pub damfc2: f64,
    pub beta: f64,
    pub e1p: f64,
    pub e2p: f64,
    pub e3p: f64,
    pub e1n: f64,
    pub e2n: f64,
    pub e3n: f64,
    pub eup: f64,
    pub eun: f64,
    pub energy_a: f64,
    pub rot_max: f64,
    pub rot_min: f64,
    pub rot_pu: f64,
    pub rot_nu: f64,
    pub energy_d: f64,
    pub load: HystereticLoad,
    pub strain: f64,
    pub stress: f64,
}

/// `HystereticMaterial`'s fixed (non-history) backbone parameters, bundled
/// so `hysteretic_positive_increment`/`hysteretic_negative_increment` don't
/// need a 20-parameter signature each. Mirrors `posEnvlpStress`/
/// `negEnvlpStress`/`posEnvlpTangent`/`negEnvlpTangent`/`posEnvlpRotlim`/
/// `negEnvlpRotlim` from `HystereticMaterial.cpp` exactly.
struct HystereticEnvelope {
    mom1p: f64,
    rot1p: f64,
    mom2p: f64,
    rot2p: f64,
    mom3p: f64,
    rot3p: f64,
    e1p: f64,
    e2p: f64,
    e3p: f64,
    mom1n: f64,
    rot1n: f64,
    mom2n: f64,
    rot2n: f64,
    mom3n: f64,
    rot3n: f64,
    e1n: f64,
    e2n: f64,
    e3n: f64,
}

/// `HystereticMaterial.cpp`'s `POS_INF_STRAIN`/`NEG_INF_STRAIN` sentinels
/// (from `UniaxialMaterial.h`) — a strain limit this large means "no
/// limit," never an actual physical bound.
const HYSTERETIC_POS_INF_STRAIN: f64 = 1.0e16;
const HYSTERETIC_NEG_INF_STRAIN: f64 = -1.0e16;

impl HystereticEnvelope {
    fn pos_stress(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            0.0
        } else if strain <= self.rot1p {
            self.e1p * strain
        } else if strain <= self.rot2p {
            self.mom1p + self.e2p * (strain - self.rot1p)
        } else if strain <= self.rot3p || self.e3p > 0.0 {
            self.mom2p + self.e3p * (strain - self.rot2p)
        } else {
            self.mom3p
        }
    }

    fn neg_stress(&self, strain: f64) -> f64 {
        if strain >= 0.0 {
            0.0
        } else if strain >= self.rot1n {
            self.e1n * strain
        } else if strain >= self.rot2n {
            self.mom1n + self.e2n * (strain - self.rot1n)
        } else if strain >= self.rot3n || self.e3n > 0.0 {
            self.mom2n + self.e3n * (strain - self.rot2n)
        } else {
            self.mom3n
        }
    }

    fn pos_tangent(&self, strain: f64) -> f64 {
        if strain < 0.0 {
            self.e1p * 1.0e-9
        } else if strain <= self.rot1p {
            self.e1p
        } else if strain <= self.rot2p {
            self.e2p
        } else if strain <= self.rot3p || self.e3p > 0.0 {
            self.e3p
        } else {
            self.e1p * 1.0e-9
        }
    }

    fn neg_tangent(&self, strain: f64) -> f64 {
        if strain > 0.0 {
            self.e1n * 1.0e-9
        } else if strain >= self.rot1n {
            self.e1n
        } else if strain >= self.rot2n {
            self.e2n
        } else if strain >= self.rot3n || self.e3n > 0.0 {
            self.e3n
        } else {
            self.e1n * 1.0e-9
        }
    }

    /// The strain at which the positive backbone, extended past `strain`,
    /// would cross zero stress on a softening segment — `POS_INF_STRAIN`
    /// if the backbone never softens (or hasn't yet, from `strain`).
    ///
    /// The two branches below returning the same `HYSTERETIC_POS_INF_STRAIN`
    /// literal are checking genuinely different conditions from
    /// `HystereticMaterial.cpp`'s `posEnvlpRotlim` (no strain limit was
    /// ever found vs. the one found doesn't correspond to a real softening
    /// crossing) — collapsing them would lose that distinction.
    #[allow(clippy::if_same_then_else)]
    fn pos_rotlim(&self, strain: f64) -> f64 {
        let mut strain_limit = HYSTERETIC_POS_INF_STRAIN;

        if strain <= self.rot1p {
            return HYSTERETIC_POS_INF_STRAIN;
        }
        if strain > self.rot1p && strain <= self.rot2p && self.e2p < 0.0 {
            strain_limit = self.rot1p - self.mom1p / self.e2p;
        }
        if strain > self.rot2p && self.e3p < 0.0 {
            strain_limit = self.rot2p - self.mom2p / self.e3p;
        }

        if strain_limit == HYSTERETIC_POS_INF_STRAIN {
            HYSTERETIC_POS_INF_STRAIN
        } else if self.pos_stress(strain_limit) > 0.0 {
            HYSTERETIC_POS_INF_STRAIN
        } else {
            strain_limit
        }
    }

    /// See `pos_rotlim` above for why the same-looking branches aren't
    /// collapsed.
    #[allow(clippy::if_same_then_else)]
    fn neg_rotlim(&self, strain: f64) -> f64 {
        let mut strain_limit = HYSTERETIC_NEG_INF_STRAIN;

        if strain >= self.rot1n {
            return HYSTERETIC_NEG_INF_STRAIN;
        }
        if strain < self.rot1n && strain >= self.rot2n && self.e2n < 0.0 {
            strain_limit = self.rot1n - self.mom1n / self.e2n;
        }
        if strain < self.rot2n && self.e3n < 0.0 {
            strain_limit = self.rot2n - self.mom2n / self.e3n;
        }

        if strain_limit == HYSTERETIC_NEG_INF_STRAIN {
            HYSTERETIC_NEG_INF_STRAIN
        } else if self.neg_stress(strain_limit) < 0.0 {
            HYSTERETIC_NEG_INF_STRAIN
        } else {
            strain_limit
        }
    }
}

impl Material {
    /// The 3-point-per-side backbone constructor (OpenSees also offers a
    /// 2-point constructor that derives the midpoint as the average of the
    /// two given points — not ported, since a caller can just pass that
    /// average as `mom2p`/`rot2p`/`mom2n`/`rot2n` directly).
    ///
    /// Debug builds assert the backbone is monotonic and one-to-one
    /// (`0 < rot1p < rot2p < rot3p`, `rot1n < rot2n < rot3n < 0`) —
    /// OpenSees hard-exits the process on this same precondition failure,
    /// which isn't a library-appropriate response here.
    #[allow(clippy::too_many_arguments)]
    pub fn hysteretic(
        mom1p: f64,
        rot1p: f64,
        mom2p: f64,
        rot2p: f64,
        mom3p: f64,
        rot3p: f64,
        mom1n: f64,
        rot1n: f64,
        mom2n: f64,
        rot2n: f64,
        mom3n: f64,
        rot3n: f64,
        pinch_x: f64,
        pinch_y: f64,
        damfc1: f64,
        damfc2: f64,
        beta: f64,
    ) -> Self {
        debug_assert!(rot1p > 0.0 && rot2p > rot1p && rot3p > rot2p, "positive backbone must be increasing");
        debug_assert!(rot1n < 0.0 && rot2n < rot1n && rot3n < rot2n, "negative backbone must be decreasing");

        let energy_a = 0.5
            * (rot1p * mom1p
                + (rot2p - rot1p) * (mom2p + mom1p)
                + (rot3p - rot2p) * (mom3p + mom2p)
                + rot1n * mom1n
                + (rot2n - rot1n) * (mom2n + mom1n)
                + (rot3n - rot2n) * (mom3n + mom2n));

        let e1p = mom1p / rot1p;
        let e2p = (mom2p - mom1p) / (rot2p - rot1p);
        let e3p = (mom3p - mom2p) / (rot3p - rot2p);
        let e1n = mom1n / rot1n;
        let e2n = (mom2n - mom1n) / (rot2n - rot1n);
        let e3n = (mom3n - mom2n) / (rot3n - rot2n);
        let eup = e1p.max(e2p).max(e3p);
        let eun = e1n.max(e2n).max(e3n);

        Material::Hysteretic(Box::new(HystereticFields {
            mom1p,
            rot1p,
            mom2p,
            rot2p,
            mom3p,
            rot3p,
            mom1n,
            rot1n,
            mom2n,
            rot2n,
            mom3n,
            rot3n,
            pinch_x,
            pinch_y,
            damfc1,
            damfc2,
            beta,
            e1p,
            e2p,
            e3p,
            e1n,
            e2n,
            e3n,
            eup,
            eun,
            energy_a,
            rot_max: 0.0,
            rot_min: 0.0,
            rot_pu: 0.0,
            rot_nu: 0.0,
            energy_d: 0.0,
            load: HystereticLoad::None,
            strain: 0.0,
            stress: 0.0,
        }))
    }
}

/// Ported from `HystereticMaterial::setTrialStrain`
/// (`HystereticMaterial.cpp`). Two of OpenSees' early-return guards
/// aren't ported, for reasons consistent with the other materials
/// above:
///
/// - `TloadIndicator == 0 && strain == 0.0`: returns the untouched
///   initial `(Tstress, Ttangent) = (0.0, E1p)` — the general formula
///   below reproduces this exactly anyway (checked by
///   `hysteretic_reproduces_initial_tangent_at_zero_strain`).
/// - `|dStrain| < DBL_EPSILON`: exists to skip re-deriving a stress/
///   tangent OpenSees would just recompute identically; the one place
///   this isn't quite idempotent (unlike every other material here) is
///   the interior branch's `dStrain < 0.0` / `dStrain > 0.0` split,
///   which has no case at all for `dStrain == 0.0` — see the comment
///   below for how that's resolved.
pub(super) fn evaluate(m: &Material, strain: f64) -> (f64, f64, Material) {
    let Material::Hysteretic(fields) = m else {
        unreachable!()
    };
    let HystereticFields {
        mom1p,
        rot1p,
        mom2p,
        rot2p,
        mom3p,
        rot3p,
        mom1n,
        rot1n,
        mom2n,
        rot2n,
        mom3n,
        rot3n,
        pinch_x,
        pinch_y,
        damfc1,
        damfc2,
        beta,
        e1p,
        e2p,
        e3p,
        e1n,
        e2n,
        e3n,
        eup,
        eun,
        energy_a,
        mut rot_max,
        mut rot_min,
        mut rot_pu,
        mut rot_nu,
        energy_d: energy_d_committed,
        mut load,
        strain: cstrain,
        stress: cstress,
    } = **fields;

    let envelope = HystereticEnvelope {
        mom1p,
        rot1p,
        mom2p,
        rot2p,
        mom3p,
        rot3p,
        e1p,
        e2p,
        e3p,
        mom1n,
        rot1n,
        mom2n,
        rot2n,
        mom3n,
        rot3n,
        e1n,
        e2n,
        e3n,
    };

    let dstrain = strain - cstrain;

    if load == HystereticLoad::None {
        load = if dstrain < 0.0 { HystereticLoad::Negative } else { HystereticLoad::Positive };
    }

    let (stress, tangent) = if strain >= rot_max {
        rot_max = strain;
        load = HystereticLoad::Positive;
        (envelope.pos_stress(strain), envelope.pos_tangent(strain))
    } else if strain <= rot_min {
        rot_min = strain;
        load = HystereticLoad::Negative;
        (envelope.neg_stress(strain), envelope.neg_tangent(strain))
    } else if dstrain < 0.0 {
        let (stress, tangent, new_rot_min, new_rot_pu, new_load) = hysteretic_negative_increment(
            dstrain, strain, cstrain, cstress, rot_max, rot_min, rot_pu, load, rot1p, rot1n, beta, eup, eun,
            pinch_x, pinch_y, damfc1, damfc2, energy_a, energy_d_committed, &envelope,
        );
        rot_min = new_rot_min;
        rot_pu = new_rot_pu;
        load = new_load;
        (stress, tangent)
    } else {
        // `dstrain >= 0.0`: OpenSees only calls `positiveIncrement` for
        // `dStrain > 0.0` (the `dStrain == 0.0` case falls through with
        // whatever `Tstress`/`Ttangent` happened to hold from before —
        // meaningless for a pure function with no "before"). Treating
        // `dstrain == 0.0` as belonging to this branch instead gives a
        // well-defined answer that's consistent with the tie-break
        // this same file already uses for the *first ever* increment
        // (`dStrain < 0.0 ? Negative : Positive`, just above) — and
        // `hysteretic_positive_increment` handles `dstrain == 0.0`
        // gracefully (`tmpmo1` reduces to `cstress` exactly).
        let (stress, tangent, new_rot_max, new_rot_nu, new_load) = hysteretic_positive_increment(
            dstrain, strain, cstrain, cstress, rot_max, rot_min, rot_nu, load, rot1p, rot1n, beta, eup, eun,
            pinch_x, pinch_y, damfc1, damfc2, energy_a, energy_d_committed, &envelope,
        );
        rot_max = new_rot_max;
        rot_nu = new_rot_nu;
        load = new_load;
        (stress, tangent)
    };

    let energy_d = energy_d_committed + 0.5 * (cstress + stress) * dstrain;

    (
        stress,
        tangent,
        Material::Hysteretic(Box::new(HystereticFields {
            mom1p,
            rot1p,
            mom2p,
            rot2p,
            mom3p,
            rot3p,
            mom1n,
            rot1n,
            mom2n,
            rot2n,
            mom3n,
            rot3n,
            pinch_x,
            pinch_y,
            damfc1,
            damfc2,
            beta,
            e1p,
            e2p,
            e3p,
            e1n,
            e2n,
            e3n,
            eup,
            eun,
            energy_a,
            rot_max,
            rot_min,
            rot_pu,
            rot_nu,
            energy_d,
            load,
            strain,
            stress,
        })),
    )
}

/// Ported from `HystereticMaterial::positiveIncrement`
/// (`HystereticMaterial.cpp`). Called only for the "interior" branch —
/// `rot_min < strain < rot_max` hasn't been hit yet by the caller's
/// envelope-extreme checks. `rot_max`/`rot_min`/`rot_nu`/`load` are the
/// values *entering* this call (`rot_min` and `rot_max` are read but only
/// `rot_max` is ever updated here); returns
/// `(stress, tangent, new_rot_max, new_rot_nu, new_load)`.
#[allow(clippy::too_many_arguments)]
fn hysteretic_positive_increment(
    dstrain: f64,
    strain: f64,
    cstrain: f64,
    cstress: f64,
    rot_max: f64,
    rot_min: f64,
    rot_nu: f64,
    mut load: HystereticLoad,
    rot1p: f64,
    rot1n: f64,
    beta: f64,
    eup: f64,
    eun: f64,
    pinch_x: f64,
    pinch_y: f64,
    damfc1: f64,
    damfc2: f64,
    energy_a: f64,
    energy_d: f64,
    envelope: &HystereticEnvelope,
) -> (f64, f64, f64, f64, HystereticLoad) {
    let kn_raw = (rot_min / rot1n).powf(beta);
    let kn = if kn_raw < 1.0 { 1.0 } else { 1.0 / kn_raw };
    let kp_raw = (rot_max / rot1p).powf(beta);
    let kp = if kp_raw < 1.0 { 1.0 } else { 1.0 / kp_raw };

    let mut rot_max = rot_max;
    let mut rot_nu = rot_nu;

    // Reversal from negative to positive loading: capture the pinching
    // reference point `rot_nu` and apply ductility-based damage to
    // `rot_max`, but only if the material had actually crossed into
    // tension-side stress (`cstress <= 0.0`) before reversing.
    if load == HystereticLoad::Negative && cstress <= 0.0 {
        rot_nu = cstrain - cstress / (eun * kn);
        let energy = energy_d - 0.5 * cstress / (eun * kn) * cstress;
        let mut damfc = 0.0;
        if rot_min < rot1n {
            damfc = damfc2 * energy / energy_a;
            damfc += damfc1 * (rot_min - rot1n) / rot1n;
        }
        rot_max *= 1.0 + damfc;
    }
    load = HystereticLoad::Positive;

    if rot_max > HYSTERETIC_POS_INF_STRAIN {
        rot_max = HYSTERETIC_POS_INF_STRAIN;
    }
    rot_max = if rot_max > rot1p { rot_max } else { rot1p };

    let maxmom = envelope.pos_stress(rot_max);
    let rotlim = envelope.neg_rotlim(rot_min);
    let rotrel = if rotlim > rot_nu { rotlim } else { rot_nu };

    let rotmp2 = rot_max - (1.0 - pinch_y) * maxmom / (eup * kp);
    let rotch = rotrel + (rotmp2 - rotrel) * pinch_x;

    let (stress, tangent) = if strain < rot_nu {
        let t = eun * kn;
        let s = cstress + t * dstrain;
        if s >= 0.0 {
            (0.0, eun * 1.0e-9)
        } else {
            (s, t)
        }
    } else if strain < rotch {
        if strain <= rotrel {
            (0.0, eup * 1.0e-9)
        } else {
            let t = maxmom * pinch_y / (rotch - rotrel);
            let tmpmo1 = cstress + eup * kp * dstrain;
            let tmpmo2 = (strain - rotrel) * t;
            if tmpmo1 < tmpmo2 {
                (tmpmo1, eup * kp)
            } else {
                (tmpmo2, t)
            }
        }
    } else {
        let t = (1.0 - pinch_y) * maxmom / (rot_max - rotch);
        let tmpmo1 = cstress + eup * kp * dstrain;
        let tmpmo2 = pinch_y * maxmom + (strain - rotch) * t;
        if tmpmo1 < tmpmo2 {
            (tmpmo1, eup * kp)
        } else {
            (tmpmo2, t)
        }
    };

    (stress, tangent, rot_max, rot_nu, load)
}

/// Ported from `HystereticMaterial::negativeIncrement`
/// (`HystereticMaterial.cpp`) — the mirror image of
/// `hysteretic_positive_increment` above; see its doc comment.
#[allow(clippy::too_many_arguments)]
fn hysteretic_negative_increment(
    dstrain: f64,
    strain: f64,
    cstrain: f64,
    cstress: f64,
    rot_max: f64,
    rot_min: f64,
    rot_pu: f64,
    mut load: HystereticLoad,
    rot1p: f64,
    rot1n: f64,
    beta: f64,
    eup: f64,
    eun: f64,
    pinch_x: f64,
    pinch_y: f64,
    damfc1: f64,
    damfc2: f64,
    energy_a: f64,
    energy_d: f64,
    envelope: &HystereticEnvelope,
) -> (f64, f64, f64, f64, HystereticLoad) {
    let kn_raw = (rot_min / rot1n).powf(beta);
    let kn = if kn_raw < 1.0 { 1.0 } else { 1.0 / kn_raw };
    let kp_raw = (rot_max / rot1p).powf(beta);
    let kp = if kp_raw < 1.0 { 1.0 } else { 1.0 / kp_raw };

    let mut rot_min = rot_min;
    let mut rot_pu = rot_pu;

    if load == HystereticLoad::Positive && cstress >= 0.0 {
        rot_pu = cstrain - cstress / (eup * kp);
        let energy = energy_d - 0.5 * cstress / (eup * kp) * cstress;
        let mut damfc = 0.0;
        if rot_max > rot1p {
            damfc = damfc2 * energy / energy_a;
            damfc += damfc1 * (rot_max - rot1p) / rot1p;
        }
        rot_min *= 1.0 + damfc;
    }
    load = HystereticLoad::Negative;

    if rot_min < HYSTERETIC_NEG_INF_STRAIN {
        rot_min = HYSTERETIC_NEG_INF_STRAIN;
    }
    rot_min = if rot_min < rot1n { rot_min } else { rot1n };

    let minmom = envelope.neg_stress(rot_min);
    let rotlim = envelope.pos_rotlim(rot_max);
    let rotrel = if rotlim < rot_pu { rotlim } else { rot_pu };

    let rotmp2 = rot_min - (1.0 - pinch_y) * minmom / (eun * kn);
    let rotch = rotrel + (rotmp2 - rotrel) * pinch_x;

    let (stress, tangent) = if strain > rot_pu {
        let t = eup * kp;
        let s = cstress + t * dstrain;
        if s <= 0.0 {
            (0.0, eup * 1.0e-9)
        } else {
            (s, t)
        }
    } else if strain > rotch {
        if strain >= rotrel {
            (0.0, eun * 1.0e-9)
        } else {
            let t = minmom * pinch_y / (rotch - rotrel);
            let tmpmo1 = cstress + eun * kn * dstrain;
            let tmpmo2 = (strain - rotrel) * t;
            if tmpmo1 > tmpmo2 {
                (tmpmo1, eun * kn)
            } else {
                (tmpmo2, t)
            }
        }
    } else {
        let t = (1.0 - pinch_y) * minmom / (rot_min - rotch);
        let tmpmo1 = cstress + eun * kn * dstrain;
        let tmpmo2 = pinch_y * minmom + (strain - rotch) * t;
        if tmpmo1 > tmpmo2 {
            (tmpmo1, eun * kn)
        } else {
            (tmpmo2, t)
        }
    };

    (stress, tangent, rot_min, rot_pu, load)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A symmetric, 3-point-per-side, bilinear-then-flat backbone shared by
    /// every `Hysteretic` test below: `mom1p/rot1p=10/0.01` (E1p=1000),
    /// `mom2p/rot2p=15/0.02` (E2p=500), `mom3p/rot3p=15/0.03` (E3p=0, flat
    /// plateau) — negative side the exact mirror. `beta=0` makes the
    /// ductility-degradation factors `kn`/`kp` identically `1.0`
    /// throughout (`x^0 == 1`, and `1.0 < 1.0` is false so the `1.0/kn`
    /// branch never engages either) — deliberately excluded from these
    /// tests to keep the hand computation tractable; `damfc1=damfc2=0`
    /// disables the deformation/energy damage factor for the same reason.
    fn hysteretic_test_material() -> Material {
        Material::hysteretic(
            10.0, 0.01, 15.0, 0.02, 15.0, 0.03, -10.0, -0.01, -15.0, -0.02, -15.0, -0.03, 0.5, 0.5, 0.0, 0.0, 0.0,
        )
    }

    #[test]
    fn hysteretic_gives_exact_initial_tangent_at_zero_strain() {
        let m = hysteretic_test_material();
        assert_eq!(m.trial_stress_tangent(0.0), (0.0, 1000.0));
    }

    #[test]
    fn hysteretic_matches_bilinear_envelope_on_first_loading() {
        let m = hysteretic_test_material();
        // Within the first segment (strain <= rot1p=0.01): E1p=1000.
        assert_eq!(m.trial_stress_tangent(0.005), (5.0, 1000.0));
        // Within the second segment (rot1p < strain <= rot2p=0.02):
        // mom1p + E2p*(strain-rot1p) = 10 + 500*0.005 = 12.5.
        assert_eq!(m.trial_stress_tangent(0.015), (12.5, 500.0));
    }

    #[test]
    fn hysteretic_pinches_on_reload_after_a_full_reversal() {
        let m = hysteretic_test_material();
        // Push to the positive envelope's second breakpoint (rot_max
        // becomes exactly rot2p=0.02, stress=mom2p=15), then reverse all
        // the way past the negative breakpoint too (rot_min becomes
        // exactly rot2n=-0.02, stress=mom2n=-15) — both reach the
        // envelope-extreme branch directly, since rot_min/rot_max both
        // start at 0.0 (any first excursion in either direction is "new").
        let m = m.commit(0.02).commit(-0.02);

        // Now reload from -0.02 toward 0.02 — strictly interior to
        // [rot_min, rot_max], so this hits `hysteretic_positive_increment`
        // for the first time, which captures the pinching reference point
        // `rot_nu = cstrain - cstress/(Eun*kn) = -0.02 - (-15)/1000 =
        // -0.005` off the just-completed reversal.
        //
        // At strain=-0.01 (still below rot_nu=-0.005): pure elastic
        // reloading at slope Eun*kn=1000 from the reversal point:
        // stress = cstress + 1000*dstrain = -15 + 1000*0.01 = -5.
        let (stress, tangent) = m.trial_stress_tangent(-0.01);
        assert!((stress - -5.0).abs() < 1e-9, "got {stress}");
        assert!((tangent - 1000.0).abs() < 1e-9, "got {tangent}");

        // At strain=0.0 (past rot_nu, in the pinched-plateau region below
        // rotch): the pinched-target line `(strain-rotrel)*maxmom*pinchY/
        // (rotch-rotrel)` governs rather than the steep elastic-unload
        // line — the hallmark of pinching (a much softer path back through
        // zero than straight elastic reload would give).
        let (stress, tangent) = m.trial_stress_tangent(0.0);
        assert!((stress - 4.285714285714286).abs() < 1e-9, "got {stress}");
        assert!((tangent - 857.1428571428571).abs() < 1e-9, "got {tangent}");
        assert!(tangent < 1000.0, "pinched stiffness should be softer than the elastic Eun, got {tangent}");
    }

    #[test]
    fn hysteretic_energy_accumulates_as_half_trapezoid_of_stress_over_strain() {
        let m = hysteretic_test_material();
        let m = m.commit(0.01); // exactly rot1p: stress=mom1p=10.
        let Material::Hysteretic(fields) = m else { panic!() };
        let energy_d = fields.energy_d;
        // Triangular area under the elastic segment: 0.5*10*0.01 = 0.05.
        assert!((energy_d - 0.05).abs() < 1e-9, "got {energy_d}");
    }
}
