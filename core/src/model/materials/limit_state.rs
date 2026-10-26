//! `LimitKind`/`LimitStateConfig`: the material-limit-state vocabulary CSFM's
//! concrete and rebar laws report through their own `limit_utilization`
//! queries (plan `docs/csfm-plan.md` §2.5 point 3), and that
//! `ElementOps::limit_utilization`/`Domain::controlling_limit_utilization`
//! (`model/elements/ops.rs`, `model/domain.rs`) reduce over an element and a
//! whole model. CSFM's answer to "how far loaded" is a utilization, not a
//! converged-or-not flag: a value >= 1 means the limit is reached, and
//! values are never clamped, so a load path can be checked just below and
//! just above a landing to distinguish exact landing (`== 1`) from a
//! strict-`>` rupture test (`Material::csfm_rebar_ruptured`).

/// Which material limit a `limit_utilization` query reports (plan §2.5
/// point 3's four kinds).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LimitKind {
    /// Concrete crushing: `max(0, -eps_2) / limit_c`, where `limit_c` is
    /// `|eps_cu2|` (plan §6 decision 5: `-0.0035` for both compression laws)
    /// unless `LimitStateConfig::concrete_crushing_limit` overrides it.
    ConcreteCrushing,
    /// Concrete tensile strain: `max(0, eps_1) / limit_t`. Unlike crushing,
    /// `limit_t` has no material-intrinsic default — it is purely
    /// problem-specific (e.g. `findCriticalLoad`'s own
    /// `maxConcreteTensileStrain` stop criterion) — so this kind is never
    /// reported unless `LimitStateConfig::concrete_tensile_limit` is given.
    ConcreteTensileStrain,
    /// Rebar rupture, same formula split as `Material::csfm_rebar_ruptured`
    /// (plan §1): tension with tension stiffening, `eps / eps_m(ft)`;
    /// compression or bare steel (same formula on both, since `-eps` and
    /// `|eps|` agree when `eps <= 0`), `|eps| / eps_u`.
    RebarRupture,
    /// `|eps| / limit`, a user-supplied rebar strain limit. Like
    /// `ConcreteTensileStrain`, there is no material-intrinsic default, so
    /// this is never reported unless `LimitStateConfig::rebar_strain_limit`
    /// is given.
    RebarStrainLimit,
}

/// Parameters a `limit_utilization` query needs beyond what the material
/// itself resolved at construction — exactly the ones with no sensible
/// material-intrinsic default (plan §2.5 point 3). `Default` (every field
/// `None`) evaluates only `ConcreteCrushing` (at the plan §6 decision-5
/// default `|eps_cu2|`) and `RebarRupture` — the two kinds with such a
/// default — and never `ConcreteTensileStrain`/`RebarStrainLimit`, which are
/// skipped entirely rather than falling back to some arbitrary value.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LimitStateConfig {
    /// Overrides the default compressive strain limit `|eps_cu2| = 0.0035`
    /// (plan §6 decision 5) for `ConcreteCrushing`. Must be finite and
    /// positive when given; an invalid value is treated as absent (falls
    /// back to the default) by the evaluating function, which never panics
    /// on caller input.
    pub concrete_crushing_limit: Option<f64>,
    /// Enables `ConcreteTensileStrain` at this strain limit. Must be finite
    /// and positive to take effect.
    pub concrete_tensile_limit: Option<f64>,
    /// Enables `RebarStrainLimit` at this strain limit. Must be finite and
    /// positive to take effect.
    pub rebar_strain_limit: Option<f64>,
}

impl LimitStateConfig {
    /// `value` if finite and positive, else `None` — every optional limit
    /// in this struct is evaluated through this filter so a caller-supplied
    /// `NaN`/non-positive override never propagates into a utilization.
    pub(super) fn usable(value: Option<f64>) -> Option<f64> {
        value.filter(|v| v.is_finite() && *v > 0.0)
    }
}
