/// Rayleigh (mass- and stiffness-proportional) damping: `C = alpha_m*M + beta_k*K`.
/// A concrete struct, not an enum — Rayleigh damping is the one form §4.4
/// asks for; other forms (e.g. explicit modal damping ratios) would be a
/// distinct, separate mechanism, not another variant of this one.
#[derive(Debug, Clone, Copy)]
pub struct RayleighDamping {
    pub alpha_m: f64,
    pub beta_k: f64,
}

impl RayleighDamping {
    pub const NONE: RayleighDamping = RayleighDamping {
        alpha_m: 0.0,
        beta_k: 0.0,
    };

    pub fn new(alpha_m: f64, beta_k: f64) -> Self {
        RayleighDamping { alpha_m, beta_k }
    }
}
