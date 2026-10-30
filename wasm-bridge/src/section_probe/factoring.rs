//! Resistance factors applied to section resultants, for factored PMM surfaces.
//!
//! - [`Factoring::Groups`]: one factor per fiber group (CSA A23.3: concrete and steel resistance factors).
//! - [`Factoring::Aci`]: a single `phi` from the largest tensile strain among the fibers of `bar_group`
//!   (ACI 318: `phi_compression` at or below `eps_y`, `phi_tension` at or above `eps_tl`, linear between). The
//!   resultants are scaled by that `phi`.
//!
//! The factored state is a function of the strain state alone, like the nominal one, so the isoline solver can
//! work in factored space directly.

use serde::{Deserialize, Serialize};

use super::{SectionProbe, SectionResponse};

#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum Factoring {
    /// `factors[g]` multiplies the force of every fiber in group `g`.
    Groups { factors: Vec<f64> },
    Aci {
        bar_group: u32,
        phi_tension: f64,
        phi_compression: f64,
        /// Yield strain of the reinforcement.
        eps_y: f64,
        /// Tension-controlled strain limit.
        #[serde(default = "default_eps_tl")]
        eps_tl: f64,
    },
}

fn default_eps_tl() -> f64 {
    0.005
}

/// Factored resultants and the factors behind them.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct FactoredResponse {
    pub response: SectionResponse,
    /// `phi` for [`Factoring::Aci`].
    pub phi: Option<f64>,
    /// Largest tensile strain of the bar group for [`Factoring::Aci`] (zero if none is in tension).
    pub eps_t: Option<f64>,
}

impl Factoring {
    pub(super) fn validate(&self, group_count: usize) -> Result<(), String> {
        match self {
            Factoring::Groups { factors } if factors.len() < group_count => Err(format!(
                "factoring needs one factor per fiber group ({group_count}), got {}",
                factors.len()
            )),
            Factoring::Aci { bar_group, .. } if *bar_group as usize >= group_count => {
                Err(format!("barGroup {bar_group} is not a fiber group"))
            }
            _ => Ok(()),
        }
    }

    fn aci_phi(&self, eps_t: f64) -> f64 {
        let Factoring::Aci {
            phi_tension,
            phi_compression,
            eps_y,
            eps_tl,
            ..
        } = self
        else {
            return 1.0;
        };
        let eps_y = eps_y.max(0.0);
        if eps_t >= *eps_tl {
            *phi_tension
        } else if eps_t <= eps_y {
            *phi_compression
        } else {
            phi_compression + (eps_t - eps_y) / (eps_tl - eps_y) * (phi_tension - phi_compression)
        }
    }
}

impl SectionProbe {
    /// Factored resultants at `(eps0, kz, ky)`.
    pub fn evaluate_factored(
        &self,
        factoring: &Factoring,
        eps0: f64,
        kz: f64,
        ky: f64,
    ) -> FactoredResponse {
        match factoring {
            Factoring::Groups { factors } => {
                let (mut n, mut mz, mut my) = (0.0, 0.0, 0.0);
                for i in 0..self.y.len() {
                    let force = factors[self.group[i] as usize] * self.fiber_force(i, eps0, kz, ky);
                    n += force;
                    mz -= self.y[i] * force;
                    my += self.z[i] * force;
                }
                FactoredResponse {
                    response: SectionResponse::new(n, mz, my),
                    phi: None,
                    eps_t: None,
                }
            }
            Factoring::Aci { bar_group, .. } => {
                let (mut n, mut mz, mut my, mut eps_t) = (0.0, 0.0, 0.0, 0.0f64);
                for i in 0..self.y.len() {
                    let force = self.fiber_force(i, eps0, kz, ky);
                    n += force;
                    mz -= self.y[i] * force;
                    my += self.z[i] * force;
                    if self.group[i] == *bar_group {
                        eps_t = eps_t.max(eps0 - self.y[i] * kz + self.z[i] * ky);
                    }
                }
                let phi = factoring.aci_phi(eps_t);
                FactoredResponse {
                    response: SectionResponse::new(phi * n, phi * mz, phi * my),
                    phi: Some(phi),
                    eps_t: Some(eps_t),
                }
            }
        }
    }
}
