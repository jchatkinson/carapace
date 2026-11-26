//! Sheathing shear laws shared by stud-wall piers and diaphragm bays.

use crate::model::Material;

/// Per-pier sheathing law: shear flow `q` (force per length) against sheathing shear strain `γ_sh`. It is the
/// whole sheathing assembly's response (panel shear and fastener slip, both sides); plate and stud axial
/// flexibility are added by the element. Every pier of the wall uses the same law with its own state.
#[derive(Debug, Clone)]
pub enum WallShearLaw {
    /// `q = ga·γ_sh`. Linear: load cases superpose and `Algorithm::Linear` is exact.
    Linear {
        /// Apparent shear stiffness (force per length), e.g. SDPWS `Ga` in consistent units.
        ga: f64,
    },
    /// CSA O86:19 11.7.1.2 shear and nail-slip terms per unit height, nonlinear elastic and odd in `q`:
    /// `γ_sh = q/Bv + 0.0025·e_n`, `e_n = (0.013·q·s/d_F²)²` (Annex A.11.7: nails in wood-based sheathing on
    /// dry lumber). **N and mm only.** Needs an iterative algorithm (Newton).
    O86 {
        /// Shear-through-thickness rigidity `Bv` (N/mm).
        bv: f64,
        /// Panel-edge nail spacing `s` (mm). Use `s·√2` for green lumber (doubles `e_n`).
        nail_spacing: f64,
        /// Nail diameter `d_F` (mm).
        nail_diameter: f64,
    },
    /// `q = σ(γ_sh)` from a uniaxial material with history, committed per pier at each converged step.
    /// The initial tangent must be positive. Needs an iterative algorithm (Newton).
    Material(Material),
}

/// Evaluate a shared wall or diaphragm bay shear law at one generalized strain.
pub(crate) fn shear_law_response(
    law: &WallShearLaw,
    material: Option<&Material>,
    gamma: f64,
) -> (f64, f64) {
    match law {
        WallShearLaw::Linear { ga } => (*ga * gamma, *ga),
        WallShearLaw::O86 {
            bv,
            nail_spacing,
            nail_diameter,
        } => {
            let a = 0.0025 * (0.013 * nail_spacing / nail_diameter.powi(2)).powi(2);
            let den = 1. / bv + (1. / (bv * bv) + 4. * a * gamma.abs()).sqrt();
            let q = 2. * gamma / den;
            (q, 1. / (1. / bv + 2. * a * q.abs()))
        }
        WallShearLaw::Material(_) => material
            .expect("material law has a state")
            .trial_stress_tangent(gamma),
    }
}
