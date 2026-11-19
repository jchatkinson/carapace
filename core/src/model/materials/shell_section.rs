//! `ShellSection`: the section seam for 3D shell elements, the shell counterpart of
//! `PlaneMaterial`. A value holds only *committed* state; `trial_resultant` is pure and
//! `commit` returns the new committed value. An element owns one independent copy per
//! Gauss point.
//!
//! The generalized strain and resultant vectors follow OpenSees' `ElasticMembranePlateSection`
//! so recorded values compare directly: `[eps_x, eps_y, gamma_xy, kappa_x, kappa_y,
//! 2 kappa_xy, gamma_xz, gamma_yz]` and `[Nx, Ny, Nxy, Mx, My, Mxy, Qx, Qy]` in the element's
//! local axes. OpenSees' bending resultants carry the opposite sign of the usual
//! `M = D kappa` (its `kappa` is `+w,xx`), so `D`'s bending block is negative; the element
//! compensates in its `B` matrix.

use std::fmt;

use nalgebra::{SMatrix, SVector};

/// Generalized shell strain or resultant: membrane (3), bending (3), transverse shear (2).
pub type ShellVector = SVector<f64, 8>;
/// The `8x8` section tangent.
pub type ShellMatrix = SMatrix<f64, 8, 8>;

/// Transverse shear correction factor `5/6`.
const SHEAR_CORRECTION: f64 = 5.0 / 6.0;

/// ETABS-style stiffness modifiers: factors on individual terms of the section tangent, all `1`
/// by default. `f11 f22 f12` scale the membrane block (`f12` the in-plane shear term), `m11 m22 m12`
/// the bending block (`m12` the twisting term) and `v13 v23` the transverse shear. The Poisson
/// coupling terms scale by the geometric mean of the two diagonal factors, which keeps the tangent
/// symmetric positive definite. Mass is not modified. `v13` and `v23` have no effect on `Shell3`,
/// which has no transverse shear.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellModifiers {
    pub f11: f64,
    pub f22: f64,
    pub f12: f64,
    pub m11: f64,
    pub m22: f64,
    pub m12: f64,
    pub v13: f64,
    pub v23: f64,
}

impl Default for ShellModifiers {
    fn default() -> Self {
        ShellModifiers {
            f11: 1.0,
            f22: 1.0,
            f12: 1.0,
            m11: 1.0,
            m22: 1.0,
            m12: 1.0,
            v13: 1.0,
            v23: 1.0,
        }
    }
}

impl ShellModifiers {
    fn all(&self) -> [f64; 8] {
        [
            self.f11, self.f22, self.f12, self.m11, self.m22, self.m12, self.v13, self.v23,
        ]
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ShellSection {
    /// Isotropic linear elastic membrane, plate bending and transverse shear
    /// (`ElasticMembranePlateSection`): modulus `e`, Poisson ratio `nu`, thickness `h`,
    /// mass density `rho`, and optional stiffness modifiers.
    ElasticMembranePlate {
        e: f64,
        nu: f64,
        h: f64,
        rho: f64,
        modifiers: ShellModifiers,
    },
}

/// Why a `ShellSection`'s constants are unusable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ShellSectionError {
    InvalidModulus,
    InvalidPoissonRatio,
    InvalidThickness,
    InvalidDensity,
    InvalidModifier,
}

impl ShellSectionError {
    pub fn message(&self) -> &'static str {
        match self {
            ShellSectionError::InvalidModulus => "modulus is not finite and positive",
            ShellSectionError::InvalidPoissonRatio => "Poisson ratio must lie in (-1, 0.5]",
            ShellSectionError::InvalidThickness => "thickness is not finite and positive",
            ShellSectionError::InvalidDensity => "density is not finite and non-negative",
            ShellSectionError::InvalidModifier => "stiffness modifier is not finite and positive",
        }
    }
}

impl fmt::Display for ShellSectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.message())
    }
}

impl std::error::Error for ShellSectionError {}

impl ShellSection {
    pub fn elastic_membrane_plate(
        e: f64,
        nu: f64,
        h: f64,
        rho: f64,
    ) -> Result<Self, ShellSectionError> {
        Self::elastic_membrane_plate_modified(e, nu, h, rho, ShellModifiers::default())
    }

    /// As [`ShellSection::elastic_membrane_plate`], with stiffness modifiers.
    pub fn elastic_membrane_plate_modified(
        e: f64,
        nu: f64,
        h: f64,
        rho: f64,
        modifiers: ShellModifiers,
    ) -> Result<Self, ShellSectionError> {
        let section = ShellSection::ElasticMembranePlate {
            e,
            nu,
            h,
            rho,
            modifiers,
        };
        section.validate()?;
        Ok(section)
    }

    /// Checks the constants of a value built from the variants directly.
    pub fn validate(&self) -> Result<(), ShellSectionError> {
        match *self {
            ShellSection::ElasticMembranePlate {
                e,
                nu,
                h,
                rho,
                modifiers,
            } => {
                if !(e.is_finite() && e > 0.0) {
                    return Err(ShellSectionError::InvalidModulus);
                }
                if !(nu.is_finite() && nu > -1.0 && nu <= 0.5) {
                    return Err(ShellSectionError::InvalidPoissonRatio);
                }
                if !(h.is_finite() && h > 0.0) {
                    return Err(ShellSectionError::InvalidThickness);
                }
                if !(rho.is_finite() && rho >= 0.0) {
                    return Err(ShellSectionError::InvalidDensity);
                }
                if !modifiers.all().iter().all(|f| f.is_finite() && *f > 0.0) {
                    return Err(ShellSectionError::InvalidModifier);
                }
                Ok(())
            }
        }
    }

    /// Whether the response is linear (independent of history).
    pub fn is_linear(&self) -> bool {
        true
    }

    /// Mass per unit area, `rho * h`.
    pub fn rho_h(&self) -> f64 {
        match *self {
            ShellSection::ElasticMembranePlate { h, rho, .. } => rho * h,
        }
    }

    /// The section tangent (OpenSees sign convention, see the module doc comment).
    pub fn tangent(&self) -> ShellMatrix {
        match *self {
            ShellSection::ElasticMembranePlate {
                e,
                nu,
                h,
                modifiers: k,
                ..
            } => {
                let m = e / (1.0 - nu * nu) * h;
                let g = 0.5 * e / (1.0 + nu) * h;
                let d = e * h * h * h / 12.0 / (1.0 - nu * nu);
                let gs = g * SHEAR_CORRECTION;
                // Each square root is taken on its own so extreme factors cannot underflow a product.
                let coupling_f = k.f11.sqrt() * k.f22.sqrt();
                let coupling_m = k.m11.sqrt() * k.m22.sqrt();
                let mut t = ShellMatrix::zeros();
                t[(0, 0)] = m * k.f11;
                t[(1, 1)] = m * k.f22;
                t[(0, 1)] = nu * m * coupling_f;
                t[(1, 0)] = nu * m * coupling_f;
                t[(2, 2)] = g * k.f12;
                t[(3, 3)] = -d * k.m11;
                t[(4, 4)] = -d * k.m22;
                t[(3, 4)] = -nu * d * coupling_m;
                t[(4, 3)] = -nu * d * coupling_m;
                t[(5, 5)] = -0.5 * d * (1.0 - nu) * k.m12;
                t[(6, 6)] = gs * k.v13;
                t[(7, 7)] = gs * k.v23;
                t
            }
        }
    }

    /// Resultants and tangent at a trial strain.
    pub fn trial_resultant(&self, strain: &ShellVector) -> (ShellVector, ShellMatrix) {
        let tangent = self.tangent();
        (tangent * strain, tangent)
    }

    /// The drilling penalty stiffness OpenSees uses: the smallest eigenvalue of the
    /// membrane block of the initial tangent (with the modifiers applied).
    pub fn drilling_stiffness(&self) -> f64 {
        match *self {
            ShellSection::ElasticMembranePlate {
                e,
                nu,
                h,
                modifiers: k,
                ..
            } => {
                let m = e / (1.0 - nu * nu) * h;
                if k.f11 == 1.0 && k.f22 == 1.0 && k.f12 == 1.0 {
                    // Eigenvalues of m * [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]].
                    return (m * (1.0 - nu) / 2.0)
                        .min(m * (1.0 - nu))
                        .min(m * (1.0 + nu));
                }
                // The 2x2 block [[a, c], [c, b]] has c^2 = nu^2 a b, so det = a b (1 - nu^2).
                // The smaller eigenvalue is det / larger, which does not cancel.
                let a = m * k.f11;
                let b = m * k.f22;
                let c = nu * m * k.f11.sqrt() * k.f22.sqrt();
                let larger = 0.5 * (a + b) + (0.5 * (a - b)).hypot(c);
                let smaller = (a / larger) * b * (1.0 - nu * nu);
                let shear = 0.5 * e / (1.0 + nu) * h * k.f12;
                smaller.min(shear)
            }
        }
    }

    /// The new committed value at `strain`. Elastic sections carry no history.
    pub fn commit(&self, _strain: &ShellVector) -> Self {
        self.clone()
    }
}
