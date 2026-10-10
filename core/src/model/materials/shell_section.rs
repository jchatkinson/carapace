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

#[derive(Debug, Clone, PartialEq)]
pub enum ShellSection {
    /// Isotropic linear elastic membrane, plate bending and transverse shear
    /// (`ElasticMembranePlateSection`): modulus `e`, Poisson ratio `nu`, thickness `h`,
    /// mass density `rho`.
    ElasticMembranePlate { e: f64, nu: f64, h: f64, rho: f64 },
}

/// Why a `ShellSection`'s constants are unusable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ShellSectionError {
    InvalidModulus,
    InvalidPoissonRatio,
    InvalidThickness,
    InvalidDensity,
}

impl ShellSectionError {
    pub fn message(&self) -> &'static str {
        match self {
            ShellSectionError::InvalidModulus => "modulus is not finite and positive",
            ShellSectionError::InvalidPoissonRatio => "Poisson ratio must lie in (-1, 0.5]",
            ShellSectionError::InvalidThickness => "thickness is not finite and positive",
            ShellSectionError::InvalidDensity => "density is not finite and non-negative",
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
        let section = ShellSection::ElasticMembranePlate { e, nu, h, rho };
        section.validate()?;
        Ok(section)
    }

    /// Checks the constants of a value built from the variants directly.
    pub fn validate(&self) -> Result<(), ShellSectionError> {
        match *self {
            ShellSection::ElasticMembranePlate { e, nu, h, rho } => {
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
            ShellSection::ElasticMembranePlate { e, nu, h, .. } => {
                let m = e / (1.0 - nu * nu) * h;
                let g = 0.5 * e / (1.0 + nu) * h;
                let d = e * h * h * h / 12.0 / (1.0 - nu * nu);
                let gs = g * SHEAR_CORRECTION;
                let mut t = ShellMatrix::zeros();
                t[(0, 0)] = m;
                t[(1, 1)] = m;
                t[(0, 1)] = nu * m;
                t[(1, 0)] = nu * m;
                t[(2, 2)] = g;
                t[(3, 3)] = -d;
                t[(4, 4)] = -d;
                t[(3, 4)] = -nu * d;
                t[(4, 3)] = -nu * d;
                t[(5, 5)] = -0.5 * d * (1.0 - nu);
                t[(6, 6)] = gs;
                t[(7, 7)] = gs;
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
    /// membrane block of the initial tangent.
    pub fn drilling_stiffness(&self) -> f64 {
        match *self {
            ShellSection::ElasticMembranePlate { e, nu, h, .. } => {
                let m = e / (1.0 - nu * nu) * h;
                // Eigenvalues of m * [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]].
                (m * (1.0 - nu) / 2.0)
                    .min(m * (1.0 - nu))
                    .min(m * (1.0 + nu))
            }
        }
    }

    /// The new committed value at `strain`. Elastic sections carry no history.
    pub fn commit(&self, _strain: &ShellVector) -> Self {
        self.clone()
    }
}
