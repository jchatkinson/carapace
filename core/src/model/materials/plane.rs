//! `PlaneMaterial`: the n-D material seam for 2D continuum elements.
//!
//! Same trial/commit design as the uniaxial `Material` (see its doc
//! comment): a value holds only *committed* state, `trial_stress_tangent` is
//! pure, and `commit` returns the new committed value. An element owns one
//! independent `PlaneMaterial` per Gauss point, cloned from the prototype
//! given to its constructor and committed individually, exactly as the fiber
//! elements own one material per fiber. That is the seam where nonlinear
//! plane materials (CSFM, damage, plasticity) plug in with per-point history.
//!
//! Conventions, fixed now so nonlinear materials later do not change them:
//! the strain vector is `[eps_x, eps_y, gamma_xy]` with **engineering
//! shear** `gamma_xy = 2 eps_xy`; the stress vector is `[sigma_x, sigma_y,
//! tau_xy]`; `D` maps one to the other.

use std::fmt;

use nalgebra::{SMatrix, SVector};

/// A plane strain or stress vector `[x, y, xy]` (engineering shear for strain).
pub type PlaneVector = SVector<f64, 3>;
/// The `3x3` plane constitutive tangent.
pub type PlaneMatrix = SMatrix<f64, 3, 3>;

#[derive(Debug, Clone, PartialEq)]
pub enum PlaneMaterial {
    /// A general symmetric `D`, given as the packed upper triangle row by
    /// row: `[d11, d12, d13, d22, d23, d33]`, with respect to the
    /// engineering-shear strain vector.
    ElasticMatrix { d: [f64; 6] },
    /// Orthotropic plane stress. `angle` is the counter-clockwise angle in
    /// radians from global x to material axis 1; `nu_xy` is the major
    /// Poisson ratio (strain in axis 2 per strain in axis 1, loading along 1).
    Orthotropic {
        ex: f64,
        ey: f64,
        nu_xy: f64,
        g_xy: f64,
        angle: f64,
    },
    /// Isotropic plane stress (`sigma_z = 0`).
    PlaneStress { e: f64, nu: f64 },
    /// Isotropic plane strain (`eps_z = 0`).
    PlaneStrain { e: f64, nu: f64 },
    /// Test-only stateful material: `D` of unit-modulus plane stress
    /// (`nu = 0`) acting on the strain measured from the last committed
    /// strain, so per-point history is observable.
    #[cfg(test)]
    Probe { last: [f64; 3] },
}

/// Why a `PlaneMaterial`'s constants are unusable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlaneMaterialError {
    /// A modulus is not finite and positive.
    InvalidModulus,
    /// A Poisson ratio is outside the range for which `D` is positive definite.
    InvalidPoissonRatio,
    /// The constitutive matrix is not finite, symmetric positive definite.
    NotPositiveDefinite,
}

impl fmt::Display for PlaneMaterialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            PlaneMaterialError::InvalidModulus => "a modulus is not finite and positive",
            PlaneMaterialError::InvalidPoissonRatio => "Poisson ratio is out of range",
            PlaneMaterialError::NotPositiveDefinite => {
                "the constitutive matrix is not symmetric positive definite"
            }
        })
    }
}

impl std::error::Error for PlaneMaterialError {}

impl PlaneMaterial {
    /// Isotropic plane stress; `e > 0` and `-1 < nu < 1`.
    pub fn plane_stress(e: f64, nu: f64) -> Result<Self, PlaneMaterialError> {
        Self::checked(PlaneMaterial::PlaneStress { e, nu })
    }

    /// Isotropic plane strain; `e > 0` and `-1 < nu < 0.5`.
    pub fn plane_strain(e: f64, nu: f64) -> Result<Self, PlaneMaterialError> {
        Self::checked(PlaneMaterial::PlaneStrain { e, nu })
    }

    /// Orthotropic plane stress; positive moduli and `nu_xy^2 < ex / ey`.
    pub fn orthotropic(
        ex: f64,
        ey: f64,
        nu_xy: f64,
        g_xy: f64,
        angle: f64,
    ) -> Result<Self, PlaneMaterialError> {
        Self::checked(PlaneMaterial::Orthotropic {
            ex,
            ey,
            nu_xy,
            g_xy,
            angle,
        })
    }

    /// A general symmetric positive-definite `D` (packed upper triangle).
    pub fn elastic_matrix(d: [f64; 6]) -> Result<Self, PlaneMaterialError> {
        Self::checked(PlaneMaterial::ElasticMatrix { d })
    }

    fn checked(material: PlaneMaterial) -> Result<Self, PlaneMaterialError> {
        material.validate()?;
        Ok(material)
    }

    /// Checks the constants of a value built from the variants directly.
    pub fn validate(&self) -> Result<(), PlaneMaterialError> {
        let positive = |v: f64| {
            if v.is_finite() && v > 0.0 {
                Ok(())
            } else {
                Err(PlaneMaterialError::InvalidModulus)
            }
        };
        let poisson = |nu: f64, upper: f64| {
            if nu.is_finite() && nu > -1.0 && nu < upper {
                Ok(())
            } else {
                Err(PlaneMaterialError::InvalidPoissonRatio)
            }
        };
        match *self {
            PlaneMaterial::PlaneStress { e, nu } => {
                positive(e)?;
                poisson(nu, 1.0)?;
            }
            PlaneMaterial::PlaneStrain { e, nu } => {
                positive(e)?;
                poisson(nu, 0.5)?;
            }
            PlaneMaterial::Orthotropic {
                ex,
                ey,
                nu_xy,
                g_xy,
                angle,
            } => {
                positive(ex)?;
                positive(ey)?;
                positive(g_xy)?;
                if !nu_xy.is_finite() || !angle.is_finite() || nu_xy * nu_xy >= ex / ey {
                    return Err(PlaneMaterialError::InvalidPoissonRatio);
                }
            }
            PlaneMaterial::ElasticMatrix { d } => {
                let m = self.d_matrix();
                if d.iter().any(|v| !v.is_finite()) || m.cholesky().is_none() {
                    return Err(PlaneMaterialError::NotPositiveDefinite);
                }
            }
            #[cfg(test)]
            PlaneMaterial::Probe { .. } => {}
        }
        Ok(())
    }

    /// The constitutive matrix at the committed state (constant for every
    /// elastic variant).
    pub fn d_matrix(&self) -> PlaneMatrix {
        match *self {
            PlaneMaterial::ElasticMatrix { d } => PlaneMatrix::new(
                d[0], d[1], d[2], //
                d[1], d[3], d[4], //
                d[2], d[4], d[5],
            ),
            PlaneMaterial::PlaneStress { e, nu } => {
                let c = e / (1.0 - nu * nu);
                PlaneMatrix::new(
                    c,
                    c * nu,
                    0.0, //
                    c * nu,
                    c,
                    0.0, //
                    0.0,
                    0.0,
                    c * (1.0 - nu) / 2.0,
                )
            }
            PlaneMaterial::PlaneStrain { e, nu } => {
                let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
                PlaneMatrix::new(
                    c * (1.0 - nu),
                    c * nu,
                    0.0, //
                    c * nu,
                    c * (1.0 - nu),
                    0.0, //
                    0.0,
                    0.0,
                    c * (1.0 - 2.0 * nu) / 2.0,
                )
            }
            PlaneMaterial::Orthotropic {
                ex,
                ey,
                nu_xy,
                g_xy,
                angle,
            } => {
                let nu_yx = nu_xy * ey / ex;
                let k = 1.0 / (1.0 - nu_xy * nu_yx);
                let local = PlaneMatrix::new(
                    k * ex,
                    k * nu_xy * ey,
                    0.0, //
                    k * nu_xy * ey,
                    k * ey,
                    0.0, //
                    0.0,
                    0.0,
                    g_xy,
                );
                // Strain transformation (global -> material axes), engineering shear;
                // stress transforms by its inverse transpose, so D = T^T D' T.
                let (s, c) = angle.sin_cos();
                let t = PlaneMatrix::new(
                    c * c,
                    s * s,
                    c * s, //
                    s * s,
                    c * c,
                    -c * s, //
                    -2.0 * c * s,
                    2.0 * c * s,
                    c * c - s * s,
                );
                t.transpose() * local * t
            }
            #[cfg(test)]
            PlaneMaterial::Probe { .. } => PlaneMatrix::identity(),
        }
    }

    /// Stress and tangent at `strain`, always relative to the committed state.
    pub fn trial_stress_tangent(&self, strain: &PlaneVector) -> (PlaneVector, PlaneMatrix) {
        let d = self.d_matrix();
        #[cfg(test)]
        if let PlaneMaterial::Probe { last } = *self {
            return (d * (strain - PlaneVector::from(last)), d);
        }
        (d * strain, d)
    }

    /// The new committed material after `strain` converged.
    pub fn commit(&self, strain: &PlaneVector) -> Self {
        #[cfg(test)]
        if let PlaneMaterial::Probe { .. } = self {
            return PlaneMaterial::Probe {
                last: [strain[0], strain[1], strain[2]],
            };
        }
        let _ = strain;
        self.clone()
    }

    /// Whether the response is linear in strain and history-free (what the
    /// enhanced `Quad4` requires).
    pub fn is_linear(&self) -> bool {
        match self {
            PlaneMaterial::ElasticMatrix { .. }
            | PlaneMaterial::Orthotropic { .. }
            | PlaneMaterial::PlaneStress { .. }
            | PlaneMaterial::PlaneStrain { .. } => true,
            #[cfg(test)]
            PlaneMaterial::Probe { .. } => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: &PlaneMatrix, b: &PlaneMatrix, tol: f64) -> bool {
        (a - b).abs().max() <= tol * a.abs().max().max(1.0)
    }

    #[test]
    fn d_is_symmetric_and_positive_definite_for_every_variant() {
        let materials = [
            PlaneMaterial::plane_stress(200e3, 0.3).unwrap(),
            PlaneMaterial::plane_strain(200e3, 0.3).unwrap(),
            PlaneMaterial::orthotropic(30e3, 10e3, 0.25, 4e3, 0.6).unwrap(),
            PlaneMaterial::elastic_matrix([4.0, 1.0, 0.5, 3.0, 0.2, 2.0]).unwrap(),
        ];
        for m in materials {
            let d = m.d_matrix();
            assert!(close(&d, &d.transpose(), 1e-14), "{m:?}");
            assert!(d.cholesky().is_some(), "{m:?}");
        }
    }

    #[test]
    fn plane_stress_with_effective_constants_equals_plane_strain() {
        let (e, nu) = (200e3, 0.3);
        let strain = PlaneMaterial::plane_strain(e, nu).unwrap();
        let (e_eff, nu_eff) = (e / (1.0 - nu * nu), nu / (1.0 - nu));
        let stress = PlaneMaterial::plane_stress(e_eff, nu_eff).unwrap();
        assert!(close(&strain.d_matrix(), &stress.d_matrix(), 1e-13));
    }

    #[test]
    fn orthotropic_with_equal_moduli_is_isotropic_for_any_angle() {
        let (e, nu) = (100.0, 0.25);
        let iso = PlaneMaterial::plane_stress(e, nu).unwrap().d_matrix();
        for angle in [0.0, 0.3, 1.1, -2.0] {
            let ortho =
                PlaneMaterial::orthotropic(e, e, nu, e / (2.0 * (1.0 + nu)), angle).unwrap();
            assert!(close(&ortho.d_matrix(), &iso, 1e-13), "angle {angle}");
        }
    }

    #[test]
    fn rotating_by_ninety_degrees_swaps_the_axes() {
        let (ex, ey, nu_xy, g) = (30.0, 10.0, 0.2, 4.0);
        let d0 = PlaneMaterial::orthotropic(ex, ey, nu_xy, g, 0.0)
            .unwrap()
            .d_matrix();
        let d90 = PlaneMaterial::orthotropic(ex, ey, nu_xy, g, std::f64::consts::FRAC_PI_2)
            .unwrap()
            .d_matrix();
        assert!((d0[(0, 0)] - d90[(1, 1)]).abs() < 1e-12);
        assert!((d0[(1, 1)] - d90[(0, 0)]).abs() < 1e-12);
        assert!((d0[(0, 1)] - d90[(0, 1)]).abs() < 1e-12);
        assert!((d0[(2, 2)] - d90[(2, 2)]).abs() < 1e-12);
        assert!(d90[(0, 2)].abs() < 1e-12 && d90[(1, 2)].abs() < 1e-12);
    }

    #[test]
    fn rotation_by_pi_over_four_couples_normal_and_shear() {
        // Uniaxial strain along the fibre direction at 45 degrees: sigma_x - sigma_y is zero by symmetry,
        // and the rotated D has a nonzero normal-shear coupling for a strongly orthotropic ply.
        let d = PlaneMaterial::orthotropic(100.0, 5.0, 0.3, 4.0, std::f64::consts::FRAC_PI_4)
            .unwrap()
            .d_matrix();
        assert!((d[(0, 0)] - d[(1, 1)]).abs() < 1e-12);
        assert!(d[(0, 2)].abs() > 1.0);
    }

    #[test]
    fn elastic_matrix_packing_matches_the_equivalent_isotropic_d() {
        let iso = PlaneMaterial::plane_stress(100.0, 0.25).unwrap().d_matrix();
        let packed = [
            iso[(0, 0)],
            iso[(0, 1)],
            iso[(0, 2)],
            iso[(1, 1)],
            iso[(1, 2)],
            iso[(2, 2)],
        ];
        let m = PlaneMaterial::elastic_matrix(packed).unwrap();
        assert!(close(&m.d_matrix(), &iso, 1e-15));
        let strain = PlaneVector::new(1e-3, -2e-3, 3e-3);
        assert!((m.trial_stress_tangent(&strain).0 - iso * strain).norm() < 1e-12);
    }

    #[test]
    fn invalid_constants_are_errors_not_panics() {
        assert_eq!(
            PlaneMaterial::plane_stress(0.0, 0.3),
            Err(PlaneMaterialError::InvalidModulus)
        );
        assert_eq!(
            PlaneMaterial::plane_stress(-1.0, 0.3),
            Err(PlaneMaterialError::InvalidModulus)
        );
        assert_eq!(
            PlaneMaterial::plane_stress(f64::NAN, 0.3),
            Err(PlaneMaterialError::InvalidModulus)
        );
        assert_eq!(
            PlaneMaterial::plane_stress(1.0, 1.0),
            Err(PlaneMaterialError::InvalidPoissonRatio)
        );
        assert_eq!(
            PlaneMaterial::plane_strain(1.0, 0.5),
            Err(PlaneMaterialError::InvalidPoissonRatio)
        );
        assert_eq!(
            PlaneMaterial::plane_strain(1.0, -1.0),
            Err(PlaneMaterialError::InvalidPoissonRatio)
        );
        assert_eq!(
            PlaneMaterial::orthotropic(1.0, 10.0, 0.5, 1.0, 0.0),
            Err(PlaneMaterialError::InvalidPoissonRatio)
        );
        assert_eq!(
            PlaneMaterial::orthotropic(1.0, 1.0, 0.0, 0.0, 0.0),
            Err(PlaneMaterialError::InvalidModulus)
        );
        assert_eq!(
            PlaneMaterial::elastic_matrix([1.0, 2.0, 0.0, 1.0, 0.0, 1.0]),
            Err(PlaneMaterialError::NotPositiveDefinite)
        );
        assert_eq!(
            PlaneMaterial::elastic_matrix([f64::NAN, 0.0, 0.0, 1.0, 0.0, 1.0]),
            Err(PlaneMaterialError::NotPositiveDefinite)
        );
        // A struct literal that bypasses the constructors is caught by `validate`.
        assert!(PlaneMaterial::PlaneStress { e: -1.0, nu: 0.0 }
            .validate()
            .is_err());
    }

    #[test]
    fn elastic_variants_are_linear_and_commit_changes_nothing() {
        let m = PlaneMaterial::plane_stress(10.0, 0.2).unwrap();
        assert!(m.is_linear());
        assert_eq!(m.commit(&PlaneVector::new(1.0, 2.0, 3.0)), m);
    }

    #[test]
    fn a_stateful_material_keeps_independent_history_per_copy() {
        let proto = PlaneMaterial::Probe { last: [0.0; 3] };
        let (mut a, mut b) = (proto.clone(), proto);
        a = a.commit(&PlaneVector::new(1.0, 0.0, 0.0));
        b = b.commit(&PlaneVector::new(0.0, 2.0, 0.0));
        let probe = PlaneVector::new(1.0, 2.0, 0.0);
        assert_eq!(
            a.trial_stress_tangent(&probe).0,
            PlaneVector::new(0.0, 2.0, 0.0)
        );
        assert_eq!(
            b.trial_stress_tangent(&probe).0,
            PlaneVector::new(1.0, 0.0, 0.0)
        );
        assert!(!a.is_linear());
    }
}
