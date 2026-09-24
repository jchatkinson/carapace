use super::Material;

/// One fiber: a small patch of cross-section area at distance `y` from the
/// section centroid, carrying a `Material`. 2D scope only (§3.1's `NDM`) —
/// `y` alone, no `z`, since bending is about the section's local z-axis
/// only.
#[derive(Debug, Clone)]
pub struct Fiber {
    pub y: f64,
    pub area: f64,
    pub material: Material,
}

impl Fiber {
    pub fn new(y: f64, area: f64, material: Material) -> Self {
        Fiber { y, area, material }
    }
}

/// A fiber-discretized cross-section: given the section's axial strain
/// `eps0` (at the centroid) and curvature `kappa`, each fiber's strain is
/// `eps0 - y*kappa`. Work-conjugate to `[eps0, kappa]` is `[N, M]`
/// (`N = sum(stress*area)`, `M = -sum(stress*area*y)`) — derived directly
/// from virtual work, not assumed: `d(stress power) = integral(stress *
/// d(strain)) dA = integral(stress*(d(eps0) - y*d(kappa))) dA =
/// d(eps0)*N + d(kappa)*M` with `M` as defined above.
///
/// Each `DispBeamColumn` integration point owns its own `FiberSection`
/// (independent material history — see `Element`'s commit/trial split),
/// even when every point starts from an identical fiber layout (the
/// prismatic-member assumption this element makes: same cross-section
/// everywhere along the length — a tapered/non-prismatic member is a
/// natural extension, not built until needed).
#[derive(Debug, Clone)]
pub struct FiberSection {
    fibers: Vec<Fiber>,
}

/// `(N, M)` and the section tangent `[[dN/deps0, dN/dkappa], [dM/deps0, dM/dkappa]]`.
pub type SectionResponse = (f64, f64, [[f64; 2]; 2]);

impl FiberSection {
    pub fn new(fibers: Vec<Fiber>) -> Self {
        FiberSection { fibers }
    }

    /// Trial evaluation at `(eps0, kappa)` — read-only, safe every Newton
    /// iteration; see `Material`'s doc comment for why this must not
    /// mutate any fiber's material.
    pub fn trial(&self, eps0: f64, kappa: f64) -> SectionResponse {
        let (mut n, mut m) = (0.0, 0.0);
        let (mut ea, mut eq, mut ei) = (0.0, 0.0, 0.0);
        for fiber in &self.fibers {
            let strain = eps0 - fiber.y * kappa;
            let (stress, tangent) = fiber.material.trial_stress_tangent(strain);
            n += stress * fiber.area;
            m += -stress * fiber.area * fiber.y;
            ea += tangent * fiber.area;
            eq += tangent * fiber.area * fiber.y;
            ei += tangent * fiber.area * fiber.y * fiber.y;
        }
        (n, m, [[ea, -eq], [-eq, ei]])
    }

    /// Commit every fiber's material at `(eps0, kappa)` — called once per
    /// converged step, by `DispBeamColumn::commit`.
    pub fn commit(&mut self, eps0: f64, kappa: f64) {
        for fiber in &mut self.fibers {
            let strain = eps0 - fiber.y * kappa;
            fiber.material = fiber.material.commit(strain);
        }
    }

    pub fn total_area(&self) -> f64 {
        self.fibers.iter().map(|f| f.area).sum()
    }
}

/// `FiberSection`'s spatial (biaxial) counterpart: each fiber carries both
/// a `y` and a `z` coordinate, since bending can now occur about either
/// local section axis. Torsion is deliberately *not* modeled here — a
/// fiber's uniaxial material only sees axial strain, so torsional
/// stiffness needs a separate (elastic, `G*J`) term, kept on the element
/// (`DispBeamColumn3`/`ForceBeamColumn3`) rather than the section, the same
/// decoupled treatment `ElasticBeamColumn3` already uses. This mirrors
/// Xara's `FiberSection3d`, which likewise keeps torsion out of the fiber
/// loop and folds it in afterward via a separate elastic term.
#[derive(Debug, Clone)]
pub struct Fiber3 {
    pub y: f64,
    pub z: f64,
    pub area: f64,
    pub material: Material,
}

impl Fiber3 {
    pub fn new(y: f64, z: f64, area: f64, material: Material) -> Self {
        Fiber3 { y, z, area, material }
    }
}

/// Biaxial fiber-discretized cross-section: given axial strain `eps0`,
/// curvature about local z `kappa_z` (bending that displaces local y, the
/// plane `FiberSection`/`ElasticBeamColumn3`'s `Iz` governs), and curvature
/// about local y `kappa_y` (displaces local z, `Iy`'s plane), each fiber's
/// strain is `eps0 - y*kappa_z + z*kappa_y` — the natural biaxial extension
/// of `FiberSection`'s `eps0 - y*kappa`, verified against Xara's
/// `FiberSection3d::addFiber`/`getStressResultant` (`matData`/`sData`
/// layout: `strain = e0 - y*e1 + z*e2`, `N = sum(stress*A)`,
/// `Mz = sum(-y*stress*A)`, `My = sum(z*stress*A)`) rather than re-derived
/// from scratch.
///
/// Work-conjugate to `[eps0, kappa_z, kappa_y]` is `[N, Mz, My]`, and the
/// tangent's off-diagonal terms follow directly from differentiating those
/// three sums (see [`FiberSection3::trial`]'s body) — independently
/// reproducing Xara's `kData` layout (`kData[1]=-y*EA` summed → `-EQz`,
/// `kData[2]=z*EA` summed → `EQy`, `kData[6]=-y*z*EA` summed → `-EIyz`)
/// gives useful cross-confirmation that the sign convention here matches
/// the reference exactly, not just plausibly.
#[derive(Debug, Clone)]
pub struct FiberSection3 {
    fibers: Vec<Fiber3>,
}

/// `(N, Mz, My)` and the section tangent, ordered
/// `[[dN/deps0, dN/dkz, dN/dky], [dMz/deps0, dMz/dkz, dMz/dky], [dMy/deps0, dMy/dkz, dMy/dky]]`.
pub type SectionResponse3 = (f64, f64, f64, [[f64; 3]; 3]);

impl FiberSection3 {
    pub fn new(fibers: Vec<Fiber3>) -> Self {
        FiberSection3 { fibers }
    }

    /// Trial evaluation at `(eps0, kappa_z, kappa_y)` — read-only, safe
    /// every Newton iteration (see `FiberSection::trial`'s doc comment for
    /// why).
    pub fn trial(&self, eps0: f64, kappa_z: f64, kappa_y: f64) -> SectionResponse3 {
        let (mut n, mut mz, mut my) = (0.0, 0.0, 0.0);
        let (mut ea, mut eqz, mut eqy, mut eizz, mut eiyy, mut eiyz) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for fiber in &self.fibers {
            let strain = eps0 - fiber.y * kappa_z + fiber.z * kappa_y;
            let (stress, tangent) = fiber.material.trial_stress_tangent(strain);
            let fs = stress * fiber.area;
            n += fs;
            mz += -fiber.y * fs;
            my += fiber.z * fs;

            let ta = tangent * fiber.area;
            ea += ta;
            eqz += ta * fiber.y;
            eqy += ta * fiber.z;
            eizz += ta * fiber.y * fiber.y;
            eiyy += ta * fiber.z * fiber.z;
            eiyz += ta * fiber.y * fiber.z;
        }
        (n, mz, my, [[ea, -eqz, eqy], [-eqz, eizz, -eiyz], [eqy, -eiyz, eiyy]])
    }

    /// Commit every fiber's material at `(eps0, kappa_z, kappa_y)` — called
    /// once per converged step.
    pub fn commit(&mut self, eps0: f64, kappa_z: f64, kappa_y: f64) {
        for fiber in &mut self.fibers {
            let strain = eps0 - fiber.y * kappa_z + fiber.z * kappa_y;
            fiber.material = fiber.material.commit(strain);
        }
    }

    pub fn total_area(&self) -> f64 {
        self.fibers.iter().map(|f| f.area).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_symmetric_elastic_fibers_reproduce_ea_ei_exactly() {
        let (e, area, iz): (f64, f64, f64) = (30000.0, 2.0, 1000.0);
        let h = (iz / area).sqrt();
        let section = FiberSection::new(vec![
            Fiber::new(h, area / 2.0, Material::Elastic { e }),
            Fiber::new(-h, area / 2.0, Material::Elastic { e }),
        ]);

        let (n, m, k) = section.trial(0.001, 0.0002);
        assert!((k[0][0] - e * area).abs() < 1e-6, "EA mismatch: {}", k[0][0]);
        assert!((k[1][1] - e * iz).abs() < 1e-6, "EI mismatch: {}", k[1][1]);
        assert!(k[0][1].abs() < 1e-9, "should be EQ=0 for a symmetric section");

        // N = EA*eps0 (kappa doesn't couple in for a symmetric section),
        // M = EI*kappa.
        assert!((n - e * area * 0.001).abs() < 1e-6);
        assert!((m - e * iz * 0.0002).abs() < 1e-6);
    }

    /// M7 stage-2 acceptance test (implementation-plan §6): a real
    /// reinforced-concrete-style section — one `Steel01` rebar layer plus
    /// two `Concrete01` layers — pushed past first yield/cracking and then
    /// partially unloaded. No closed form exists for a nonlinear composite
    /// section (unlike the elastic case above), so each state is checked
    /// against a hand computation: evaluate every fiber's `Material`
    /// directly at that fiber's own strain and sum `N = sum(stress*area)`,
    /// `M = -sum(stress*area*y)` — the same definition `FiberSection::trial`
    /// uses, but computed independently, fiber by fiber, outside it.
    #[test]
    fn steel_and_concrete_fiber_section_matches_hand_computed_n_m_through_yield_and_partial_unload() {
        let y_rebar = 5.0;
        let area_rebar = 0.4;
        let y_c1 = -2.0;
        let area_c1 = 3.0;
        let y_c2 = -4.0;
        let area_c2 = 3.0;

        let rebar = || Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);
        let concrete = || Material::concrete01(4.0, 0.002, 3.0, 0.006);

        let mut section = FiberSection::new(vec![
            Fiber::new(y_rebar, area_rebar, rebar()),
            Fiber::new(y_c1, area_c1, concrete()),
            Fiber::new(y_c2, area_c2, concrete()),
        ]);

        // Independent "by hand" copies of each fiber's material, committed
        // in lockstep with the section, to compute the expected (N, M) at
        // each state from scratch.
        let mut hand_rebar = rebar();
        let mut hand_c1 = concrete();
        let mut hand_c2 = concrete();

        let check = |section: &mut FiberSection,
                          hand_rebar: &mut Material,
                          hand_c1: &mut Material,
                          hand_c2: &mut Material,
                          eps0: f64,
                          kappa: f64,
                          label: &str| {
            let strain_rebar = eps0 - y_rebar * kappa;
            let strain_c1 = eps0 - y_c1 * kappa;
            let strain_c2 = eps0 - y_c2 * kappa;

            let (stress_rebar, _) = hand_rebar.trial_stress_tangent(strain_rebar);
            let (stress_c1, _) = hand_c1.trial_stress_tangent(strain_c1);
            let (stress_c2, _) = hand_c2.trial_stress_tangent(strain_c2);

            let expected_n = stress_rebar * area_rebar + stress_c1 * area_c1 + stress_c2 * area_c2;
            let expected_m = -(stress_rebar * area_rebar * y_rebar
                + stress_c1 * area_c1 * y_c1
                + stress_c2 * area_c2 * y_c2);

            let (n, m, _k) = section.trial(eps0, kappa);
            assert!((n - expected_n).abs() < 1e-9, "{label}: N mismatch, expected {expected_n}, got {n}");
            assert!((m - expected_m).abs() < 1e-9, "{label}: M mismatch, expected {expected_m}, got {m}");

            section.commit(eps0, kappa);
            *hand_rebar = hand_rebar.commit(strain_rebar);
            *hand_c1 = hand_c1.commit(strain_c1);
            *hand_c2 = hand_c2.commit(strain_c2);
        };

        // State 1: curvature alone (eps0=0) puts the rebar (y=5) into
        // strain +0.005 — well past yield (epsy = 60/29000 ≈ 0.00207) —
        // and the two concrete layers into -0.002 (right at the peak,
        // epsc0) and -0.004 (past the peak, on the descending branch).
        check(&mut section, &mut hand_rebar, &mut hand_c1, &mut hand_c2, 0.0, -0.001, "state 1: first yield/cracking");

        // State 2: push further — deeper into the steel hardening range and
        // the concrete's descending branch.
        check(&mut section, &mut hand_rebar, &mut hand_c1, &mut hand_c2, 0.0, -0.002, "state 2: deeper into the nonlinear range");

        // State 3: partial unload (smaller-magnitude curvature) — history
        // (steel's isotropic-hardening shift, concrete's degraded
        // unload-reload slope) must carry over correctly for the (N, M) to
        // still match the independently-committed hand computation.
        check(&mut section, &mut hand_rebar, &mut hand_c1, &mut hand_c2, 0.0, -0.0008, "state 3: partial unload");
    }

    /// `FiberSection3` biaxial analogue of `two_symmetric_elastic_fibers_
    /// reproduce_ea_ei_exactly`: four corner fibers at `(±hy, ±hz)`
    /// reproduce `EA`, `EIy`, `EIz` exactly and give zero product-of-inertia
    /// coupling (`dMz/dky` / `dMy/dkz` = 0) by symmetry — verifying the
    /// `eps0 - y*kappa_z + z*kappa_y` strain field and the `N`/`Mz`/`My`
    /// sums independently of any element that will later drive this
    /// section.
    #[test]
    fn four_corner_elastic_fibers_reproduce_ea_eiy_eiz_with_zero_cross_coupling() {
        let (e, area, iy, iz): (f64, f64, f64, f64) = (30000.0, 4.0, 500.0, 2000.0);
        let hz = (iy / area).sqrt();
        let hy = (iz / area).sqrt();
        let a4 = area / 4.0;
        let section = FiberSection3::new(vec![
            Fiber3::new(hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(hy, -hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, -hz, a4, Material::Elastic { e }),
        ]);

        let (n, mz, my, k) = section.trial(0.001, 0.0002, 0.0003);
        assert!((k[0][0] - e * area).abs() < 1e-6, "EA mismatch: {}", k[0][0]);
        assert!((k[1][1] - e * iz).abs() < 1e-6, "EIz mismatch: {}", k[1][1]);
        assert!((k[2][2] - e * iy).abs() < 1e-6, "EIy mismatch: {}", k[2][2]);
        assert!(k[0][1].abs() < 1e-9, "EQz should be 0 for a symmetric section");
        assert!(k[0][2].abs() < 1e-9, "EQy should be 0 for a symmetric section");
        assert!(k[1][2].abs() < 1e-9, "EIyz should be 0 for a symmetric section");

        assert!((n - e * area * 0.001).abs() < 1e-6);
        assert!((mz - e * iz * 0.0002).abs() < 1e-6);
        assert!((my - e * iy * 0.0003).abs() < 1e-6);
    }
}
