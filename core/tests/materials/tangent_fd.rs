//! The generic tangent-vs-central-difference checker
//! (`carapace_core::testkit::{assert_material_tangent_matches_fd,
//! assert_plane_material_tangent_matches_fd}`, `docs/csfm-plan.md` phase 0)
//! against existing nonlinear uniaxial materials and every elastic
//! `PlaneMaterial` variant. This is the exit criterion "the checker passes
//! on existing materials and elements" — CSFM's own laws (phase 1) will
//! reuse the same checker on their smooth branches.

use carapace_core::model::{Material, PlaneMaterial, PlaneVector};
use carapace_core::testkit::{
    assert_material_tangent_matches_fd, assert_plane_material_tangent_matches_fd, FdConfig,
};

/// `Steel01`'s bilinear kinematic-hardening law, checked on its elastic
/// branch and on its post-yield hardening branch, each away from the yield
/// kink at `strain = fy / e0`. Reference: central difference of the law's
/// own stress, `rel` 1e-7 (`FdConfig::default`).
#[test]
fn steel01_smooth_branches_match_central_difference() {
    let (fy, e0, b) = (250.0, 200_000.0, 0.02);
    let material = Material::steel01(fy, e0, b, 0.0, 1.0, 0.0, 1.0);
    let yield_strain = fy / e0;
    let kink = move |x: f64| (x.abs() - yield_strain).abs() < 1e-4 * yield_strain;
    let config = FdConfig::default();

    // Elastic branch.
    assert!(assert_material_tangent_matches_fd(
        &material, 0.0005, &config, kink
    ));
    // Hardening branch, well past yield.
    assert!(assert_material_tangent_matches_fd(
        &material, 0.004, &config, kink
    ));
}

/// `Concrete01`'s Kent-Scott-Park envelope, checked on its parabolic
/// ascending branch and on its linear descending branch, each away from
/// the kinks at `epsc0` (parabola/line) and `epscu` (line/flat). Reference:
/// central difference of the law's own stress, `rel` 1e-7.
#[test]
fn concrete01_smooth_branches_match_central_difference() {
    let (fpc, epsc0, fpcu, epscu) = (30.0, 0.002, 6.0, 0.006);
    let material = Material::concrete01(fpc, epsc0, fpcu, epscu);
    let kinks = [0.0, -epsc0, -epscu];
    let kink = move |x: f64| kinks.iter().any(|k| (x - k).abs() < 1e-4 * epscu);
    let config = FdConfig::default();

    // Ascending (parabolic) branch, strictly between 0 and epsc0.
    assert!(assert_material_tangent_matches_fd(
        &material, -0.0012, &config, kink
    ));
    // Descending (linear) branch, strictly between epsc0 and epscu.
    assert!(assert_material_tangent_matches_fd(
        &material, -0.004, &config, kink
    ));
}

/// Every elastic `PlaneMaterial` variant (`PlaneStress`, `PlaneStrain`,
/// `Orthotropic`, `ElasticMatrix`) is linear, so its tangent is a constant
/// `D` with no kink anywhere; the `kink` predicate is `|_| false`. This is
/// the only `PlaneMaterial` coverage possible today — `CsfmConcrete`
/// (phase 1) is the first nonlinear variant and will reuse this same
/// checker on its smooth branches, with a real `kink` predicate in
/// principal-strain space.
#[test]
fn elastic_plane_materials_match_central_difference() {
    let config = FdConfig::default();
    let strain = PlaneVector::new(0.0021, -0.0008, 0.0013);
    let no_kink = |_: &PlaneVector| false;

    let variants = [
        PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap(),
        PlaneMaterial::plane_strain(30_000.0, 0.2).unwrap(),
        PlaneMaterial::orthotropic(30_000.0, 20_000.0, 0.15, 8_000.0, 0.4).unwrap(),
        PlaneMaterial::elastic_matrix([30_000.0, 6_000.0, 0.0, 30_000.0, 0.0, 12_000.0]).unwrap(),
    ];
    for material in &variants {
        assert!(assert_plane_material_tangent_matches_fd(
            material, &strain, &config, no_kink
        ));
    }
}

/// The checker itself refuses a stencil that straddles a kink, rather than
/// silently averaging across it — the plan's named requirement. Checked
/// directly against `Steel01`'s yield kink: a point within the default
/// stencil width of `fy / e0` must make the checker return `false`.
#[test]
fn checker_refuses_a_stencil_straddling_a_kink() {
    let (fy, e0, b) = (250.0, 200_000.0, 0.02);
    let material = Material::steel01(fy, e0, b, 0.0, 1.0, 0.0, 1.0);
    let yield_strain = fy / e0;
    let kink = move |x: f64| (x.abs() - yield_strain).abs() < 1e-4 * yield_strain;
    let config = FdConfig::default();

    assert!(!assert_material_tangent_matches_fd(
        &material,
        yield_strain,
        &config,
        kink
    ));
}
