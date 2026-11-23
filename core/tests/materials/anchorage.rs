//! Pier-end anchorage stress and tangent checked against the specified
//! WoodWorks Shearwalls crushing curve and closed-form hold-down law; branch
//! values use 1e-14 relative tolerance and smooth tangents use central
//! differences at 1e-7.

use carapace_core::model::{AnchorageBearing, Material, MaterialError};
use carapace_core::testkit::{assert_material_tangent_matches_fd, FdConfig};

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * a.abs().max(b.abs()), "{a} != {b}");
}

fn crush_from_r(r: f64) -> f64 {
    if r <= 0.73 {
        0.5 * r / 0.73
    } else if r <= 1.0 {
        0.5 * (1.0 + (r - 0.73) / 0.27)
    } else {
        r.powi(3)
    }
}

/// Tension, slack, compression branches, the two crush-curve kinks and
/// inverse consistency follow the equations in `docs/floors.md`, with
/// relative tolerance 1e-14.
#[test]
fn branches_kinks_and_inverse_match_closed_forms() {
    let qr = 20_000.0;
    let m = Material::anchorage(12_000.0, 0.0063, AnchorageBearing::O86Crush { qr }).unwrap();
    assert_eq!(m.trial_stress_tangent(0.0), (-0.0, 1.46 * qr));
    assert_eq!(m.trial_stress_tangent(0.003), (0.0, 0.0));
    assert_eq!(m.trial_stress_tangent(0.0063), (0.0, 0.0));
    assert_eq!(m.trial_stress_tangent(0.0083), (24.0, 12_000.0));
    for r in [0.2, 0.73, 0.85, 1.0, 1.4, 3.0] {
        let c = crush_from_r(r);
        let (stress, tangent) = m.trial_stress_tangent(-c);
        near(-stress / qr, r, 1e-14);
        let expected = if c <= 0.5 {
            1.46 * qr
        } else if c <= 1.0 {
            0.54 * qr
        } else {
            qr / (3.0 * c.powf(2.0 / 3.0))
        };
        near(tangent, expected, 1e-14);
        let recovered_r = if c <= 0.5 {
            0.73 * c / 0.5
        } else if c <= 1.0 {
            0.73 + 0.27 * (c / 0.5 - 1.0)
        } else {
            c.cbrt()
        };
        near(qr * recovered_r, -stress, 1e-14);
    }
    // Both piecewise boundaries are continuous in force.
    near(m.trial_stress_tangent(-0.5).0, -0.73 * qr, 1e-14);
    near(m.trial_stress_tangent(-0.5 - 1e-12).0, -0.73 * qr, 1e-11);
    near(m.trial_stress_tangent(-1.0).0, -qr, 1e-14);
    near(m.trial_stress_tangent(-1.0 - 1e-12).0, -qr, 1e-11);
}

/// Smooth compression and engaged-tension tangents match central
/// differences of the material response, `rel` 1e-7.
#[test]
fn smooth_tangents_match_finite_difference() {
    let m = Material::anchorage(
        12_000.0,
        0.0063,
        AnchorageBearing::O86Crush { qr: 20_000.0 },
    )
    .unwrap();
    let no_kink = |x: f64| x == 0.0 || x == 0.5 || x == 1.0 || x == 0.0063;
    let config = FdConfig::default();
    for strain in [-0.2, -0.8, -2.0, 0.01] {
        assert!(assert_material_tangent_matches_fd(
            &m, strain, &config, no_kink
        ));
    }
}

/// Constructor validation rejects negative, non-positive and non-finite
/// inputs; the slack split reports shrinkage separately from elastic uplift.
#[test]
fn validation_and_slack_split() {
    let bearing = AnchorageBearing::Linear { k: 100.0 };
    for (tension, slack, b) in [
        (-1.0, 0.0, bearing),
        (1.0, -1.0, bearing),
        (f64::NAN, 0.0, bearing),
        (f64::INFINITY, 0.0, bearing),
        (1.0, f64::INFINITY, bearing),
        (1.0, 0.0, AnchorageBearing::Linear { k: 0.0 }),
        (1.0, 0.0, AnchorageBearing::Linear { k: f64::INFINITY }),
        (1.0, 0.0, AnchorageBearing::O86Crush { qr: -1.0 }),
        (1.0, 0.0, AnchorageBearing::O86Crush { qr: f64::NAN }),
        (1.0, 0.0, AnchorageBearing::O86Crush { qr: f64::INFINITY }),
    ] {
        assert_eq!(
            Material::anchorage(tension, slack, b).unwrap_err(),
            MaterialError::InvalidAnchorageProperties
        );
    }
    let m = Material::anchorage(1.0, 2.0, bearing).unwrap();
    assert_eq!(m.anchorage_split(3.0), Some((2.0, 1.0)));
    assert_eq!(m.anchorage_split(1.0), Some((1.0, 0.0)));
    assert_eq!(m.anchorage_split(-1.0), Some((0.0, -1.0)));
    assert_eq!(Material::Elastic { e: 1.0 }.anchorage_split(1.0), None);
}

/// Reversing a pier-end displacement closes the tension slack and switches
/// the opposite end into compression; the stateless law returns the same
/// response on a repeated cycle. Reference: the specified tension and
/// bearing equations, exact within floating-point roundoff.
#[test]
fn reversal_reopens_slack_without_material_memory() {
    let m = Material::anchorage(5000.0, 0.2, AnchorageBearing::Linear { k: 10_000.0 }).unwrap();
    let displacements = [0.3, -0.3, 0.1, 0.3];
    let responses: Vec<_> = displacements
        .iter()
        .map(|&u| m.trial_stress_tangent(u))
        .collect();
    near(responses[0].0, 500.0, 1e-14);
    assert_eq!(responses[0], responses[3]);
    assert_eq!(responses[1], (-3000.0, 10_000.0));
    assert_eq!(responses[2], (0.0, 0.0));
    // Opposite ends of a symmetric pair swap tension and compression under
    // reversal; returning to the first displacement restores the response.
    let ends_forward = (m.trial_stress_tangent(0.3), m.trial_stress_tangent(-0.3));
    let ends_reversed = (m.trial_stress_tangent(-0.3), m.trial_stress_tangent(0.3));
    assert_eq!(ends_forward.0, ends_reversed.1);
    assert_eq!(ends_forward.1, ends_reversed.0);
}
