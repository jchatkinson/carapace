//! Plan `docs/csfm-plan.md` §4 C: `CsfmRebar`'s TCM branch and `CsfmConcrete` against the
//! independent references written in phase 0 (`core/tests/csfm/reference/tension_chord.{py,md}`,
//! `cmm_panel.{py,md}`) — the Tension Chord Model and Cracked Membrane Model themselves, not
//! calcs.app's TypeScript. The closed-form formulas below are transcribed from those documents
//! independently of `core/src/model/materials/csfm_rebar.rs`/`csfm_concrete.rs`'s own
//! implementation, so agreement is a real check, not a tautology.
//!
//! This uses the bare `Material`/`PlaneMaterial::trial_stress_tangent` calls directly rather than
//! building a `Domain`/`EmbeddedBar`: for a *uniform* strain state (what both references assume),
//! `EmbeddedBar::local_response` reduces to exactly this call (its own strain, under an affine host
//! displacement field, equals the panel's imposed strain exactly) - a full FE model here would
//! exercise `embed_polyline`'s geometry, which §4 D already covers, not anything new about §4 C.

use carapace_core::model::{CsfmConcreteLaw, Material, PlaneMaterial, PlaneVector};

fn ecm_of(fck: f64) -> f64 {
    22000.0 * ((fck + 8.0) / 10.0).powf(0.3)
}

fn fctm_of(fck: f64) -> f64 {
    0.3 * fck.powf(2.0 / 3.0)
}

/// `|a - b| <= abs + rel * max(|a|, |b|)`.
fn close(a: f64, b: f64, abs: f64, rel: f64) -> bool {
    (a - b).abs() <= abs + rel * a.abs().max(b.abs())
}

/// `tension_chord.md`'s closed-form mean strain, transcribed independently of
/// `csfm_rebar.rs::forward_strain`'s TCM branch (same formula, different source, by design — see
/// module doc). The intercept `sigma_min = tau_b0 * sr / d` uses the *uncracked* bond `tau_b0`
/// (the crack-formation stage), not `tau_b1` — see `tension_chord.md`'s "Mean strain" section for
/// why an earlier version of this reference used the wrong one.
fn reference_eps_m(es: f64, fctm: f64, d: f64, rho: f64, lambda: f64, sigma0: f64) -> f64 {
    let tau_b0 = 2.0 * fctm;
    let sr = lambda * d * fctm * (1.0 - rho) / (2.0 * tau_b0 * rho);
    let sigma_min = tau_b0 * sr / d;
    (sigma0 - sigma_min) / es
}

/// §4 C, tension chord: `CsfmRebar`'s TCM inverse must reproduce `tension_chord.md`'s closed-form
/// mean strain at `sigma0 = sigma_cr`, `sigma_cr + 50`, `400`, and `fy` (the same points
/// `tension_chord.py::main` prints), to `rel 1e-8` on the round trip.
#[test]
fn csfm_rebar_tcm_branch_matches_the_tension_chord_reference() {
    let (es, fy, ft, eps_u) = (200_000.0, 500.0, 550.0, 0.05);
    let fck = 30.0;
    let (ecm, fctm) = (ecm_of(fck), fctm_of(fck));
    let (d, rho, lambda) = (16.0, 0.01, 0.67);

    let rho_crit = fctm / (fy - (es / ecm - 1.0) * fctm);
    assert!(
        rho > rho_crit,
        "this example must sit on the TCM branch, not POM"
    );

    let rebar = Material::csfm_rebar(es, fy, ft, eps_u, fctm, ecm, d, rho, lambda, true)
        .expect("valid TCM rebar parameters");

    let sigma_cr = fctm / rho + fctm * (es / ecm - 1.0);
    for &sigma0 in &[sigma_cr, sigma_cr + 50.0, 400.0, fy] {
        let eps = reference_eps_m(es, fctm, d, rho, lambda, sigma0);
        let (stress, _tangent) = rebar.trial_stress_tangent(eps);
        assert!(
            close(stress, sigma0, 0.0, 1e-8),
            "sigma0={sigma0}: CsfmRebar gives {stress} at eps_m={eps} (reference round trip)"
        );
    }
}

struct PanelCase {
    rho_x: f64,
    d_x: f64,
    rho_y: f64,
    d_y: f64,
    eps_x: f64,
    eps_y: f64,
    gamma_xy: f64,
    // Expected (sigma_c_x, sigma_c_y, tau_c_xy): concrete-only, to rounding.
    concrete: (f64, f64, f64),
    // Expected (sigma_sx, sigma_sy): bar stress, within the crack-spacing band (cmm_panel.md).
    bar: (f64, f64),
}

/// `cmm_panel.md`'s cases A and B (full precision, not the `.md` table's 6-decimal display
/// values), re-derived by running `python3 -I core/tests/csfm/reference/cmm_panel.py` — not
/// hand-copied from any Carapace or calcs.app output.
fn panel_cases() -> [PanelCase; 2] {
    [
        PanelCase {
            rho_x: 0.01,
            d_x: 16.0,
            rho_y: 0.01,
            d_y: 16.0,
            eps_x: 0.0016,
            eps_y: 0.0016,
            gamma_xy: 0.0045,
            concrete: (-6.098300970873784, -6.098300970873784, 6.098300970873784),
            bar: (416.0613663196506, 416.0613663196506),
        },
        PanelCase {
            rho_x: 0.012,
            d_x: 16.0,
            rho_y: 0.008,
            d_y: 12.0,
            // 0.0020, not 0.0022: at 0.0022 it falls outside bar_x's TCM branch under the
            // corrected sigma_min-anchored formula (cmm_panel.md's case B note).
            eps_x: 0.0020,
            eps_y: 0.0017,
            gamma_xy: 0.0045,
            concrete: (-3.707470477687189, -4.235852459756758, 3.962864865521766),
            bar: (479.8894191293998, 460.31928710826674),
        },
    ]
}

/// §4 C, CMM panel: `CsfmConcrete`'s Cartesian stress at the panel's prescribed uniform strain
/// matches the reference's concrete-only contribution to rounding (`abs 1e-9` MPa, `rel 1e-10` -
/// shared `kc2`/zero-tension formulas, no modeling slack), and `CsfmRebar`'s stress at each
/// direction's own normal strain matches within the `rel 10%` crack-spacing band `cmm_panel.md`
/// measures (CSFM evaluates each bar's tension stiffening from its own 1D crack spacing, not the
/// angle-corrected 2D spacing CMM nominally uses).
#[test]
fn csfm_concrete_and_rebar_match_the_cmm_panel_reference() {
    let fck = 30.0;
    let (ecm, fctm) = (ecm_of(fck), fctm_of(fck));
    let (es, fy, ft, eps_u, lambda) = (200_000.0, 500.0, 550.0, 0.05, 0.67);

    let concrete =
        PlaneMaterial::csfm_concrete(fck, None, CsfmConcreteLaw::ParabolaRectangle, true, None)
            .expect("valid concrete parameters");

    for case in panel_cases() {
        let strain = PlaneVector::new(case.eps_x, case.eps_y, case.gamma_xy);
        let (stress, _tangent) = concrete.trial_stress_tangent(&strain);
        let (want_x, want_y, want_xy) = case.concrete;
        assert!(
            close(stress[0], want_x, 1e-9, 1e-10),
            "sigma_c_x: got {}, reference {want_x}",
            stress[0]
        );
        assert!(
            close(stress[1], want_y, 1e-9, 1e-10),
            "sigma_c_y: got {}, reference {want_y}",
            stress[1]
        );
        assert!(
            close(stress[2], want_xy, 1e-9, 1e-10),
            "tau_c_xy: got {}, reference {want_xy}",
            stress[2]
        );

        let rebar_x = Material::csfm_rebar(
            es, fy, ft, eps_u, fctm, ecm, case.d_x, case.rho_x, lambda, true,
        )
        .expect("valid rebar parameters");
        let rebar_y = Material::csfm_rebar(
            es, fy, ft, eps_u, fctm, ecm, case.d_y, case.rho_y, lambda, true,
        )
        .expect("valid rebar parameters");
        let (sigma_sx, _) = rebar_x.trial_stress_tangent(case.eps_x);
        let (sigma_sy, _) = rebar_y.trial_stress_tangent(case.eps_y);
        let (want_sx, want_sy) = case.bar;
        assert!(
            close(sigma_sx, want_sx, 0.0, 0.10),
            "sigma_sx: got {sigma_sx}, reference {want_sx} (rel 10% crack-spacing band)"
        );
        assert!(
            close(sigma_sy, want_sy, 0.0, 0.10),
            "sigma_sy: got {sigma_sy}, reference {want_sy} (rel 10% crack-spacing band)"
        );
    }
}
