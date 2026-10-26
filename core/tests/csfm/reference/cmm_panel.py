#!/usr/bin/env python3
"""Independent Cracked Membrane Model (CMM, Kaufmann & Marti 1998) reference for an
orthogonally reinforced panel, for CSFM's concrete + rebar laws together (plan
docs/csfm-plan.md S4 C). Written in phase 0, before CsfmConcrete/CsfmRebar/
EmbeddedBar exist in Carapace: this is the a priori target, not a fit to any
Carapace output.

Model: a single uniform membrane element, orthogonal reinforcement in x and y,
displacements prescribed so the strain [eps_x, eps_y, gamma_xy] is uniform.
Rotating-angle equilibrium (concrete principal stresses act at the strain's own
principal direction -- the standard CMM/compression-field assumption, and also
CSFM's: both use one isotropic rotating-principal-strain concrete law):

    sigma_x  = rho_x * sigma_sx + sigma_c_x
    sigma_y  = rho_y * sigma_sy + sigma_c_y
    tau_xy   =                    tau_c_xy

(sigma_c_x, sigma_c_y, tau_c_xy) is (sigma1, sigma2) at angle theta rotated back
to x-y by the ordinary stress-transformation (Mohr's circle), where
(eps1, eps2, theta) is the Mohr decomposition of the prescribed strain. Bar
strain in each direction is taken equal to the membrane's own normal strain
there (the standard smeared-reinforcement compatibility assumption shared by
CMM and CSFM), and bar stress is read off the Tension Chord Model relation in
tension_chord.py (shared, since TCM is exactly what CSFM's rebar law also
implements per the plan S1).

The concrete law (parabola-rectangle compression, Kaufmann kc2 softening, zero
tension) is stated directly from the plan's S1/S2.1 formulas, which are
Eurocode/Kaufmann 1998 constants, not calcs.app internals.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tension_chord import TcmParams, eps_m, sigma_of_eps_m  # noqa: E402


# --------------------------------------------------------------------- concrete

def ecm_of(fck: float) -> float:
    return 22000.0 * ((fck + 8.0) / 10.0) ** 0.3


def fctm_of(fck: float) -> float:
    return 0.3 * fck ** (2.0 / 3.0)


def eta_fc_of(fck: float) -> float:
    return min(1.0, (30.0 / fck) ** (1.0 / 3.0))


EPS_C2 = -0.002
EPS_CU2 = -0.0035


def parabola_rectangle(fck: float, eps: float) -> float:
    """Base (unsoftened) compression law: parabola to the peak -eta_fc*fck at
    eps_c2, then a plateau to eps_cu2. Zero for eps > 0 (this function is only
    ever evaluated on the compressive principal direction)."""
    fc_eff = eta_fc_of(fck) * fck
    if eps > 0.0:
        return 0.0
    if eps < EPS_CU2:
        raise ValueError(f"eps={eps} exceeds the crushing strain eps_cu2={EPS_CU2}")
    xi = eps / EPS_C2  # in [0, 1] on the parabola, >= 1 on the plateau
    if xi <= 1.0:
        return -fc_eff * (2.0 * xi - xi * xi)
    return -fc_eff


def kc2(eps1: float) -> float:
    """Kaufmann's compression-softening coefficient for the minor principal
    direction, active only when the coexisting major principal strain eps1 is
    tensile. This is the original Kaufmann & Marti 1998 formula -- it is shared
    between CMM and CSFM by construction (CSFM borrows it verbatim per the
    plan), so it is NOT part of the "banded" comparison below: it is exact
    agreement to rounding."""
    if eps1 <= 0.0:
        return 1.0
    return min(1.0, 1.0 / (0.8 + 140.0 * eps1))


def concrete_principal_stress(fck: float, eps1: float, eps2: float) -> tuple[float, float]:
    """(sigma1, sigma2): zero tension on the major direction, Kaufmann-softened
    compression on the minor direction. eps1 >= eps2 is assumed (eps1 the major)."""
    sigma1 = 0.0 if eps1 > 0.0 else parabola_rectangle(fck, eps1)
    sigma2 = kc2(eps1) * parabola_rectangle(fck, eps2)
    return sigma1, sigma2


# ---------------------------------------------------------------- Mohr circle

def principal(ex: float, ey: float, gxy: float) -> tuple[float, float, float]:
    """(eps1, eps2, theta): eps1 >= eps2, theta the CCW angle from x to the
    eps1 direction. Engineering shear convention (gxy = 2 * eps_xy)."""
    avg = (ex + ey) / 2.0
    radius = math.hypot((ex - ey) / 2.0, gxy / 2.0)
    theta = 0.5 * math.atan2(gxy, ex - ey) if (ex != ey or gxy != 0.0) else 0.0
    return avg + radius, avg - radius, theta


def from_principal(s1: float, s2: float, theta: float) -> tuple[float, float, float]:
    """Inverse of `principal` for a symmetric 2-tensor: (x, y, xy) at angle
    theta from principal values (s1, s2). Same transformation for stress and
    (engineering-shear) strain -- ordinary tensor rotation, not a CSFM-specific
    formula."""
    avg = (s1 + s2) / 2.0
    half_diff = (s1 - s2) / 2.0
    c2, s2a = math.cos(2.0 * theta), math.sin(2.0 * theta)
    return avg + half_diff * c2, avg - half_diff * c2, 2.0 * half_diff * s2a


def assert_round_trip(ex: float, ey: float, gxy: float, tol: float = 1e-12) -> None:
    """Mohr decomposition then its own inverse must reproduce the input
    exactly (to rounding) -- the 'equilibrium to rounding' identity bucket."""
    e1, e2, theta = principal(ex, ey, gxy)
    bx, by, bgxy = from_principal(e1, e2, theta)
    scale = max(1.0, abs(ex), abs(ey), abs(gxy))
    assert abs(bx - ex) < tol * scale, (bx, ex)
    assert abs(by - ey) < tol * scale, (by, ey)
    assert abs(bgxy - gxy) < tol * scale, (bgxy, gxy)


# --------------------------------------------------------------------- panel

@dataclass(frozen=True)
class PanelCase:
    name: str
    fck: float
    ex: float
    ey: float
    gxy: float
    bar_x: TcmParams
    bar_y: TcmParams


def bar_stress(bar: TcmParams, eps: float) -> float:
    """Stress in one reinforcement direction at membrane strain `eps`, equal
    to the bar's own average strain by the smeared compatibility assumption.
    TCM branch for eps on [eps_m(sigma_cr), eps_m(fy)]; bare elastic steel
    below that (not reached by either case below -- both are chosen to sit on
    the TCM branch, away from the first-crack transition the plan's "startup
    bridge" discussion flags as a separate, non-closed-form concern)."""
    lo, hi = eps_m(bar, bar.sigma_cr), eps_m(bar, bar.fy)
    if not (lo <= eps <= hi):
        raise ValueError(f"eps={eps} outside the TCM branch [{lo}, {hi}] for this case; "
                         f"choose a case strain that avoids the first-crack transition")
    return sigma_of_eps_m(bar, eps)


def solve_panel(case: PanelCase) -> dict:
    assert_round_trip(case.ex, case.ey, case.gxy)
    eps1, eps2, theta = principal(case.ex, case.ey, case.gxy)
    assert eps1 >= eps2 - 1e-15

    sigma1, sigma2 = concrete_principal_stress(case.fck, eps1, eps2)
    sigma_c_x, sigma_c_y, tau_c_xy = from_principal(sigma1, sigma2, theta)

    sigma_sx = bar_stress(case.bar_x, case.ex)
    sigma_sy = bar_stress(case.bar_y, case.ey)

    rho_x, rho_y = case.bar_x.rho, case.bar_y.rho
    sigma_x = rho_x * sigma_sx + sigma_c_x
    sigma_y = rho_y * sigma_sy + sigma_c_y
    tau_xy = tau_c_xy

    return dict(eps1=eps1, eps2=eps2, theta=theta, kc2=kc2(eps1),
                sigma1=sigma1, sigma2=sigma2,
                sigma_c_x=sigma_c_x, sigma_c_y=sigma_c_y, tau_c_xy=tau_c_xy,
                sigma_sx=sigma_sx, sigma_sy=sigma_sy,
                sigma_x=sigma_x, sigma_y=sigma_y, tau_xy=tau_xy)


# ------------------------------------------------- the "band": crack-spacing simplification
#
# CSFM evaluates tension stiffening for each embedded bar in isolation, using
# that bar's own 1D crack spacing sr_x = sr(rho_x) (tension_chord.py). A fuller
# 2D CMM treatment accounts for the crack running obliquely to both
# reinforcement directions by combining the two 1D spacings along the crack
# normal (the standard Vecchio-Collins-style combination, also used in
# Kaufmann's own biaxial crack-spacing work):
#
#     1 / sr_full(theta) = |cos(theta)| / sr_x + |sin(theta)| / sr_y
#
# This is the one place CSFM's simplification (plan S4 C: "its crack spacing
# ... [is] simplified") enters quantitatively. We bound the resulting bar
# stress band by re-evaluating each bar's tension-stiffening term with
# sr_full(theta) in place of its own 1D sr, and reporting the stress shift.

def sr_full(bar_x: TcmParams, bar_y: TcmParams, theta: float) -> float:
    c, s = abs(math.cos(theta)), abs(math.sin(theta))
    return 1.0 / (c / bar_x.sr + s / bar_y.sr)


def bar_stress_with_sr(bar: TcmParams, eps: float, sr_override: float) -> float:
    """Same TCM inversion as `bar_stress`, but with the crack spacing replaced
    by `sr_override` (used only to measure the simplification's effect, not as
    part of the CMM reference solution itself)."""
    patched = TcmParams(bar.es, bar.ecm, bar.fy, bar.ft, bar.fctm, bar.d, bar.rho, bar.lam)
    # eps_m(sigma0) = sigma0/Es - tau_b1*sr/(Es*d); solve for sigma0 at the
    # overridden sr by the same closed form, bypassing TcmParams.sr.
    lo, hi = patched.sigma_cr, patched.fy

    def eps_m_override(sigma0: float) -> float:
        return sigma0 / patched.es - patched.tau_b1 * sr_override / (patched.es * patched.d)

    flo, fhi = eps_m_override(lo) - eps, eps_m_override(hi) - eps
    if flo > 0 or fhi < 0:
        raise ValueError("eps outside the sr-overridden TCM branch")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = eps_m_override(mid) - eps
        if abs(f) < 1e-14 or (hi - lo) < 1e-14:
            return mid
        if f > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def crack_spacing_band(case: PanelCase) -> dict:
    sr_f = sr_full(case.bar_x, case.bar_y, abs(principal(case.ex, case.ey, case.gxy)[2]))
    sx_full = bar_stress_with_sr(case.bar_x, case.ex, sr_f)
    sy_full = bar_stress_with_sr(case.bar_y, case.ey, sr_f)
    sx_csfm = bar_stress(case.bar_x, case.ex)
    sy_csfm = bar_stress(case.bar_y, case.ey)
    band_x = abs(sx_full - sx_csfm) / max(1.0, abs(sx_csfm))
    band_y = abs(sy_full - sy_csfm) / max(1.0, abs(sy_csfm))
    return dict(sr_full=sr_f, sr_x=case.bar_x.sr, sr_y=case.bar_y.sr,
                sigma_sx_full=sx_full, sigma_sy_full=sy_full, band_x=band_x, band_y=band_y)


def main() -> None:
    fck = 30.0
    ecm, fctm = ecm_of(fck), fctm_of(fck)

    bar_a = TcmParams(es=200000.0, ecm=ecm, fy=500.0, ft=550.0, fctm=fctm, d=16.0, rho=0.01, lam=0.67)

    case_a = PanelCase(
        name="A: pure shear, symmetric reinforcement",
        fck=fck, ex=0.0016, ey=0.0016, gxy=0.0045,
        bar_x=bar_a, bar_y=bar_a,
    )

    bar_bx = TcmParams(es=200000.0, ecm=ecm, fy=500.0, ft=550.0, fctm=fctm, d=16.0, rho=0.012, lam=0.67)
    bar_by = TcmParams(es=200000.0, ecm=ecm, fy=500.0, ft=550.0, fctm=fctm, d=12.0, rho=0.008, lam=0.67)

    case_b = PanelCase(
        name="B: biaxial, asymmetric reinforcement",
        fck=fck, ex=0.0022, ey=0.0017, gxy=0.0045,
        bar_x=bar_bx, bar_y=bar_by,
    )

    for case in (case_a, case_b):
        print(f"\n--- {case.name} ---")
        print(f"eps_x={case.ex}, eps_y={case.ey}, gamma_xy={case.gxy}")
        r = solve_panel(case)
        print(f"eps1={r['eps1']:.6e}  eps2={r['eps2']:.6e}  theta={math.degrees(r['theta']):.3f} deg  "
              f"kc2={r['kc2']:.6f}")
        print(f"sigma1={r['sigma1']:.6f} MPa  sigma2={r['sigma2']:.6f} MPa")
        print(f"sigma_c: x={r['sigma_c_x']:.6f}  y={r['sigma_c_y']:.6f}  xy={r['tau_c_xy']:.6f} MPa")
        print(f"sigma_sx={r['sigma_sx']:.4f} MPa  sigma_sy={r['sigma_sy']:.4f} MPa")
        print(f"panel: sigma_x={r['sigma_x']:.6f}  sigma_y={r['sigma_y']:.6f}  tau_xy={r['tau_xy']:.6f} MPa")

        band = crack_spacing_band(case)
        print(f"crack spacing: sr_x={band['sr_x']:.4f}  sr_y={band['sr_y']:.4f}  "
              f"sr_full(theta)={band['sr_full']:.4f} mm")
        print(f"bar stress with full-CMM spacing: sx={band['sigma_sx_full']:.4f}  "
              f"sy={band['sigma_sy_full']:.4f} MPa")
        print(f"band (bar stress, CSFM 1D spacing vs full-CMM combined spacing): "
              f"x={band['band_x']*100:.3f}%  y={band['band_y']*100:.3f}%")


if __name__ == "__main__":
    main()
