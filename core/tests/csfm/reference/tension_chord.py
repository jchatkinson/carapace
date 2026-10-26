#!/usr/bin/env python3
"""Independent closed-form Tension Chord Model (TCM) reference for CSFM's rebar law.

Derivation (see tension_chord.md for the full writeup): a bar of diameter `d`
embedded in a concrete tie of effective reinforcement ratio `rho`, pulled so that
the stress at a crack is `sigma0`. Between cracks, rigid-plastic bond (constant
shear stress) gives a linear steel-stress profile; averaging it over half the
crack spacing gives the mean strain

    eps_m(sigma0) = sigma0 / Es - tau_b1 * sr / (Es * d)

valid for sigma0 in [sigma_cr, fy], where sigma_cr = fctm * (1 - rho) / rho is the
stress at first cracking and fy is yield. This is run BEFORE any Carapace result
exists: it is the a priori reference for CsfmRebar's inverse relation (plan
docs/csfm-plan.md, S4 C), not something fit to it.

No numpy (not installed in this environment); everything here is scalar closed form.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class TcmParams:
    es: float       # steel modulus, MPa
    ecm: float      # concrete modulus, MPa
    fy: float       # yield stress, MPa
    ft: float       # ultimate stress, MPa
    fctm: float     # concrete mean tensile strength, MPa
    d: float        # bar diameter, mm
    rho: float      # effective reinforcement ratio As / Ac,eff (0 < rho < 1)
    lam: float = 0.67  # crack-spacing calibration factor (CSFM UI range 0.5-1.0)

    @property
    def tau_b0(self) -> float:
        return 2.0 * self.fctm  # uncracked bond strength

    @property
    def tau_b1(self) -> float:
        return self.fctm  # cracked (stabilized) bond strength

    @property
    def sr(self) -> float:
        """Stabilized crack spacing from first-crack bond equilibrium (tau_b0),
        with calibration factor lambda. Derivation: at first cracking, bond
        transfers force from steel to concrete at rate tau_b0 * (perimeter/d-area)
        = tau_b0 * 4/d per unit stress; a new crack forms sr/2 away once concrete
        regains fctm, i.e. fctm = rho/(1-rho) * (2 * tau_b0 * sr0 / d), giving
        sr0 = d * fctm * (1 - rho) / (2 * rho * tau_b0). lambda >= 1 is a
        calibration margin over this theoretical minimum spacing."""
        return self.lam * self.d * self.fctm * (1.0 - self.rho) / (2.0 * self.tau_b0 * self.rho)

    @property
    def sigma_cr(self) -> float:
        """Steel stress at a crack the instant it forms. Just before cracking the
        bar and the concrete share one strain eps_cr = fctm/Ecm (uncracked,
        bonded, both elastic), so the bar carries Es*fctm/Ecm and the concrete
        carries fctm over its net area (1-rho)*A. Total axial force is preserved
        across cracking, but immediately after the crack the concrete there
        carries zero, so it all transfers to the bar:

            A*[rho*Es*(fctm/Ecm) + (1-rho)*fctm] = A*rho*sigma_cr
            sigma_cr = fctm/rho + fctm*(Es/Ecm - 1)

        This is the sigma_cr(rho) whose fy-level inverse is exactly the
        TCM/POM routing boundary rho_crit below (checked in self_checks)."""
        return self.fctm / self.rho + self.fctm * (self.es / self.ecm - 1.0)

    def rho_crit(self) -> float:
        """TCM-vs-bare-steel routing boundary (plan S1): below this reinforcement
        ratio the bar yields before the concrete can crack again (stress at first
        crack already exceeds yield under TCM's own bookkeeping), so routing falls
        back to the Pull-Out Model. fy - (Es/Ecm - 1) * fctm is the stress the bar
        would carry at yield-strain compatibility with the concrete just before
        cracking; fctm / that quantity is the rho at which sigma_cr == fy."""
        return self.fctm / (self.fy - (self.es / self.ecm - 1.0) * self.fctm)


def eps_m(p: TcmParams, sigma0: float) -> float:
    """Mean (average) steel strain along the tie at crack stress `sigma0`,
    valid for sigma0 in [p.sigma_cr, p.fy] (the 'stabilized cracking' branch)."""
    if not (p.sigma_cr - 1e-9 <= sigma0 <= p.fy + 1e-9):
        raise ValueError(f"sigma0={sigma0} outside stabilized-cracking branch "
                          f"[{p.sigma_cr}, {p.fy}]")
    return sigma0 / p.es - p.tau_b1 * p.sr / (p.es * p.d)


def sigma_of_eps_m(p: TcmParams, eps: float, tol: float = 1e-14) -> float:
    """Inverse of eps_m by bisection (eps_m is affine/monotone increasing in
    sigma0 on the stabilized branch, so the inverse is unique there)."""
    lo, hi = p.sigma_cr, p.fy
    flo, fhi = eps_m(p, lo) - eps, eps_m(p, hi) - eps
    if flo > 0 or fhi < 0:
        raise ValueError(f"eps={eps} outside the range reachable on the stabilized branch "
                         f"[{eps_m(p, lo)}, {eps_m(p, hi)}]")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fmid = eps_m(p, mid) - eps
        if abs(fmid) < tol or (hi - lo) < tol:
            return mid
        if fmid > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def self_checks(p: TcmParams) -> None:
    """Internal consistency of the closed form, independent of any FE result."""
    # 1. eps_m is affine in sigma0 with slope exactly 1/Es (bond term is constant).
    s1, s2 = p.sigma_cr + 5.0, p.sigma_cr + 55.0
    slope = (eps_m(p, s2) - eps_m(p, s1)) / (s2 - s1)
    assert abs(slope - 1.0 / p.es) < 1e-12 * (1.0 / p.es), "eps_m slope must be exactly 1/Es"

    # 2. At sigma0 = sigma_cr, eps_m must be LESS than the bare-steel strain
    # sigma_cr/Es (tension stiffening stiffens the response: same crack stress,
    # smaller mean strain, because the concrete between cracks is still
    # carrying some of the load).
    bare = p.sigma_cr / p.es
    assert eps_m(p, p.sigma_cr) < bare, "tension stiffening must reduce mean strain below bare steel"

    # 3. Forward-then-inverse round trip.
    for sigma0 in (p.sigma_cr, p.sigma_cr + 30.0, p.fy):
        e = eps_m(p, sigma0)
        back = sigma_of_eps_m(p, e)
        assert abs(back - sigma0) < 1e-8 * max(1.0, abs(sigma0)), (
            f"round trip failed at sigma0={sigma0}: got {back}")

    # 4. rho_crit is exactly the rho at which sigma_cr(rho) == fy, by construction
    # of both formulas (both express the same "bond-limited yield" condition).
    rc = p.rho_crit()
    p_at_rc = TcmParams(p.es, p.ecm, p.fy, p.ft, p.fctm, p.d, rc, p.lam)
    assert abs(p_at_rc.sigma_cr - p.fy) < 1e-9 * p.fy, "rho_crit definition is inconsistent with sigma_cr"


def main() -> None:
    # Representative parameters in the CSFM UI's range (CSFM_CONTEXT.md, calcs.app).
    fck = 30.0
    ecm = 22000.0 * ((fck + 8.0) / 10.0) ** 0.3
    fctm = 0.3 * fck ** (2.0 / 3.0)
    p = TcmParams(es=200000.0, ecm=ecm, fy=500.0, ft=550.0, fctm=fctm, d=16.0, rho=0.01, lam=0.67)

    print(f"Ecm = {ecm:.4f} MPa, fctm = {fctm:.4f} MPa")
    print(f"tau_b0 = {p.tau_b0:.4f} MPa, tau_b1 = {p.tau_b1:.4f} MPa, sr = {p.sr:.4f} mm")
    print(f"sigma_cr = {p.sigma_cr:.4f} MPa, rho_crit = {p.rho_crit():.6f} (rho = {p.rho})")
    assert p.rho > p.rho_crit(), "this example must sit on the TCM branch, not POM"

    self_checks(p)
    print("self-checks passed")

    for sigma0 in (p.sigma_cr, p.sigma_cr + 50.0, 400.0, p.fy):
        e = eps_m(p, sigma0)
        back = sigma_of_eps_m(p, e)
        print(f"sigma0={sigma0:9.4f} MPa  ->  eps_m={e:.8e}  ->  sigma(eps_m)={back:9.4f} MPa  "
              f"(residual {back - sigma0:.3e} MPa)")


if __name__ == "__main__":
    main()
