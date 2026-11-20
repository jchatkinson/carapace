# CSFM independent references

Phase 0 of the CSFM port (`docs/csfm-plan.md` S4 C). These check the
*method*, against independent closed-form derivations fixed before any
Carapace CSFM result exists — not the port, which is checked instead against
the golden fixtures in `../fixtures/` (see the distinction the plan draws
between "parity" and "independent" references).

| File | Checks |
|---|---|
| `tension_chord.py` / `.md` | The Tension Chord Model average-strain relation `CsfmRebar`'s TCM branch must invert, `rel 1e-8`. |
| `cmm_panel.py` / `.md` | An orthogonally reinforced panel's equilibrium/compatibility solution (Kaufmann & Marti 1998), for the eventual `CsfmConcrete` + `CsfmRebar` + `EmbeddedBar` combination: equilibrium identities and the shared `kc2` formula to rounding, bar stress to a measured-and-margined band. |

Both are plain-stdlib Python (no numpy in this environment); run with
`python3 -I <file>.py`. Neither depends on Carapace or on calcs.app's
TypeScript — they derive from the Tension Chord Model and the Cracked
Membrane Model themselves, using only the law *parameters* (`fck`, `rho`,
`lambda`, ...) the plan's S1 describes.
