"""Overlays carapace's and openseespy's cyclic stress-strain traces for one
material and reports the peak-stress relative error between them.

Usage: plot_compare.py <material>
"""

import csv
import sys
import pathlib

import matplotlib.pyplot as plt

OUT_DIR = pathlib.Path(__file__).parent / "out"


def load(path):
    with open(path) as f:
        rows = list(csv.reader(f))[1:]
    return [float(r[0]) for r in rows], [float(r[1]) for r in rows]


def main():
    material = sys.argv[1]
    strain_c, stress_c = load(OUT_DIR / f"{material}_carapace.csv")
    strain_o, stress_o = load(OUT_DIR / f"{material}_ospy.csv")

    fig, ax = plt.subplots(figsize=(6, 5))
    # carapace: a thick, translucent line underneath so it reads as a band
    # even where openseespy's line sits exactly on top of it.
    ax.plot(strain_c, stress_c, "-", color="tab:blue", linewidth=4.5, alpha=0.4, solid_capstyle="round", label="carapace")
    # openseespy: a thin, high-contrast dashed line with sparse markers on
    # top, so a perfect overlap still shows two distinct series rather than
    # one line swallowing the other.
    ax.plot(
        strain_o, stress_o, "--", color="black", linewidth=1.0,
        marker="o", markersize=3.5, markevery=max(1, len(strain_o) // 60),
        markerfacecolor="tab:orange", markeredgecolor="black", markeredgewidth=0.5,
        label="openseespy",
    )
    ax.axhline(0, color="0.8", linewidth=0.8, zorder=0)
    ax.axvline(0, color="0.8", linewidth=0.8, zorder=0)
    ax.set_xlabel("strain")
    ax.set_ylabel("stress")
    ax.set_title(material)
    ax.legend()
    fig.tight_layout()

    out_path = OUT_DIR / f"{material}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    peak_c = max(abs(s) for s in stress_c)
    peak_o = max(abs(s) for s in stress_o)
    rel_err = abs(peak_c - peak_o) / max(abs(peak_o), 1e-12)
    print(f"{material}: wrote {out_path}  peak|stress| carapace={peak_c:.4f} openseespy={peak_o:.4f} rel_err={rel_err:.4%}")


if __name__ == "__main__":
    main()
