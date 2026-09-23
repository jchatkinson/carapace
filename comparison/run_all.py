"""Runs the full carapace-vs-openseespy material comparison: for every
material in materials.MATERIALS, drives the shared cyclic protocol
(protocol.txt) through both engines and overlays the result.

Usage: run_all.py  (from the comparison/ directory, using comparison/.venv)
"""

import subprocess
import sys
import pathlib

ROOT = pathlib.Path(__file__).parent.parent
COMPARISON = ROOT / "comparison"
PYTHON = COMPARISON / ".venv" / "bin" / "python"

sys.path.insert(0, str(COMPARISON))
from materials import MATERIALS  # noqa: E402


def run(cmd, **kwargs):
    print(f"$ {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, check=True, **kwargs)


def main():
    for material in MATERIALS:
        run(
            [
                "cargo", "run", "-q", "-p", "carapace-core", "--example", "cyclic_material",
                "--", material, str(COMPARISON / "out" / f"{material}_carapace.csv"),
            ],
            cwd=ROOT,
        )
        run([str(PYTHON), "run_openseespy.py", material, f"out/{material}_ospy.csv"], cwd=COMPARISON)
        run([str(PYTHON), "plot_compare.py", material], cwd=COMPARISON)


if __name__ == "__main__":
    main()
