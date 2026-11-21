"""A small builder over the columnar ``CarapaceInputV1`` tables.

Nodes and elements are added one at a time and referred to by the integer ids the ``add_*`` calls return,
so the parallel-array bookkeeping of the wire format stays out of sight. Anything the builder does not cover
can be set on ``model.raw`` (a plain dict merged into the output last).
"""

from __future__ import annotations

from typing import Any, Sequence

_FREE, _FIXED = 0, 1


class Model:
    """``Model(ndm=2)`` for a 2D frame (dofs ux, uy, rz) or ``Model(ndm=3)`` (ux, uy, uz, rx, ry, rz)."""

    def __init__(self, ndm: int = 2):
        if ndm not in (2, 3):
            raise ValueError("ndm must be 2 or 3")
        self.ndm = ndm
        self.ndf = 3 if ndm == 2 else 6
        self._coords: list[float] = []
        self._fixed: list[int] = []
        self._materials: list[dict] = []
        self._trusses: dict[str, list] = {k: [] for k in ("nodeI", "nodeJ", "area", "material", "density")}
        self._beams: dict[str, list] = {}
        self._loads: dict[str, list] = {k: [] for k in ("pattern", "node", "dof", "value", "stage")}
        self._stages: list[dict] = []
        self._recorders: list[dict] = []
        self.raw: dict[str, Any] = {}

    # -- geometry -----------------------------------------------------------------------------------
    def add_node(self, *coords: float, fixed: Sequence[int] | str = ()) -> int:
        """Add a node at ``(x, y[, z])``. ``fixed`` is a list of fixed dof indices (0-based) or ``"all"``."""
        if len(coords) != self.ndm:
            raise ValueError(f"expected {self.ndm} coordinates, got {len(coords)}")
        dofs = range(self.ndf) if fixed == "all" else fixed
        mask = 0
        for d in dofs:
            if not 0 <= d < self.ndf:
                raise ValueError(f"dof {d} out of range for ndm={self.ndm}")
            mask |= 1 << d
        self._coords.extend(float(c) for c in coords)
        self._fixed.append(mask)
        return len(self._fixed) - 1

    def add_elastic_material(self, e: float) -> int:
        self._materials.append({"kind": "elastic", "e": e})
        return len(self._materials) - 1

    def add_truss(self, i: int, j: int, area: float, material: int, density: float = 0.0) -> int:
        t = self._trusses
        for key, value in zip(("nodeI", "nodeJ", "area", "material", "density"), (i, j, area, material, density)):
            t[key].append(value)
        return len(t["nodeI"]) - 1

    def add_elastic_beam(self, i: int, j: int, e: float, a: float, iz: float, *, g: float | None = None,
                         j_torsion: float | None = None, iy: float | None = None,
                         vec_xz: Sequence[float] = (0.0, 0.0, 1.0), density: float = 0.0,
                         p_delta: bool = False) -> int:
        """Elastic beam-column. 2D needs ``e, a, iz``; 3D also ``g, j_torsion, iy`` and the local-plane ``vec_xz``."""
        if self.ndm == 2:
            row = {"nodeI": i, "nodeJ": j, "e": e, "a": a, "iz": iz,
                   "transform": "pDelta" if p_delta else "linear", "density": density}
        else:
            if g is None or j_torsion is None or iy is None:
                raise ValueError("3D beams need g, j_torsion and iy")
            kind = "pDelta3" if p_delta else "linear3"
            row = {"nodeI": i, "nodeJ": j, "e": e, "g": g, "a": a, "j": j_torsion, "iy": iy, "iz": iz,
                   "transform": {"kind": kind, "vecXz": list(vec_xz)}, "density": density}
        for key, value in row.items():
            self._beams.setdefault(key, []).append(value)
        return len(self._beams["nodeI"]) - 1

    # -- loads and analysis -----------------------------------------------------------------------
    def add_nodal_load(self, node: int, dof: int, value: float, stage: int = 0) -> None:
        """A reference load, ramped by the load factor of the static ``stage`` that registers it."""
        for key, v in zip(("pattern", "node", "dof", "value", "stage"), (0, node, dof, value, stage)):
            self._loads[key].append(v)

    def add_static_stage(self, stage_id: str = "static", steps: int = 1, increment: float = 1.0,
                         algorithm: str = "linear", tol: float = 1e-8, max_iter: int = 20) -> int:
        """A load-controlled static stage; ``algorithm`` is ``"linear"`` or ``"newton"``."""
        self._stages.append({
            "kind": "static", "id": stage_id, "steps": steps,
            "integrator": {"kind": "loadControl", "increment": increment},
            "algorithm": {"kind": algorithm},
            "convergence": {"kind": "normUnbalance", "tol": tol, "maxIter": max_iter},
            "holdPatternsAfter": [],
        })
        return len(self._stages) - 1

    def add_modal_stage(self, stage_id: str = "modal", modes: int = 3) -> int:
        self._stages.append({"kind": "modal", "id": stage_id, "modes": modes})
        return len(self._stages) - 1

    def record_node_disp(self, node: int, dof: int) -> int:
        """Record a node displacement; returns the recorder index to read results with."""
        self._recorders.append({"response": "nodeDisp", "node": node, "dof": dof})
        return len(self._recorders) - 1

    # -- output -----------------------------------------------------------------------------------
    def to_dict(self) -> dict:
        d: dict[str, Any] = {
            "header": {"schemaVersion": 1, "ndm": self.ndm, "engineVersion": "python"},
            "nodes": {"coords": self._coords, "fixed": self._fixed, "massNodeIndex": [], "mass": []},
            "sequence": {"stages": self._stages, "recorders": self._recorders},
        }
        if self._materials:
            d["materials"] = self._materials
        if self._trusses["nodeI"]:
            d["trusses"] = self._trusses
        if self._beams:
            d["elasticBeamColumns2d" if self.ndm == 2 else "elasticBeamColumns3d"] = self._beams
        if self._loads["node"]:
            d["loadPatterns"] = {"series": [{"kind": "linear", "slope": 1.0}], "scaleFactor": [1.0]}
            d["nodalLoads"] = self._loads
        d.update(self.raw)
        return d
