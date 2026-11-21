"""An OpenSeesPy-style model: the same command names, argument order and 1-based tags and DOFs.

``Model`` collects commands (``node``, ``element``, ``load``, ``analyze``, ...) into the columnar
``CarapaceInputV1`` tables and steps the engine for ``analyze``. Where it differs from OpenSeesPy it says so
loudly: an unsupported command or option raises ``NotImplementedError`` naming it, it is never ignored. The
deviations are listed in ``docs/python-api.md``.
"""

from __future__ import annotations

import json
from typing import Any

from . import _carapace
from .errors import ModelError

# Accepted and ignored: Carapace has one numbering and one sparse direct solver (see docs/python-api.md).
_IGNORED = {
    "numberer": {"RCM", "Plain", "AMD"},
    "system": {"BandGeneral", "BandSPD", "ProfileSPD", "SparseGeneral", "SparseSPD", "SparseSYM", "UmfPack",
               "FullGeneral", "SuperLU", "Mumps"},
}
_TESTS = {"NormUnbalance": "normUnbalance", "NormDispIncr": "normDispIncr", "EnergyIncr": "energyIncr"}
_BIG = 2**31  # a stage is stepped by `analyze`, so its nominal length is just "long enough"

# OpenSeesPy commands Carapace does not implement. Calling one raises instead of silently doing nothing.
_UNSUPPORTED = {
    "beamIntegration", "block2D", "block3D", "damping", "database", "domainChange", "equalDOF", "fixX", "fixY",
    "fixZ", "frictionModel", "getEleClassTags", "groundMotion", "imposedMotion", "layer", "load_const",
    "loadConst", "mp", "nDMaterial", "patch", "printModel", "randomVariable", "recorder", "region", "reset",
    "restore", "rayleigh", "remove", "rigidDiaphragm", "rigidLink", "save", "section", "sdfResponse", "sp",
    "wipeAnalysis", "setTime", "sensitivityAlgorithm", "modalDamping", "eleResponse", "eleForce", "basicForce",
    "nodeVel", "nodeAccel", "nodeBounds", "nodeMass", "updateElementDomain", "setNodeDisp", "setNodeVel",
    "setNodeAccel", "sectionForce", "sectionDeformation", "sectionStiffness", "sectionFlexibility", "analysisNames",
}


def _unsupported(what: str) -> NotImplementedError:
    return NotImplementedError(f"{what} is not supported by carapace (see docs/python-api.md)")


class Model:
    """One model plus its analysis, ``Model(ndm=2, ndf=3)`` or ``Model(ndm=3, ndf=6)``.

    Unlike OpenSeesPy there is no hidden global: create as many as you like, and use ``carapace.opensees`` for
    a module-level instance with the plain ``ops.node(...)`` spelling.
    """

    def __init__(self, ndm: int = 2, ndf: int | None = None, *, reactions: bool = True):
        self._configure(ndm, ndf)
        self._record_reactions = reactions

    # -- state ---------------------------------------------------------------------------------------
    def _configure(self, ndm: int, ndf: int | None) -> None:
        if ndm not in (2, 3) or ndf not in (None, 3 if ndm == 2 else 6):
            raise _unsupported(f"model with ndm={ndm}, ndf={ndf} (only ndm=2/ndf=3 and ndm=3/ndf=6)")
        self.ndm, self.ndf = ndm, 3 if ndm == 2 else 6
        self._clear()

    def _clear(self) -> None:
        self._tags: dict[int, int] = {}  # node tag -> table index
        self._coords: list[float] = []
        self._fixed: list[int] = []
        self._mass: dict[int, list[float]] = {}
        self._materials: list[dict] = []
        self._material_tags: dict[int, int] = {}
        self._transforms: dict[int, dict] = {}
        self._tables: dict[str, dict[str, list]] = {}
        self._element_tags: dict[int, tuple[str, int]] = {}
        self._series: dict[int, tuple[dict, float]] = {}
        self._patterns: list[tuple[dict, float]] = []
        self._pattern_tags: dict[int, int] = {}
        self._current_pattern: int | None = None
        self._loads: dict[str, list] = {k: [] for k in ("pattern", "node", "dof", "value", "stage")}
        self._test = ("NormUnbalance", 1e-8, 20)
        self._algorithm: dict = {"kind": "linear"}
        self._integrator: dict | None = None
        self._analysis: str | None = None
        self._session: Any = None
        self._disp: list[float] = []
        self._react: list[float] = []
        self._time = 0.0
        self._modes: list[dict] | None = None
        self.raw: dict[str, Any] = {}
        """Escape hatch: a dict merged over the generated ``CarapaceInputV1`` document."""
        self.last_error: dict | None = None
        """The analysis error behind the last negative ``analyze`` return, as a tagged dict."""

    def _open(self, command: str) -> None:
        if self._session is not None:
            raise RuntimeError(
                f"{command}() after analyze(): a Carapace model is fixed once analysis starts "
                "(call wipe() to start over)")

    def __getattr__(self, name: str):
        if name in _UNSUPPORTED:
            raise _unsupported(f"{name}()")
        raise AttributeError(name)

    # -- model ---------------------------------------------------------------------------------------
    def model(self, *args) -> None:
        """``model('basic', '-ndm', n, '-ndf', f)``; discards everything defined so far."""
        if not args or args[0] != "basic":
            raise _unsupported(f"model type {args[:1]}")
        opts = dict(zip(args[1::2], args[2::2]))
        if set(opts) - {"-ndm", "-ndf"} or "-ndm" not in opts:
            raise ValueError("model('basic', '-ndm', n [, '-ndf', f])")
        self._configure(int(opts["-ndm"]), opts.get("-ndf") and int(opts["-ndf"]))

    def wipe(self) -> None:
        """Remove the model and analysis (``ndm``/``ndf`` stay as constructed)."""
        self._clear()

    def node(self, tag: int, *args) -> None:
        """``node(tag, x, y[, z] [, '-mass', m1, ...])``."""
        self._open("node")
        args = list(args)
        mass = None
        if "-mass" in args:
            at = args.index("-mass")
            mass, args = args[at + 1:], args[:at]
        if any(isinstance(a, str) for a in args):
            raise _unsupported(f"node option {[a for a in args if isinstance(a, str)][0]!r}")
        if len(args) != self.ndm:
            raise ValueError(f"node {tag}: expected {self.ndm} coordinates, got {len(args)}")
        if tag in self._tags:
            raise ValueError(f"node {tag} is already defined")
        self._tags[tag] = len(self._fixed)
        self._coords.extend(float(c) for c in args)
        self._fixed.append(0)
        if mass is not None:
            self.mass(tag, *mass)

    def mass(self, tag: int, *values: float) -> None:
        """``mass(tag, m1, ..., m_ndf)``."""
        self._open("mass")
        if len(values) != self.ndf:
            raise ValueError(f"mass at node {tag}: expected {self.ndf} values, got {len(values)}")
        self._mass[self._node(tag)] = [float(v) for v in values]

    def fix(self, tag: int, *flags: int) -> None:
        """``fix(tag, f1, ..., f_ndf)``; 1 fixes a DOF, 0 leaves it free (a 0 does not free a fixed DOF)."""
        self._open("fix")
        if len(flags) != self.ndf:
            raise ValueError(f"fix at node {tag}: expected {self.ndf} flags, got {len(flags)}")
        at = self._node(tag)
        for dof, flag in enumerate(flags):
            if flag:
                self._fixed[at] |= 1 << dof

    def _node(self, tag: int) -> int:
        try:
            return self._tags[tag]
        except KeyError:
            raise ValueError(f"node {tag} is not defined") from None

    def uniaxialMaterial(self, kind: str, tag: int, *args: float) -> None:
        """``uniaxialMaterial('Elastic', tag, E)``."""
        self._open("uniaxialMaterial")
        if kind != "Elastic":
            raise _unsupported(f"uniaxialMaterial {kind!r}")
        if len(args) != 1:
            raise _unsupported("Elastic with eta or Eneg (only uniaxialMaterial('Elastic', tag, E))")
        self._material_tags[tag] = len(self._materials)
        self._materials.append({"kind": "elastic", "e": float(args[0])})

    def geomTransf(self, kind: str, tag: int, *vecxz: float) -> None:
        """2D: ``geomTransf('Linear'|'PDelta'|'Corotational', tag)``; 3D: ``geomTransf('Linear'|'PDelta', tag, *vecxz)``."""
        self._open("geomTransf")
        if self.ndm == 2:
            names = {"Linear": "linear", "PDelta": "pDelta", "Corotational": "corotational"}
            if kind not in names or vecxz:
                raise _unsupported(f"2D geomTransf {kind!r} with arguments {vecxz}")
            self._transforms[tag] = names[kind]
        else:
            names = {"Linear": "linear3", "PDelta": "pDelta3"}
            if kind not in names or len(vecxz) != 3 or any(isinstance(v, str) for v in vecxz):
                raise _unsupported(f"3D geomTransf {kind!r} (only Linear/PDelta with vecxz x, y, z)")
            self._transforms[tag] = {"kind": names[kind], "vecXz": [float(v) for v in vecxz]}

    def element(self, kind: str, tag: int, *args) -> None:
        """``element('truss', tag, i, j, A, matTag [, '-rho', r])`` and ``element('elasticBeamColumn', tag, i, j,
        A, E, Iz, transfTag [, '-mass', m])`` (3D: ``A, E, G, J, Iy, Iz, transfTag``)."""
        self._open("element")
        if tag in self._element_tags:
            raise ValueError(f"element {tag} is already defined")
        args = list(args)
        per_length = 0.0  # OpenSees' -rho / -mass: mass per unit length; the wire wants mass per volume
        for flag in ("-rho", "-mass"):
            if flag in args:
                at = args.index(flag)
                per_length = float(args[at + 1])
                del args[at:at + 2]
        if any(isinstance(a, str) for a in args):
            raise _unsupported(f"{kind} option {[a for a in args if isinstance(a, str)][0]!r}")
        if kind == "truss":
            if len(args) != 4:
                raise ValueError("element('truss', tag, i, j, A, matTag)")
            i, j, area, mat = args
            density = self._volumetric(tag, per_length, area)
            if mat not in self._material_tags:
                raise ValueError(f"element {tag}: material {mat} is not defined")
            row = {"nodeI": self._node(i), "nodeJ": self._node(j), "area": float(area),
                   "material": self._material_tags[mat], "density": density}
            self._add_row("trusses", tag, row)
        elif kind == "elasticBeamColumn":
            n = 6 if self.ndm == 2 else 9
            if len(args) != n:
                raise ValueError("element('elasticBeamColumn', tag, i, j, " +
                                 ("A, E, Iz, transfTag)" if self.ndm == 2 else "A, E, G, J, Iy, Iz, transfTag)"))
            i, j, *props, transf = args
            density = self._volumetric(tag, per_length, props[0])
            if transf not in self._transforms:
                raise ValueError(f"element {tag}: geomTransf {transf} is not defined")
            if self.ndm == 2:
                a, e, iz = props
                row = {"nodeI": self._node(i), "nodeJ": self._node(j), "e": float(e), "a": float(a),
                       "iz": float(iz), "transform": self._transforms[transf], "density": density}
            else:
                a, e, g, jt, iy, iz = props
                row = {"nodeI": self._node(i), "nodeJ": self._node(j), "e": float(e), "g": float(g),
                       "a": float(a), "j": float(jt), "iy": float(iy), "iz": float(iz),
                       "transform": self._transforms[transf], "density": density}
            self._add_row("elasticBeamColumns2d" if self.ndm == 2 else "elasticBeamColumns3d", tag, row)
        else:
            raise _unsupported(f"element {kind!r}")

    @staticmethod
    def _volumetric(tag: int, per_length: float, area: float) -> float:
        if per_length and not area > 0:
            raise ValueError(f"element {tag}: a mass per length needs a positive area")
        return per_length / area if per_length else 0.0

    def _add_row(self, table: str, tag: int, row: dict) -> None:
        columns = self._tables.setdefault(table, {k: [] for k in row})
        for key, value in row.items():
            columns[key].append(value)
        self._element_tags[tag] = (table, len(columns["nodeI"]) - 1)

    # -- loads ---------------------------------------------------------------------------------------
    def timeSeries(self, kind: str, tag: int, *args) -> None:
        """``timeSeries('Linear'|'Constant', tag [, '-factor', f])``."""
        self._open("timeSeries")
        if kind not in ("Linear", "Constant"):
            raise _unsupported(f"timeSeries {kind!r}")
        factor = 1.0
        if args:
            if len(args) != 2 or args[0] != "-factor":
                raise _unsupported(f"timeSeries option {args[0]!r}")
            factor = float(args[1])
        self._series[tag] = ({"kind": "linear", "slope": 1.0} if kind == "Linear" else {"kind": "constant"}, factor)

    def pattern(self, kind: str, tag: int, series: int, *args) -> None:
        """``pattern('Plain', tag, timeSeriesTag [, '-fact', f])``. Later ``load`` calls go to this pattern."""
        self._open("pattern")
        if kind != "Plain":
            raise _unsupported(f"pattern {kind!r}")
        if series not in self._series:
            raise ValueError(f"pattern {tag}: timeSeries {series} is not defined")
        factor = 1.0
        if args:
            if len(args) != 2 or args[0] != "-fact":
                raise _unsupported(f"pattern option {args[0]!r}")
            factor = float(args[1])
        spec, series_factor = self._series[series]
        self._pattern_tags[tag] = len(self._patterns)
        self._patterns.append((spec, series_factor * factor))
        self._current_pattern = self._pattern_tags[tag]

    def load(self, tag: int, *values: float) -> None:
        """``load(nodeTag, f1, ..., f_ndf)`` into the current pattern. Zero components are skipped."""
        self._open("load")
        if self._current_pattern is None:
            raise RuntimeError("load() needs a pattern(); there is no current load pattern")
        if len(values) != self.ndf:
            raise ValueError(f"load at node {tag}: expected {self.ndf} values, got {len(values)}")
        at = self._node(tag)
        for dof, value in enumerate(values):
            if value:
                for key, v in zip(("pattern", "node", "dof", "value", "stage"),
                                  (self._current_pattern, at, dof, float(value), 0)):
                    self._loads[key].append(v)

    # -- analysis ------------------------------------------------------------------------------------
    def constraints(self, kind: str, *args) -> None:
        if kind not in ("Plain", "Transformation"):
            raise _unsupported(f"constraints {kind!r}")

    def numberer(self, kind: str, *args) -> None:
        self._ignored("numberer", kind)

    def system(self, kind: str, *args) -> None:
        self._ignored("system", kind)

    def _ignored(self, command: str, kind: str) -> None:
        if kind not in _IGNORED[command]:
            raise ValueError(f"unknown {command} {kind!r}")

    def test(self, kind: str, tol: float, iterations: int, *flags) -> None:
        """``test('NormUnbalance'|'NormDispIncr'|'EnergyIncr', tol, maxIter [, printFlag])``."""
        self._open_analysis("test")
        if kind not in _TESTS:
            raise _unsupported(f"test {kind!r}")
        self._test = (kind, float(tol), int(iterations))

    def algorithm(self, kind: str, *flags) -> None:
        """``algorithm('Linear'|'Newton'|'ModifiedNewton' [, '-initial'])``."""
        self._open_analysis("algorithm")
        initial = "-initial" in flags
        if set(flags) - {"-initial"}:
            raise _unsupported(f"algorithm option {[f for f in flags if f != '-initial'][0]!r}")
        if kind == "Linear":
            self._algorithm = {"kind": "linear"}
        elif kind in ("Newton", "ModifiedNewton"):
            tangent = "initial" if initial else "current" if kind == "Newton" else "reuseAtStepStart"
            self._algorithm = {"kind": "newton", "tangent": tangent}
        else:
            raise _unsupported(f"algorithm {kind!r}")

    def integrator(self, kind: str, *args) -> None:
        """``integrator('LoadControl', dLambda)`` or ``integrator('DisplacementControl', node, dof, dU)``."""
        self._open_analysis("integrator")
        if kind == "LoadControl" and len(args) == 1:
            self._integrator = {"kind": "loadControl", "increment": float(args[0])}
        elif kind == "DisplacementControl" and len(args) == 3:
            node, dof, inc = args
            self._integrator = {"kind": "displacementControl", "node": self._node(node),
                                "dof": self._dof(dof), "increment": float(inc)}
        else:
            raise _unsupported(f"integrator {kind!r} with {len(args)} arguments")

    def analysis(self, kind: str, *args) -> None:
        """``analysis('Static')``."""
        self._open_analysis("analysis")
        if kind != "Static" or args:
            raise _unsupported(f"analysis {kind!r}")
        self._analysis = kind

    def _open_analysis(self, command: str) -> None:
        if self._session is not None:
            raise _unsupported(f"{command}() after analyze() (changing the analysis between analyze calls)")

    def _dof(self, dof: int) -> int:
        if not 1 <= dof <= self.ndf:
            raise ValueError(f"dof {dof} is out of range 1..{self.ndf}")
        return dof - 1

    def analyze(self, steps: int = 1, *rest) -> int:
        """Take ``steps`` steps. Returns 0, or -1 on failure (details in ``last_error``)."""
        if rest:
            raise _unsupported("analyze with a time step (static analysis only)")
        if self._session is None:
            self._start()
        outcome = json.loads(self._session.advance(int(steps)))
        for batch in outcome["bulkBatches"]:
            size = batch["elements"] * batch["points"] * batch["components"]
            if batch["values"]:
                (self._disp if batch["recorderIndex"] == 0 else self._react)[:] = batch["values"][-size:]
        self._time = outcome["loadFactor"]
        self.last_error = outcome["error"]
        return -1 if outcome["error"] else 0

    def _start(self) -> None:
        if self._integrator is None or self._analysis is None:
            raise RuntimeError("analyze() needs integrator() and analysis() first")
        self._session = _session(self.to_dict(steps=_BIG))
        n = len(self._fixed)
        self._disp, self._react = [0.0] * (n * self.ndf), [0.0] * (n * self.ndf)

    # -- results -------------------------------------------------------------------------------------
    def nodeDisp(self, tag: int, dof: int | None = None):
        """Displacement of one DOF (1-based) or, without ``dof``, a list of all of them."""
        return self._nodal(self._disp, tag, dof)

    def nodeReaction(self, tag: int, dof: int | None = None):
        """Reaction at a node (OpenSees needs ``reactions()`` first; see ``reactions``)."""
        if not self._record_reactions:
            raise RuntimeError("Model(reactions=False) does not record reactions")
        return self._nodal(self._react, tag, dof)

    def reactions(self) -> None:
        """Accepted for compatibility: reactions are recorded at every step, so there is nothing to compute."""

    def _nodal(self, values: list[float], tag: int, dof: int | None):
        at = self._node(tag) * self.ndf
        if not values:
            values = [0.0] * (len(self._fixed) * self.ndf)
        if dof is None:
            return values[at:at + self.ndf]
        return values[at + self._dof(dof)]

    def getTime(self) -> float:
        """The load factor reached by the last step (the pseudo-time of a static analysis)."""
        return self._time

    def nodeCoord(self, tag: int, dim: int | None = None):
        at = self._node(tag) * self.ndm
        coords = self._coords[at:at + self.ndm]
        return coords if dim is None else coords[dim - 1]

    def getNodeTags(self) -> list[int]:
        return list(self._tags)

    def getEleTags(self) -> list[int]:
        return list(self._element_tags)

    def eigen(self, count: int, *flags) -> list[float]:
        """The ``count`` lowest eigenvalues (omega squared) of the undeformed model. Solver flags are accepted
        and ignored."""
        from . import run

        results = run(self.to_dict(modes=count))
        self._modes = results.modes(0)
        return [m["frequency"] ** 2 for m in self._modes]

    def nodeEigenvector(self, tag: int, mode: int, dof: int | None = None):
        """Component(s) of mode ``mode`` (1-based) at a node, from the last ``eigen`` call."""
        if self._modes is None:
            raise RuntimeError("nodeEigenvector() needs eigen() first")
        shape = self._modes[mode - 1]["shape"]
        at = self._node(tag) * self.ndf
        return shape[at:at + self.ndf] if dof is None else shape[at + self._dof(dof)]

    # -- wire format ---------------------------------------------------------------------------------
    def to_dict(self, steps: int = 1, modes: int | None = None) -> dict:
        """The ``CarapaceInputV1`` document: a static stage of ``steps`` steps (or, with ``modes``, a modal
        stage) using the integrator, algorithm and test set so far. Recorder 0 holds every node displacement
        per step, recorder 1 every reaction (``Results.bulk``)."""
        d: dict[str, Any] = {
            "header": {"schemaVersion": 1, "ndm": self.ndm, "engineVersion": "python"},
            "nodes": {"coords": self._coords, "fixed": self._fixed,
                      "massNodeIndex": list(self._mass), "mass": [v for m in self._mass.values() for v in m]},
        }
        if self._materials:
            d["materials"] = self._materials
        d.update({table: columns for table, columns in self._tables.items()})
        if self._patterns:
            d["loadPatterns"] = {"series": [s for s, _ in self._patterns], "scaleFactor": [f for _, f in self._patterns]}
        if self._loads["node"]:
            d["nodalLoads"] = self._loads
        if modes is not None:
            stages = [{"kind": "modal", "id": "modal", "modes": modes}]
            recorders: list[dict] = []
        else:
            if self._integrator is None:
                raise RuntimeError("no integrator(): set one before building the analysis")
            kind, tol, iterations = self._test
            stages = [{"kind": "static", "id": "static", "steps": steps, "integrator": self._integrator,
                       "algorithm": self._algorithm,
                       "convergence": {"kind": _TESTS[kind], "tol": tol, "maxIter": iterations},
                       "holdPatternsAfter": []}]
            recorders = [{"response": "bulk", "table": "nodes", "quantity": "displacement"}]
            if self._record_reactions:
                recorders.append({"response": "bulk", "table": "nodes", "quantity": "reaction"})
        d["sequence"] = {"stages": stages, "recorders": recorders}
        d.update(self.raw)
        return d


def _session(doc: dict):
    try:
        return _carapace.Session(json.dumps(doc))
    except ValueError as e:
        raise ModelError(json.loads(str(e))["errors"]) from None
