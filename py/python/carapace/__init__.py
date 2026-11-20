"""Python interface to the Carapace finite-element engine.

Build a model with :class:`Model` (or pass any ``CarapaceInputV1`` dict / JSON text), then
``carapace.validate(model)`` and ``carapace.run(model)``. Errors come back as structured data, never
as a bare exception, so a script (or an agent) can read the path of the bad value, fix it and retry.
"""

from __future__ import annotations

import json
from typing import Any, Union

import numpy as np

from . import _carapace
from .model import Model

__all__ = ["Model", "Results", "ModelError", "validate", "run", "run_many"]

ModelLike = Union["Model", dict, str]


class ModelError(Exception):
    """The model is invalid. ``errors`` is a list of ``{stage, path?, detail}`` dicts."""

    def __init__(self, errors: list[dict]):
        self.errors = errors
        super().__init__(json.dumps(errors, indent=2))


def _text(model: ModelLike) -> str:
    if isinstance(model, Model):
        return json.dumps(model.to_dict())
    if isinstance(model, dict):
        return json.dumps(model)
    return model


def validate(model: ModelLike) -> dict:
    """Decode without running. Returns ``{"valid": True, "ndm": n}`` or ``{"valid": False, "errors": [...]}``."""
    return json.loads(_carapace.validate(_text(model)))


class Results:
    """The report of a run. ``raw`` is the full JSON document; helpers expose recorders as numpy arrays."""

    def __init__(self, raw: dict):
        self.raw = raw

    @property
    def completed(self) -> bool:
        return self.raw["completed"]

    @property
    def error(self) -> dict | None:
        """The analysis error that stopped the run (tagged dict), or ``None``."""
        return self.raw["error"]

    def recorder(self, index: int, stage: int | None = None) -> np.ndarray:
        """Samples of recorder ``index`` as an ``(n, 2)`` array of ``(pseudo_time, value)``."""
        rows: list[Any] = []
        for series in self.raw["recorders"]:
            if series["recorderIndex"] == index and (stage is None or series["stageIndex"] == stage):
                rows.extend(series["samples"])
        return np.asarray(rows, dtype=float).reshape(-1, 2)

    def bulk(self, index: int, stage: int | None = None) -> tuple[np.ndarray, np.ndarray]:
        """``(progress, values)`` of a bulk recorder; ``values`` is ``(frames, elements, points, components)``."""
        for batch in self.raw["bulk"]:
            if batch["recorderIndex"] == index and (stage is None or batch["stageIndex"] == stage):
                shape = (-1, batch["elements"], batch["points"], batch["components"])
                return np.asarray(batch["progress"]), np.asarray(batch["values"]).reshape(shape)
        raise KeyError(f"no bulk recorder {index}")

    def modes(self, stage: int = 0) -> list[dict]:
        """Modal results of a modal stage: frequency (rad/time), shape, participation, mass ratio."""
        for result in self.raw["modal"]["stages"]:
            if result["stageIndex"] == stage:
                return result["modes"]
        raise KeyError(f"no modal stage {stage}")


def run(model: ModelLike, step_budget: int = 64, raise_on_invalid: bool = True, threads: int = 1) -> Results:
    """Decode and run to completion. An invalid model raises :class:`ModelError` (or is returned as the
    ``{"valid": False, ...}`` dict when ``raise_on_invalid=False``). An analysis failure (non-convergence)
    does not raise: check ``results.completed`` / ``results.error``; results before the failure are kept.
    ``threads`` sets the solver/assembly thread count (``0`` = all cores); it never changes the results.
    The setting is process-wide."""
    out = json.loads(_carapace.run(_text(model), step_budget, threads))
    if not out["valid"]:
        if raise_on_invalid:
            raise ModelError(out["errors"])
        return out
    return Results(out["report"])


def run_many(models: list[ModelLike], threads: int = 0) -> list[Results | dict]:
    """Run independent models concurrently (``threads=0`` uses every core), in input order. An invalid model
    yields its ``{"valid": False, "errors": [...]}`` dict in its slot instead of raising, so one bad model
    does not lose the others."""
    out = []
    for text in _carapace.run_many([_text(m) for m in models], threads):
        doc = json.loads(text)
        out.append(Results(doc["report"]) if doc["valid"] else doc)
    return out
