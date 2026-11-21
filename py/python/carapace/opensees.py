"""``import carapace.opensees as ops``: the OpenSeesPy spelling over one module-level :class:`~carapace.Model`.

Every command is the corresponding ``Model`` method (``ops.node``, ``ops.analyze``, ...). ``ops.model(...)`` and
``ops.wipe()`` reset it. For several models at once (``run_many``, threads) use ``Model`` objects instead.
"""

from __future__ import annotations

from .model import Model

_model = Model()


def __getattr__(name: str):
    if name == "model":
        return _model.model
    return getattr(_model, name)


def __dir__():
    return sorted(n for n in dir(Model) if not n.startswith("_"))
