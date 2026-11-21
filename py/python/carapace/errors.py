"""Exceptions shared by the package."""

from __future__ import annotations

import json


class ModelError(Exception):
    """The model is invalid. ``errors`` is a list of ``{stage, path?, detail}`` dicts."""

    def __init__(self, errors: list[dict]):
        self.errors = errors
        super().__init__(json.dumps(errors, indent=2))
