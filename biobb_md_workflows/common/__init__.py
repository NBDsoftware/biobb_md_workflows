#!/usr/bin/env python3

"""Shared helpers for the biobb_md_workflows CLI pipelines."""

from typing import Any


def to_yaml(value: Any) -> str:
    """Render a Python value as a YAML scalar for injection into a config template.

    Ensures ``None`` becomes ``null`` (not the string ``"None"``), booleans become
    lowercase ``true``/``false``, and lists become YAML flow sequences. Everything
    else is rendered with ``str()``.
    """
    if value is None:
        return "null"
    # bool must be checked before int/float (bool is a subclass of int)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(to_yaml(v) for v in value) + "]"
    return str(value)
