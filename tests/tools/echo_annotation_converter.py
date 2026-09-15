"""Compatibility shim — use ``samovar.baselines.identity_converter``."""

from samovar.baselines.identity_converter import dump, load

__all__ = ["dump", "load"]
