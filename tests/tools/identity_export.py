"""Compatibility shim — use ``samovar.baselines.identity_export``."""

from samovar.baselines.identity_export import export, _cli

__all__ = ["export"]

if __name__ == "__main__":
    raise SystemExit(_cli())
