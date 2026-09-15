"""Compatibility shim — use ``samovar.baselines.constant_taxid``."""

from samovar.baselines.constant_taxid import main, parse_output

__all__ = ["main", "parse_output"]

if __name__ == "__main__":
    raise SystemExit(main())
