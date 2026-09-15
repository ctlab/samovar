"""Compatibility shim — use ``samovar.baselines.identity_binner``."""

from samovar.baselines.identity_binner import main

if __name__ == "__main__":
    raise SystemExit(main())
