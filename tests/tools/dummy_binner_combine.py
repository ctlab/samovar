"""Compatibility shim — use ``samovar.baselines.identity_binner_combine``."""

from samovar.baselines.identity_binner_combine import main

if __name__ == "__main__":
    raise SystemExit(main())
