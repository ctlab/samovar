"""Compatibility shim — use ``samovar.baselines.identity_read_assigner``."""

from samovar.baselines.identity_read_assigner import main

if __name__ == "__main__":
    raise SystemExit(main())
