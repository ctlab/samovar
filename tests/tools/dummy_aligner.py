"""Compatibility shim — use ``samovar.baselines.identity_aligner``."""

from samovar.baselines.identity_aligner import main

if __name__ == "__main__":
    raise SystemExit(main())
