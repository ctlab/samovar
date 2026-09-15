"""Compatibility shim — use ``samovar.baselines.constant_mag_taxonomy``."""

from samovar.baselines.constant_mag_taxonomy import main

if __name__ == "__main__":
    raise SystemExit(main())
