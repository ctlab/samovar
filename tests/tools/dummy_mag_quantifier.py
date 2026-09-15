"""Compatibility shim — use ``samovar.baselines.constant_mag_quantifier``."""

from samovar.baselines.constant_mag_quantifier import main

if __name__ == "__main__":
    raise SystemExit(main())
