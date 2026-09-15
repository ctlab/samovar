"""Compatibility shim — use ``samovar.baselines.translate_orfs``."""

from samovar.baselines.translate_orfs import main

if __name__ == "__main__":
    raise SystemExit(main())
