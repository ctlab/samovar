"""Compatibility shim — use ``samovar.baselines.identity_assembler``."""

from samovar.baselines.identity_assembler import main

if __name__ == "__main__":
    raise SystemExit(main())
