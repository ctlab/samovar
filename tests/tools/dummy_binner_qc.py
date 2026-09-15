"""Compatibility shim — use ``samovar.baselines.passthrough_binner_qc``."""

from samovar.baselines.passthrough_binner_qc import main

if __name__ == "__main__":
    raise SystemExit(main())
