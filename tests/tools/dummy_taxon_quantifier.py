"""Compatibility shim — use ``samovar.baselines.identity_taxon_quantifier``."""

from samovar.baselines.identity_taxon_quantifier import main

if __name__ == "__main__":
    raise SystemExit(main())
