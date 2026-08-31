#!/usr/bin/env python3
"""Logistic abundance correction — thin wrapper around the builtin.

Prefer ``samovar prepare --export logistic`` (the preprocess default).
Import this file only if you need a custom name/flags overlay::

    samovar tools import -n logistic_correction \\
        --exec-path examples/logistic-correction/logistic_corrector.py \\
        --type export --pytest --flags "--min-efficiency 0.05"
"""

from __future__ import annotations

from samovar.abundance_correctors import export_logistic as export

__all__ = ["export"]


def _cli(argv=None) -> int:
    from samovar.abundance_correctors import main

    return main(argv)


if __name__ == "__main__":
    raise SystemExit(_cli())
