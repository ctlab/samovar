#!/usr/bin/env python3
"""Backward-compatible CLI for ML reprofiling (see ``samovar.reprofilers``)."""

from samovar.reprofilers import main

if __name__ == "__main__":
    raise SystemExit(main())
