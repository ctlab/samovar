#!/usr/bin/env python3
"""Python replacement for workflow/compare_annotations.R.

Prepare calls ``samovar score`` instead. This entry point keeps the previous
combined invocation for existing callers.
"""

from samovar.viz_annotation import legacy_compare_main


def main(argv=None) -> int:
    return legacy_compare_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
