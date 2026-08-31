"""Identity Annotation → abundance export (``samovar tools import --type export``)."""

from __future__ import annotations

from pathlib import Path

from samovar.abundance_correctors import export_identity
from samovar.parse_annotators import Annotation


def export(annotation: Annotation, dest, config):
    return export_identity(annotation, dest, config)


def _cli(argv=None) -> int:
    import argparse
    from samovar.abundance_correctors import as_annotation

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True)
    parser.add_argument("-r", dest="reference", default="")
    args, extra = parser.parse_known_args(argv)
    cfg = {"to": "abundance", "extra_argv": extra}
    if args.reference:
        cfg["reference"] = args.reference
    export(as_annotation(args.i), Path(args.o), cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
