"""Convert ``Annotation`` between named formats.

Built-ins: ``annotation`` (per-read long table), ``abundance``
(``taxid`` × ``N_<sample>``), ``kraken2`` (kreport), and ``cami`` (bioboxes
``.profile``). Other names are ``annotation_converter`` tools imported with
``samovar tools import --type annotation-converter``.

A converter is the hub around ``Annotation``:

* ``load(path, config) -> Annotation`` when ``--from`` is a custom format
* ``dump(annotation, dest, config) -> path`` when ``--to`` is a custom format
* or a single ``convert(src, dest, config)`` (``src`` is ``Annotation`` or Path)

``samovar convert -i IN -o OUT --to abundance`` uses the builtin export.
``samovar convert -i IN -o reports --to kraken2`` writes Kraken 2 kreport files
(needs taxdump; ``--taxonomy ncbi`` or ``gtdb``).
``samovar convert -i IN -o profiles --to cami`` writes CAMI bioboxes ``.profile``
files (RFC 0.10.0; ``@TaxonomyID`` follows ``--taxonomy``).
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Union

import pandas as pd

from samovar.abundance import (
    input_to_abundance_tables,
    is_abundance_table,
    load_table_input,
    n_sample_columns,
    write_abundance_dir,
)
from samovar.main_config import (
    iter_tools,
    parse_tool_entry,
    tool_flags,
    tool_path,
)
from samovar.parse_annotators import Annotation
from samovar.paths import load_config
from samovar.table_regenerators import extra_flags_argv

PathLike = Union[str, Path]

BUILTIN_FORMATS = ("annotation", "abundance", "kraken2", "kraken2_mpa", "cami")
FORMAT_ALIASES = {
    "annotation": "annotation",
    "annot": "annotation",
    "long": "annotation",
    "reads": "annotation",
    "abundance": "abundance",
    "abund": "abundance",
    "table": "abundance",
    "otu": "abundance",
    "counts": "abundance",
    "kraken2": "kraken2",
    "kreport": "kraken2",
    "k2report": "kraken2",
    "kraken2_report": "kraken2",
    "kraken2_style": "kraken2",
    "kraken2style": "kraken2",
    "kraken2_mpa": "kraken2_mpa",
    "kreport_mpa": "kraken2_mpa",
    "mpa": "kraken2_mpa",
    "cami": "cami",
    "cami_profile": "cami",
    "camiprofile": "cami",
    "bioboxes": "cami",
    "profiling": "cami",
}
CONVERTER_GROUP = "annotation_converter"


class UnknownFormatError(ValueError):
    """``--from`` / ``--to`` is neither builtin nor an imported converter."""


def normalize_format(name: str) -> str:
    key = str(name or "").strip().lower().replace("-", "_")
    if not key:
        return ""
    return FORMAT_ALIASES.get(key, key)


def is_builtin_format(name: str) -> bool:
    return normalize_format(name) in BUILTIN_FORMATS


def rescale_abundance_tables(
    tables: Dict[str, pd.DataFrame],
    n_reads: Optional[int],
) -> Dict[str, pd.DataFrame]:
    if not n_reads:
        return tables
    out: Dict[str, pd.DataFrame] = {}
    for name, table in tables.items():
        work = table.copy()
        cols = n_sample_columns(work)
        if not cols:
            out[name] = work
            continue
        totals = work[cols].sum(axis=0).replace(0, 1)
        work[cols] = (work[cols].div(totals, axis=1) * int(n_reads)).round()
        out[name] = work
    return out


def annotation_to_abundance(
    annotation: Annotation,
    *,
    n_reads: Optional[int] = None,
) -> Dict[str, pd.DataFrame]:
    tables = input_to_abundance_tables(annotation)
    tables = rescale_abundance_tables(tables, n_reads)
    annotation.abundance_tables = tables
    return tables


def annotation_long_table(annotation: Annotation) -> pd.DataFrame:
    """Per-read table (``taxID_*``, optional ``seq`` / ``sample`` / ``true``)."""
    df = annotation.DataFrame.copy()
    if df.empty:
        return df
    if "seq" not in df.columns:
        df = df.reset_index()
        if "index" in df.columns:
            df = df.rename(columns={"index": "seq"})
    if getattr(annotation, "true_annotation", None) and "true" not in df.columns:
        n = min(len(annotation.true_annotation), len(df))
        if n:
            df = df.copy()
            df["true"] = list(annotation.true_annotation)[:n] + [""] * (len(df) - n)
    return df


def write_abundance(dest: PathLike, tables: Dict[str, pd.DataFrame]) -> Path:
    dest_path = Path(dest)
    if dest_path.suffix.lower() in {".csv", ".tsv", ".txt"}:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        if len(tables) == 1:
            next(iter(tables.values())).to_csv(dest_path, index=False)
            return dest_path
        dest_path = dest_path.with_suffix("")
    dest_path.mkdir(parents=True, exist_ok=True)
    return write_abundance_dir(dest_path, tables)


def write_annotation_table(dest: PathLike, frame: pd.DataFrame) -> Path:
    dest_path = Path(dest)
    if dest_path.suffix.lower() in {".csv", ".tsv", ".txt"} or dest_path.suffix:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(dest_path, index=False)
        return dest_path
    dest_path.mkdir(parents=True, exist_ok=True)
    out = dest_path / "combined_annotation_table.csv"
    frame.to_csv(out, index=False)
    return out


def annotation_to_cami_map(
    annotation: Annotation,
    *,
    n_reads: Optional[int] = None,
    taxdump: Optional[PathLike] = None,
    taxonomy: str = "ncbi",
) -> Dict[str, Dict[str, str]]:
    """``{annotator: {sample: profile_text}}`` via abundance tables + taxdump."""
    from samovar.cami_profile import abundance_tables_to_cami, try_taxonomy

    tables = annotation_to_abundance(annotation, n_reads=n_reads)
    if not tables:
        return {}
    return abundance_tables_to_cami(tables, try_taxonomy(taxdump, taxonomy=taxonomy))


def annotation_to_kreport_map(
    annotation: Annotation,
    *,
    n_reads: Optional[int] = None,
    taxdump: Optional[PathLike] = None,
    mpa: bool = False,
    taxonomy: str = "ncbi",
) -> Dict[str, Dict[str, str]]:
    """``{annotator: {sample: report_text}}`` via abundance tables + taxdump."""
    from samovar.kreport import KReportTaxonomy, abundance_tables_to_reports

    tables = annotation_to_abundance(annotation, n_reads=n_reads)
    if not tables:
        return {}
    tax = KReportTaxonomy.from_taxdump(taxdump, taxonomy=taxonomy)
    return abundance_tables_to_reports(tables, tax, mpa=mpa)


def _kreport_options(fmt: str, kwargs: Dict[str, Any]) -> tuple:
    """``(mpa, taxdump, taxonomy)`` from format name, kwargs, and ``--flags``."""
    from samovar.taxonomy import normalize_taxonomy

    mpa = bool(kwargs.get("mpa")) or normalize_format(fmt) == "kraken2_mpa"
    taxdump = kwargs.get("taxdump") or kwargs.get("taxdump_dir")
    taxonomy = kwargs.get("taxonomy") or kwargs.get("taxonomy_id") or "ncbi"
    argv = list(kwargs.get("extra_argv") or extra_flags_argv(kwargs.get("extra_flags")))
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in {"--mpa", "--use-mpa-style"}:
            mpa = True
        elif tok.startswith("--taxdump="):
            taxdump = tok.split("=", 1)[1]
        elif tok == "--taxdump" and i + 1 < len(argv):
            taxdump = argv[i + 1]
            i += 1
        elif tok.startswith("--taxonomy="):
            taxonomy = tok.split("=", 1)[1]
        elif tok == "--taxonomy" and i + 1 < len(argv):
            taxonomy = argv[i + 1]
            i += 1
        i += 1
    return mpa, taxdump or None, normalize_taxonomy(taxonomy)


def dump_builtin(annotation: Annotation, dest: PathLike, fmt: str, **kwargs) -> Path:
    name = normalize_format(fmt)
    if name == "abundance":
        tables = annotation_to_abundance(annotation, n_reads=kwargs.get("n_reads"))
        if not tables:
            raise ValueError("no abundance tables to write (empty annotation)")
        return write_abundance(dest, tables)
    if name == "annotation":
        frame = annotation_long_table(annotation)
        if frame.empty:
            raise ValueError(
                "no per-read annotation table to write; use --to abundance"
            )
        return write_annotation_table(dest, frame)
    if name in {"kraken2", "kraken2_mpa"}:
        from samovar.kreport import dump_kreport

        tables = annotation_to_abundance(annotation, n_reads=kwargs.get("n_reads"))
        if not tables:
            raise ValueError("no abundance tables to write (empty annotation)")
        mpa, taxdump, taxonomy = _kreport_options(name, kwargs)
        return dump_kreport(
            tables, dest, taxdump=taxdump, mpa=mpa, taxonomy=taxonomy
        )
    if name == "cami":
        from samovar.cami_profile import dump_cami

        tables = annotation_to_abundance(annotation, n_reads=kwargs.get("n_reads"))
        if not tables:
            raise ValueError("no abundance tables to write (empty annotation)")
        _mpa, taxdump, taxonomy = _kreport_options(name, kwargs)
        return dump_cami(tables, dest, taxdump=taxdump, taxonomy=taxonomy)
    raise UnknownFormatError(fmt)


def detect_format(path: PathLike) -> str:
    source = Path(path)
    if source.is_file():
        try:
            peek = pd.read_csv(source, nrows=8)
        except Exception:
            return "annotation"
        if is_abundance_table(peek):
            return "abundance"
        return "annotation"
    from samovar.abundance import dir_looks_like_annotation

    if dir_looks_like_annotation(source):
        return "annotation"
    return "abundance"


def load_builtin(path: PathLike, fmt: str = "") -> Annotation:
    name = normalize_format(fmt) or detect_format(path)
    if name in {"kraken2", "kraken2_mpa", "cami"}:
        raise ValueError(
            f"{name}-style reports are export-only; use --from annotation or abundance"
        )
    data = load_table_input(path)
    if name == "abundance" and not getattr(data, "abundance_tables", None):
        tables = input_to_abundance_tables(data)
        if tables:
            return Annotation.from_abundance_tables(tables)
    return data


def imported_converters() -> Dict[str, Any]:
    cfg = load_config()
    out: Dict[str, Any] = {}
    for name, spec in iter_tools(cfg).items():
        rec = parse_tool_entry(spec, name)
        group = (rec[3] if len(rec) > 3 else "").strip()
        if group != CONVERTER_GROUP:
            continue
        out[name.split(":")[0]] = spec
        out[name] = spec
    return out


def lookup_converter(name: str) -> Optional[List[str]]:
    token = str(name or "").strip()
    if not token:
        return None
    converters = imported_converters()
    if token in converters:
        return parse_tool_entry(converters[token], token)
    bare = token.split(":")[0]
    if bare in converters:
        return parse_tool_entry(converters[bare], bare)
    return None


def _load_python_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(f"samovar_annconv_{name}", path)
    if spec is None or spec.loader is None:
        return None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _python_hooks(path: Path, name: str):
    if path.suffix.lower() != ".py" or not path.is_file():
        return None, None, None
    try:
        module = _load_python_module(path, name)
    except Exception:
        return None, None, None
    if module is None:
        return None, None, None
    load_fn = getattr(module, "load", None)
    dump_fn = getattr(module, "dump", None)
    convert_fn = getattr(module, "convert", None)
    cls = getattr(module, "AnnotationConverter", None)
    if cls is not None and not (callable(load_fn) or callable(dump_fn) or callable(convert_fn)):
        try:
            inst = cls()
        except Exception:
            inst = None
        if inst is not None:
            load_fn = getattr(inst, "load", None) or load_fn
            dump_fn = getattr(inst, "dump", None) or dump_fn
            convert_fn = getattr(inst, "convert", None) or convert_fn
    return (
        load_fn if callable(load_fn) else None,
        dump_fn if callable(dump_fn) else None,
        convert_fn if callable(convert_fn) else None,
    )


def _run_cli_converter(exe: Path, src: PathLike, dest: PathLike, config: Dict[str, Any]) -> Path:
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(exe), "-i", str(src), "-o", str(dest_path)]
    src_fmt = str(config.get("from") or "").strip()
    dst_fmt = str(config.get("to") or "").strip()
    if src_fmt:
        cmd.extend(["--from", src_fmt])
    if dst_fmt:
        cmd.extend(["--to", dst_fmt])
    extra = list(config.get("extra_argv") or extra_flags_argv(config.get("extra_flags")))
    cmd.extend(extra)
    if exe.suffix.lower() == ".py":
        cmd = [sys.executable, *cmd]
    subprocess.run(cmd, check=True)
    return dest_path


def _spec_path_and_name(spec: List[str], name: str) -> tuple:
    path = Path(tool_path(spec, name)).expanduser()
    return path, name.split(":")[0]


def load_custom(path: PathLike, spec: List[str], name: str, config: Dict[str, Any]) -> Annotation:
    exe, bare = _spec_path_and_name(spec, name)
    load_fn, _dump, convert_fn = _python_hooks(exe, bare)
    cfg = dict(config or {})
    if load_fn is not None:
        loaded = load_fn(path, cfg)
        if isinstance(loaded, Annotation):
            return loaded
        if isinstance(loaded, pd.DataFrame):
            return Annotation.from_long_table(loaded)
        tables = input_to_abundance_tables(loaded)
        if tables:
            return Annotation.from_abundance_tables(tables)
        raise TypeError(f"{exe} load() must return Annotation or a table")
    if convert_fn is not None:
        tmp = Path(cfg.get("output") or path).parent / f".{bare}_loaded"
        loaded = convert_fn(path, tmp, cfg)
        if isinstance(loaded, Annotation):
            return loaded
        return load_builtin(tmp if loaded is None else loaded, cfg.get("to") or "annotation")
    staging = Path(cfg.get("output") or Path(path).parent) / f".{bare}_in"
    _run_cli_converter(exe, path, staging, cfg)
    return load_builtin(staging, cfg.get("to") or "")


def dump_custom(
    annotation: Annotation,
    dest: PathLike,
    spec: List[str],
    name: str,
    config: Dict[str, Any],
) -> Path:
    exe, bare = _spec_path_and_name(spec, name)
    _load, dump_fn, convert_fn = _python_hooks(exe, bare)
    cfg = dict(config or {})
    dest_path = Path(dest)
    if dump_fn is not None:
        written = dump_fn(annotation, dest_path, cfg)
        return Path(written) if written else dest_path
    if convert_fn is not None:
        written = convert_fn(annotation, dest_path, cfg)
        return Path(written) if written else dest_path
    staging = dest_path.parent / f".{bare}_annotation.csv"
    if annotation_long_table(annotation).empty and annotation_to_abundance(annotation):
        src = write_abundance(dest_path.parent / f".{bare}_abundance", annotation_to_abundance(annotation))
    else:
        src = write_annotation_table(staging, annotation_long_table(annotation))
    return _run_cli_converter(exe, src, dest_path, cfg)


def available_formats() -> List[str]:
    names = list(BUILTIN_FORMATS)
    for name in imported_converters():
        if ":" not in name and name not in names:
            names.append(name)
    return names


def resolve_converter_name(
    *,
    dest_format: str,
    source_format: str = "",
    converter: str = "",
) -> str:
    if converter.strip():
        return converter.strip()
    dest = normalize_format(dest_format)
    src = normalize_format(source_format)
    if lookup_converter(dest):
        return dest
    if lookup_converter(src):
        return src
    return ""


def parse_export_formats(*values: Any) -> List[str]:
    """Normalize a mix of ``--to`` tokens / lists / comma-separated names.

    ``mpa`` is Kraken 2 ``--use-mpa-style``. Unknown names raise ``UnknownFormatError``.
    """
    out: List[str] = []
    for raw in values:
        if raw in (None, False, ""):
            continue
        if isinstance(raw, (list, tuple, set)):
            out.extend(parse_export_formats(*raw))
            continue
        if isinstance(raw, bool):
            continue
        for piece in str(raw).replace(",", " ").split():
            name = normalize_format(piece)
            if not name:
                continue
            if not is_builtin_format(name) and not lookup_converter(name):
                raise UnknownFormatError(
                    f"unknown export format {piece!r}. "
                    f"Builtins: {', '.join(BUILTIN_FORMATS)} (mpa → kraken2-mpa)"
                )
            if name not in out:
                out.append(name)
    return out


def formats_from_prepare_args(args: Any) -> List[str]:
    """Collect ``--to`` / ``--to-mpa`` / ``--to-cami`` / … from prepare argparse."""
    pieces: List[Any] = list(getattr(args, "export_to", None) or [])
    flag_map = (
        ("to_mpa", "mpa"),
        ("to_cami", "cami"),
        ("to_kraken2", "kraken2"),
        ("to_abundance", "abundance"),
        ("to_annotation", "annotation"),
    )
    for attr, fmt in flag_map:
        if getattr(args, attr, False):
            pieces.append(fmt)
    return parse_export_formats(*pieces)


def export_annotation_formats(
    source: PathLike,
    dest_root: PathLike,
    formats: Iterable[str],
    *,
    source_format: str = "",
    taxdump: str = "",
    n_reads: Optional[int] = None,
    extra_flags: str = "",
    taxonomy: str = "ncbi",
) -> Dict[str, Path]:
    """Write each format under ``dest_root/<format>/`` (or ``dest_root`` if one format)."""
    names = parse_export_formats(*list(formats))
    if not names:
        raise ValueError("--to is required")
    dest_path = Path(dest_root)
    written: Dict[str, Path] = {}
    if len(names) == 1:
        written[names[0]] = convert_annotation(
            source,
            dest_path,
            source_format=source_format,
            dest_format=names[0],
            taxdump=taxdump,
            n_reads=n_reads,
            extra_flags=extra_flags,
            taxonomy=taxonomy,
        )
        return written
    dest_path.mkdir(parents=True, exist_ok=True)
    for fmt in names:
        written[fmt] = convert_annotation(
            source,
            dest_path / fmt,
            source_format=source_format,
            dest_format=fmt,
            taxdump=taxdump,
            n_reads=n_reads,
            extra_flags=extra_flags,
            taxonomy=taxonomy,
        )
    return written


def convert_annotation(
    source: PathLike,
    dest: PathLike,
    *,
    source_format: str = "",
    dest_format: str,
    converter: str = "",
    n_reads: Optional[int] = None,
    extra_flags: str = "",
    taxdump: str = "",
    mpa: bool = False,
    taxonomy: str = "ncbi",
) -> Path:
    src_fmt = normalize_format(source_format) or detect_format(source)
    dst_fmt = normalize_format(dest_format)
    if not dst_fmt:
        raise ValueError("--to is required")
    conv_name = resolve_converter_name(
        dest_format=dst_fmt, source_format=src_fmt, converter=converter
    )
    spec = lookup_converter(conv_name) if conv_name else None
    cfg = {
        "from": src_fmt,
        "to": dst_fmt,
        "input": str(source),
        "output": str(dest),
        "n_reads": n_reads,
        "extra_flags": extra_flags,
        "extra_argv": extra_flags_argv(extra_flags),
        "taxdump": taxdump or "",
        "mpa": mpa,
        "taxonomy": taxonomy or "ncbi",
    }
    if is_builtin_format(src_fmt):
        annotation = load_builtin(source, src_fmt)
    elif spec is not None:
        annotation = load_custom(source, spec, conv_name, cfg)
    else:
        raise UnknownFormatError(
            f"unknown --from {source_format or src_fmt!r}. "
            f"Builtins: {', '.join(BUILTIN_FORMATS)}. "
            f"Imported converters: {', '.join(sorted(imported_converters()) or ['(none)'])}"
        )
    if is_builtin_format(dst_fmt) and (not conv_name or is_builtin_format(conv_name)):
        return dump_builtin(
            annotation,
            dest,
            dst_fmt,
            n_reads=n_reads,
            taxdump=taxdump or None,
            mpa=mpa,
            taxonomy=taxonomy or "ncbi",
            extra_flags=extra_flags,
            extra_argv=extra_flags_argv(extra_flags),
        )
    if spec is None:
        spec = lookup_converter(dst_fmt)
        conv_name = dst_fmt
    if spec is None:
        raise UnknownFormatError(
            f"unknown --to {dest_format!r}. "
            f"Builtins: {', '.join(BUILTIN_FORMATS)}. "
            f"Import a converter: samovar tools import -n {dst_fmt} "
            f"--type annotation-converter --exec-path …"
        )
    flags = merge_converter_flags(spec, conv_name, extra_flags)
    cfg["extra_flags"] = flags
    cfg["extra_argv"] = extra_flags_argv(flags)
    return dump_custom(annotation, dest, spec, conv_name, cfg)


def merge_converter_flags(spec: List[str], name: str, extra: str) -> str:
    imported = tool_flags(spec, name)
    return " ".join(p for p in (imported, extra or "") if p).strip()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="samovar convert",
        description=(
            "Convert an Annotation between formats. Built-ins: annotation, "
            "abundance, kraken2 (kreport), kraken2-mpa, cami. "
            "Lineage formats use --taxonomy ncbi (default) or gtdb. "
            "Other --to/--from names are annotation-converter tools from "
            "`samovar tools import --type annotation-converter`."
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="File or directory")
    parser.add_argument("-o", "--output", required=True, help="File or directory")
    parser.add_argument(
        "--from",
        dest="source_format",
        default="",
        help="Input format (default: detect). Built-ins: annotation, abundance.",
    )
    parser.add_argument(
        "--to",
        dest="dest_format",
        action="append",
        required=True,
        help=(
            "Output format (repeatable): abundance, annotation, kraken2, mpa "
            "(kraken2 --use-mpa-style), cami, or an imported converter. "
            "Several --to values write dest/<format>/"
        ),
    )
    parser.add_argument(
        "--converter",
        default="",
        help="annotation-converter tool name (default: --to / --from if imported)",
    )
    parser.add_argument(
        "--n-reads",
        dest="n_reads",
        type=int,
        default=None,
        help="Rescale abundance so each sample sums to this many reads",
    )
    parser.add_argument(
        "--flags",
        default="",
        help="Extra CLI flags forwarded to a custom converter (or --use-mpa-style / --taxdump)",
    )
    parser.add_argument(
        "--taxdump",
        default="",
        help="Taxonomy dump directory for --to kraken2 / cami (NCBI or GTDB)",
    )
    parser.add_argument(
        "--taxonomy",
        default="ncbi",
        help="Taxonomy system: ncbi (default) or gtdb",
    )
    parser.add_argument(
        "--mpa",
        action="store_true",
        help="With --to kraken2, write Kraken 2 --use-mpa-style reports",
    )
    return parser


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    try:
        formats = parse_export_formats(*(args.dest_format or []))
        if len(formats) > 1:
            written = export_annotation_formats(
                args.input,
                args.output,
                formats,
                source_format=args.source_format,
                taxdump=args.taxdump,
                n_reads=args.n_reads,
                extra_flags=args.flags,
                taxonomy=args.taxonomy,
            )
            dest = Path(args.output)
            print(" ".join(str(p) for p in written.values()))
            return 0
        dest = convert_annotation(
            args.input,
            args.output,
            source_format=args.source_format,
            dest_format=formats[0] if formats else "",
            converter=args.converter,
            n_reads=args.n_reads,
            extra_flags=args.flags,
            taxdump=args.taxdump,
            mpa=args.mpa,
            taxonomy=args.taxonomy,
        )
    except (UnknownFormatError, ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
