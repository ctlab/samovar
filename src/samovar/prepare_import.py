"""Import programs and databases named on ``samovar prepare``.

* no flag: the same dialogue as ``./install.sh ask``
* ``--autoimport 1``: record whatever can be resolved, then ask about the rest
* ``--autoimport 0``: import nothing
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Tuple

from samovar.abundance_correctors import resolve_export
from samovar.feature_importance import BUILTIN_NAMES as FEATURE_IMPORTANCE_BUILTINS
from samovar.feature_importance import NONE_TOKENS as FEATURE_IMPORTANCE_NONE
from samovar.install_ask import (
    ASK_TOOLS,
    AskTool,
    ask_tool,
    ask_yes_no,
    conda_layout,
    prompts_enabled,
    resolve_conda,
    sibling_executable,
    which_all,
)
from samovar.metagenome_generators import resolve_metagenome_generator
from samovar.qc import canonical_qc_name, is_identity_qc, parse_qc_postfix
from samovar.reads_generators import resolve_reads_generator
from samovar.regenerate import resolve_regeneration_mode
from samovar.reprofilers import resolve_reprofiler
from samovar.sample_filters import BUILTIN_FILTERS
from samovar.sample_scorers import BUILTIN_SAMPLE_SCORERS
from samovar.table_scorers import BUILTIN_TABLE_SCORERS
from samovar.tools_import import import_database, import_tool

# Annotators implemented as ``python -m samovar.*``, not external binaries.
_PYTHON_ANNOTATORS = frozenset(
    {
        "dummy",
        "dummy9606",
        "constant9606",
        "constant_taxid",
        "constant",
        "random",
        "kmer2",
        "kmer_counter",
        "dinuc",
        "dinucleotide",
        "gc",
        "gc_content",
        "gcfrac",
        "gc_fraction",
        "length",
        "read_length",
        "readlen",
        "kmer_encoder",
        "kmerenc",
        "encoder",
        "kmerencoder",
    }
)

_SKIP = frozenset({"", "none", "off", "false", "0", "all"})


def autoimport_mode(value: object) -> str:
    """``ask`` when the flag is omitted, ``auto`` for 1, ``off`` for 0."""
    if value is None or value == "":
        return "ask"
    if int(value) == 0:
        return "off"
    return "auto"


def offer_prepare_imports(
    args: argparse.Namespace,
    *,
    ask: Callable[[str], str] = input,
    interactive: Optional[bool] = None,
    ask_tool_fn: Callable = ask_tool,
    import_tool_fn: Callable = import_tool,
    import_db_fn: Callable = import_database,
    configured_fn: Optional[Callable[..., object]] = None,
    db_known_fn: Optional[Callable[[str, str], bool]] = None,
    which_fn: Callable[[Sequence[str]], Sequence[str]] = which_all,
    resolve_conda_fn: Callable[[str, str], Tuple[str, str]] = resolve_conda,
    sibling_fn: Callable[[str, str], str] = sibling_executable,
) -> None:
    """Import programs and databases this prepare command names."""
    mode = autoimport_mode(getattr(args, "autoimport", None))
    if mode == "off":
        return
    if interactive is None:
        interactive = prompts_enabled()
    if mode == "ask" and not interactive:
        return
    if configured_fn is None:
        from samovar.install_ask import configured_import

        configured_fn = configured_import
    if db_known_fn is None:
        db_known_fn = _database_known
    seen = set()
    for program, database in iter_prepare_refs(args):
        if program is not None and program[0] not in seen:
            seen.add(program[0])
            _offer_program(
                program,
                mode=mode,
                ask=ask,
                interactive=interactive,
                ask_tool_fn=ask_tool_fn,
                import_tool_fn=import_tool_fn,
                configured_fn=configured_fn,
                which_fn=which_fn,
                resolve_conda_fn=resolve_conda_fn,
                sibling_fn=sibling_fn,
            )
        if database is not None:
            _offer_database(
                database[0],
                database[1],
                mode=mode,
                ask=ask,
                interactive=interactive,
                import_db_fn=import_db_fn,
                db_known_fn=db_known_fn,
            )


def iter_prepare_refs(
    args: argparse.Namespace,
) -> List[Tuple[Optional[Tuple[str, str, str]], Optional[Tuple[str, str]]]]:
    """``((name, group, explicit_path)|None, (tool, db_name)|None)`` rows."""
    rows: List[Tuple[Optional[Tuple[str, str, str]], Optional[Tuple[str, str]]]] = []
    for attr in dir(args):
        if not attr.startswith("cmd_") or attr == "cmd_args":
            continue
        raw = getattr(args, attr, None)
        if not raw:
            continue
        for item in raw:
            text = item[0] if isinstance(item, (list, tuple)) else item
            program, database = _annotator_ref(attr, str(text or ""))
            if program or database:
                rows.append((program, database))
    for name, group in _named_programs(args):
        if _should_offer(name, group):
            rows.append(((name, group, ""), None))
    return rows


def _annotator_ref(
    attr: str, text: str
) -> Tuple[Optional[Tuple[str, str, str]], Optional[Tuple[str, str]]]:
    parts = text.split()
    if not parts:
        return None, None
    if "samovar." in text and (
        "-m" in parts or any(token.endswith("python") or token.endswith("python3") for token in parts[:1])
    ):
        return None, None
    cmd = parts[0]
    base = Path(cmd).name
    run_name = attr[4:] if attr.startswith("cmd_") else attr
    if _norm(run_name) in _PYTHON_ANNOTATORS or _norm(base) in _PYTHON_ANNOTATORS:
        return None, None
    explicit = ""
    candidate = Path(cmd).expanduser()
    try:
        if candidate.is_file():
            explicit = str(candidate)
    except OSError:
        explicit = ""
    program = None
    if _should_offer(base, "annotator") or explicit:
        program = (base, "annotator", explicit)
    db_token = parts[1] if len(parts) > 1 and not parts[1].startswith("-") else ""
    database = (base, db_token) if _database_token(db_token) else None
    return program, database


def _named_programs(args: argparse.Namespace) -> List[Tuple[str, str]]:
    found: List[Tuple[str, str]] = []

    def add(value: object, group: str, *, colon: bool = False) -> None:
        for token in _tokens(value):
            name = _after_colon(token) if colon else token
            if name:
                found.append((name, group))

    add(getattr(args, "reads_generator", None), "reads_generator")
    add(getattr(args, "metagenome_generator", None), "metagenome_generator")
    add(getattr(args, "table_reads_generator", None), "table_reads_generator")
    add(getattr(args, "reprofiler", None), "reprofiler")
    add(getattr(args, "export_corrector", None), "export")
    add(getattr(args, "feature_importance", None), "feature_importance")
    add(getattr(args, "table_score", None), "table_scoring")
    add(getattr(args, "scoring", None), "scoring")
    for attr in ("qc", "qc_initial", "qc_generated"):
        add(getattr(args, attr, None), "qc")
    for tool in parse_qc_postfix(getattr(args, "qc_postfix", None)).values():
        if tool:
            found.append((tool, "qc"))
    add(getattr(args, "sample_score", None), "sample_scoring", colon=True)
    add(getattr(args, "sample_filter", None), "sample_filtering", colon=True)
    add(getattr(args, "assembly_assembler", None), "assembler")
    add(getattr(args, "assembly_gene_caller", None), "gene_caller")
    add(getattr(args, "assembly_binners", None), "binner")
    add(getattr(args, "assembly_binner_qc", None), "binner_qc")
    add(getattr(args, "assembly_binner_combine", None), "binner_combine")
    add(getattr(args, "assembly_mag_taxonomy", None), "mag_taxonomy")
    add(getattr(args, "assembly_aligner", None), "aligner")
    add(getattr(args, "assembly_mag_quantifier", None), "mag_quantifier")
    return found


def _should_offer(name: str, group: str) -> bool:
    token = str(name or "").strip()
    if not token or token.lower() in _SKIP:
        return False
    if _ask_spec(token) is not None:
        return True
    key = _norm(token)
    if group == "annotator":
        return key not in _PYTHON_ANNOTATORS
    if group == "reads_generator":
        return resolve_reads_generator(token)[0] == "custom"
    if group == "metagenome_generator":
        return resolve_metagenome_generator(token)[0] == "custom"
    if group == "table_reads_generator":
        return resolve_regeneration_mode(token)[0] == "custom"
    if group == "reprofiler":
        return resolve_reprofiler(token)[0] == "custom"
    if group == "export":
        return resolve_export(token)[0] == "custom"
    if group == "qc":
        canon = canonical_qc_name(token)
        return bool(canon) and not is_identity_qc(canon)
    if group == "table_scoring":
        return key not in {_norm(item) for item in BUILTIN_TABLE_SCORERS}
    if group == "sample_scoring":
        return key not in {_norm(item) for item in BUILTIN_SAMPLE_SCORERS}
    if group == "sample_filtering":
        return key not in {_norm(item) for item in BUILTIN_FILTERS}
    if group == "feature_importance":
        return token.lower() not in FEATURE_IMPORTANCE_NONE and key not in {
            _norm(item) for item in FEATURE_IMPORTANCE_BUILTINS
        }
    return True


def _offer_program(
    program: Tuple[str, str, str],
    *,
    mode: str,
    ask: Callable[[str], str],
    interactive: bool,
    ask_tool_fn: Callable,
    import_tool_fn: Callable,
    configured_fn: Callable[..., object],
    which_fn: Callable[[Sequence[str]], Sequence[str]],
    resolve_conda_fn: Callable[[str, str], Tuple[str, str]],
    sibling_fn: Callable[[str, str], str],
) -> None:
    name, group, explicit = program
    spec = _ask_spec(name) or AskTool(name, (name,), group, "")
    if configured_fn(spec.name, spec.which):
        return
    if mode == "auto":
        found = _auto_program_path(spec, explicit, which_fn, resolve_conda_fn)
        if found:
            path, env = found
            _record_program(
                spec,
                path,
                env=env,
                import_tool_fn=import_tool_fn,
                sibling_fn=sibling_fn,
            )
            return
        if not interactive:
            print(f"  {spec.name}: not imported (no executable found to record)")
            return
    elif explicit and interactive:
        if _confirm_path(ask, spec.name, explicit):
            _record_program(
                spec,
                explicit,
                env="",
                import_tool_fn=import_tool_fn,
                sibling_fn=sibling_fn,
            )
            return
    ask_tool_fn(spec, ask=ask, prefer_config=False, check=False)


def _auto_program_path(
    spec: AskTool,
    explicit: str,
    which_fn: Callable[[Sequence[str]], Sequence[str]],
    resolve_conda_fn: Callable[[str, str], Tuple[str, str]],
) -> Optional[Tuple[str, str]]:
    """One executable we can record without a question: ``(path, env)``."""
    if explicit:
        env, _prefix = conda_layout(explicit)
        return explicit, env
    hits = [str(path) for path in which_fn(spec.which) if str(path).strip()]
    if len(hits) == 1:
        env, _prefix = conda_layout(hits[0])
        return hits[0], env
    if len(hits) > 1:
        return None
    exe = spec.which[0] if spec.which else spec.name
    try:
        path, _prefix = resolve_conda_fn(spec.name, exe)
    except (OSError, FileNotFoundError, ValueError):
        return None
    if path:
        return path, "conda"
    return None


def _record_program(
    spec: AskTool,
    path: str,
    *,
    env: str,
    import_tool_fn: Callable,
    sibling_fn: Callable[[str, str], str],
) -> None:
    import_tool_fn(
        name=spec.name,
        tool_type=spec.type,
        env=env,
        exec_name=Path(path).name,
        exec_path=path,
    )
    print(f"  recorded tools.{spec.name} -> {path}")
    for companion in spec.companions:
        sibling = sibling_fn(path, companion)
        if not sibling:
            continue
        import_tool_fn(
            name=companion,
            tool_type=spec.type,
            env=env,
            exec_name=companion,
            exec_path=sibling,
        )
        print(f"  recorded tools.{companion} -> {sibling}")


def _confirm_path(ask: Callable[[str], str], name: str, path: str) -> bool:
    return ask_yes_no(
        ask,
        f"Import {name} from {path}? [Y/n] ",
        default=True,
        detail=(
            f"Import {name} from {path}? Answer y or n.\n"
            "  y     record this executable in the install config\n"
            "  n     leave it unregistered\n"
            "> "
        ),
    )


def _offer_database(
    tool: str,
    name: str,
    *,
    mode: str,
    ask: Callable[[str], str],
    interactive: bool,
    import_db_fn: Callable,
    db_known_fn: Callable[[str, str], bool],
) -> None:
    if not _database_token(name):
        return
    catalog = _catalog_name(name)
    if db_known_fn(tool, catalog) or (catalog != name and db_known_fn(tool, name)):
        return
    existing = _existing_path(name)
    if mode == "auto" and existing:
        _record_database(import_db_fn, tool, catalog, existing)
        return
    if mode == "ask" and existing:
        if _confirm_path(ask, f"database {catalog} for {tool}", existing):
            _record_database(import_db_fn, tool, catalog, existing)
        return
    if not interactive:
        print(f"  {tool} database {catalog}: not imported (no path found to record)")
        return
    short = f"Import database {catalog} for {tool}? PATH | no\n> "
    long = (
        f"Import of database {catalog} for {tool} is not done yet.\n"
        "Reply with one line:\n"
        "  /path/to/db\n"
        "      directory or index prefix to record\n"
        "  no\n"
        "      leave it unregistered\n"
        "> "
    )
    prompt = short
    while True:
        reply = str(ask(prompt) or "").strip()
        if reply.lower() in {"n", "no"}:
            return
        candidate = Path(reply).expanduser() if reply else None
        try:
            exists = bool(candidate and candidate.exists())
        except OSError:
            exists = False
        if reply and exists:
            _record_database(import_db_fn, tool, catalog, str(candidate))
            return
        prompt = long


def _record_database(import_db_fn: Callable, tool: str, name: str, path: str) -> None:
    import_db_fn(name=name, tool=tool, exec_path=path)
    print(f"  recorded databases.{tool}.{name} -> {Path(path).expanduser().resolve()}")


def _database_known(tool: str, name: str) -> bool:
    from samovar.db_spec import lookup_database_record
    from samovar.paths import load_config

    return lookup_database_record(load_config(), tool, name) is not None


def _database_token(token: str) -> bool:
    text = str(token or "").strip()
    return bool(text) and text != "." and not text.startswith("-")


def _existing_path(token: str) -> str:
    candidate = Path(str(token or "")).expanduser()
    try:
        if candidate.exists():
            return str(candidate)
    except OSError:
        return ""
    return ""


def _catalog_name(token: str) -> str:
    existing = _existing_path(token)
    if existing:
        return Path(existing).name
    return str(token or "").strip()


def _ask_spec(name: str) -> Optional[AskTool]:
    key = _norm(name)
    if not key:
        return None
    for spec in ASK_TOOLS:
        names = {_norm(spec.name), *(_norm(item) for item in spec.which)}
        if key in names:
            return spec
    return None


def _tokens(value: object) -> Sequence[str]:
    if value in (None, False, ""):
        return []
    if isinstance(value, (list, tuple)):
        out: List[str] = []
        for item in value:
            out.extend(_tokens(item))
        return out
    return [str(value).strip()] if str(value).strip() else []


def _after_colon(token: str) -> str:
    text = str(token or "").strip()
    if ":" in text:
        return text.split(":", 1)[1].strip()
    return text


def _norm(value: str) -> str:
    return str(value or "").strip().lower().replace("-", "_")
