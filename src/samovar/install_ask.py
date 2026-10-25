"""``./install.sh ask`` — import each important tool into the main config.

One row in :data:`ASK_TOOLS` is one prompt. Columns:

* ``name`` — key written under ``tools.<name>``
* ``which`` — executables to look up on ``PATH``
* ``type`` — tool group (``annotator``, ``scoring``, …)
* ``install`` — what “install into the samovar environment” runs:
  ``conda:<package>`` (bioconda into the active env) or
  ``optional:<install.sh name>`` (the existing installer, including sidecars)
"""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Tuple, Union

from samovar.paths import repo_root, sidecar_envs_dir
from samovar.sidecar import conda_executable
from samovar.tools_import import import_tool

EnvAnswer = Union[str, Tuple[str, str, str]]


@dataclass(frozen=True)
class ConfiguredImport:
    """A tool already stored in the main config, including its environment."""

    path: str
    env: str
    prefix: str
    exec_name: str
    answer: str


@dataclass(frozen=True)
class AskTool:
    name: str
    which: Tuple[str, ...]
    type: str
    install: str
    # Executables that ship next to ``which[0]`` (kraken2-build beside kraken2).
    companions: Tuple[str, ...] = ()
    # ``core`` is ``./install.sh ask``. ``full`` is only ``./install.sh full ask``.
    scope: str = "full"
    # Argument lists tried before the generic ``--help`` / ``-h`` launch check.
    probes: Tuple[Tuple[str, ...], ...] = ()


# One row per program the pipeline can launch. ``companions`` are not separate
# dialogues: after this row is imported they are offered yes/no from the same
# directory and the same conda/base answer.
# ``scope="core"`` is the obligatory ``./install.sh ask`` set: Kraken2 and Kaiju.
# Kraken, KrakenUniq, and every other row are asked only by ``./install.sh full ask``.
# ``install=""`` is import-only.
ASK_TOOLS: Tuple[AskTool, ...] = (
    AskTool(
        "kraken2", ("kraken2",), "annotator", "optional:kraken2",
        companions=("kraken2-build", "kraken2-inspect"),
        scope="core",
    ),
    AskTool(
        "kaiju", ("kaiju",), "annotator", "optional:kaiju",
        companions=("kaiju-mkbwt", "kaiju-mkfmi", "kaiju-multi"),
        scope="core",
    ),
    AskTool(
        "kraken", ("kraken",), "annotator", "conda:kraken",
        companions=("kraken-build",),
    ),
    AskTool(
        "krakenuniq", ("krakenuniq",), "annotator", "conda:krakenuniq",
        companions=("krakenuniq-build", "jellyfish"),
    ),
    AskTool("metaphlan", ("metaphlan", "metaphlan4"), "annotator", ""),
    AskTool("centrifuge", ("centrifuge",), "annotator", "conda:centrifuge"),
    AskTool(
        "metax", ("metax",), "annotator", "optional:metax",
        probes=(("profile", "-h"), ("-h",)),
    ),
    AskTool(
        "bowtie2", ("bowtie2",), "aligner", "conda:bowtie2",
        companions=("bowtie2-build",),
    ),
    AskTool("fastp", ("fastp",), "qc", "conda:fastp"),
    AskTool("cutadapt", ("cutadapt",), "qc", "conda:cutadapt"),
    AskTool("trimmomatic", ("trimmomatic",), "qc", "conda:trimmomatic"),
    AskTool("chopper", ("chopper",), "qc", "conda:chopper"),
    AskTool("nanofilt", ("NanoFilt",), "qc", "conda:nanofilt"),
    AskTool("opal", ("opal.py",), "scoring", "optional:opal"),
    AskTool("multiqc", ("multiqc",), "scoring", "optional:multiqc"),
    AskTool("camisim", ("camisim",), "metagenome_generator", "optional:camisim"),
    AskTool("nanosim", ("simulator.py",), "metagenome_generator", "optional:nanosim"),
    AskTool("art", ("art_illumina",), "reads_generator", "optional:art"),
    AskTool("seqtk", ("seqtk",), "reads_generator", "optional:seqtk"),
    AskTool("nextflow", ("nextflow",), "workflow", "optional:nextflow"),
    AskTool("wgsim", ("wgsim",), "reads_generator", "conda:wgsim"),
    AskTool("sparsedossa2", ("sparsedossa2",), "table_reads_generator", "optional:sparsedossa2"),
    AskTool("R", ("R",), "runtime", "optional:r-package"),
    AskTool("megahit", ("megahit",), "assembler", "optional:megahit"),
    AskTool("spades", ("spades.py",), "assembler", "optional:spades"),
    AskTool("flye", ("flye",), "assembler", "optional:flye"),
    AskTool("prodigal", ("prodigal",), "gene_caller", "optional:prodigal"),
    AskTool("minimap2", ("minimap2",), "aligner", "optional:minimap2"),
    AskTool("samtools", ("samtools",), "aligner", "conda:samtools"),
    AskTool("coverm", ("coverm",), "mag_quantifier", "optional:coverm"),
    AskTool("dastool", ("DAS_Tool",), "binner_combine", "optional:dastool"),
    AskTool("metabat2", ("metabat2",), "binner", "conda:metabat2"),
    AskTool(
        "anvio", ("anvi-gen-contigs-database",), "binner", "optional:anvio",
        companions=(
            "anvi-profile",
            "anvi-merge",
            "anvi-cluster-contigs",
            "anvi-summarize",
        ),
    ),
    AskTool("checkm2", ("checkm2",), "binner_qc", "optional:checkm2"),
    AskTool("gtdbtk", ("gtdbtk",), "mag_taxonomy", "optional:gtdbtk"),
    AskTool("cuml", ("cuml",), "reprofiler", "optional:cuml"),
)


def tools_for_ask(*, full: bool) -> Tuple[AskTool, ...]:
    """Obligatory rows, or every supported program when ``full`` is set."""
    if full:
        return ASK_TOOLS
    return tuple(tool for tool in ASK_TOOLS if tool.scope == "core")


def ask_full_requested(argv: Optional[Sequence[str]] = None) -> bool:
    """True for ``./install.sh full ask`` (env or ``--full``)."""
    if argv and "--full" in argv:
        return True
    flag = os.environ.get("SAMOVAR_INSTALL_ASK_FULL", "").strip().lower()
    return flag in {"1", "true", "yes"}


def launch_check_enabled(argv: Optional[Sequence[str]] = None, *, full: bool) -> bool:
    """Full ask runs each imported binary unless ``--no-check`` is set."""
    if not full:
        return False
    args = list(argv if argv is not None else sys.argv[1:])
    if "--no-check" in args:
        return False
    flag = os.environ.get("SAMOVAR_INSTALL_ASK_NO_CHECK", "").strip().lower()
    return flag not in {"1", "true", "yes"}


def parse_env_answer(text: str) -> Optional[EnvAnswer]:
    """``conda ENV EXEC``, ``base EXEC``, or ``n``/``no``. Anything else is invalid."""
    parts = str(text or "").split()
    if not parts:
        return None
    head = parts[0].lower()
    if len(parts) == 1 and head in {"n", "no"}:
        return "no"
    if len(parts) == 2 and head == "base" and parts[1]:
        return ("base", "", parts[1])
    if len(parts) == 3 and head == "conda" and parts[1] and parts[2]:
        return ("conda", parts[1], parts[2])
    return None


def which_all(names: Sequence[str]) -> List[str]:
    """Every executable on ``PATH`` for these names, in PATH order."""
    found: List[str] = []
    for name in names:
        binary = str(name or "").strip()
        if not binary:
            continue
        for directory in os.environ.get("PATH", "").split(os.pathsep):
            if not directory:
                continue
            candidate = Path(directory) / binary
            try:
                ok = candidate.is_file() and os.access(candidate, os.X_OK)
            except OSError:
                ok = False
            if not ok:
                continue
            try:
                resolved = str(candidate.resolve())
            except OSError:
                resolved = str(candidate)
            if resolved not in found:
                found.append(resolved)
    return found


def pick_which(answer: str, hits: Sequence[str]) -> Optional[str]:
    """``Y``/empty/``1`` picks the first hit, a number picks that hit, ``n`` declines.

    Returns ``""`` when the answer is not one of those.
    """
    text = str(answer or "").strip().lower()
    if not hits:
        return None
    if text in {"", "y", "yes", "1"}:
        return hits[0]
    if text in {"n", "no"}:
        return None
    if text.isdigit():
        index = int(text) - 1
        if 0 <= index < len(hits):
            return hits[index]
    return ""


def is_yes(answer: str, *, default: bool) -> Optional[bool]:
    """``True`` / ``False`` for y/n. Empty uses ``default``. Anything else is ``None``."""
    text = str(answer or "").strip().lower()
    if not text:
        return default
    if text in {"y", "yes"}:
        return True
    if text in {"n", "no"}:
        return False
    return None


def ask_yes_no(
    ask: Callable[[str], str],
    prompt: str,
    *,
    default: bool,
    detail: str = "",
) -> bool:
    """Repeat until the answer is y or n. An empty line uses ``default``.

    ``detail`` is the prompt used after the first answer is neither y nor n.
    """
    current = prompt
    while True:
        answer = ask_yes_no_once(ask(current), default=default)
        if answer is not None:
            return answer
        current = detail or prompt


def ask_yes_no_once(answer: str, *, default: bool) -> Optional[bool]:
    return is_yes(answer, default=default)


def which_prompt(spec: AskTool, hits: Sequence[str], *, detailed: bool = False) -> str:
    lines = [f"{spec.name}: found on PATH:"]
    for index, path in enumerate(hits, start=1):
        lines.append(f"  {index}) {path}")
    if not detailed:
        if len(hits) == 1:
            lines.append("Import 1? [Y/n] ")
        else:
            lines.append("Number, or n: ")
        return "\n".join(lines)
    lines.append("The previous answer did not import a copy that runs.")
    lines.append("  1 or y    import the first path")
    lines.append("  2         import that number, when it is listed")
    lines.append("  n         do not import any of these copies")
    lines.append("> ")
    return "\n".join(lines)


def _exe_name(spec: AskTool) -> str:
    return spec.which[0] if spec.which else spec.name


def env_prompt_short(spec: AskTool, default: str = "") -> str:
    exe = _exe_name(spec)
    if default:
        return f"Import {spec.name}? [{default}] "
    return f"Import {spec.name}? conda ENV {exe} | base {exe} | no\n> "


def env_prompt(spec: AskTool, default: str = "") -> str:
    """Long form, shown after an import attempt did not finish."""
    exe = _exe_name(spec)
    env = spec.name
    shown = f"\nEmpty keeps [{default}]." if default else ""
    return (
        f"Import of {spec.name} is not done yet.{shown}\n"
        "Reply with one line:\n"
        f"  conda {env} {exe}\n"
        f"      conda environment `{env}`, executable `{exe}`\n"
        f"  base {exe}\n"
        f"      the active environment, executable `{exe}`\n"
        f"  no\n"
        f"      do not import; the next question is whether to install it\n"
        f"> "
    )


def install_prompt(spec: AskTool, *, detailed: bool = False) -> str:
    if not detailed:
        return f"Install {spec.name} using samovar? [y/N] "
    return (
        f"Install {spec.name} using samovar? Answer y or n.\n"
        f"  y     install it into the samovar environment and record it\n"
        f"  n     skip {spec.name}\n"
        f"> "
    )


def config_keep_prompt(spec: AskTool, current: ConfiguredImport) -> str:
    label = current.answer or current.path
    return (
        f"{spec.name}: already in config as `{label}`\n"
        f"  {current.path}\n"
        f"Keep this? [Y/n] "
    )


def program_launches(
    path: str,
    *,
    conda_env: str = "",
    timeout: float = 12,
    probes: Sequence[Sequence[str]] = (),
) -> Tuple[bool, str]:
    """True when the binary prints something or exits 0.

    Tries tool-specific ``probes``, then ``--help``, ``-h``, ``--version``, then no
    arguments. A hang is a failure. ``conda_env`` runs the binary inside that
    environment so its libraries resolve.
    """
    attempts = [list(item) for item in probes]
    for generic in (["--help"], ["-h"], ["--version"], []):
        if generic not in attempts:
            attempts.append(generic)
    notes: List[str] = []
    conda = conda_executable() if conda_env else ""
    for args in attempts:
        if conda_env and conda:
            cmd = [conda, "run", "-n", conda_env, path, *args]
        else:
            cmd = [path, *args]
        label = " ".join(args) if args else "(no args)"
        try:
            proc = subprocess.run(
                cmd,
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
                stdin=subprocess.DEVNULL,
            )
        except subprocess.TimeoutExpired:
            notes.append(f"{label} timed out")
            continue
        except OSError as exc:
            return False, str(exc)
        text = "\n".join(part for part in (proc.stdout, proc.stderr) if part).strip()
        if text or proc.returncode == 0:
            line = text.splitlines()[0].strip() if text else f"exit {proc.returncode}"
            return True, line[:180]
        notes.append(f"{label} exited {proc.returncode} with no output")
    return False, "; ".join(notes) or "no output"


def record_tool(
    spec: AskTool,
    path: str,
    *,
    env: str = "",
    prefix: str = "",
    exec_name: str = "",
) -> None:
    """Write ``tools.<name>`` in the main config.

    ``prefix`` is the conda env that contains ``path``. A path under ``<prefix>/bin``
    is what the config exposes as ``tool_envs.<name>``.
    """
    del prefix  # recorded via the executable path (…/bin/<exe>)
    names = [spec.name]
    for alias in spec.which:
        if alias and alias not in names:
            names.append(alias)
    for name in names:
        import_tool(
            name=name,
            tool_type=spec.type,
            env=env,
            exec_name=exec_name or name,
            exec_path=path,
        )
        print(f"  recorded tools.{name} -> {path}")


def conda_env_prefix(env_name: str) -> str:
    conda = conda_executable()
    if not conda:
        raise FileNotFoundError("conda/mamba was not found; cannot resolve a conda env")
    proc = subprocess.run(
        [conda, "run", "-n", env_name, "printenv", "CONDA_PREFIX"],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        raise FileNotFoundError(f"conda env {env_name!r} was not found. {detail}".rstrip())
    lines = [line.strip() for line in (proc.stdout or "").splitlines() if line.strip()]
    if not lines or not Path(lines[-1]).is_dir():
        raise FileNotFoundError(f"conda env {env_name!r} has no prefix")
    return lines[-1]


def executable_in_prefix(prefix: str, exe: str) -> str:
    candidate = Path(prefix) / "bin" / exe
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate.resolve())
    raise FileNotFoundError(f"{exe} is not in {prefix}/bin")


def resolve_base(exe: str) -> Tuple[str, str]:
    prefix = os.environ.get("CONDA_PREFIX", "").strip()
    if prefix:
        try:
            return executable_in_prefix(prefix, exe), prefix
        except FileNotFoundError:
            pass
    for path in which_all((exe,)):
        return path, ""
    raise FileNotFoundError(f"{exe} is not in the active environment or on PATH")


def resolve_conda(env_name: str, exe: str) -> Tuple[str, str]:
    prefix = conda_env_prefix(env_name)
    return executable_in_prefix(prefix, exe), prefix


def install_action(spec: AskTool) -> None:
    """Run the row's install action (conda package or ``./install.sh <name>``)."""
    kind, _, arg = spec.install.partition(":")
    if not arg:
        raise ValueError(f"{spec.name}: install action must be conda:PKG or optional:NAME")
    if kind == "conda":
        conda = conda_executable()
        if not conda:
            raise FileNotFoundError("conda/mamba was not found")
        cmd = [conda, "install", "-y", "-c", "conda-forge", "-c", "bioconda", arg]
        if os.environ.get("CONDA_PREFIX"):
            cmd[2:2] = ["-p", os.environ["CONDA_PREFIX"]]
        print(" ", " ".join(cmd))
        subprocess.check_call(cmd)
        return
    if kind == "optional":
        script = repo_root() / "install.sh"
        print(f"  bash {script} {arg}")
        subprocess.check_call(["bash", str(script), arg])
        return
    raise ValueError(f"{spec.name}: unknown install action {spec.install!r}")


def locate_installed(spec: AskTool) -> Tuple[str, str]:
    """Find the installed sidecar first, then PATH and the active environment."""
    kind, _, arg = spec.install.partition(":")
    if kind == "optional" and arg:
        installed = configured_import(spec.name, spec.which)
        if installed is not None and Path(installed.path).is_file():
            return installed.path, installed.prefix
        prefix = sidecar_envs_dir() / arg
        for name in spec.which:
            candidate = prefix / "bin" / name
            if candidate.is_file():
                return str(candidate.resolve()), str(prefix)
    hits = which_all(spec.which)
    active = os.environ.get("CONDA_PREFIX", "").strip()
    if hits:
        prefix = ""
        parent = Path(hits[0]).parent
        if parent.name == "bin":
            prefix = str(parent.parent)
        return hits[0], prefix or active
    if active:
        for name in spec.which:
            try:
                return executable_in_prefix(active, name), active
            except FileNotFoundError:
                continue
    return "", ""


def conda_env_name(prefix: str) -> str:
    path = Path(prefix)
    if path.name == "bin":
        path = path.parent
    if path.parent.name in {"envs", "env"}:
        return path.name
    return path.name


def describe_configured(path: str, prefix: str, env_field: str, exe: str) -> Tuple[str, str, str]:
    """Return ``(answer, env, prefix)`` for a config record.

    ``answer`` is ``conda ENV EXEC``, ``base EXEC``, or empty when the binary
    is not inside a conda prefix. The active conda env is ``base``.
    """
    if not prefix and path:
        parent = Path(path).parent
        if parent.name == "bin":
            prefix = str(parent.parent)
    active = os.environ.get("CONDA_PREFIX", "").strip()
    if prefix and active and os.path.normpath(prefix) == os.path.normpath(active):
        return f"base {exe}", "", prefix
    condaish = str(env_field or "").lower() == "conda" or (
        bool(prefix) and (Path(prefix) / "conda-meta").is_dir()
    )
    if prefix and condaish:
        return f"conda {conda_env_name(prefix)} {exe}", "conda", prefix
    return "", "", prefix


def configured_import(name: str, which_names: Sequence[str] = ()) -> Optional[ConfiguredImport]:
    """Read ``tools.<name>`` from the main config, including its conda env."""
    from samovar.main_config import iter_tools, parse_tool_entry, tool_env_prefix, tool_path
    from samovar.paths import load_config

    tools = iter_tools(load_config())
    spec = None
    key = ""
    for candidate in (name, *(which_names or ())):
        if candidate in tools:
            spec = tools[candidate]
            key = candidate
            break
    if spec is None:
        return None
    entry = parse_tool_entry(spec, key or name)
    env_field = entry[0]
    raw_path = entry[2]
    resolved = str(tool_path(spec, key or name) or raw_path or "").strip()
    if not resolved:
        return None
    exe = Path(resolved).name
    if not exe or exe in {".", ""}:
        exe = (which_names or (name,))[0]
    prefix = tool_env_prefix(spec, key or name)
    answer, env, prefix = describe_configured(resolved, prefix, env_field, exe)
    return ConfiguredImport(
        path=resolved,
        env=env,
        prefix=prefix,
        exec_name=exe,
        answer=answer,
    )


def conda_layout(path: str) -> Tuple[str, str]:
    """``(<env>, <prefix>)`` when ``path`` is ``<conda-prefix>/bin/<exe>``."""
    parent = Path(path).parent
    if parent.name != "bin":
        return "", ""
    prefix = parent.parent
    if (prefix / "conda-meta").is_dir():
        return "conda", str(prefix)
    return "", ""


def sibling_executable(primary_path: str, name: str) -> str:
    candidate = Path(primary_path).parent / name
    try:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    except OSError:
        return ""
    return ""


def offer_companions(
    spec: AskTool,
    primary_path: str,
    *,
    env: str,
    prefix: str,
    ask: Callable[[str], str],
    record,
    sibling_fn: Callable[[str, str], str] = sibling_executable,
    prefer_config: bool = False,
    configured_fn: Callable[[str, Sequence[str]], Optional[ConfiguredImport]] = configured_import,
    check: bool = False,
    launch_fn: Callable[..., Tuple[bool, str]] = program_launches,
    conda_env: str = "",
) -> None:
    """Yes/no for helper binaries that sit next to the one just imported."""
    for name in spec.companions:
        path = sibling_fn(primary_path, name)
        if not path:
            continue
        note = ""
        if prefer_config:
            current = configured_fn(name, (name,))
            if current is not None and os.path.normpath(current.path) == os.path.normpath(path):
                note = ", already in config"
        prompt = (
            f"Import {name} from the same place as {spec.name}{note} ({path})? [Y/n] "
        )
        if not ask_yes_no(ask, prompt, default=True):
            print(f"  {name}: skipped")
            continue
        if check:
            if launch_fn is program_launches:
                ok, detail = program_launches(path, conda_env=conda_env)
            else:
                ok, detail = launch_fn(path, conda_env=conda_env)
            if not ok:
                print(f"  {name} did not run ({detail}). Not recorded.")
                continue
            print(f"  {name} runs: {detail}")
        companion = AskTool(name, (name,), spec.type, "")
        record(companion, path, env=env, prefix=prefix, exec_name=name)


def ask_tool(
    spec: AskTool,
    *,
    ask: Callable[[str], str],
    which_fn: Callable[[Sequence[str]], Sequence[str]] = which_all,
    resolve_base_fn: Callable[[str], Tuple[str, str]] = resolve_base,
    resolve_conda_fn: Callable[[str, str], Tuple[str, str]] = resolve_conda,
    record=record_tool,
    install=install_action,
    locate=locate_installed,
    sibling_fn: Callable[[str, str], str] = sibling_executable,
    prefer_config: bool = False,
    configured_fn: Callable[[str, Sequence[str]], Optional[ConfiguredImport]] = configured_import,
    check: bool = False,
    launch_fn: Callable[..., Tuple[bool, str]] = program_launches,
) -> str:
    """Prompt for one tool. Returns ``imported-which``, ``imported-env``, ``installed``, or ``skipped``."""

    def _launch(path: str, conda_env: str = "") -> Tuple[bool, str]:
        if launch_fn is program_launches:
            return program_launches(path, conda_env=conda_env, probes=spec.probes)
        return launch_fn(path, conda_env=conda_env)

    def commit(
        path: str,
        *,
        env: str,
        prefix: str,
        exec_name: str,
        detect_conda: bool = True,
        conda_env: str = "",
    ) -> bool:
        if check:
            ok, detail = _launch(path, conda_env)
            if not ok:
                print(f"  {spec.name} did not run ({detail}).")
                return False
            print(f"  {spec.name} runs: {detail}")
        if detect_conda and not env:
            detected_env, detected_prefix = conda_layout(path)
            if detected_env:
                env, prefix = detected_env, detected_prefix
        record(spec, path, env=env, prefix=prefix, exec_name=exec_name)
        offer_companions(
            spec,
            path,
            env=env,
            prefix=prefix,
            ask=ask,
            record=record,
            sibling_fn=sibling_fn,
            prefer_config=prefer_config,
            configured_fn=configured_fn,
            check=check,
            launch_fn=launch_fn,
            conda_env=conda_env,
        )
        return True

    current = configured_fn(spec.name, spec.which) if prefer_config else None
    default_env = ""
    if current is not None and current.answer.startswith(("conda ", "base ")):
        default_env = current.answer
    detailed = False
    if current is not None and current.path and Path(current.path).is_file():
        while True:
            reply = ask(config_keep_prompt(spec, current)).strip().lower()
            if reply in {"", "y", "yes"}:
                conda_env = ""
                if current.answer.startswith("conda "):
                    parts = current.answer.split()
                    conda_env = parts[1] if len(parts) > 1 else ""
                if commit(
                    current.path,
                    env=current.env,
                    prefix=current.prefix,
                    exec_name=current.exec_name,
                    detect_conda=False,
                    conda_env=conda_env,
                ):
                    return "imported-config"
                detailed = True
                break
            if reply in {"n", "no"}:
                break
            print("Answer Y or n.")

    hits = list(which_fn(spec.which))
    which_detailed = detailed
    if hits:
        while True:
            picked = pick_which(ask(which_prompt(spec, hits, detailed=which_detailed)), hits)
            if picked is None:
                break
            if picked == "":
                which_detailed = True
                continue
            if commit(picked, env="", prefix="", exec_name=Path(picked).name):
                return "imported-which"
            which_detailed = True
            detailed = True

    while True:
        prompt = env_prompt(spec, default_env) if detailed else env_prompt_short(spec, default_env)
        reply = ask(prompt)
        if not str(reply or "").strip() and default_env:
            reply = default_env
        parsed = parse_env_answer(reply)
        if parsed is None:
            detailed = True
            continue
        if parsed == "no":
            break
        kind, env_name, exe = parsed
        try:
            if kind == "base":
                path, prefix = resolve_base_fn(exe)
                recorded = commit(path, env="", prefix=prefix, exec_name=exe)
            else:
                path, prefix = resolve_conda_fn(env_name, exe)
                recorded = commit(
                    path, env="conda", prefix=prefix, exec_name=exe, conda_env=env_name
                )
        except (OSError, FileNotFoundError, ValueError) as exc:
            print(f"  {exc}")
            detailed = True
            continue
        if recorded:
            return "imported-env"
        detailed = True

    if not spec.install:
        print(f"  {spec.name}: skipped (optional; not installed automatically)")
        return "skipped"
    install_detailed = False
    while True:
        if not ask_yes_no(
            ask,
            install_prompt(spec, detailed=install_detailed),
            default=False,
            detail=install_prompt(spec, detailed=True),
        ):
            print(f"  {spec.name}: skipped")
            return "skipped"
        try:
            install(spec)
        except (OSError, subprocess.CalledProcessError, ValueError) as exc:
            print(f"  {spec.name}: install failed: {exc}")
            print(f"  {spec.name} is not imported yet. Answer the install question again, or n to skip.")
            continue
        path, prefix = locate(spec)
        if not path:
            print(f"  {spec.name}: installer finished, but the executable was not found.")
            print(f"  {spec.name} is not imported yet. Answer the install question again, or n to skip.")
            continue
        env = "conda" if prefix else ""
        if commit(path, env=env, prefix=prefix, exec_name=Path(path).name):
            return "installed"
        install_detailed = True
        print(f"  {spec.name} is not imported yet. Answer the install question again, or n to skip.")


def prompts_enabled() -> bool:
    if os.environ.get("CI", "").strip():
        return False
    try:
        return sys.stdin.isatty()
    except Exception:
        return False


def format_smoketest_report(
    code: int,
    output: str,
    log_path: Path,
    out_dir: Path,
) -> str:
    """Short success line, or the error lines from a failed smoketest."""
    lines = [line for line in str(output or "").splitlines() if line.strip()]
    if code == 0:
        return (
            "Smoketest passed.\n"
            f"  output: {out_dir}\n"
            f"  log: {log_path}"
        )
    errors = [
        line
        for line in lines
        if any(
            token in line.lower()
            for token in ("error", "traceback", "failed", "exception", "not found")
        )
    ]
    shown = errors[-20:] or lines[-40:]
    body = "\n".join(f"  {line}" for line in shown) or "  (no output)"
    return (
        f"Smoketest failed (exit {code}).\n"
        f"  output: {out_dir}\n"
        f"  log: {log_path}\n"
        f"{body}"
    )


def run_smoketest(root: Optional[Path] = None) -> Tuple[int, str]:
    """Run ``examples/smoketest/pipeline.sh`` and return ``(exit_code, report)``."""
    root = Path(root) if root is not None else repo_root()
    script = root / "examples" / "smoketest" / "pipeline.sh"
    out_dir = root / "examples_outdir" / "smoketest"
    log_path = out_dir / "ask-smoketest.log"
    if not script.is_file():
        report = format_smoketest_report(
            1, f"error: smoketest script not found: {script}", log_path, out_dir
        )
        print(report)
        return 1, report
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Running {script}")
    env = os.environ.copy()
    env["SAMOVAR"] = str(root)
    env.setdefault("SAMOVAR_ROOT", str(root))
    chunks: List[str] = []
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.Popen(
            ["bash", str(script)],
            cwd=str(root),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            handle.write(line)
            chunks.append(line)
        code = proc.wait()
    report = format_smoketest_report(code, "".join(chunks), log_path, out_dir)
    print("")
    print(report)
    return code, report


def default_smoketest_threads(cpu_count: Optional[int] = None) -> int:
    """Three quarters of the CPUs, at least one. Suggested smoketest default."""
    n = cpu_count if cpu_count is not None else (os.cpu_count() or 1)
    return max(1, (int(n) * 3) // 4)


def ask_smoketest_threads(ask: Callable[[str], str]) -> int:
    """Ask for a thread count. An empty line keeps 3/4 of the CPUs."""
    default = default_smoketest_threads()
    prompt = f"examples/smoketest threads [{default}]: "
    detail = (
        f"Enter a positive integer. An empty line keeps {default} "
        f"(3/4 of the CPUs).\n{prompt}"
    )
    current = prompt
    while True:
        text = str(ask(current) or "").strip()
        if not text:
            print(f"Using {default} threads.")
            return default
        if text.isdigit() and int(text) >= 1:
            chosen = int(text)
            print(f"Using {chosen} threads.")
            return chosen
        current = detail


def ask_smoketest(
    *,
    ask: Callable[[str], str],
    run: Callable[[], Tuple[int, str]] = run_smoketest,
) -> str:
    """Ask to run ``examples/smoketest``. Yes asks for threads, then runs it."""
    if not ask_yes_no(ask, "run examples/smoketest? [y/N] ", default=False):
        print("Smoketest skipped.")
        return "skipped"
    threads = ask_smoketest_threads(ask)
    previous = os.environ.get("SAMOVAR_CORES")
    os.environ["SAMOVAR_CORES"] = str(threads)
    try:
        code, _report = run()
    finally:
        if previous is None:
            os.environ.pop("SAMOVAR_CORES", None)
        else:
            os.environ["SAMOVAR_CORES"] = previous
    return "passed" if code == 0 else "failed"


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    full = ask_full_requested(args)
    check = launch_check_enabled(args, full=full)
    if not prompts_enabled():
        print("install.sh ask: not a terminal; skipping tool prompts.")
        return 0
    tools = tools_for_ask(full=full)
    if full:
        print("Full ask: mandatory programs first, then optional ones.")
        sections = (
            ("Mandatory", tuple(tool for tool in tools if tool.scope == "core")),
            ("Optional", tuple(tool for tool in tools if tool.scope != "core")),
        )
    else:
        print("Ask: obligatory annotators. ./install.sh full ask also asks about every supported program.")
        sections = (("", tools),)
    for title, group in sections:
        if title:
            print("")
            print(f"{title}:")
        for spec in group:
            print("")
            try:
                ask_tool(spec, ask=input, prefer_config=full, check=check)
            except KeyboardInterrupt:
                print("")
                return 130
            except EOFError:
                print("")
                return 0
    print("")
    try:
        status = ask_smoketest(ask=input)
    except KeyboardInterrupt:
        print("")
        return 130
    except EOFError:
        print("")
        return 0
    return 0 if status in {"skipped", "passed"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
