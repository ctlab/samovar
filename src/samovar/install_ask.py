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
class AskTool:
    name: str
    which: Tuple[str, ...]
    type: str
    install: str
    # Executables that ship next to ``which[0]`` (kraken2-build beside kraken2).
    companions: Tuple[str, ...] = ()
    # ``core`` is ``./install.sh ask``. ``full`` is only ``./install.sh full ask``.
    scope: str = "full"


# One row per program the pipeline can launch. ``companions`` are not separate
# dialogues: after this row is imported they are offered yes/no from the same
# directory and the same conda/base answer.
# ``scope="core"`` is the obligatory ``./install.sh ask`` set. Every other row
# is asked only by ``./install.sh full ask``. ``install=""`` is import-only.
ASK_TOOLS: Tuple[AskTool, ...] = (
    AskTool(
        "kraken2", ("kraken2",), "annotator", "conda:kraken2",
        companions=("kraken2-build", "kraken2-inspect"),
        scope="core",
    ),
    AskTool(
        "kaiju", ("kaiju",), "annotator", "conda:kaiju",
        companions=("kaiju-mkbwt", "kaiju-mkfmi"),
        scope="core",
    ),
    AskTool(
        "kraken", ("kraken",), "annotator", "conda:kraken",
        companions=("kraken-build",),
        scope="core",
    ),
    AskTool(
        "krakenuniq", ("krakenuniq",), "annotator", "conda:krakenuniq",
        companions=("krakenuniq-build", "jellyfish"),
        scope="core",
    ),
    AskTool("metaphlan", ("metaphlan", "metaphlan4"), "annotator", ""),
    AskTool("centrifuge", ("centrifuge",), "annotator", "conda:centrifuge"),
    AskTool("metax", ("metax",), "annotator", "optional:metax"),
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


def is_yes(answer: str, *, default: bool) -> bool:
    text = str(answer or "").strip().lower()
    if not text:
        return default
    if text in {"y", "yes"}:
        return True
    if text in {"n", "no"}:
        return False
    return default


def which_prompt(spec: AskTool, hits: Sequence[str]) -> str:
    lines = [f"{spec.name}: found on PATH:"]
    for index, path in enumerate(hits, start=1):
        lines.append(f"  {index}) {path}")
    if len(hits) == 1:
        lines.append(f"Import {spec.name} from 1? [Y/n] ")
    else:
        lines.append("Import which number? [1/n] ")
    return "\n".join(lines)


def env_prompt(spec: AskTool) -> str:
    exe = spec.which[0] if spec.which else spec.name
    return (
        f"Import {spec.name} from an environment?\n"
        f"  conda ENV {exe}\n"
        f"  base {exe}\n"
        f"  n\n"
        f"> "
    )


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
    import_tool(
        name=spec.name,
        tool_type=spec.type,
        env=env,
        exec_name=exec_name or spec.name,
        exec_path=path,
    )
    print(f"  recorded tools.{spec.name} -> {path}")


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
    """Find the binary after an install: PATH, active env, then the sidecar prefix."""
    hits = which_all(spec.which)
    active = os.environ.get("CONDA_PREFIX", "").strip()
    if hits:
        prefix = ""
        parent = Path(hits[0]).parent
        if parent.name == "bin":
            prefix = str(parent.parent)
        return hits[0], prefix or active
    kind, _, arg = spec.install.partition(":")
    if kind == "optional" and arg:
        prefix = sidecar_envs_dir() / arg
        for name in spec.which:
            candidate = prefix / "bin" / name
            if candidate.is_file():
                return str(candidate.resolve()), str(prefix)
    if active:
        for name in spec.which:
            try:
                return executable_in_prefix(active, name), active
            except FileNotFoundError:
                continue
    return "", ""


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
) -> None:
    """Yes/no for helper binaries that sit next to the one just imported."""
    for name in spec.companions:
        path = sibling_fn(primary_path, name)
        if not path:
            continue
        prompt = (
            f"Import {name} from the same place as {spec.name} ({path})? [Y/n] "
        )
        if not is_yes(ask(prompt), default=True):
            print(f"  {name}: skipped")
            continue
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
) -> str:
    """Prompt for one tool. Returns ``imported-which``, ``imported-env``, ``installed``, or ``skipped``."""

    def commit(path: str, *, env: str, prefix: str, exec_name: str) -> None:
        if not env:
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
        )

    hits = list(which_fn(spec.which))
    if hits:
        while True:
            picked = pick_which(ask(which_prompt(spec, hits)), hits)
            if picked is None:
                break
            if picked == "":
                print("Answer Y, a number from the list, or n.")
                continue
            commit(picked, env="", prefix="", exec_name=Path(picked).name)
            return "imported-which"

    while True:
        parsed = parse_env_answer(ask(env_prompt(spec)))
        if parsed is None:
            print("Answer `conda ENV EXEC`, `base EXEC`, or n/no.")
            continue
        if parsed == "no":
            break
        kind, env_name, exe = parsed
        try:
            if kind == "base":
                path, prefix = resolve_base_fn(exe)
                commit(path, env="", prefix=prefix, exec_name=exe)
            else:
                path, prefix = resolve_conda_fn(env_name, exe)
                commit(path, env="conda", prefix=prefix, exec_name=exe)
        except (OSError, FileNotFoundError, ValueError) as exc:
            print(f"  {exc}")
            continue
        return "imported-env"

    if not spec.install:
        print(f"  {spec.name}: skipped (optional; not installed automatically)")
        return "skipped"
    if is_yes(ask(f"Install {spec.name} into the samovar environment? [y/N] "), default=False):
        try:
            install(spec)
        except (OSError, subprocess.CalledProcessError, ValueError) as exc:
            print(f"  {spec.name}: install failed: {exc}")
            return "skipped"
        path, prefix = locate(spec)
        if path:
            env = "conda" if prefix else ""
            commit(path, env=env, prefix=prefix, exec_name=Path(path).name)
            return "installed"
        print(f"  {spec.name}: installer finished; executable not found to record.")
        return "installed"
    print(f"  {spec.name}: skipped")
    return "skipped"


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
    """Short success line, or the error lines from a failed toy example."""
    lines = [line for line in str(output or "").splitlines() if line.strip()]
    if code == 0:
        return (
            "Toy example passed.\n"
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
        f"Toy example failed (exit {code}).\n"
        f"  output: {out_dir}\n"
        f"  log: {log_path}\n"
        f"{body}"
    )


def run_toy_example(root: Optional[Path] = None) -> Tuple[int, str]:
    """Run ``examples/toy/pipeline.sh`` and return ``(exit_code, report)``."""
    root = Path(root) if root is not None else repo_root()
    script = root / "examples" / "toy" / "pipeline.sh"
    out_dir = root / "examples_outdir" / "toy"
    log_path = out_dir / "ask-smoketest.log"
    if not script.is_file():
        report = format_smoketest_report(
            1, f"error: toy example script not found: {script}", log_path, out_dir
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


def ask_smoketest(
    *,
    ask: Callable[[str], str],
    run: Callable[[], Tuple[int, str]] = run_toy_example,
) -> str:
    """Ask ``run smoketest?``. Yes runs the toy example and prints its report."""
    if not is_yes(ask("run smoketest? [y/N] "), default=False):
        print("Smoketest skipped.")
        return "skipped"
    code, _report = run()
    return "passed" if code == 0 else "failed"


def main(argv: Optional[Sequence[str]] = None) -> int:
    full = ask_full_requested(argv)
    if not prompts_enabled():
        print("install.sh ask: not a terminal; skipping tool prompts.")
        return 0
    tools = tools_for_ask(full=full)
    if full:
        print("Full ask: every supported program (anvio, centrifuge, gtdbtk, metax, and the rest).")
    else:
        print("Ask: obligatory annotators. ./install.sh full ask also asks about every supported program.")
    for spec in tools:
        print("")
        try:
            ask_tool(spec, ask=input)
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
