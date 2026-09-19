"""Bundled BibTeX for built-in tools, plus harvest on ``./install.sh``.

``cite/citations.json`` maps each ``TOOL_GROUP_BY_NAME`` key to ``*.bib`` files
under ``cite/``. ``python -m samovar.citations rebuild`` refreshes those files
from doi.org, R ``citation()``, and CLI ``--citation`` when available.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from samovar.main_config import TOOL_GROUP_BY_NAME
from samovar.paths import _code_repo_root, user_config_dir
from samovar.tool_spec import parse_citation_refs

DOI_ACCEPT = "application/x-bibtex"
USER_AGENT = "samovar-citations/1.0 (https://github.com/ctlab/samovar)"

# filename (no directory) -> DOI. Rebuild overwrites the matching cite/*.bib.
DOI_RECORDS: Dict[str, str] = {
    "kraken2_wood2019.bib": "10.1186/s13059-019-1891-0",
    "kraken_wood2014.bib": "10.1186/gb-2014-15-3-r46",
    "krakenuniq_breitwieser2018.bib": "10.1186/s13059-018-1568-0",
    "kaiju_menzel2016.bib": "10.1038/ncomms11257",
    "metaphlan_blanco2023.bib": "10.1038/s41587-023-01688-w",
    "metaphlan_beghini2021.bib": "10.7554/eLife.65088",
    "centrifuge_kim2016.bib": "10.1101/gr.210641.116",
    "iss_gourle2019.bib": "10.1093/bioinformatics/bty630",
    "art_huang2012.bib": "10.1093/bioinformatics/btr708",
    "samtools_danecek2021.bib": "10.1093/gigascience/giab008",
    "samtools_li2009.bib": "10.1093/bioinformatics/btp352",
    "camisim_fritz2019.bib": "10.1186/s40168-019-0633-6",
    "nanosim_yang2017.bib": "10.1093/gigascience/gix010",
    "sparsedossa2_ma2021.bib": "10.1371/journal.pcbi.1008913",
    "opal_meyer2019.bib": "10.1186/s13059-019-1646-y",
    "cami_sczyrba2017.bib": "10.1038/nmeth.4458",
    "cami2_meyer2022.bib": "10.1038/s41592-022-01431-4",
    "multiqc_ewels2016.bib": "10.1093/bioinformatics/btw354",
    "fastp_chen2018.bib": "10.1093/bioinformatics/bty560",
    "cutadapt_martin2011.bib": "10.14806/ej.17.1.200",
    "trimmomatic_bolger2014.bib": "10.1093/bioinformatics/btu170",
    "nanopack2_decoster2023.bib": "10.1093/bioinformatics/btad311",
    "nanopack_decoster2018.bib": "10.1093/bioinformatics/bty149",
    "megahit_li2015.bib": "10.1093/bioinformatics/btv033",
    "prodigal_hyatt2010.bib": "10.1186/1471-2105-11-119",
    "anvio_eren2015.bib": "10.7717/peerj.1319",
    "anvio_eren2021.bib": "10.1038/s41564-020-00834-3",
    "metabat2_kang2019.bib": "10.7717/peerj.7359",
    "checkm2_chklovski2023.bib": "10.1038/s41592-023-01940-w",
    "dastool_sieber2018.bib": "10.1038/s41564-018-0171-1",
    "gtdbtk_chaumeil2020.bib": "10.1093/bioinformatics/btz848",
    "gtdbtk2_chaumeil2022.bib": "10.1038/s41587-022-01347-6",
    "minimap2_li2018.bib": "10.1093/bioinformatics/bty191",
    "coverm_aroney2025.bib": "10.1093/bioinformatics/btaf147",
    "snakemake_koster2012.bib": "10.1093/bioinformatics/bts480",
    "snakemake_molder2021.bib": "10.12688/f1000research.29032.2",
    "nextflow_ditommaso2017.bib": "10.1038/nbt.3820",
    "magscot_ruhlemann2022.bib": "10.1093/bioinformatics/btac694",
    "bray_curtis1957.bib": "10.2307/1942268",
}

# Records doi.org does not serve (or that we keep as the project citation).
STATIC_RECORDS: Dict[str, str] = {
    "sklearn_pedregosa2011.bib": """@article{JMLR:v12:pedregosa11a,
  author  = {Fabian Pedregosa and Ga{\\"e}l Varoquaux and Alexandre Gramfort and Vincent Michel and Bertrand Thirion and Olivier Grisel and Mathieu Blondel and Peter Prettenhofer and Ron Weiss and Vincent Dubourg and Jake Vanderplas and Alexandre Passos and David Cournapeau and Matthieu Brucher and Matthieu Perrot and {\\'E}douard Duchesnay},
  title   = {Scikit-learn: Machine Learning in Python},
  journal = {Journal of Machine Learning Research},
  year    = {2011},
  volume  = {12},
  number  = {85},
  pages   = {2825--2830},
  url     = {http://jmlr.org/papers/v12/pedregosa11a.html}
}
""",
    "samovar_chechenina2023.bib": """@article{samovar2023,
  title={Development of in-silico models of metagenomic communities with given properties and a pipeline for their generation},
  author={Chechenina, A. and Vaulin, N. and Ivanov, A. and Ulyantsev, V.},
  journal={Bioinformatics institute 2022/23},
  year={2023},
  pages={22--24}
}
""",
    "r_core_team.bib": """@Manual{R-core,
  title = {R: A Language and Environment for Statistical Computing},
  author = {{R Core Team}},
  organization = {R Foundation for Statistical Computing},
  address = {Vienna, Austria},
  year = {2024},
  url = {https://www.R-project.org/}
}
""",
}

# Canonical tool name -> bundled bib filenames. Aliases share the same files.
TOOL_FILES: Dict[str, List[str]] = {
    "kraken2": ["kraken2_wood2019.bib"],
    "kraken": ["kraken_wood2014.bib"],
    "krakenuniq": ["krakenuniq_breitwieser2018.bib"],
    "kaiju": ["kaiju_menzel2016.bib"],
    "metaphlan": ["metaphlan_blanco2023.bib", "metaphlan_beghini2021.bib"],
    "metaphlan4": ["metaphlan_blanco2023.bib", "metaphlan_beghini2021.bib"],
    "centrifuge": ["centrifuge_kim2016.bib"],
    "iss": ["iss_gourle2019.bib"],
    "art": ["art_huang2012.bib"],
    "art_illumina": ["art_huang2012.bib"],
    "wgsim": ["samtools_li2009.bib", "samtools_danecek2021.bib"],
    "samtools": ["samtools_danecek2021.bib", "samtools_li2009.bib"],
    "camisim": ["camisim_fritz2019.bib", "cami_sczyrba2017.bib", "cami2_meyer2022.bib"],
    "nanosim": ["nanosim_yang2017.bib"],
    "nanosim3": ["nanosim_yang2017.bib"],
    "simulator.py": ["nanosim_yang2017.bib"],
    "sparsedossa2-fit": ["sparsedossa2_ma2021.bib"],
    "sparsedossa2-stool": ["sparsedossa2_ma2021.bib"],
    "sparsedossa2-vaginal": ["sparsedossa2_ma2021.bib"],
    "sparsedossa2-ibd": ["sparsedossa2_ma2021.bib"],
    "sparsedossa2-cv": ["sparsedossa2_ma2021.bib"],
    "bray_curtis": ["bray_curtis1957.bib"],
    "bray-curtis": ["bray_curtis1957.bib"],
    "opal": ["opal_meyer2019.bib", "cami_sczyrba2017.bib", "cami2_meyer2022.bib"],
    "opal.py": ["opal_meyer2019.bib", "cami_sczyrba2017.bib", "cami2_meyer2022.bib"],
    "multiqc": ["multiqc_ewels2016.bib"],
    "fastp": ["fastp_chen2018.bib"],
    "cutadapt": ["cutadapt_martin2011.bib"],
    "trimmomatic": ["trimmomatic_bolger2014.bib"],
    "chopper": ["nanopack2_decoster2023.bib"],
    "nanofilt": ["nanopack_decoster2018.bib"],
    "megahit": ["megahit_li2015.bib"],
    "prodigal": ["prodigal_hyatt2010.bib"],
    "anvio": ["anvio_eren2015.bib", "anvio_eren2021.bib"],
    "metabat2": ["metabat2_kang2019.bib"],
    "checkm2": ["checkm2_chklovski2023.bib"],
    "dastool": ["dastool_sieber2018.bib"],
    "DAS_Tool": ["dastool_sieber2018.bib"],
    "magscot": ["magscot_ruhlemann2022.bib"],
    "gtdbtk": ["gtdbtk_chaumeil2020.bib", "gtdbtk2_chaumeil2022.bib"],
    "minimap2": ["minimap2_li2018.bib"],
    "coverm": ["coverm_aroney2025.bib"],
    "random_forest": ["sklearn_pedregosa2011.bib"],
    "adaboost": ["sklearn_pedregosa2011.bib"],
    "snakemake": ["snakemake_koster2012.bib", "snakemake_molder2021.bib"],
    "nextflow": ["nextflow_ditommaso2017.bib"],
    "assembly": ["samovar_chechenina2023.bib"],
    "assembly_profiling": ["samovar_chechenina2023.bib"],
    "R": ["r_core_team.bib"],
    "Rscript": ["r_core_team.bib"],
}

# Extra names that are built-in but not in TOOL_GROUP_BY_NAME.
EXTRA_TOOLS: Tuple[str, ...] = ("samovar", "sparsedossa2")

R_PACKAGE_TOOLS: Dict[str, str] = {
    "SparseDOSSA2": "sparsedossa2",
    "samovaR": "samovar",
}

CLI_CITATION_FLAGS: Dict[str, Sequence[str]] = {
    "kraken2": ("--citation",),
    "kraken": ("--citation",),
    "krakenuniq": ("--citation",),
    "centrifuge": ("--citation",),
    "multiqc": ("--citation",),
    "snakemake": ("--citation",),
    "cutadapt": ("--citation",),
    "fastp": ("--citation",),
}


def cite_dir(root: Optional[Path] = None) -> Path:
    base = Path(root) if root is not None else _code_repo_root()
    return base / "cite"


def user_cite_dir() -> Path:
    """Writable cite overlay next to the active install config.

    ``$SAMOVAR_CITE`` overrides. Imported tools merge here so they do not
    rewrite the bundled ``<repo>/cite/citations.json``.
    """
    override = os.environ.get("SAMOVAR_CITE", "").strip()
    if override:
        return Path(override).expanduser()
    return user_config_dir() / "cite"


_BIB_ENTRY_START = re.compile(r"@(?!comment\b|string\b|preamble\b)\w+\s*\{", re.IGNORECASE)
_BIB_KEY = re.compile(r"@\w+\s*\{\s*([^,\s}]+)", re.IGNORECASE)
_SAFE_TOKEN = re.compile(r"[^A-Za-z0-9._-]+")


def parse_bibtex_entries(text: str) -> List[str]:
    blob = str(text or "")
    starts = [m.start() for m in _BIB_ENTRY_START.finditer(blob)]
    if not starts:
        stripped = blob.strip()
        return [stripped] if stripped else []
    entries: List[str] = []
    for i, start in enumerate(starts):
        end = starts[i + 1] if i + 1 < len(starts) else len(blob)
        body = blob[start:end].strip()
        if body:
            entries.append(body if body.endswith("\n") else body + "\n")
    return entries


def bibtex_citekey(entry: str) -> str:
    match = _BIB_KEY.search(entry or "")
    if not match:
        return ""
    return _SAFE_TOKEN.sub("_", match.group(1).strip())


def _content_tag(text: str) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:10]


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _read_registry(path: Path) -> Dict[str, List[str]]:
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(data, dict):
        return {}
    out: Dict[str, List[str]] = {}
    for key, value in data.items():
        files = parse_citation_refs(value)
        if files:
            out[str(key)] = files
        else:
            out[str(key)] = []
    return out


def merge_citation_registry(
    updates: Mapping[str, Sequence[str]],
    *,
    dest: Optional[Path] = None,
) -> Dict[str, List[str]]:
    """Merge ``{tool: [file.bib]}`` into the user registry without dropping other tools."""
    path = dest if dest is not None else user_cite_dir() / "citations.json"
    current = _read_registry(path)
    for tool, files in updates.items():
        name = str(tool or "").strip()
        if not name:
            continue
        merged = parse_citation_refs(list(current.get(name) or []) + list(files))
        current[name] = merged
    _atomic_write_text(path, json.dumps(current, indent=2, sort_keys=True) + "\n")
    return current


def _write_bib_file(dest: Path, text: str) -> None:
    body = (text or "").strip()
    if body and not body.endswith("\n"):
        body += "\n"
    if dest.is_file():
        existing = dest.read_text(encoding="utf-8")
        if existing.strip() == body.strip():
            return
    _atomic_write_text(dest, body)


def _link_citation_file(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    src = src.expanduser().resolve()
    if dest.exists() or dest.is_symlink():
        try:
            if dest.resolve() == src:
                return
        except OSError:
            pass
        if dest.is_file() and not dest.is_symlink():
            try:
                if dest.read_bytes() == src.read_bytes():
                    return
            except OSError:
                pass
        dest.unlink()
    try:
        os.symlink(str(src), dest)
    except OSError:
        shutil.copy2(src, dest)


def _store_inline_entry(cite_root: Path, tool: str, entry: str) -> str:
    key = bibtex_citekey(entry) or _content_tag(entry)
    filename = f"{tool}_{key}.bib"
    path = cite_root / filename
    if path.is_file() and path.read_text(encoding="utf-8").strip() != entry.strip():
        filename = f"{tool}_{key}_{_content_tag(entry)}.bib"
        path = cite_root / filename
    _write_bib_file(path, entry)
    return filename


def register_tool_citations(
    tool: str,
    *,
    inline: Sequence[str] = (),
    files: Sequence[str] = (),
    dest_dir: Optional[Path] = None,
) -> List[str]:
    """Write/link BibTeX into the user cite dir and merge ``citations.json``.

    Inline strings are split into one ``.bib`` per ``@entry``. Paths are linked
    (symlink, copy fallback) as a single file each. Re-registering the same
    content reuses the same filenames.
    """
    name = str(tool or "").strip()
    if not name:
        raise ValueError("tool name is required for citations")
    cite_root = Path(dest_dir) if dest_dir is not None else user_cite_dir()
    cite_root.mkdir(parents=True, exist_ok=True)
    stored: List[str] = []
    extra_files = [str(p) for p in files]

    for blob in inline:
        text = str(blob or "").strip()
        if not text:
            continue
        if text == "-":
            text = sys.stdin.read()
        elif text.startswith("@") and "{" not in text:
            extra_files.append(text[1:])
            continue
        for entry in parse_bibtex_entries(text):
            stored.append(_store_inline_entry(cite_root, name, entry))

    for raw in extra_files:
        src = Path(str(raw or "").strip()).expanduser()
        if not str(raw or "").strip():
            continue
        if not src.is_file():
            raise FileNotFoundError(f"--bibtex-file not found: {raw}")
        stem = _SAFE_TOKEN.sub("_", src.stem) or _content_tag(src.read_text(encoding="utf-8"))
        suffix = src.suffix if src.suffix in {".bib", ".txt"} else ".bib"
        filename = f"{name}_{stem}{suffix}"
        _link_citation_file(src, cite_root / filename)
        stored.append(filename)

    stored = parse_citation_refs(stored)
    if stored:
        merge_citation_registry({name: stored}, dest=cite_root / "citations.json")
    return stored


def builtin_tool_names() -> List[str]:
    names = list(TOOL_GROUP_BY_NAME)
    for extra in EXTRA_TOOLS:
        if extra not in names:
            names.append(extra)
    return names


def mapping_for_tools() -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for name in builtin_tool_names():
        files = list(TOOL_FILES.get(name, []))
        if name == "samovar":
            files = ["samovar_chechenina2023.bib"]
        if name == "sparsedossa2":
            files = ["sparsedossa2_ma2021.bib"]
        out[name] = files
    return out


def missing_tools(mapping: Optional[Mapping[str, Sequence[str]]] = None) -> List[str]:
    data = mapping if mapping is not None else mapping_for_tools()
    return [name for name, files in data.items() if not files]


def _normalize_bib(text: str) -> str:
    body = (text or "").strip()
    if not body:
        return ""
    return body + "\n"


def fetch_doi_bibtex(doi: str, timeout: float = 30.0) -> str:
    url = "https://doi.org/" + str(doi).strip()
    req = urllib.request.Request(
        url,
        headers={"Accept": DOI_ACCEPT, "User-Agent": USER_AGENT},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return _normalize_bib(resp.read().decode("utf-8", errors="replace"))


def _which(name: str) -> Optional[str]:
    path = shutil.which(name)
    return path


def harvest_r_citation(package: str = "") -> str:
    rscript = _which("Rscript") or _which("R")
    if not rscript:
        return ""
    if Path(rscript).name == "R":
        cmd = [rscript, "--vanilla", "-s", "-e"]
    else:
        cmd = [rscript, "--vanilla", "-e"]
    if package:
        expr = (
            f'if (requireNamespace("{package}", quietly=TRUE)) '
            f'print(citation("{package}"), style="Bibtex")'
        )
    else:
        expr = 'print(citation(), style="Bibtex")'
    try:
        proc = subprocess.run(
            cmd + [expr],
            check=False,
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    if proc.returncode != 0:
        return ""
    lines = [ln for ln in (proc.stdout or "").splitlines() if not ln.startswith(">")]
    body = "\n".join(lines).strip()
    if "@" not in body:
        return ""
    return _normalize_bib(body)


def harvest_cli_citation(binary: str, flags: Sequence[str]) -> str:
    exe = _which(binary)
    if not exe:
        return ""
    try:
        proc = subprocess.run(
            [exe, *flags],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    text = ((proc.stdout or "") + "\n" + (proc.stderr or "")).strip()
    if "@" not in text:
        return ""
    return _normalize_bib(text)


def _write_bib(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_normalize_bib(text), encoding="utf-8")


def _append_unique(mapping: Dict[str, List[str]], tool: str, filename: str) -> None:
    files = mapping.setdefault(tool, [])
    if filename not in files:
        files.append(filename)


def rebuild(
    root: Optional[Path] = None,
    *,
    online: Optional[bool] = None,
    harvest: bool = True,
    timeout: float = 30.0,
) -> Dict[str, List[str]]:
    """Write ``cite/*.bib``, ``cite/citations.json``, and ``cite/missing.json``."""
    dest = cite_dir(root)
    dest.mkdir(parents=True, exist_ok=True)
    if online is None:
        online = os.environ.get("SAMOVAR_OFFLINE", "0") == "0"

    mapping = mapping_for_tools()

    for filename, body in STATIC_RECORDS.items():
        _write_bib(dest / filename, body)

    for filename, doi in DOI_RECORDS.items():
        path = dest / filename
        if online:
            try:
                body = fetch_doi_bibtex(doi, timeout=timeout)
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
                print(f"Warning: doi.org failed for {doi} ({exc}); keeping local {filename}", file=sys.stderr)
                body = path.read_text(encoding="utf-8") if path.is_file() else ""
            if body:
                _write_bib(path, body)
            elif not path.is_file():
                _write_bib(path, f"% failed to fetch https://doi.org/{doi}\n")
        elif not path.is_file():
            _write_bib(path, f"% placeholder; fetch https://doi.org/{doi}\n")

    if not harvest:
        (dest / "citations.json").write_text(
            json.dumps(mapping, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        missing = missing_tools(mapping)
        (dest / "missing.json").write_text(
            json.dumps(missing, indent=2) + "\n",
            encoding="utf-8",
        )
        return mapping

    r_core = harvest_r_citation("")
    if r_core:
        name = "r_core_harvest.bib"
        _write_bib(dest / name, r_core)
        _append_unique(mapping, "R", name)
        _append_unique(mapping, "Rscript", name)

    for pkg, tool in R_PACKAGE_TOOLS.items():
        harvested = harvest_r_citation(pkg)
        if not harvested:
            continue
        name = f"{tool}_r_harvest.bib"
        _write_bib(dest / name, harvested)
        _append_unique(mapping, tool, name)
        if tool == "sparsedossa2":
            for alias in (
                "sparsedossa2-fit",
                "sparsedossa2-stool",
                "sparsedossa2-vaginal",
                "sparsedossa2-ibd",
                "sparsedossa2-cv",
            ):
                _append_unique(mapping, alias, name)

    for binary, flags in CLI_CITATION_FLAGS.items():
        harvested = harvest_cli_citation(binary, flags)
        if not harvested:
            continue
        name = f"{binary}_cli_harvest.bib"
        _write_bib(dest / name, harvested)
        if binary in mapping:
            _append_unique(mapping, binary, name)

    (dest / "citations.json").write_text(
        json.dumps(mapping, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    missing = missing_tools(mapping)
    (dest / "missing.json").write_text(
        json.dumps(missing, indent=2) + "\n",
        encoding="utf-8",
    )
    return mapping


def load_citations(root: Optional[Path] = None) -> Dict[str, List[str]]:
    path = cite_dir(root) / "citations.json"
    if path.is_file():
        bundled = _read_registry(path)
    else:
        bundled = mapping_for_tools()
    overlay = _read_registry(user_cite_dir() / "citations.json")
    merged = dict(bundled)
    for tool, files in overlay.items():
        merged[tool] = parse_citation_refs(list(merged.get(tool) or []) + list(files))
    return merged


USED_CITATIONS_NAME = "used_citations.bib"

_SKIP_PREPARE_TOOLS = frozenset(
    {
        "",
        "none",
        "off",
        "false",
        "0",
        "skip",
        "no",
        "identity",
        "raw",
        "python",
        "python3",
        "bash",
        "R",
        "Rscript",
    }
)

# Extra names to try when looking up BibTeX (canonical tool files).
_CITE_LOOKUP = {
    "ensemble": ("random_forest", "adaboost"),
    "linear": ("random_forest",),
    "builtin": ("random_forest",),
    "native": ("random_forest",),
    "sklearn": ("random_forest",),
    "permutation": ("random_forest",),
    "default": ("random_forest",),
    "auto": ("random_forest",),
    "illumina": ("fastp",),
    "bgi": ("fastp",),
    "mgi": ("fastp",),
    "ont": ("chopper",),
    "nanopore": ("chopper",),
    "nanosim": ("nanosim",),
    "opal.py": ("opal",),
    "art_illumina": ("art",),
    "metaphlan4": ("metaphlan",),
    "simulator.py": ("nanosim",),
    "sparsedossa2": ("sparsedossa2-fit",),
    "sparsedossa2_cv": ("sparsedossa2-cv",),
    "sd2_cv": ("sparsedossa2-cv",),
    "bray_ks": ("bray_curtis",),
    "bray": ("bray_curtis",),
    "bray_curtis_ks": ("bray_curtis",),
    "bray-curtis": ("bray_curtis",),
    "DAS_Tool": ("dastool",),
    "samovar": ("samovar",),
}


def _unique_names(names: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for raw in names:
        name = str(raw or "").strip()
        if not name:
            continue
        key = name.lower()
        if key in _SKIP_PREPARE_TOOLS or key in seen:
            continue
        seen.add(key)
        out.append(name)
    return out


def _add_cmd_token(out: List[str], cmd: str) -> None:
    text = str(cmd or "").strip()
    if not text:
        return
    base = Path(text.split()[0]).name
    if base.endswith(".py"):
        base = base[:-3]
    if base and base.lower() not in _SKIP_PREPARE_TOOLS:
        out.append(base)


def _annotator_tool_names(ann: Any) -> List[str]:
    names: List[str] = []
    typ = str(getattr(ann, "type", "") or "").strip()
    run = str(getattr(ann, "run_name", "") or "").strip()
    if run.endswith("-test"):
        run = run[: -len("-test")]
    elif run.endswith("_test"):
        run = run[: -len("_test")]
    if typ:
        names.append(typ)
    dummy_alias = {
        "dummy",
        "dummy9606",
        "constant9606",
        "constant",
        "random",
        "constant_taxid",
    }
    if run and run.lower() not in dummy_alias and run.lower() != typ.lower():
        names.append(run)
    _add_cmd_token(names, str(getattr(ann, "cmd", "") or ""))
    extra = str(getattr(ann, "extra", "") or "")
    if extra and typ.lower() in {"assembly", "assembly_profiling"}:
        from samovar.assembly_profiling import parse_slot_extra

        slots = parse_slot_extra(extra)
        names.append(str(slots.get("assembler") or ""))
        names.append(str(slots.get("gene_caller") or ""))
        names.extend(str(x) for x in slots.get("binners") or [])
        names.append(str(slots.get("binner_qc") or ""))
        names.append(str(slots.get("binner_combine") or ""))
        names.append(str(slots.get("mag_taxonomy") or ""))
        names.append(str(slots.get("aligner") or ""))
        names.append(str(slots.get("mag_quantifier") or ""))
    return names


def selected_prepare_tools(config: Any) -> List[str]:
    """Tool names actually wired by this ``samovar prepare`` config."""
    names: List[str] = []
    for ann in getattr(config, "annotators", None) or []:
        names.extend(_annotator_tool_names(ann))

    try:
        from samovar.qc import canonical_qc_name

        for raw in (
            getattr(config, "qc", ""),
            getattr(config, "qc_initial", ""),
            getattr(config, "qc_generated", ""),
        ):
            names.append(canonical_qc_name(raw) or str(raw or ""))
        for raw in (getattr(config, "qc_postfix", None) or {}).values():
            names.append(canonical_qc_name(raw) or str(raw or ""))
        for raw in (getattr(config, "qc_tool_flags", None) or {}):
            names.append(canonical_qc_name(raw) or str(raw or ""))
    except Exception:
        names.extend(
            [
                str(getattr(config, "qc", "") or ""),
                str(getattr(config, "qc_initial", "") or ""),
                str(getattr(config, "qc_generated", "") or ""),
            ]
        )

    names.extend(getattr(config, "regeneration_modes", None) or [getattr(config, "regeneration_mode", "")])
    names.append(getattr(config, "table_score", "") or "")
    names.append(getattr(config, "sample_score", "") or "")
    names.extend((getattr(config, "sample_score_by_annotator", None) or {}).values())
    names.extend((getattr(config, "sample_score_by_method", None) or {}).values())
    names.extend((getattr(config, "sample_filter_by_method", None) or {}).values())
    names.append(getattr(config, "sample_filter", "") or "")
    names.append(getattr(config, "reads_generator", "") or "iss")
    names.append(getattr(config, "metagenome_generator", "") or "")

    repro = str(getattr(config, "reprofiler", "") or "ensemble")
    try:
        from samovar.reprofilers import resolve_reprofiler

        _kind, repro = resolve_reprofiler(repro)
    except Exception:
        pass
    names.append(repro)

    fi = str(getattr(config, "feature_importance", "") or "builtin")
    if fi.strip().lower() not in _SKIP_PREPARE_TOOLS:
        names.append(fi)

    export_name = str(getattr(config, "export_corrector", "") or "logistic")
    try:
        from samovar.abundance_correctors import is_skipped_export, require_known_export

        if not is_skipped_export(export_name):
            try:
                names.append(require_known_export(export_name))
            except Exception:
                names.append(export_name)
    except Exception:
        names.append(export_name)

    scoring = getattr(config, "scoring_tools", None)
    if scoring is not None:
        if isinstance(scoring, (list, tuple)):
            names.extend(str(x) for x in scoring)
        else:
            names.extend(str(scoring).replace(",", " ").split())
    else:
        try:
            from samovar.scorers import iter_custom_scoring_names

            names.extend(iter_custom_scoring_names())
        except Exception:
            pass
    names.extend(getattr(config, "scoring_tool_flags", None) or {})

    names.append("snakemake")
    meta = str(getattr(config, "metagenome_generator", "") or "").lower()
    reads = str(getattr(config, "reads_generator", "") or "").lower()
    if meta in {"camisim"} or reads in {"camisim", "art", "wgsim", "hybrid"}:
        names.append("camisim")
    if meta in {"camisim", "nanosim", "hybrid"}:
        names.append("nextflow")

    run_mqc = getattr(config, "run_multiqc", None)
    if run_mqc is True:
        names.append("multiqc")
    elif run_mqc is None:
        try:
            from samovar.paths import discover_multiqc

            if discover_multiqc():
                names.append("multiqc")
        except Exception:
            pass

    return _unique_names(names)


def resolve_citation_file(ref: str) -> Optional[Path]:
    text = str(ref or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    if path.is_file():
        return path
    name = path.name
    for root in (user_cite_dir(), cite_dir()):
        cand = root / name
        try:
            if cand.is_file() or cand.is_symlink():
                return cand
        except OSError:
            continue
    return None


def citation_refs_for_tool(name: str) -> List[str]:
    """Filenames/paths for one tool: record field, then registry, then bundled map."""
    bare = str(name or "").strip()
    if not bare:
        return []
    refs: List[str] = []
    try:
        from samovar.main_config import lookup_tool_record
        from samovar.paths import load_config

        rec = lookup_tool_record(load_config(), bare) or {}
        refs.extend(parse_citation_refs(rec.get("citation")))
    except Exception:
        pass
    mapping = load_citations()
    aliases = [bare, *list(_CITE_LOOKUP.get(bare, ())), *list(_CITE_LOOKUP.get(bare.lower(), ()))]
    if bare.lower() != bare:
        aliases.append(bare.lower())
    seen = set()
    for key in aliases:
        if not key or key in seen:
            continue
        seen.add(key)
        refs.extend(mapping.get(key) or [])
        refs.extend(TOOL_FILES.get(key) or [])
        refs.extend(mapping_for_tools().get(key) or [])
    return parse_citation_refs(refs)


def bibtex_entries_for_refs(refs: Sequence[str]) -> List[str]:
    chunks: List[str] = []
    seen_keys = set()
    seen_body = set()
    for ref in refs:
        path = resolve_citation_file(ref)
        if path is None:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        for entry in parse_bibtex_entries(text):
            body = entry.strip()
            if not body:
                continue
            key = bibtex_citekey(body) or _content_tag(body)
            if key in seen_keys or body in seen_body:
                continue
            seen_keys.add(key)
            seen_body.add(body)
            chunks.append(body if body.endswith("\n") else body + "\n")
    return chunks


def write_used_citations(config: Any, outdir: Optional[Path] = None) -> Path:
    """Write ``used_citations.bib`` for tools selected in this prepare.

    Missing papers warn (``UserWarning``) and do not abort. Rewrite is atomic
    and stable across identical prepares.
    """
    import warnings

    dest_root = Path(outdir) if outdir is not None else Path(getattr(config, "output_dir"))
    dest_root = dest_root.expanduser()
    dest = dest_root / USED_CITATIONS_NAME
    tools = selected_prepare_tools(config)
    missing: List[str] = []
    blocks: List[str] = []
    found_tools: List[str] = []
    for name in tools:
        entries = bibtex_entries_for_refs(citation_refs_for_tool(name))
        if not entries:
            missing.append(name)
            continue
        found_tools.append(name)
        blocks.append(f"% --- {name} ---\n")
        blocks.extend(entries)
        if not blocks[-1].endswith("\n"):
            blocks[-1] += "\n"
        if not blocks[-1].endswith("\n\n"):
            blocks.append("\n")
    header = [
        "% SamovaR used_citations.bib (tools selected at prepare)\n",
        "% Tools: " + (", ".join(tools) if tools else "(none)") + "\n",
    ]
    if missing:
        header.append("% Missing citations: " + ", ".join(missing) + "\n")
        warnings.warn(
            "No BibTeX citation for prepare tool(s): " + ", ".join(missing),
            UserWarning,
            stacklevel=2,
        )
    header.append("\n")
    text = "".join(header + blocks)
    if not text.endswith("\n"):
        text += "\n"
    _atomic_write_text(dest, text)
    configs_copy = dest_root / ".log" / "configs" / USED_CITATIONS_NAME
    try:
        configs_copy.parent.mkdir(parents=True, exist_ok=True)
        if configs_copy.resolve() != dest.resolve():
            _atomic_write_text(configs_copy, text)
    except OSError:
        pass
    return dest


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="SamovaR tool BibTeX harvest")
    sub = parser.add_subparsers(dest="cmd")
    rebuild_p = sub.add_parser("rebuild", help="Refresh cite/*.bib from DOIs and installed tools")
    rebuild_p.add_argument("--root", type=Path, default=None, help="Checkout root (default: package repo)")
    rebuild_p.add_argument(
        "--offline",
        action="store_true",
        help="Do not call doi.org (still try R citation() and CLI --citation)",
    )
    rebuild_p.add_argument(
        "--no-harvest",
        action="store_true",
        help="Only write bundled DOI/static records (skip R citation() and CLI --citation)",
    )
    list_p = sub.add_parser("missing", help="Print tools with no BibTeX")
    list_p.add_argument("--root", type=Path, default=None)
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.cmd == "rebuild":
        mapping = rebuild(args.root, online=not args.offline, harvest=not args.no_harvest)
        dest = cite_dir(args.root)
        missing = missing_tools(mapping)
        print(f"Wrote {dest / 'citations.json'} ({len(mapping)} tools)")
        if missing:
            print("No citation for: " + ", ".join(missing))
        return 0
    if args.cmd == "missing":
        mapping = load_citations(args.root)
        for name in missing_tools(mapping):
            print(name)
        return 0
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
