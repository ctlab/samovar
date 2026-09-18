"""Bundled BibTeX for built-in tools, plus harvest on ``./install.sh``.

``cite/citations.json`` maps each ``TOOL_GROUP_BY_NAME`` key to ``*.bib`` files
under ``cite/``. ``python -m samovar.citations rebuild`` refreshes those files
from doi.org, R ``citation()``, and CLI ``--citation`` when available.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from samovar.main_config import TOOL_GROUP_BY_NAME
from samovar.paths import _code_repo_root

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
        data = json.loads(path.read_text(encoding="utf-8"))
        return {str(k): list(v) for k, v in data.items()}
    return mapping_for_tools()


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
