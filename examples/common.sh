#!/usr/bin/env bash
# Shared helpers for example pipelines. Source after setting SCRIPT_DIR.
#
# Paths are resolved from this file so examples work from any cwd:
#   bash /absolute/path/to/examples/toy/pipeline.sh
#
# Optional environment:
#   SAMOVAR                 repo root (inferred from this file)
#   CONDA_PREFIX            prepended to PATH
#   NCBI_EMAIL / ENTREZ_EMAIL / SAMOVAR_EMAIL
#   SAMOVAR_KRAKEN2_DB_ROOT  unused by examples; indexes resolve from the
#                            install catalog, then official lazy-download
#   SAMOVAR_SLURM=1         wrap `samovar exec` in sbatch (not used inside example
#                           scripts themselves — set this when launching)
#   SAMOVAR_SLURM_CPUS / SAMOVAR_SLURM_MEM / SAMOVAR_SLURM_TIME / SLURM_PARTITION
#   SAMOVAR_SLURM_WAIT=1    sbatch --wait (default when SAMOVAR_SLURM=1)
#   SAMOVAR_CORES           passed to generate/prepare

if [[ -z "${SAMOVAR:-}" || ! -d "${SAMOVAR}/src/samovar" ]]; then
  SAMOVAR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

samovar_setup_env() {
  local conda_bin=""
  if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/bin" ]]; then
    conda_bin="${CONDA_PREFIX}/bin"
  fi
  export PATH="${SAMOVAR}/bin${conda_bin:+:${conda_bin}}:${PATH}"
  export PYTHONPATH="${SAMOVAR}/src${PYTHONPATH:+:${PYTHONPATH}}"
  export NCBI_EMAIL="${NCBI_EMAIL:-${ENTREZ_EMAIL:-${SAMOVAR_EMAIL:-anonymous@example.com}}}"
}

# Full run trees live under examples_outdir/<name>/, not under examples/<name>/.
# examples/ keeps scripts + harvested figures/multiqc_report.html only.
# Override with SAMOVAR_OUTDIR (absolute path to one run).
samovar_example_outdir() {
  local name="${1:-}"
  if [[ -z "$name" ]]; then
    name="$(basename "${SCRIPT_DIR:-.}")"
  fi
  if [[ -n "${SAMOVAR_OUTDIR:-}" ]]; then
    echo "$SAMOVAR_OUTDIR"
    return 0
  fi
  echo "${SAMOVAR}/examples_outdir/${name}"
}

samovar_toy_database_dir() {
  echo "${SAMOVAR}/examples_outdir/toy/.database"
}

# CI public examples: catalog name phage_test, imported by examples/phage/pipeline.sh.
# Set SAMOVAR_CI_LIGHT_INDEXES=1 (full-integration examples-public job).
samovar_light_public_indexes() {
  [[ "${SAMOVAR_CI_LIGHT_INDEXES:-0}" == "1" ]]
}

# Catalog names + official URLs. Paths are not set here: a preinstalled index
# comes from the SamovaR database catalog; otherwise lazy-download uses the URL.
# Light mode → phage_test (imported by examples/phage/pipeline.sh).
# Sets: K2_NAME KAIJU_NAME K2_URL KAIJU_URL
samovar_public_index_vars() {
  if samovar_light_public_indexes; then
    K2_NAME="phage_test"
    K2_URL=""
    KAIJU_NAME="phage_test"
    KAIJU_URL=""
  else
    K2_NAME="standard_8GB"
    K2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz"
    KAIJU_NAME="refseq"
    KAIJU_URL="https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_2024-08-14.tgz"
  fi
}

# Ensure catalog names for the full public indexes (S3 lazy-download if missing).
# Light mode only selects the name phage_test. examples/phage/pipeline.sh imports it.
samovar_ensure_public_indexes() {
  samovar_public_index_vars
  if samovar_light_public_indexes; then
    return 0
  fi
  K2_NAME="$(samovar_ensure_named_database kraken2 "$K2_NAME" "hash.k2d" "$K2_URL")"
  KAIJU_NAME="$(samovar_ensure_named_database kaiju "$KAIJU_NAME" "*.fmi" "$KAIJU_URL")"
}

# Resolve a catalog database without embedding a machine path.
# Uses an installed catalog copy when its marker is on disk (another version of
# the same tool counts, so a local run can use a preinstalled index).
# Otherwise lazy-downloads URL into examples_outdir/databases/<tool>/<name>.
# Stdout is the catalog name to pass to prepare. Logs go to stderr.
# Usage: samovar_ensure_named_database TOOL NAME MARKER [URL] [DOWNLOAD]
# DOWNLOAD=0 returns 1 instead of fetching when nothing is installed.
samovar_ensure_named_database() {
  local tool="$1"
  local name="$2"
  local marker="$3"
  local url="${4:-}"
  local download="${5:-1}"
  local line
  line="$(python - "$tool" "$name" "$marker" "$url" "$download" "${SAMOVAR}" <<'PY'
import sys
from pathlib import Path

from samovar.db_spec import iter_database_records, lookup_database_record, official_url_for
from samovar.paths import load_config

tool, name, marker, url, download, samovar = sys.argv[1:]
cfg = load_config()

def has_marker(path: str) -> bool:
    root = Path(path).expanduser()
    if not root.is_dir():
        return False
    # Marker files sit at the index root. Do not walk the whole tree.
    return any(root.glob(marker)) or any(root.glob(f"*/{marker}"))

def emit(status: str, resolved: str, path: str, href: str) -> None:
    print("\t".join((status, resolved, path, href)))

grouped = iter_database_records(cfg).get(tool) or {}
rec = lookup_database_record(cfg, tool, name) or {}
href = url or str(rec.get("url") or "") or official_url_for(
    tool, name, str(rec.get("_version") or ""), url
)
path = str(rec.get("path") or "")
if path and has_marker(path):
    emit("installed", str(rec.get("name") or name), path, href)
    raise SystemExit(0)
if download != "1":
    emit("missing", name, "", href)
    raise SystemExit(0)
for other in grouped.values():
    opath = str(other.get("path") or "")
    oname = str(other.get("name") or "")
    if oname == name or not opath or not has_marker(opath):
        continue
    # Local stand-in: another installed version of this tool. Keep its own URL.
    emit("installed", oname, opath, str(other.get("url") or ""))
    raise SystemExit(0)
dest = str(Path(samovar) / "examples_outdir" / "databases" / tool / name)
emit("download", name, dest, href)
PY
)"
  local status resolved dest href
  IFS=$'\t' read -r status resolved dest href <<<"$line"
  if [[ "$status" == "missing" ]]; then
    echo "REPORT: ${tool}/${name} is not installed and download is off" >&2
    return 1
  fi
  if [[ "$status" == "installed" ]]; then
    echo "Using installed ${tool} index ${resolved}" >&2
    echo "$resolved"
    return 0
  fi
  echo "lazy-download ${tool}/${resolved} (not installed)" >&2
  samovar_ensure_database "$tool" "$resolved" "$dest" "$marker" "$href" >&2
  echo "$resolved"
}

# Print REPORT if a registered database has a path but no lazy-download recipe.
samovar_report_missing_lazy() {
  python - <<'PY'
from samovar.db_spec import iter_database_records
from samovar.paths import load_config
from pathlib import Path
cfg = load_config()
for tool, grouped in iter_database_records(cfg).items():
    for key, rec in grouped.items():
        path = str(rec.get("path") or "").strip()
        lazy = str(rec.get("lazy-download") or "").strip()
        url = str(rec.get("url") or "").strip()
        if path and not lazy:
            exists = Path(path).expanduser().exists()
            print(
                f"REPORT: databases.{tool}.{key} has no lazy-download "
                f"(path exists={exists}: {path})"
            )
PY
}

# Import a database record. URL fills lazy-download via db_spec defaults.
# Usage: samovar_import_database TOOL NAME PATH [URL]
samovar_import_database() {
  local tool="$1"
  local name="$2"
  local path="$3"
  local url="${4:-}"
  local extra=()
  if [[ -n "$url" ]]; then
    extra+=(--url "$url")
  fi
  samovar tools import -n "$name" --type database --tool "$tool" \
    --exec-path "$path" "${extra[@]}"
}

# If PATH already has MARKER, import and return. Else run stored/default
# lazy-download with PREFIX=PATH. If the index has no recipe, print REPORT.
# Usage: samovar_ensure_database TOOL NAME PATH MARKER [URL]
samovar_ensure_database() {
  local tool="$1"
  local name="$2"
  local dest="$3"
  local marker="${4:-}"
  local url="${5:-}"
  mkdir -p "$dest"
  if [[ -n "$marker" ]] && find -L "$dest" -name "$marker" 2>/dev/null | grep -q .; then
    echo "Database already on disk: ${dest} (${marker})"
    samovar_import_database "$tool" "$name" "$dest" "$url"
    return 0
  fi
  local recipe
  recipe="$(python - <<PY
from samovar.db_spec import lazy_download_for, lookup_database_record
from samovar.paths import load_config
tool, name, url = "$tool", "$name", "$url"
rec = lookup_database_record(load_config(), tool, name) or {}
text = str(rec.get("lazy-download") or "").strip()
if not text:
    text = lazy_download_for(tool, name, str(rec.get("_version") or ""), url or rec.get("url") or "")
print(text)
PY
)"
  if [[ -z "$recipe" ]]; then
    echo "REPORT: databases.${tool}.${name} is not on disk at ${dest} and has no lazy-download" >&2
    return 1
  fi
  echo "lazy-download ${tool}/${name} -> ${dest}"
  PREFIX="$dest" bash -c "$recipe"
  samovar_import_database "$tool" "$name" "$dest" "$url"
}

# Run samovar exec. Examples call this instead of embedding sbatch.
samovar_run_exec() {
  local out="$1"
  shift || true
  if [[ "${SAMOVAR_SLURM:-0}" == "1" ]] && command -v sbatch >/dev/null 2>&1; then
    local cpus="${SAMOVAR_SLURM_CPUS:-${SAMOVAR_CORES:-16}}"
    local mem="${SAMOVAR_SLURM_MEM:-32G}"
    local time="${SAMOVAR_SLURM_TIME:-12:00:00}"
    local part="${SLURM_PARTITION:-main}"
    mkdir -p "${out}/.log"
    local script="${out}/.log/slurm_exec.sh"
    cat > "$script" << EOF
#!/bin/bash
#SBATCH --job-name=samovar_example
#SBATCH --partition=${part}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --ntasks=1
#SBATCH --mem=${mem}
#SBATCH --time=${time}
#SBATCH --output=${out}/.log/slurm_exec_%j.out
#SBATCH --error=${out}/.log/slurm_exec_%j.err
set -euo pipefail
export PATH="${SAMOVAR}/bin:\${CONDA_PREFIX:+\${CONDA_PREFIX}/bin:}\${PATH}"
export PYTHONPATH="${SAMOVAR}/src:\${PYTHONPATH:-}"
export NCBI_EMAIL="${NCBI_EMAIL:-}"
export SAMOVAR_CORES=${cpus}
samovar exec --output_dir "${out}"
EOF
    chmod +x "$script"
    local wait_flag=()
    if [[ "${SAMOVAR_SLURM_WAIT:-1}" == "1" ]]; then
      wait_flag=(--wait)
    fi
    echo "Submitting samovar exec via sbatch (${cpus} CPUs)"
    sbatch "${wait_flag[@]}" --parsable "$script" | tee "${out}/.log/slurm_jobid.txt"
  else
    samovar exec --output_dir "$out" "$@"
  fi
}

# Build toy Kraken2/Kaiju from bundled test genomes if not already imported.
samovar_ensure_toy_annotators() {
  local dest="${1:-}"
  if [[ -z "$dest" ]]; then
    dest="$(samovar_toy_database_dir)"
  fi
  mkdir -p "$dest"
  local yaml="$dest/config.yaml"
  cat > "$yaml" << EOF
input_dir:
  - "${SAMOVAR}/data/test_genomes/meta"
  - "${SAMOVAR}/data/test_genomes/host"
output_dir: "${dest}/database_prep"
mutation_rate: 0.02
include_percent: 70.0
EOF
  if [[ ! -e "$dest/kraken2_db/hash.k2d" && ! -e "$dest/kraken2_db/taxo.k2d" ]]; then
    samovar build_database --type kraken2 --config_path "$yaml" \
      --db_path "$dest/kraken2_db" --example-omit --index toy --flags ""
  fi
  if ! find -L "$dest/kaiju_db" -name '*.fmi' 2>/dev/null | grep -q .; then
    samovar build_database --type kaiju --config_path "$yaml" \
      --db_path "$dest/kaiju_db" --example-omit --index toy --flags ""
  fi
  local toy_lazy
  toy_lazy="$(cat <<EOF
#!/bin/bash
set -euo pipefail
DEST="\${PREFIX:-.}"
SAMOVAR="\${SAMOVAR:-}"
if [[ -z "\$SAMOVAR" || ! -d "\$SAMOVAR/src/samovar" ]]; then
  echo "Set SAMOVAR to the SamovaR repo root to rebuild the toy index." >&2
  exit 1
fi
echo "Rebuild toy indexes with: bash \$SAMOVAR/examples/toy/pipeline.sh"
EOF
)"
  samovar tools import -n toy --type database --tool kraken2 \
    --exec-path "$dest/kraken2_db" --lazy-download "$toy_lazy"
  samovar tools import -n toy --type database --tool kaiju \
    --exec-path "$dest/kaiju_db" --lazy-download "$toy_lazy"
}

# Copy a few key plots + multiqc_report.html into examples/<name>/{figures,multiqc}.
# Full run stays under examples_outdir/; do not copy multiqc_data.
samovar_harvest_example() {
  local run="$1"
  local example="${2:-${SCRIPT_DIR}}"
  python - "$run" "$example" <<'PY'
from __future__ import annotations

import shutil
import sys
from pathlib import Path

STAGES = (
    "initial_annotations_plots",
    "regenerated_annotations_plots",
    "reprofiled_annotations_plots",
)

def copy_one(src: Path, dest_dir: Path, stem: str) -> Path | None:
    if not src.is_file():
        return None
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{stem}{src.suffix}"
    shutil.copy2(src, dest)
    return dest

def first_existing(plots: Path, names: list[str]) -> Path | None:
    for name in names:
        p = plots / name
        if p.is_file():
            return p
    return None

run = Path(sys.argv[1])
example = Path(sys.argv[2])
figures = example / "figures"
multiqc_out = example / "multiqc"
figures.mkdir(parents=True, exist_ok=True)
copied: list[Path] = []

for folder in STAGES:
    plots = run / folder
    if not plots.is_dir():
        continue
    short = folder.replace("_annotations_plots", "")
    f1 = first_existing(plots, ["F1.png", "F1_SAMOVAR.png", "F1_kraken2.png", "F1_kaiju.png"])
    if f1 is None:
        pngs = sorted(plots.glob("F1_*.png"))
        f1 = pngs[0] if pngs else None
    scores = first_existing(plots, ["scores.png", "opal_scores.png"])
    roc = plots / "roc_comparison.png" if folder.startswith("reprofiled") else None
    for src, label in (
        (scores, "scores"),
        (f1, "F1"),
        (roc if roc is not None and roc.is_file() else None, "roc"),
    ):
        if src is None:
            continue
        dest = copy_one(src, figures, f"{short}_{label}")
        if dest:
            copied.append(dest)

table_pngs = sorted((run / "regenerated_annotations_plots").glob("TableScore_*.png"))
if not table_pngs:
    table_pngs = sorted(run.glob("**/TableScore_*.png"))
if table_pngs:
    dest = copy_one(table_pngs[0], figures, "table_score")
    if dest:
        copied.append(dest)

html_candidates = [
    run / "multiqc_samovar" / "multiqc_report.html",
    run / "multiqc" / "SAMOVAR_multiqc_report.html",
    run / "multiqc" / "multiqc_report.html",
    run / "SAMOVAR_multiqc_report.html",
    run / "multiqc_report.html",
]
mqc_html = next((p for p in html_candidates if p.is_file()), None)
if mqc_html is None:
    reports = list(run.glob("**/SAMOVAR_multiqc_report.html")) + list(run.glob("**/multiqc_report.html"))
    mqc_html = reports[0] if reports else None
if mqc_html:
    multiqc_out.mkdir(parents=True, exist_ok=True)
    # Drop leftover multiqc_data from older harvests.
    data_dir = multiqc_out / "multiqc_data"
    if data_dir.is_dir():
        shutil.rmtree(data_dir)
    shutil.copy2(mqc_html, multiqc_out / "multiqc_report.html")
    print("multiqc", multiqc_out / "multiqc_report.html")

print("figures", len(copied))
for p in copied:
    print(" ", p.relative_to(example))
PY
}
