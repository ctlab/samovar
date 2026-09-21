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
#   SAMOVAR_KRAKEN2_DB_ROOT default cluster Kraken2 store (override on other hosts)
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
  export SAMOVAR_KRAKEN2_DB_ROOT="${SAMOVAR_KRAKEN2_DB_ROOT:-/mnt/tank/scratch/partition-metagenomics/databases/kraken2}"
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

samovar_phage_database_root() {
  echo "${SAMOVAR_PUBLIC_DB_ROOT:-${SAMOVAR}/examples_outdir/phage/database1}"
}

# CI public examples: use locally built phage_test indexes (same as examples/phage).
# Set SAMOVAR_CI_LIGHT_INDEXES=1 (full-integration examples-public job).
samovar_light_public_indexes() {
  [[ "${SAMOVAR_CI_LIGHT_INDEXES:-0}" == "1" ]]
}

# Resolve annotator DB names/paths for realistic / assembly / db comparison.
# Light mode → phage_test (built under DEST/.database or reused from examples/phage).
# Sets: K2_NAME K2_DIR KAIJU_NAME KAIJU_DIR  (K2_URL/KAIJU_URL empty for phage)
samovar_public_index_vars() {
  if samovar_light_public_indexes; then
    local root="${SAMOVAR_PUBLIC_DB_ROOT:-$(samovar_phage_database_root)}"
    K2_NAME="phage_test"
    K2_DIR="${root}/.database/kraken2_db"
    K2_URL=""
    KAIJU_NAME="phage_test"
    KAIJU_DIR="${root}/.database/kaiju_db"
    KAIJU_URL=""
  else
    K2_NAME="standard_8GB"
    K2_DIR="${SAMOVAR_KRAKEN2_DB_ROOT}/standard_8GB_2025oct"
    K2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz"
    KAIJU_NAME="refseq"
    KAIJU_DIR="${SAMOVAR_KAIJU_DB:-/mnt/tank/scratch/partition-metagenomics/databases/kaiju/refseq_2024aug}"
    KAIJU_URL="https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_2024-08-14.tgz"
  fi
}

# Build or reuse phage_test kraken2+kaiju indexes and import them.
samovar_ensure_phage_indexes() {
  local root="${1:-$(samovar_phage_database_root)}"
  mkdir -p "$root/.database"
  local have_k2=0 have_kj=0
  [[ -e "$root/.database/kraken2_db/hash.k2d" ]] && have_k2=1
  find -L "$root/.database/kaiju_db" -name '*.fmi' 2>/dev/null | grep -q . && have_kj=1
  if [[ "$have_k2" == 1 && "$have_kj" == 1 && "${SAMOVAR_REBUILD_DB:-0}" != "1" ]]; then
    echo "Reusing phage_test indexes in ${root}/.database"
  else
    local kaiju_acc=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000867865.1)
    local kraken_acc=(GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)
    echo "Building phage_test indexes under ${root}/.database"
    python -m samovar.genome_fetcher \
      --output-dir "$root/.database/kaiju_src" \
      --accessions "${kaiju_acc[@]}" \
      --reindex 0 \
      --email "$NCBI_EMAIL"
    python -m samovar.genome_fetcher \
      --output-dir "$root/.database/kraken2_src" \
      --accessions "${kraken_acc[@]}" \
      --reindex 0 \
      --email "$NCBI_EMAIL"
    cat > "$root/.database/kaiju.yaml" << EOF
input_dir:
  - "${root}/.database/kaiju_src/.genomes/processed"
output_dir: "${root}/database_prep_kaiju"
mutation_rate: 0.0
include_percent: 100.0
EOF
    cat > "$root/.database/kraken2.yaml" << EOF
input_dir:
  - "${root}/.database/kraken2_src/.genomes/processed"
output_dir: "${root}/database_prep_kraken2"
mutation_rate: 0.0
include_percent: 100.0
EOF
    samovar build --type kaiju \
      --config_path "$root/.database/kaiju.yaml" \
      --db_path "$root/.database/kaiju_db" \
      --index phage_test --flags ""
    samovar build --type kraken2 \
      --config_path "$root/.database/kraken2.yaml" \
      --db_path "$root/.database/kraken2_db" \
      --index phage_test --flags ""
  fi
  local lazy
  lazy="$(cat <<EOF
#!/bin/bash
set -euo pipefail
echo "Rebuild phage_test with: SAMOVAR_REBUILD_DB=1 and examples/phage or samovar_ensure_phage_indexes"
EOF
)"
  samovar tools import -n phage_test --type database --tool kaiju \
    --exec-path "$root/.database/kaiju_db" --lazy-download "$lazy"
  samovar tools import -n phage_test --type database --tool kraken2 \
    --exec-path "$root/.database/kraken2_db" --lazy-download "$lazy"
}

# Ensure public-example indexes (phage_test in light mode, else S3 lazy-download).
samovar_ensure_public_indexes() {
  samovar_public_index_vars
  if samovar_light_public_indexes; then
    local root="${SAMOVAR_PUBLIC_DB_ROOT:-$(samovar_phage_database_root)}"
    samovar_ensure_phage_indexes "$root"
    samovar_public_index_vars
  else
    samovar_ensure_database kraken2 "$K2_NAME" "$K2_DIR" "hash.k2d" "$K2_URL"
    samovar_ensure_database kaiju "$KAIJU_NAME" "$KAIJU_DIR" "*.fmi" "$KAIJU_URL"
  fi
}

# Fill DEST with phage community genomes when SAMOVAR_CI_LIGHT_INDEXES=1.
# Returns 0 if seeded / already enough files, 1 if caller should do the full NCBI fetch.
samovar_seed_public_genomes() {
  local dest="$1"
  local min_count="${2:-3}"
  mkdir -p "$dest"
  local existing
  existing="$(find "$dest" -maxdepth 1 \( -name '*.fa' -o -name '*.fa.gz' -o -name '*.fna' -o -name '*.fna.gz' -o -name '*.fasta' -o -name '*.fasta.gz' \) 2>/dev/null | wc -l)"
  if [[ "$existing" -ge "$min_count" ]]; then
    echo "Found ${existing} genomes under ${dest}; skipping fetch"
    return 0
  fi
  if ! samovar_light_public_indexes; then
    return 1
  fi
  echo "CI phage indexes: fetching phage accessions into ${dest}"
  local phage_acc=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)
  local tmp="${dest}/_tmp_phage"
  rm -rf "$tmp"
  python -m samovar.genome_fetcher \
    --output-dir "$tmp" \
    --accessions "${phage_acc[@]}" \
    --reindex 0 \
    --email "$NCBI_EMAIL"
  shopt -s nullglob
  local f base
  for f in "${tmp}"/*-processed.fasta "${tmp}"/*-processed.fasta.gz \
           "${tmp}/.genomes/processed"/* \
           "${tmp}"/*.fa.gz "${tmp}"/*.fna.gz "${tmp}"/*.fasta.gz \
           "${tmp}"/*.fa "${tmp}"/*.fna "${tmp}"/*.fasta; do
    [[ -f "$f" ]] || continue
    base="$(basename "$f")"
    if [[ ! -e "${dest}/${base}" ]]; then
      cp -a "$f" "${dest}/${base}"
    fi
  done
  shopt -u nullglob
  rm -rf "$tmp"
  existing="$(find "$dest" -maxdepth 1 \( -name '*.fa' -o -name '*.fa.gz' -o -name '*.fna' -o -name '*.fna.gz' -o -name '*.fasta' -o -name '*.fasta.gz' \) 2>/dev/null | wc -l)"
  [[ "$existing" -ge 1 ]]
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
