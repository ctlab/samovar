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
    dest="$(cd "$(dirname "${BASH_SOURCE[0]}")/toy/run/.database" && pwd)"
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
