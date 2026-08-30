#!/usr/bin/env bash
# Named phage_test indexes + generate --reindex 0/1 and samovar reindex.
#
# Phases (SAMOVAR_PHASE, default all):
#   database1  download generate genomes, build kaiju/kraken2 as phage_test,
#               samovar generate --reindex 1, prepare, exec
#   database2   samovar generate --reindex 0 with phage_test DBs
#   reindex     samovar reindex on database2
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

PHASE="${SAMOVAR_PHASE:-all}"
out1="${SAMOVAR_OUTDIR_1:-${SCRIPT_DIR}/run/database1}"
out2="${SAMOVAR_OUTDIR_2:-${SCRIPT_DIR}/run/database2}"
n_samples="${N_SAMPLES:-2}"
total_reads="${TOTAL_READS:-2000}"
host_fraction="${HOST_FRACTION:-0.15}"
host="${SAMOVAR}/data/test_genomes/host/9606.fna"

GEN1=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)
KAIJU=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000867865.1)
KRAKEN=(GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)
GEN2=(GCF_000819615.1 GCA_035764635.1 GCF_000836945.1 GCF_000844825.1 GCF_000867865.1)

python_download() {
  local dest="$1"
  shift
  python -m samovar.genome_fetcher \
    --output-dir "$dest" \
    --accessions "$@" \
    --reindex 0 \
    --email "$NCBI_EMAIL"
}

drop_acc() {
  python - <<'PY' "$1"
import sys
from pathlib import Path
from samovar.genome_index import drop_indexed, resolve_indexed_path
acc = sys.argv[1]
path = resolve_indexed_path(acc)
drop_indexed(acc)
if path:
    try:
        Path(path).unlink()
        print(f"dropped {acc} file {path}")
    except OSError as exc:
        print(f"dropped {acc} index; file: {exc}")
else:
    print(f"dropped {acc} (not on disk)")
PY
}

write_db_yaml() {
  local yaml="$1"
  local src="$2"
  local prep="$3"
  cat > "$yaml" << EOF
input_dir:
  - "${src}"
output_dir: "${prep}"
mutation_rate: 0.0
include_percent: 100.0
EOF
}

import_phage_dbs() {
  local root="$1"
  local kaiju_lazy kraken_lazy
  kaiju_lazy="$(cat <<EOF
#!/bin/bash
set -euo pipefail
echo "Rebuild kaiju phage_test with examples/phage/pipeline.sh (SAMOVAR_PHASE=database1)"
EOF
)"
  kraken_lazy="$(cat <<EOF
#!/bin/bash
set -euo pipefail
echo "Rebuild kraken2 phage_test with examples/phage/pipeline.sh (SAMOVAR_PHASE=database1)"
EOF
)"
  samovar tools import -n phage_test --type database --tool kaiju \
    --exec-path "$root/.database/kaiju_db" --lazy-download "$kaiju_lazy"
  samovar tools import -n phage_test --type database --tool kraken2 \
    --exec-path "$root/.database/kraken2_db" --lazy-download "$kraken_lazy"
}

phase_database1() {
  mkdir -p "$out1/.database" "$out1/.genomes"
  local have_k2 have_kj
  have_k2=0
  have_kj=0
  [[ -e "$out1/.database/kraken2_db/hash.k2d" ]] && have_k2=1
  find -L "$out1/.database/kaiju_db" -name '*.fmi' 2>/dev/null | grep -q . && have_kj=1

  if [[ "$have_k2" == 1 && "$have_kj" == 1 && "${SAMOVAR_REBUILD_DB:-0}" != "1" ]]; then
    echo "==> Reusing phage_test indexes in ${out1}/.database"
    import_phage_dbs "$out1"
  else
    echo "==> Download generate accessions into ${out1}"
  python_download "$out1" "${GEN1[@]}"

  echo "==> Download kaiju genomes (include GCF_000867865.1)"
  python_download "$out1/.database/kaiju_src" "${KAIJU[@]}"
  echo "==> Download kraken2 genomes"
  python_download "$out1/.database/kraken2_src" "${KRAKEN[@]}"

  write_db_yaml "$out1/.database/kaiju.yaml" \
    "$out1/.database/kaiju_src/.genomes/processed" \
    "$out1/database_prep_kaiju"
  write_db_yaml "$out1/.database/kraken2.yaml" \
    "$out1/.database/kraken2_src/.genomes/processed" \
    "$out1/database_prep_kraken2"

  samovar build --type kaiju \
    --config_path "$out1/.database/kaiju.yaml" \
    --db_path "$out1/.database/kaiju_db" \
    --index phage_test --flags ""
  samovar build --type kraken2 \
    --config_path "$out1/.database/kraken2.yaml" \
    --db_path "$out1/.database/kraken2_db" \
    --index phage_test --flags ""
  import_phage_dbs "$out1"
  fi

  rm -rf "$out1/initial" "$out1/.generate" "$out1/.log/checkpoints"
  echo "==> generate --reindex 1"
  samovar generate \
    --accessions "${GEN1[@]}" \
    --reindex 1 \
    --raw-genomes 0 \
    --host_genome "$host" \
    --output_dir "$out1" \
    --n_samples "$n_samples" \
    --total_reads "$total_reads" \
    --host_fraction "$host_fraction" \
    --cores 1

  samovar prepare \
    --output_dir "$out1" \
    --kraken2-test "kraken2 phage_test" \
    --kaiju-test "kaiju phage_test"

  samovar_run_exec "$out1"

  echo "==> Check index + outdir after database1"
  python - <<PY
from pathlib import Path
from samovar.genome_index import processed_storage_dir, resolve_indexed_path, run_processed_dir
accs = "GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1".split()
store = processed_storage_dir()
run = run_processed_dir("$out1")
missing = []
for acc in accs:
    idx = resolve_indexed_path(acc)
    staged = run / f"{acc}.fa.gz"
    stored = store / f"{acc}.fa.gz"
    print(f"{acc} indexed={idx} staged={staged.is_file()} stored={stored.is_file()}")
    if idx is None or not staged.is_file() or not stored.is_file():
        missing.append(acc)
if missing:
    raise SystemExit("database1 missing index/files: " + ",".join(missing))
print("genome_dir for ISS is", run)
PY

  drop_acc GCF_000867865.1
}

phase_database2() {
  rm -rf "$out2"
  mkdir -p "$out2"

  echo "==> generate --reindex 0 with phage_test DBs"
  samovar generate \
    --accessions "${GEN2[@]}" \
    --reindex 0 \
    --raw-genomes 0 \
    --host_genome "$host" \
    --output_dir "$out2" \
    --n_samples "$n_samples" \
    --total_reads "$total_reads" \
    --host_fraction "$host_fraction" \
    --cores 1

  python - <<PY
from pathlib import Path
from samovar.genome_index import genome_data_map, run_processed_dir
run = run_processed_dir("$out2")
names = sorted(p.name for p in run.glob("*.fa.gz"))
print("database2 processed:", names)
data = genome_data_map()
new = [n for n in names if n.startswith("GCF_000867865") or n.startswith("GCA_035764635")]
print("new-or-unstored names in run:", new)
if "GCF_000867865.1" in data:
    raise SystemExit("GCF_000867865.1 should not be in the main index after --reindex 0")
PY

  samovar prepare \
    --output_dir "$out2" \
    --kraken2-test "kraken2 phage_test" \
    --kaiju-test "kaiju phage_test"

  samovar_run_exec "$out2"
}

phase_reindex() {
  echo "==> samovar reindex ${out2}"
  before="$(python - <<'PY'
from samovar.genome_index import genome_data_map
print(len(genome_data_map()))
PY
)"
  samovar reindex "$out2"
  python - <<PY
from pathlib import Path
from samovar.genome_index import genome_data_map, processed_storage_dir, resolve_indexed_path, run_processed_dir
acc = "GCF_000867865.1"
idx = resolve_indexed_path(acc)
store = processed_storage_dir() / f"{acc}.fa.gz"
run = run_processed_dir("$out2") / f"{acc}.fa.gz"
print("after reindex", acc, "indexed=", idx, "stored=", store.is_file(), "run_left=", run.is_file())
if idx is None or not store.is_file():
    raise SystemExit("reindex did not move/index GCF_000867865.1")
print("index size", len(genome_data_map()))
PY
  echo "index size before reindex: ${before}"
}

case "$PHASE" in
  database1) phase_database1 ;;
  database2) phase_database2 ;;
  reindex) phase_reindex ;;
  all)
    phase_database1
    phase_database2
    phase_reindex
    ;;
  *)
    echo "Unknown SAMOVAR_PHASE=$PHASE (database1|database2|reindex|all)"
    exit 1
    ;;
esac

if [[ -d "$out1" ]]; then
  samovar multiqc --output_dir "$out1" -- --export --interactive || true
fi
