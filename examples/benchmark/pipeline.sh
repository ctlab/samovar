#!/usr/bin/env bash
# Benchmark package: four phages in the reads, five in every de novo index.
# Kraken2, Kaiju, Centrifuge, and Kraken2+Bracken. MetaX has no samovar build type.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

out="$(samovar_example_outdir)"
db="${SAMOVAR}/examples_outdir/benchmark_db"
total_reads="${TOTAL_READS:-2000}"
threads="${SAMOVAR_CORES:-1}"
community=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)
library=("${community[@]}" GCF_000867865.1)

mkdir -p "$db/src"
ready=1
for tool in kraken2 kaiju centrifuge bracken; do
  if ! samovar import --status | grep -Eq "^${tool}"$'\t'"phage_bench(:[^[:space:]]*)?"$'\t'"ok"$'\t'; then
    ready=0
  fi
done
if [[ "$ready" -eq 0 ]]; then
  # Already indexed genomes are reused from the shared store and are not copied
  # into --output-dir. Link them so build sees a FASTA directory.
  while IFS= read -r genome; do
    [[ -f "$genome" ]] || continue
    ln -sfn "$genome" "$db/src/$(basename "$genome")"
  done < <(samovar genome download --output-dir "$db/src" "${library[@]}")
  if ! compgen -G "$db/src/*" > /dev/null; then
    echo "genome download produced no FASTA under $db/src" >&2
    exit 1
  fi
  for type in kraken2 kaiju centrifuge; do
    cat > "$db/${type}.yaml" << EOF
input_dir:
  - ${db}/src
output_dir: ${db}/${type}_prep
mutation_rate: 0.0
include_percent: 100.0
EOF
    samovar build --type "$type" \
      --config_path "$db/${type}.yaml" \
      --db_path "$db/${type}_db" \
      --threads "$threads" \
      --index phage_bench --flags "" --no-example-omit
    samovar import -n phage_bench --type database --tool "$type" --exec-path "$db/${type}_db"
  done
  cat > "$db/bracken.yaml" << EOF
input_dir:
  - ${db}/src
output_dir: ${db}/bracken_prep
mutation_rate: 0.0
include_percent: 100.0
kmer_len: 35
read_lengths:
  - 150
EOF
  samovar build --type bracken \
    --config_path "$db/bracken.yaml" \
    --db_path "$db/bracken_db" \
    --threads "$threads" \
    --index phage_bench --flags "" --no-example-omit
  samovar import -n phage_bench --type database --tool bracken --exec-path "$db/bracken_db"
fi

samovar generate \
  --accessions "${community[@]}" \
  --reindex 1 \
  --host_genome "" \
  --host_fraction 0 \
  --output_dir "$out" \
  --n_samples 1 \
  --total_reads "$total_reads" \
  --repeat 3 \
  --cores "${SAMOVAR_GENERATE_CORES:-1}"
samovar prepare \
  --output_dir "$out" \
  --mode benchmark \
  --cores "$threads" \
  --kraken2 "kraken2 phage_bench" \
  --kaiju "kaiju phage_bench" \
  --centrifuge "centrifuge phage_bench" \
  --bracken "bracken phage_bench"
samovar exec --output_dir "$out"
samovar multiqc --output_dir "$out"
# Package plots live under summary/ and out_XXXX/, not the single-run folders
# samovar_harvest_example expects.
python - "$out" "$SCRIPT_DIR" <<'PY'
import shutil, sys
from pathlib import Path
run, example = Path(sys.argv[1]), Path(sys.argv[2])
figures, report = example / "figures", example / "multiqc"
figures.mkdir(parents=True, exist_ok=True)
report.mkdir(parents=True, exist_ok=True)
html = run / "summary" / "multiqc" / "SAMOVAR_multiqc_report.html"
if html.is_file():
    shutil.copy2(html, report / "multiqc_report.html")
for name in ("annotators.png", "profilers.png"):
    src = run / "summary" / "figures" / name
    if src.is_file():
        shutil.copy2(src, figures / name)
print("harvest", report / "multiqc_report.html")
PY
