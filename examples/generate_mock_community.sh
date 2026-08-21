#!/bin/bash
# Synthesize a mixed mock community (10 bacteria / archaea / viruses / eukaryotes).
# Large downloads go under SAMOVAR_MOCK_ROOT (default: repo-relative genomes/).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MOCK_ROOT="${SAMOVAR_MOCK_ROOT:-$ROOT}"
MOCK_DIR="${MOCK_ROOT}/genomes/mock_community"
OUT_DIR="${SAMOVAR_MOCK_OUT:-$ROOT/tests_outs/benchmarking/initial}"
N_READS="${SAMOVAR_MOCK_READS:-5000000}"
CPUS="${SAMOVAR_MOCK_CPUS:-16}"

mkdir -p "$MOCK_DIR" "$OUT_DIR"

BACTERIA=("NC_000913.3" "NC_007795.1" "NC_002516.2" "NC_000962.3" "NC_003197.2" "NC_004557.1" "NC_000919.1" "NC_003888.3" "NC_012489.1" "NC_000964.3")
ARCHAEA=("NC_000909.1" "NC_002607.1" "NC_002754.1" "NC_013790.1" "NC_000854.2" "NC_003551.1" "NC_003552.1" "NC_003364.1" "NC_005791.1" "NC_003413.1")
VIRUSES=("NC_001416.1" "NC_001422.1" "NC_001806.2" "NC_002549.1" "NC_001604.1" "NC_001451.1" "NC_001489.1" "NC_004162.2" "NC_001405.1" "NC_045512.2")
EUKARYOTA=("NC_001133.9" "NC_001134.8" "NC_001135.5" "NC_001136.10" "NC_004325.2" "NC_004314.2" "NC_004318.2" "NC_000932.1" "NC_013967.1" "NC_013968.1")
ALL_GENOMES=("${BACTERIA[@]}" "${ARCHAEA[@]}" "${VIRUSES[@]}" "${EUKARYOTA[@]}")

echo "==================================================="
echo "  SamovaR: Complex Mock Community Generator"
echo "  genomes: $MOCK_DIR"
echo "  reads:   $OUT_DIR  (n_reads=$N_READS)"
echo "==================================================="

for ACC in "${ALL_GENOMES[@]}"; do
    if [ ! -s "${MOCK_DIR}/${ACC}.fa" ]; then
        echo "  -> Downloading ${ACC}..."
        wget -qO "${MOCK_DIR}/${ACC}.fa" \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text"
        sleep 1
    fi
done

echo "[INFO] Generating ${N_READS} paired-end reads..."
iss generate \
    --genomes "${MOCK_DIR}"/*.fa \
    --model hiseq \
    --n_reads "$N_READS" \
    --cpus "$CPUS" \
    --output "${OUT_DIR}/complex_mock"
