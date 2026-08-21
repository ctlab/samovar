#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-samovar_realistic}"
mkdir -p "$output_dir/.genomes" "$output_dir/.database" "$output_dir/.log"

# ---------------------------------------------------------------------------
# 0. Download 21 raw genomes (3 per NCBI Organism group) into .genomes
# ---------------------------------------------------------------------------
# Protista has no RefSeq "complete genome" hits; use Alveolata as protist proxy.
# NOTE: do not name this GROUPS — bash reserves GROUPS as the user GID list
ORG_GROUPS=(Archaea Bacteria Viridiplantae Alveolata Fungi Metazoa Viruses)
existing_genomes="$(find "${output_dir}/.genomes" -maxdepth 1 -name '*-processed.fasta' 2>/dev/null | wc -l)"

if [[ "$existing_genomes" -ge 21 ]]; then
    echo "Found ${existing_genomes} genomes under ${output_dir}/.genomes; skipping fetch"
else
    for org_group in "${ORG_GROUPS[@]}"; do
        tmp_dir="${output_dir}/.genomes/_tmp_${org_group}"
        rm -rf "$tmp_dir"
        mkdir -p "$tmp_dir"

        echo "Fetching 3 genomes for ${org_group}..."
        python -m samovar.genome_fetcher \
            --output-dir "$tmp_dir" \
            --N 3 \
            --group "$org_group" \
            --max-genome-mb 100 \
            --email "$NCBI_EMAIL" \
            --silent
        sleep 2

        shopt -s nullglob
        for f in "${tmp_dir}"/*-processed.fasta; do
            base="$(basename "$f")"
            dest="${output_dir}/.genomes/${base}"
            if [[ -e "$dest" ]]; then
                dest="${output_dir}/.genomes/${org_group}_${base}"
            fi
            mv "$f" "$dest"
        done
        shopt -u nullglob
        rm -rf "$tmp_dir"

        got="$(find "${output_dir}/.genomes" -maxdepth 1 -name '*-processed.fasta' | wc -l)"
        echo "  total genomes so far: ${got}"
    done
fi

# ---------------------------------------------------------------------------
# 1. Public indexes into .database/ (or SAMOVAR_KAIJU_DB / SAMOVAR_KRAKEN2_DB)
# ---------------------------------------------------------------------------
KAIJU_DB_DIR="${output_dir}/.database/kaiju_db"
KRAKEN2_DB_DIR="${output_dir}/.database/kraken2_db"
KAIJU_INDEX_URL="${KAIJU_INDEX_URL:-https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_ref_2024-08-14.tgz}"
KRAKEN2_INDEX_URL="${KRAKEN2_INDEX_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260626.tar.gz}"

samovar_use_or_fetch_db "${SAMOVAR_KAIJU_DB:-}" "$KAIJU_DB_DIR" "$KAIJU_INDEX_URL" "*.fmi"
if [[ -n "${SAMOVAR_KRAKEN2_DB:-}" ]]; then
    samovar_use_or_fetch_db "$SAMOVAR_KRAKEN2_DB" "$KRAKEN2_DB_DIR" "$KRAKEN2_INDEX_URL" "hash.k2d"
else
    samovar_fetch_archive "$KRAKEN2_INDEX_URL" "$KRAKEN2_DB_DIR" "hash.k2d"
fi

KAIJU_FMI="$(find -L "$KAIJU_DB_DIR" -name '*.fmi' -print -quit)"
if [[ -z "$KAIJU_FMI" ]]; then
    echo "Error: no .fmi file found under ${KAIJU_DB_DIR}" >&2
    exit 1
fi
KAIJU_DB_PATH="$(dirname "$KAIJU_FMI")"

# ---------------------------------------------------------------------------
# 2. Generate 10 metagenomes
# ---------------------------------------------------------------------------
samovar generate \
    --genome_dir "${output_dir}/.genomes" \
    --host_genome "${SAMOVAR}/data/test_genomes/host/9606.fna" \
    --n_samples 10 \
    --total_reads 100000 \
    --output_dir "$output_dir" \
    --cores "${SAMOVAR_CORES:-32}"

# ---------------------------------------------------------------------------
# 3. Preprocess with kraken2 + kaiju (downloaded DBs)
# ---------------------------------------------------------------------------
samovar preprocess \
    --output_dir "$output_dir" \
    --kraken2-test "kraken2 ${KRAKEN2_DB_DIR}" \
    --kaiju-test "kaiju ${KAIJU_DB_PATH}" \
    --cores "${SAMOVAR_CORES:-32}"

# ---------------------------------------------------------------------------
# 4. Submit Slurm exec job (override partition/mem via env)
# ---------------------------------------------------------------------------
SLURM_PARTITION="${SLURM_PARTITION:-main}"
SLURM_CPUS="${SLURM_CPUS:-32}"
SLURM_MEM="${SLURM_MEM:-250G}"
SLURM_TIME="${SLURM_TIME:-48:00:00}"
SLURM_SCRIPT="${output_dir}/.log/slurm_exec.sh"
cat > "$SLURM_SCRIPT" << SBATCH_EOF
#!/bin/bash
#SBATCH --job-name=samovar_realistic
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --ntasks=1
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output=${output_dir}/.log/slurm_exec_%j.out
#SBATCH --error=${output_dir}/.log/slurm_exec_%j.err

set -euo pipefail
cd "${SAMOVAR}"
export PATH="${SAMOVAR}/bin:\${CONDA_PREFIX:+\${CONDA_PREFIX}/bin:}\${PATH}"
export SAMOVAR_CORES=${SLURM_CPUS}
export PYTHONPATH="${SAMOVAR}/src:\${PYTHONPATH:-}"
export NCBI_EMAIL="${NCBI_EMAIL}"

samovar exec --output_dir "${output_dir}"
SBATCH_EOF

chmod +x "$SLURM_SCRIPT"
JOBID="$(sbatch --parsable "$SLURM_SCRIPT")"
echo "Submitted Slurm job: ${JOBID}"
echo "$JOBID" > "${output_dir}/.log/slurm_jobid.txt"
