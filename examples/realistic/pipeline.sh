#!/usr/bin/env bash
set -euo pipefail

# Resolve SamovaR repo root from this script's location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMOVAR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$SAMOVAR"

output_dir="samovar_realistic"
mkdir -p "$output_dir/.genomes" "$output_dir/.database" "$output_dir/.log"

# Prefer conda tools (kraken2, kaiju, iss, snakemake) when available
CONDA_BIN="${CONDA_PREFIX:-/nfs/home/dsmutin/miniconda3}/bin"
export PATH="${SAMOVAR}/bin:${CONDA_BIN}:${PATH}"

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
        # Sequential fetch avoids NCBI 429; prefer genomes <=100MB for practicality
        python -m samovar.genome_fetcher \
            --output-dir "$tmp_dir" \
            --N 3 \
            --group "$org_group" \
            --max-genome-mb 100 \
            --silent
        sleep 2

        # Flatten into .genomes (avoid clobbering same taxids across groups)
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
# 1. Download small real databases into .database/
# ---------------------------------------------------------------------------
KAIJU_DB_DIR="${output_dir}/.database/kaiju_db"
KRAKEN2_DB_DIR="${output_dir}/.database/kraken2_db"
# Prefer shared RefSeq Kaiju index when present (avoids huge RVDB / re-download).
KAIJU_LOCAL_REFSEQ="/mnt/tank/scratch/partition-metagenomics/databases/kaiju/refseq_2024aug"

if [[ -f "${KAIJU_LOCAL_REFSEQ}/kaiju_db_refseq.fmi" ]]; then
    echo "Using local Kaiju RefSeq DB: ${KAIJU_LOCAL_REFSEQ}"
    mkdir -p "$KAIJU_DB_DIR"
    # Symlink so preprocess/exec paths stay under output_dir/.database
    ln -sfn "${KAIJU_LOCAL_REFSEQ}/kaiju_db_refseq.fmi" "${KAIJU_DB_DIR}/kaiju_db_refseq.fmi"
    ln -sfn "${KAIJU_LOCAL_REFSEQ}/nodes.dmp" "${KAIJU_DB_DIR}/nodes.dmp"
    ln -sfn "${KAIJU_LOCAL_REFSEQ}/names.dmp" "${KAIJU_DB_DIR}/names.dmp"
elif ! find "$KAIJU_DB_DIR" -name '*.fmi' 2>/dev/null | grep -q .; then
    ## fetch Kaiju RefSeq-ref index (toy-style download)
    mkdir -p "$KAIJU_DB_DIR"
    wget -O "$KAIJU_DB_DIR/kaiju_refseq_ref.tar.gz" \
        "https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_ref_2024-08-14.tgz"
    tar -xzf "$KAIJU_DB_DIR/kaiju_refseq_ref.tar.gz" -C "$KAIJU_DB_DIR"
    rm -f "$KAIJU_DB_DIR/kaiju_refseq_ref.tar.gz"
else
    echo "Kaiju DB already present under ${KAIJU_DB_DIR}; skipping download"
fi

if [[ ! -f "${KRAKEN2_DB_DIR}/hash.k2d" ]]; then
    mkdir -p "$KRAKEN2_DB_DIR"
    kraken_tgz="${output_dir}/.database/k2_standard_08_GB_20260626.tar.gz"
    wget -O "$kraken_tgz" \
        "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260626.tar.gz"
    tar -xzf "$kraken_tgz" -C "$KRAKEN2_DB_DIR"
    # Some archives nest contents one level deeper; hoist if needed
    if [[ ! -f "${KRAKEN2_DB_DIR}/hash.k2d" ]]; then
        nested="$(find "$KRAKEN2_DB_DIR" -name 'hash.k2d' -print -quit)"
        if [[ -n "$nested" ]]; then
            nested_dir="$(dirname "$nested")"
            shopt -s dotglob nullglob
            mv "${nested_dir}"/* "$KRAKEN2_DB_DIR"/
            shopt -u dotglob nullglob
        fi
    fi
    rm -f "$kraken_tgz"
else
    echo "Kraken2 DB already present under ${KRAKEN2_DB_DIR}; skipping download"
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
    --cores 32

# ---------------------------------------------------------------------------
# 3. Preprocess with kraken2 + kaiju (downloaded DBs)
# ---------------------------------------------------------------------------
samovar preprocess \
    --output_dir "$output_dir" \
    --kraken2-test "kraken2 ${KRAKEN2_DB_DIR}" \
    --kaiju-test "kaiju ${KAIJU_DB_PATH}" \
    --cores 32

# ---------------------------------------------------------------------------
# 4. Submit Slurm exec job (32 CPUs, 0 GPUs, partition main)
# ---------------------------------------------------------------------------
SLURM_SCRIPT="${output_dir}/.log/slurm_exec.sh"
cat > "$SLURM_SCRIPT" << SBATCH_EOF
#!/bin/bash
#SBATCH --job-name=samovar_realistic
#SBATCH --partition=main
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --output=${output_dir}/.log/slurm_exec_%j.out
#SBATCH --error=${output_dir}/.log/slurm_exec_%j.err

set -euo pipefail
cd "${SAMOVAR}"
CONDA_BIN="${CONDA_PREFIX:-/nfs/home/dsmutin/miniconda3}/bin"
export PATH="${SAMOVAR}/bin:\${CONDA_BIN}:\${PATH}"
export SAMOVAR_CORES=32
export PYTHONPATH="${SAMOVAR}/src:\${PYTHONPATH:-}"

samovar exec --output_dir "${output_dir}"
SBATCH_EOF

chmod +x "$SLURM_SCRIPT"
JOBID="$(sbatch --parsable "$SLURM_SCRIPT")"
echo "Submitted Slurm job: ${JOBID}"
