# Realistic pipeline

End-to-end SamovaR benchmark on **real NCBI genomes**, public **Kaiju** / **Kraken2** indexes, and **Slurm** (`generate` → `preprocess` → `exec`).

Entry script: [`examples/realistic/pipeline.sh`](../examples/realistic/pipeline.sh)

Example figures from a completed run live in [`examples/realistic/results/`](../examples/realistic/results/).

---

## Overview

| Step | What happens |
|------|----------------|
| 0 | Fetch **21** RefSeq genomes (3 each: Archaea, Bacteria, Viridiplantae, Alveolata, Fungi, Metazoa, Viruses), prefer ≤100 MB assemblies |
| 1 | Attach **Kaiju RefSeq** + **Kraken2 Standard 8 GB** indexes |
| 2 | `samovar generate` — **10** metagenomes × **100 000** reads, `--cores 32` |
| 3 | `samovar preprocess` — wire Kraken2 + Kaiju, `--cores 32`, `--max_genomes 50` |
| 4 | `sbatch` — 32 CPUs, 250 G RAM, 0 GPUs → `samovar exec` |

---

## Prerequisites

- SamovaR installed (`./install.sh` or `pip install -e .`) with `build/config.json`
- On `PATH`: `bin/samovar`, `kraken2`, `kaiju`, `iss`, `snakemake`, `R` (conda: activate the env so `CONDA_PREFIX` is set)
- Slurm (`sbatch`) and network access to NCBI / public index URLs
- Optional local indexes: `SAMOVAR_KAIJU_DB`, `SAMOVAR_KRAKEN2_DB`
- NCBI contact: `NCBI_EMAIL` (or `ENTREZ_EMAIL` / `SAMOVAR_EMAIL`)
- Roughly **≥120 GB** free for Kraken2 + working files

Check tools:

```bash
which samovar kraken2 kaiju iss snakemake R sbatch
cat build/config.json
```

---

## Quick start

From the repository root:

```bash
cd samovar
bash examples/realistic/pipeline.sh
```

The script prints a Slurm job id. Follow progress:

```bash
JOB=$(cat samovar_realistic/.log/slurm_jobid.txt)
squeue -j "$JOB"
tail -f samovar_realistic/.log/slurm_exec_${JOB}.err
```

Re-runs skip genome/DB download when `.genomes` already has ≥21 `*-processed.fasta` files and indexes are present.

---

## Stage-by-stage manual

Set a common environment once:

```bash
export PATH="$PWD/bin${CONDA_PREFIX:+:$CONDA_PREFIX/bin}:$PATH"
export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"
export NCBI_EMAIL="${NCBI_EMAIL:-you@example.org}"
# Optional: skip the public Kaiju download
# export SAMOVAR_KAIJU_DB=/path/to/kaiju_index_dir
OUT=samovar_realistic
mkdir -p "$OUT/.genomes" "$OUT/.database" "$OUT/.log"
```

### 1. Fetch genomes

Do **not** name the bash array `GROUPS` — bash reserves that name for user GIDs.

```bash
# Protista has no RefSeq "complete genome" hits → use Alveolata
ORG_GROUPS=(Archaea Bacteria Viridiplantae Alveolata Fungi Metazoa Viruses)

for g in "${ORG_GROUPS[@]}"; do
  tmp="$OUT/.genomes/_tmp_$g"
  rm -rf "$tmp"; mkdir -p "$tmp"
  python -m samovar.genome_fetcher \
    --output-dir "$tmp" --N 3 --group "$g" --max-genome-mb 100 \
    --email "$NCBI_EMAIL" --silent
  mv "$tmp"/*-processed.fasta "$OUT/.genomes/" 2>/dev/null || true
  rm -rf "$tmp"
  sleep 2
done

ls "$OUT/.genomes"/*-processed.fasta | wc -l   # expect 21
```

### 2. Databases

Indexes are downloaded into `$OUT/.database/` unless you point at a local copy:

```bash
# export SAMOVAR_KAIJU_DB=/path/to/existing/kaiju_dir
# export SAMOVAR_KRAKEN2_DB=/path/to/existing/kraken2_dir
KAIJU_DB="$OUT/.database/kaiju_db"
KRAKEN2_DB="$OUT/.database/kraken2_db"
mkdir -p "$KAIJU_DB" "$KRAKEN2_DB"

if [[ -n "${SAMOVAR_KAIJU_DB:-}" ]]; then
  ln -sfn "$SAMOVAR_KAIJU_DB" "$KAIJU_DB"
else
  wget -O "$KAIJU_DB/kaiju_refseq_ref.tar.gz" \
    "${KAIJU_INDEX_URL:-https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_ref_2024-08-14.tgz}"
  tar -xzf "$KAIJU_DB/kaiju_refseq_ref.tar.gz" -C "$KAIJU_DB"
  rm -f "$KAIJU_DB/kaiju_refseq_ref.tar.gz"
fi

if [[ -n "${SAMOVAR_KRAKEN2_DB:-}" ]]; then
  ln -sfn "$SAMOVAR_KRAKEN2_DB" "$KRAKEN2_DB"
elif [[ ! -f "$KRAKEN2_DB/hash.k2d" ]]; then
  wget -O "$OUT/.database/kraken.tgz" \
    "${KRAKEN2_INDEX_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260626.tar.gz}"
  tar -xzf "$OUT/.database/kraken.tgz" -C "$KRAKEN2_DB"
  rm -f "$OUT/.database/kraken.tgz"
fi
```

### 3. Generate metagenomes

```bash
samovar generate \
  --genome_dir "$OUT/.genomes" \
  --host_genome data/test_genomes/host/9606.fna \
  --n_samples 10 \
  --total_reads 100000 \
  --output_dir "$OUT" \
  --cores 32
```

### 4. Preprocess (annotators + pipeline script)

```bash
samovar preprocess \
  --output_dir "$OUT" \
  --kraken2-test "kraken2 $KRAKEN2_DB" \
  --kaiju-test "kaiju $KAIJU_DB" \
  --cores 32 \
  --max_genomes 50
```

`--max_genomes` caps how many abundant taxids annotation2iss will fetch from NCBI (large public DBs can report thousands of hits).

### 5. Execute on Slurm

`examples/realistic/pipeline.sh` writes `$OUT/.log/slurm_exec.sh`. Override with `SLURM_PARTITION`, `SLURM_CPUS`, `SLURM_MEM`, `SLURM_TIME`. Defaults:

```bash
#SBATCH --partition=main
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=48:00:00
```

Submit:

```bash
sbatch "$OUT/.log/slurm_exec.sh"
# or locally:
# samovar exec --output_dir "$OUT"
```

If a previous run left a stale lock after `scancel`:

```bash
snakemake -s workflow/iss_test/Snakefile \
  --configfile "$OUT/.generate/configs/iss_config.yaml" --unlock
rm -rf .snakemake/locks/*
```

---

## Output layout

```
samovar_realistic/
  .genomes/                      # 21 input *-processed.fasta
  .database/{kaiju_db,kraken2_db}/
  .generate/                     # ISS config + generate.sh
  .log/                          # preprocess configs, samovar.sh, Slurm logs
  initial/                       # 10 paired FASTQs
  initial_reports/               # raw Kraken2 + Kaiju
  initial_annotations/           # + plots/
  regenerated/                   # re-simulated reads
  regenerated_annotations/       # + plots/, combined_annotation_table.csv
  reprofiled_annotations/        # ML outputs + roc_comparison.png
  reprofiled_annotations_plots/  # F1 / R2 / CV
```

---

## Reading the figures

Copied example artifacts: [`examples/realistic/results/`](../examples/realistic/results/).

### F1 confusion matrices (all true taxa)

F1 panels keep **every true taxID** on both axes (empty cells = 0). Predictions outside the true set are collapsed to `other`. Unclassified stays `0`.

Stacked overview:

![F1 overview](../examples/realistic/results/F1.png)

Per annotator:

| Annotator | Figure |
|-----------|--------|
| Kaiju | ![F1 kaiju](../examples/realistic/results/F1_kaiju.png) |
| Kraken2 | ![F1 kraken2](../examples/realistic/results/F1_kraken2.png) |
| SAMOVAR (reprofiled) | ![F1 SAMOVAR](../examples/realistic/results/F1_SAMOVAR.png) |

Regenerate after a run with the same workflow script the pipeline uses:

```bash
python workflow/compare_annotations.py \
  --annotation_dir samovar_realistic/reprofiled_annotations \
  --output_dir samovar_realistic/reprofiled_annotations_plots \
  --show_top 0
```

`--show_top 0` means “do not truncate” for CV / confidence plots (also the pipeline default).

### R² and cross-validation

Abundance R² panels:

![R2](../examples/realistic/results/R2.png)

Cross-validation (all predicted taxIDs; pairwise):

| Pair | Figure |
|------|--------|
| Kraken2 vs Kaiju | ![CV k2 vs kaiju](../examples/realistic/results/CV_kraken2_vs_kaiju.png) |
| SAMOVAR vs Kaiju | ![CV SAMOVAR vs kaiju](../examples/realistic/results/CV_SAMOVAR_vs_kaiju.png) |
| SAMOVAR vs Kraken2 | ![CV SAMOVAR vs kraken2](../examples/realistic/results/CV_SAMOVAR_vs_kraken2.png) |

### ML ROC

![ROC](../examples/realistic/results/roc_comparison.png)

### Supporting tables

- [`true_taxa_counts.csv`](../examples/realistic/results/true_taxa_counts.csv) — read counts per true taxID  
- [`combined_annotation_table_head.csv`](../examples/realistic/results/combined_annotation_table_head.csv) — first 200 rows of the combined table  

---

## Pitfalls

1. **Bash `GROUPS`** is reserved → use `ORG_GROUPS`.
2. **`samovar exec`** takes `--output_dir` (not `--output-dir`).
3. **Kaiju RefSeq FMI ~98 G** → request **≥250 G** Slurm memory, or set `SAMOVAR_KAIJU_DB` to a local copy.
4. **NCBI rate limits** → fetch genome groups sequentially; use `--max-genome-mb`.
5. **Protista** → use **Alveolata**.
6. **Thousands of taxids** from public DBs → `--max_genomes 50` (top by abundance).
7. After **`scancel`**, unlock Snakemake before resubmitting.

---

## Reference run log

| Field | Value |
|-------|-------|
| Job | **800190** on `griffin` (`main`, 32 CPU, 250 G) |
| Wall time | ~1 h 23 m |
| Peak RSS | ~108 G |
| Status | `Pipeline execution completed.` |
| Inputs | 21 genomes · 10 × 100k reads |
| Annotators | `kraken2-test`, `kaiju-test` |
| `max_genomes` | 50 |
| ML | RandomForest **0.970**, AdaBoost **0.689** → `trained_model.joblib` |
| True taxa in F1 | **22** (host + microbial set) |

Harvest:

```bash
JOB=$(cat samovar_realistic/.log/slurm_jobid.txt)
sacct -j "$JOB" --format=JobID,State,ExitCode,Elapsed,MaxRSS -P
ls samovar_realistic/reprofiled_annotations/
ls examples/realistic/results/
```
