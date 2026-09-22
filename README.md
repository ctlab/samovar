# SAMOVAR <img src="data/img/logos/logo_stable.png" align="right" width="150" alt="SAMOVAR logo">
### Metagenomic ensemble taxonomic annotation, cross-validation, and ML re-profiling

[![ITMO](https://raw.githubusercontent.com/aimclub/open-source-ops/43bb283758b43d75ec1df0a6bb4ae3eb20066323/badges/ITMO_badge.svg)](https://itmo.ru/)
[![conda](https://github.com/ctlab/samovar/actions/workflows/python-package.yml/badge.svg?branch=main)](https://github.com/ctlab/samovar/actions/workflows/python-package.yml?label=conda)
[![stability testing](https://img.shields.io/github/actions/workflow/status/ctlab/samovar/full-integration.yml?label=stability)](https://github.com/ctlab/samovar/actions/workflows/full-integration.yml)
[![version](https://img.shields.io/badge/version-0.11-blue)](pyproject.toml)
[![license](https://img.shields.io/github/license/ctlab/samovar)](LICENSE.md)
[![conda](https://img.shields.io/badge/conda-environment.yml-44A833?logo=anaconda&logoColor=white)](environment.yml)

In metagenomics, we often do not know which tool to use (or, which is much worse - know because they are SOTA). SAMOVAR team try to solve this problem with the automated benchmark based on the real inputed data to include in the model selection process information about the real community properties

Metagenomic classifiers disagree. SAMOVAR treats **multiple annotators as an ensemble**: it runs them on the same reads, cross-validates calls, regenerates in-silico communities from those calls, and trains a supervised **re-profiler** (SAMOVAR) that combines the tools.

What does the tool do? It gets the metagenome input & taxonomy profiling tools and SAMOVAR that (regenerate artficial metagenomes, evaluate & combine the tools).

We strongly recommend to understand [SAMOVAR main concepts](https://github.com/ctlab/samovar/wiki) before the installation & usage, because the workflow is large & depends on a lot of other different tools

## Installation

Python 3.10+; **conda is recommended**.

```bash
git clone https://github.com/ctlab/samovar
cd samovar
conda env create -f environment.yml
conda activate samovar
chmod +x install.sh
./install.sh
```

`install.sh` may ask you some questions, like e-mail for the NCBI API.

Some tools are optional but are required for some actions. Install them all with:

```bash
./install.sh full
```

Already-installed extras are reused and written into the install config. More details: [installation wiki](https://github.com/ctlab/samovar/wiki/How-to-install).

The tool main config layout is described in [config wiki](https://github.com/ctlab/samovar/wiki/Configs-&-data).

## Ensemble annotation

Built-in ensemble members (wired through `samovar prepare`):

| Tool | Role in the ensemble |
|------|----------------------|
| **Kraken2** | k-mer LCA classifier |
| **Kaiju** | protein-level (translated) classifier |
| **Kraken / KrakenUniq** | additional k-mer votes |
| **MetaPhlAn** | marker-gene profiler |
| **Custom** | extra votes; can be easily additionally implemented for the developers |

Workflow:

1. Annotate real or ISS-simulated reads with every configured tool.
2. Cross-validate taxIDs across tools (CV heatmaps) and against known truth when available (F1 / R²).
3. Re-simulate a community from the annotation table (annotation2iss).
4. Re-annotate the synthetic reads and train an ML ensemble (`workflow/ML.py`: RandomForest / AdaBoost) that maps tool votes → corrected taxID.

```mermaid
%%{init: {'theme': 'base', 'themeVariables': { 'fontSize': '16px', 'fontFamily': 'arial', 'primaryColor': '#fff', 'primaryTextColor': '#000', 'primaryBorderColor': '#000', 'lineColor': '#000', 'secondaryColor': '#fff', 'tertiaryColor': '#fff'}}}%%
graph TD
    subgraph Input
        A1[FastQ / ISS genomes]
        A3[Annotator configs]
    end

    subgraph Ensemble
        A1 --> C[Initial annotation]
        A3 --> C
        C --> CV[Cross-validation]
        C --> E[Metagenome regeneration]
        E --> F[Re-annotation]
    end

    subgraph Re-profiling
        C --> ML[Train ensemble]
        F --> ML
        F --> R
        ML --> R[SAMOVAR results]
    end
```

Each stage have built-in & custom options. Their usage & integration approaches ar described in the [tool wiki](https://github.com/ctlab/samovar/wiki/Custom-tools-import).

## Brief Usage

```bash
# Generate metagenome (skip for running SAMOVAR on real data as ensemble)
samovar generate \
    --genome_dir $SAMOVAR/data/test_genomes/meta \
    --host_genome $SAMOVAR/data/test_genomes/host/9606.fna \ # optional
    --output_dir samovar_out
# -> output: bash script for the in silico generated metagenome

# Prepare workflow & scripts, create generation config
samovar prepare \
    --output_dir samovar_out \
    --kraken2-test "kraken2 $DB_KRAKEN2" \ # format: --NAME "type database_path"
    --kaiju-test "kaiju $DB_KAIJU"
# -> output: bash script for the downstream pipeline

# Do SAMOVARing (resumes from `.log/checkpoints`; `--redo` reruns every step)
samovar exec --output_dir samovar_out
# -> pipeline execution

# MultiQC report: optional
# samovar multiqc --output_dir samovar_out -- --export --interactive
# -> interactive multiqc report
```

## R package

The optional R generator (`samovar_boil`) lives on the **[`r-package`](https://github.com/ctlab/samovar/tree/r-package)** branch and is **not** part of this tree. Install it with `./install.sh R-package` (requires R + remotes). That writes a small R driver into `~/.config/samovar/` that only calls exported `samovaR` functions. Then set `regeneration_mode: samovar`.

## References

BibTeX for built-in tools lives in [`cite/`](cite/citations.json) (`cite/*.bib`). Refresh after install with `./install.sh --rebuild-citations 1` if needed (default) or `python -m samovar.citations rebuild`. 

Do not forget to cite every annotator used in the ensemble.

*SAMOVAR itself is now approaching the publication, but have no valid citation.*