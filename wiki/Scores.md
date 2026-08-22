# Scores

SamovaR scores a taxonomic annotator against a **known true taxid per read** (ISS labels, or a `true` column on the combined annotation table). Regenerated and reprofiled plot folders write the same figures; reprofiled plots also include **SAMOVAR**.

Two families of numbers appear on the figures:

1. **Current plot formulas** — F1 on the confusion heatmap and R² on the abundance scatter, computed exactly as in `viz_annotation` today.
2. **Standard classification metrics** — accuracy, precision, recall, and F1 with sklearn `macro` / `weighted` / `micro` averages, plus the purity scores.

Formulas below use HTML (`<sub>`, `<sup>`) so they render in GitHub, VS Code, and Cursor markdown previews (LaTeX delimiters are not required).

## Contents

- [Where the files land](#where-the-files-land)
- [Current F1 (heatmap)](#current-f1-heatmap)
- [Current R² (abundance)](#current-r-abundance)
- [Standard classification metrics](#standard-classification-metrics)
- [Accuracy purity](#accuracy-purity)
- [F1 purity](#f1-purity)
- [Barplots](#barplots)
- [Worked example](#worked-example)
- [Rank](#rank)

## Where the files land

`samovar exec` calls `workflow/compare_annotations.py` twice:

| Directory | What is scored |
|-----------|----------------|
| `regenerated_annotations_plots/` | Each annotator, plus a **consensus** majority vote |
| `reprofiled_annotations_plots/` | The same tools, consensus, and **SAMOVAR** |

| File | Contents |
|------|----------|
| `scores.png` / `scores.html` | Grouped bars: Accuracy, P-macro, R-macro, F1-macro, F1 (current), R² (current), Accuracy purity, F1 purity |
| `quality_scores.csv` | One row per annotator (all metrics, including F1-micro and F1-weighted) |
| `purity_by_taxon.csv` | Per-true-taxid majority label, recall, precision, F1 |
| `F1_<annotator>.png` | Confusion heatmap; caption is the **current F1 formula** plus standard P/R/F1 |
| `R2_<annotator>.png` | Predicted vs true abundance; caption is the **current R² formula** |

Before scoring, predictions whose taxid is not in the true set are collapsed to `other` (same as the heatmaps), except for **purity**, which uses the raw predicted taxid so a consistent wrong name can still look pure. Unclassified / `0` / `other` are not treated as a true source taxid.

## Current F1 (heatmap)

This is the number printed as **F1 (current)** on `F1_*.png` and in the barplot. It is the formula already used in `viz_annotation` / the original R helper.

Let *n*<sub>t,c</sub> be the number of reads with true taxid *t* and predicted taxid *c*, and

<p align="center"><b>N = Σ<sub>t,c</sub> n<sub>t,c</sub></b></p>

<p align="center"><b>TP = Σ<sub>t</sub> n<sub>t,t</sub></b> &nbsp;&nbsp;(reads where predicted taxid = true taxid)</p>

Every off-diagonal cell is treated as a single shared error pool:

<p align="center"><b>FP + FN = N − TP</b></p>

<p align="center"><b>P = TP / (TP + FP + FN) = TP / N</b></p>

<p align="center"><b>R = P</b> &nbsp;&nbsp;(same off-diagonal mass)</p>

<p align="center"><b>F1<sub>current</sub> = 2PR / (P + R) = TP / N</b></p>

So **current F1 equals per-read accuracy**. That is intentional in the present code, not a separate balanced F1. A consistent wrong name (every *Escherichia* read called *Shigella*) scores **0** here.

The heatmap caption also prints sklearn **P-macro**, **R-macro**, **F1-macro**, and **F1-weighted** (see below) so the figure is not limited to this diagonal identity.

## Current R² (abundance)

This is the number printed as **R² (current)** on `R2_*.png` and in the barplot.

For each taxid *t* that appears as a **true** label:

<p align="center"><b>n<sub>t</sub></b> = number of reads with true taxid *t*</p>

<p align="center"><b>n̂<sub>t</sub></b> = number of reads predicted as *t*</p>

(Predicted labels that are not in the true set were already mapped to `other`. Taxa that appear only as predictions and never as true labels are not extra rows; they only show up if that token exists in the true column.)

<p align="center"><b>R² = 1 − Σ<sub>t</sub> (n̂<sub>t</sub> − n<sub>t</sub>)² / Σ<sub>t</sub> (n<sub>t</sub> − n̄)²</b></p>

where n̄ is the mean of the true abundances *n*<sub>t</sub>. This is a profile-level coefficient of determination, not sklearn `r2_score` on per-read class labels.

- **1** — predicted counts match true counts for every taxid.
- **0** — the predictor is no better than “every taxon has the mean abundance”.
- **&lt; 0** — worse than that mean baseline (typical when a tool piles reads onto the wrong taxa).

## Standard classification metrics

Computed with scikit-learn on the same collapsed true/pred labels as the heatmap (`average=`, `zero_division=0`). Every read has exactly one label, so **F1-micro = accuracy**.

For a class *t*:

<p align="center"><b>TP<sub>t</sub> = n<sub>t,t</sub></b></p>

<p align="center"><b>FP<sub>t</sub> = (Σ<sub>u</sub> n<sub>u,t</sub>) − n<sub>t,t</sub></b> &nbsp;&nbsp;(predicted *t*, true not *t*)</p>

<p align="center"><b>FN<sub>t</sub> = (Σ<sub>c</sub> n<sub>t,c</sub>) − n<sub>t,t</sub></b> &nbsp;&nbsp;(true *t*, predicted not *t*)</p>

<p align="center"><b>P<sub>t</sub> = TP<sub>t</sub> / (TP<sub>t</sub> + FP<sub>t</sub>)</b></p>

<p align="center"><b>R<sub>t</sub> = TP<sub>t</sub> / (TP<sub>t</sub> + FN<sub>t</sub>)</b></p>

<p align="center"><b>F1<sub>t</sub> = 2 P<sub>t</sub> R<sub>t</sub> / (P<sub>t</sub> + R<sub>t</sub>)</b></p>

| Name in CSV / bars | sklearn call | Meaning |
|--------------------|--------------|---------|
| **Accuracy** | `accuracy_score` | Fraction of reads with pred = true. Same value as **F1 (current)** / F1-micro. |
| **P-macro** | `precision_score(..., average="macro")` | Unweighted mean of P<sub>t</sub> over taxa (rare taxa count as much as common ones). |
| **R-macro** | `recall_score(..., average="macro")` | Unweighted mean of R<sub>t</sub>. |
| **F1-macro** | `f1_score(..., average="macro")` | Unweighted mean of F1<sub>t</sub>. The usual “standard F1” when classes are imbalanced. |
| **F1-weighted** | `f1_score(..., average="weighted")` | Mean of F1<sub>t</sub> weighted by support (how many true reads of *t*). |
| **F1-micro** | `f1_score(..., average="micro")` | Global TP / N. Written to CSV; equal to accuracy here. |

Unlike **F1 (current)**, macro-F1 **does** drop when a rare taxon is missed even if overall accuracy stays high.

## Accuracy purity

For each true taxid *t*, look **only** at reads whose true label is *t*. Let *m*(*t*) be the most abundant *classified* prediction among those reads (ties: lexicographically first). Unclassified / `0` / `other` cannot be *m*(*t*); if a taxon is entirely unclassified, its majority count is 0.

<p align="center"><b>P<sub>acc</sub> = Σ<sub>t</sub> n<sub>t, m(t)</sub> / Σ<sub>t</sub> n<sub>t, ·</sub></b></p>

That is the read-weighted mean of per-taxon majority fractions (cluster purity of true taxa). Range 0–1.

This score does **not** require *m*(*t*) = *t*. It only asks how concentrated the calls are *inside* each true taxid.

| Situation | Accuracy purity |
|-----------|-----------------|
| All reads of *t* get the same predicted name, even if that name is wrong | high |
| Reads of *t* split across many predicted names | low |
| Tool dumps *t* into unclassified | 0 for that taxon |

## F1 purity

Accuracy purity is a **recall** of the majority mapping. It ignores reads from *other* true taxa that the tool also assigned to *m*(*t*).

<p align="center"><b>recall<sub>t</sub> = n<sub>t, m(t)</sub> / n<sub>t, ·</sub></b></p>

<p align="center"><b>precision<sub>t</sub> = n<sub>t, m(t)</sub> / n<sub>·, m(t)</sub></b></p>

The denominator of precision is every read predicted as *m*(*t*). Then

<p align="center"><b>F1<sub>t</sub> = 2 · precision<sub>t</sub> · recall<sub>t</sub> / (precision<sub>t</sub> + recall<sub>t</sub>)</b></p>

<p align="center"><b>F1<sub>purity</sub> = Σ<sub>t</sub> F1<sub>t</sub> · n<sub>t, ·</sub> / N<sub>true</sub></b></p>

Range 0–1. F1<sub>purity</sub> ≤ P<sub>acc</sub>. They match when each majority label *m*(*t*) is used only for reads from *t* (a permutation of names, including the identity). They diverge when several true taxa collapse onto the same predicted taxid.

## Barplots

`scores.png` groups these bars per annotator (then **consensus**, then **SAMOVAR**):

| Bar | What it is |
|-----|------------|
| Accuracy | Standard; = current F1 |
| P-macro / R-macro / F1-macro | Standard sklearn macro averages |
| F1 (current) | Heatmap formula TP/N |
| R² (current) | Abundance formula above |
| Accuracy purity / F1 purity | Majority-mapping scores |

A useful ensemble improves **F1-macro** (or current F1) **and** F1 purity together. Purity-only gains can be a single dump taxid.

## Worked example

Six reads, two true taxa `A` and `B`:

| true | A | A | A | B | B | B |
|------|---|---|---|---|---|---|
| pred | A | A | A | A | B | B |

**Current F1 / accuracy / F1-micro** = 5/6 ≈ 0.833 (one `B` read called `A`).

**Per-class (standard):**

- A: P = 3/4 = 0.75, R = 3/3 = 1, F1 ≈ 0.857
- B: P = 2/2 = 1, R = 2/3 ≈ 0.667, F1 = 0.800
- **F1-macro** ≈ 0.829, **F1-weighted** ≈ 0.829

**Purity:** *m*(A)=A, *m*(B)=B. Accuracy purity = 5/6 ≈ 0.833. F1 purity ≈ 0.829 (same numbers as standard F1 here, because the majority labels *are* the true labels).

If instead every `A` and every `B` read were called `A`:

- Current F1 / accuracy = 0.5
- Accuracy purity = 1 (each true taxon is internally unanimous)
- F1 purity = 2/3 (shared majority label `A` hurts precision)
- F1-macro is also low (class B has recall 0)

## Rank

By default plots remap taxIDs to **genus** (`--rank genus`) before scoring. Pass `--rank none` to score exact taxIDs. Combined CSVs on disk are never rewritten; remapping is in-memory only.

Implementation: `src/samovar/scores.py`, plots in `src/samovar/viz_annotation.py`.
