#!/usr/bin/env Rscript
# Export phyloseq::GlobalPatterns OTU counts as a SamovaR abundance CSV.
# Usage: Rscript scripts/export_globalpatterns.R [output.csv] [n_taxa]
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "tests/data/globalpatterns_abundance.csv"
n_taxa <- if (length(args) >= 2) as.integer(args[[2]]) else 80L
suppressPackageStartupMessages(library(phyloseq))
data(GlobalPatterns)
otu <- as(otu_table(GlobalPatterns), "matrix")
if (!taxa_are_rows(GlobalPatterns)) otu <- t(otu)
keep <- rowSums(otu > 0) >= 2
otu <- otu[keep, , drop = FALSE]
ord <- order(rowSums(otu), decreasing = TRUE)
otu <- otu[ord[seq_len(min(n_taxa, nrow(otu)))], , drop = FALSE]
df <- data.frame(taxid = rownames(otu), otu, check.names = FALSE)
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(df, out, row.names = FALSE)
message("Wrote ", nrow(df), " taxa x ", ncol(df) - 1L, " samples to ", out)
