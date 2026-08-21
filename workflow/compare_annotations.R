library(samovaR)
library(tidyverse)

# Prefer in-repo viz_annotation.R when present (dev tree / editable installs)
local_viz <- file.path(getwd(), "R", "viz_annotation.R")
if (file.exists(local_viz)) {
  source(local_viz)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
annotation_dir <- NULL
show_top <- 0
output_dir <- NULL
types <- c("f1", "R2", "cv")
split <- F
rank <- "genus"
# Parse arguments
i <- 1
csv_file <- NULL
while (i <= length(args)) {
  if (args[i] == "--annotation_dir") {
    annotation_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--show_top") {
    show_top <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--output_dir") {
    output_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--csv") {
    csv_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--rank") {
    rank <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--split") {
    split <- T
    i <- i + 1
  } else {
    i <- i + 1
  }
}

# Check required arguments
if (is.null(annotation_dir)) {
  stop("--annotation_dir is required")
}

# Create output directory if specified
if (!is.null(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# Restore leftovers from a previous interrupted viz, then hide combined
# tables so they are not mixed into the plot.
combined_hidden <- character()
restore_combined_plotbak <- function() {
  leftovers <- list.files(
    annotation_dir,
    pattern = "^combined_annotation_table.*\\.csv\\.plotbak$",
    full.names = TRUE
  )
  for (hidden in unique(c(combined_hidden, leftovers))) {
    orig <- sub("\\.plotbak$", "", hidden)
    if (file.exists(hidden) && !file.exists(orig)) {
      file.rename(hidden, orig)
    }
  }
}
restore_combined_plotbak()

combined_files <- list.files(
  annotation_dir,
  pattern = "^combined_annotation_table.*\\.csv$",
  full.names = TRUE
)
for (f in combined_files) {
  hidden <- paste0(f, ".plotbak")
  if (file.rename(f, hidden)) {
    combined_hidden <- c(combined_hidden, hidden)
  }
}
on.exit(restore_combined_plotbak(), add = TRUE)

# Read and process data
data <- read_annotation_dir(annotation_dir)

# `cv` plots require at least two annotators; with one annotator the internal
# plot list is empty and ggpubr raises subscript errors.
tax_cols <- names(data)[grepl("^taxID_", names(data))]
annotators <- unique(sapply(strsplit(tax_cols, "_"), function(x) if (length(x) >= 2) x[2] else NA))
annotators <- annotators[!is.na(annotators)]
if (length(annotators) < 2) {
  types <- setdiff(types, "cv")
}

# Generate visualizations (best-effort: do not break pipeline on plot-only errors)
results <- tryCatch({
  viz_annotation(
    data = data,
    type = types,
    show_top = show_top,
    output_dir = output_dir,
    plot = FALSE,
    split = split,
    rank = rank
  )
}, error = function(e) {
  warning(sprintf("Skipping visualization due to error: %s", conditionMessage(e)))
  NULL
})

if (!is.null(csv_file)) {
  write.csv(data, csv_file, row.names = FALSE)
}
