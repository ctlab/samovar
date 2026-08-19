library(samovaR)
library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
annotation_dir <- NULL
show_top <- 10
output_dir <- NULL
types <- c("f1", "R2", "cv")
split <- F
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

# Installed samovaR may concatenate every CSV in the directory, including a
# previously written combined_annotation_table.csv. Hide those files so they
# cannot be mixed back into the plot (that is how true taxID 0 reappeared).
combined_files <- list.files(
  annotation_dir,
  pattern = "^combined_annotation_table.*\\.csv$",
  full.names = TRUE
)
combined_hidden <- character()
for (f in combined_files) {
  hidden <- paste0(f, ".plotbak")
  if (file.rename(f, hidden)) {
    combined_hidden <- c(combined_hidden, hidden)
  }
}
on.exit({
  for (hidden in combined_hidden) {
    orig <- sub("\\.plotbak$", "", hidden)
    if (file.exists(hidden) && !file.exists(orig)) {
      file.rename(hidden, orig)
    } else if (file.exists(hidden)) {
      unlink(hidden)
    }
  }
}, add = TRUE)

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
    split = split
  )
}, error = function(e) {
  warning(sprintf("Skipping visualization due to error: %s", conditionMessage(e)))
  NULL
})

if (!is.null(csv_file)) {
  write.csv(data, csv_file, row.names = FALSE)
}
