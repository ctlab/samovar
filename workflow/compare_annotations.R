library(samovaR)
library(tidyverse)
library(data.table)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
annotation_dir <- NULL
show_top <- 100
output_dir <- NULL
types <- c("f1", "R2", "cv")
split <- FALSE

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
    split <- TRUE
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

# ==========================================
# SUPERVISOR'S FIX: Custom robust reading function
# Upgraded with data.table::fread for 5M reads performance
# ==========================================
read_annotation_dir_custom <- function(data_dir, sample_name_position = 0) {
  results <- tibble()
  files <- dir(data_dir, pattern = "\\.csv$|\\.tsv$", full.names = TRUE)
  
  if (length(files) == 0) {
    stop(paste("[ERROR] No annotation files found in", data_dir))
  }
  
  for (data_path in files) {
    # Extract sample name
    basename_parts <- stringr::str_split(basename(data_path), "\\.")[[1]]
    sample_name <- paste(basename_parts[1:(sample_name_position + 1)], collapse = "_")
    
    # Fast reading
    tmp <- data.table::fread(data_path, data.table = FALSE, fill = TRUE, sep = "\t")
    tmp <- tmp %>% mutate(sample = sample_name)
    
    # Standardize column names based on supervisor's logic
    colnames(tmp) <- tolower(colnames(tmp)) # Lowercase all first
    colnames(tmp) <- colnames(tmp) %>% 
      stringr::str_remove("_[0-9]*$") %>% 
      stringr::str_replace("taxid", "taxID")
    
    results <- dplyr::bind_rows(results, tmp)
  }
  return(results)
}

cat(sprintf("[INFO] Reading annotations from %s...\n", annotation_dir))
data <- read_annotation_dir_custom(annotation_dir)

# ==========================================
# SUPERVISOR'S FIX 2: Bypassing the viz_annotation edge-case bug
# Prevent NAs in 'true' column from crashing the package
# ==========================================
if ("true" %in% colnames(data)) {
  data$true <- as.character(data$true)
  data$true[is.na(data$true)] <- "0"
} else if ("seq" %in% colnames(data)) {
  # If dealing with InSilicoSeq generated reads
  data$true <- stringr::str_extract(data$seq, "(?<=taxid:)[0-9]+")
  data$true[is.na(data$true)] <- "0"
}

# Filter out NA rows if they slipped through
data <- data[!is.na(data$seq), ]

# ==========================================
# Generate visualizations
# ==========================================
cat("[INFO] Generating visualizations...\n")
results <- viz_annotation(
  data = data,
  type = types,
  show_top = show_top,
  output_dir = output_dir,
  split = split
)

# ==========================================
# RESTORED: Save final CSV for ML.py
# ==========================================
if (!is.null(csv_file)) {
  cat(sprintf("[INFO] Saving combined annotation table to %s...\n", csv_file))
  # Using \t as separator because ML.py expects a tab-separated file
  write.table(data, file = csv_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("[SUCCESS] compare_annotations.R finished successfully.\n")