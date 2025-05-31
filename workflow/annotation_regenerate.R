library(tidyverse)
library(samovar)

# Default values
annotation_dir <- NULL
output_dir <- NULL
config_samovar <- NULL

# Parse arguments
i <- 1
while (i <= length(args)) {
  if (args[i] == "--annotation_dir") {
    annotation_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--output_dir") {
    output_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--config_samovar") {
    config_samovar <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(config_samovar)) {
  config <- list()
} else {
  # Parse SAMOVAR config file
  config <- unpack_config(config_samovar)
}

# Fix config
if (!("N" %in% names(config)) ) N <- 1
if (!("N_reads" %in% names(config)) ) N_reads <- 100
if("output_dir" %in% names(config)) output_dir <- config$output_dir
if("input_dir" %in% names(config)) input_dir <- config$input_dir

# Check required arguments
if (is.null(annotation_dir)) {
  stop("--annotation_dir is required")
}

# Create output directory if specified
if (!is.null(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# Process taxonomy tables to SAMOVAR tables
samovar_data_long <- read_annotation_dir(annotation_dir)
samovar_data_list <- annotation2samovar(samovar_data_long)

res <- list()
# Process each SAMOVAR data object
for (i in 1:length(samovar_data_list)) {
  # Process SAMOVAR data object
  config$samovar_data <- samovar_data_list[[i]]$copy()
  samovar <- do.call(
    samovar_preprocess,
    config
  )
  new_data <- samovar_boil(samovar, N = config$N)
  res[[ names(samovar_data_list)[i] ]] <- new_data$data * config$N_reads
  write.csv(
    res[[ names(samovar_data_list)[i] ]],
    paste0(output_dir, "/", names(samovar_data_list)[i], ".csv")
  )
}

