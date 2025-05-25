library(samovar)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
annotation_dir <- NULL
show_top <- 10
output_dir <- NULL
types <- c("f1", "R2", "cv")

# Parse arguments
i <- 1
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

# Read and process data
data <- read_annotation_dir(annotation_dir)

# Generate visualizations
results <- viz_annotation(
  data = data,
  type = types,
  show_top = show_top,
  output_dir = output_dir
)

if (!is.null(csv_file)) {
  write.csv(data, csv_file, row.names = FALSE)
}
