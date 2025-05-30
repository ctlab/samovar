#' Convert samovar object to phyloseq object
#'
#' @param samovar_data A samovar object
#' @return A phyloseq object
#' @export
samovar2phyloseq <- function(samovar_data) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required. Please install it using BiocManager::install('phyloseq')")
  }

  # Extract abundance table
  abund_table <- samovar_data$data

  # Create phyloseq object
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(abund_table, taxa_are_rows = TRUE)
  )

  return(ps)
}

#' Convert phyloseq object to samovar object
#'
#' @param phyloseq_data A phyloseq object
#' @return A samovar object
#' @export
phyloseq2samovar <- function(phyloseq_data) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required. Please install it using BiocManager::install('phyloseq')")
  }

  # Extract abundance table
  abund_table <- as.data.frame(phyloseq::otu_table(phyloseq_data))

  # Create samovar object
  samovar_data <- table2samovar(abund_table)

  return(samovar_data)
}
