#' Build an example samovar_data object for tests and vignettes
#'
#' Creates a reproducible abundance matrix when GMrepo is unavailable.
#'
#' @param n_species Number of species to simulate.
#' @param n_samples Number of samples to simulate.
#' @param seed Random seed for reproducibility.
#' @return A \code{samovar_data} object.
#' @export
example_samovar_data <- function(n_species = 50, n_samples = 20, seed = 1) {
  set.seed(seed)
  abundance <- matrix(
    stats::runif(n_species * n_samples, min = 0.001, max = 1),
    nrow = n_species,
    dimnames = list(
      paste0("Species_", seq_len(n_species)),
      paste0("sample", seq_len(n_samples))
    )
  )
  abundance <- apply(abundance, 2, function(x) x / sum(x))

  metadata <- data.frame(
    body_site = sample(c("gut", "oral", "skin"), n_samples, replace = TRUE),
    row.names = colnames(abundance),
    stringsAsFactors = FALSE
  )

  table2samovar(abundance, metadata = metadata, min_sp = 0, min_samp = 0)
}

.packaged_abundance_path <- function(filename) {
  path <- system.file("testdata", filename, package = "samovaR")
  if (!nzchar(path) || !file.exists(path)) {
    path <- file.path("inst", "testdata", filename)
  }
  if (!file.exists(path)) {
    return(NA_character_)
  }
  path
}

#' Load a packaged abundance table as samovar_data
#'
#' @param filename File name in \code{inst/testdata}.
#' @return A \code{samovar_data} object, or NULL if the file is missing.
#' @noRd
abundance_rds_to_samovar <- function(filename) {
  path <- .packaged_abundance_path(filename)
  if (is.na(path)) {
    return(NULL)
  }
  obj <- readRDS(path)
  metadata <- obj$metadata
  if (is.null(metadata) || !is.data.frame(metadata) || nrow(metadata) == 0) {
    metadata <- FALSE
  }
  table2samovar(obj$data, metadata = metadata, min_sp = 0, min_samp = 0)
}

#' Download GMrepo data with an offline fallback
#'
#' Tries the GMrepo API, then packaged example abundance tables, then a
#' synthetic matrix from \code{\link{example_samovar_data}}.
#'
#' @param ... Arguments passed to \code{\link{GMrepo_type2data}}.
#' @param fallback Packaged RDS name in \code{inst/testdata} used when GMrepo fails.
#' @return A \code{samovar_data} object.
#' @export
GMrepo_type2data_or_example <- function(..., fallback = "teatree_abundance.rds") {
  result <- tryCatch(
    GMrepo_type2data(...),
    error = function(e) e
  )
  valid <- inherits(result, "samovar_data") &&
    is.data.frame(result@data) &&
    nrow(result@data) > 0 &&
    ncol(result@data) > 0
  if (valid) {
    return(result)
  }

  packaged <- abundance_rds_to_samovar(fallback)
  if (!is.null(packaged)) {
    message("GMrepo unavailable; using packaged example abundance table.")
    return(packaged)
  }

  message("GMrepo unavailable; using synthetic example data.")
  example_samovar_data()
}
