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

#' Download GMrepo data with an offline fallback
#'
#' @param ... Arguments passed to \code{\link{GMrepo_type2data}}.
#' @return A \code{samovar_data} object.
#' @export
GMrepo_type2data_or_example <- function(...) {
  tryCatch(
    GMrepo_type2data(...),
    error = function(e) {
      message("GMrepo unavailable; using synthetic example data.")
      example_samovar_data()
    }
  )
}
