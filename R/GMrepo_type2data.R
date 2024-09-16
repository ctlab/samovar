#' Get data from GMrepo
#'
#' @description
#' Wrapper around GMrepo_type2run and GMrepo_run2data functions
#' @param mesh_ids Character. All types of meshID to use. List of relations between meshID and phenotype could be obtained using GMrepo_meshID(). Health meshID by default
#' @param number_to_process False by default, or maximum number of runs per meshID
#' @param number_to_out False by default, maximum number of obtained data
#' @param at_level "species" by default. level to obtain classification from GMrepo
# @param keep_metadata To be implemented. Keep metadata from query
#' @param QC_filter QCStatus by default. Perform auto QC filtering based on metadata column, or False for no checking.
#' @import tidyverse
#' @example R/examples/GMrepo.R
#' @export

GMrepo_type2data <- function(
    mesh_ids = c("D006262"),
    number_to_process = F,
    number_to_out = F,
    at_level = "species",
    QC_filter = "QCStatus") {
  GMrepo_type2run(
    mesh_ids = mesh_ids,
    number_to_process = number_to_process
  ) %>%
    GMrepo_run2data(
      number_to_out = number_to_out,
      at_level = at_level,
      QC_filter = QC_filter
    ) %>%
    return
}
