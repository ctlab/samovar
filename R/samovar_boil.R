#' Generate artificial data
#'
#' @description
#' Use pre-built samovar_data with its parameters
#'
#' @param samovar_base samovar data after preprocessing and building stages
#' @param n number of artificial samples to generate
#' @param init_sp species vector for initializing data generation, or FALSE for automatic usage
#' @param init_value species amount vector (values from 0 to 1) for initializing data generation, or FALSE for automatic usage
#' @example R/examples/processing.R
#' @export

samovar_boil <- function(samovar_base, n = 1,
                         init_sp, init_value) {
  return(cat("Generation done"))
  #__ ----
  # Additional functions ----

  #__ ----

  # Main ----
  ## Initializing all ----
  res <- data.frame()

  ## make iterations ----
  for(iter in 1:n) {
    ## Initializing iter ----
    clusters_todo <- 1:samovar_base$samovar_data$cluster_n()

    ## Generation loop ----

    ### Re-initializing cluster ----

    ### Generating cluster ----

    ### Predict new cluster ----


  }
}
