#' @title Preprocess SAMOVAR data
#' @description Wrapper for SAMOVAR preprocess commands
#' @param samovar_data samovar data object
#' @inheritParams teatree_trim
#' @inheritParams tealeaves_pack
#' @inheritParams teabag_brew
#' @inheritParams concotion_pour
#' @param ... Additional arguments, passed
#' @return Build SAMOVAR object
#' @example R/examples/processing.R
#' @export
samovar_preprocess <- function(
    samovar_data,

    #teatree_trim
    metadata_filter = F,
    treshhold_amount = 10^(-5),
    treshhold_samples = 1,
    treshhold_species = 1,
    drop_species = F,
    drop_unclassified = T,

    #tealeaves_pack
    normalization_function = function(x) log10(x + 1),
    plot_log = T,

    #teabag_brew
    dist_function = function(x) dist(x),
    network = F,
    min_cluster_size = 2,
    max_cluster_size = 100,

    #concotion_pour
    inner_method = "glm",
    inter_method = "glm",
    inner_model = "gaussian",
    inter_model = "gaussian",
    probability_calculation = "oriented",
    cluster_connection = "mean",

    #other
    ...
) {
  data <- teatree_trim(
    samovar_data$copy(),
    metadata_filter,
    treshhold_amount,
    treshhold_samples,
    treshhold_species,
    drop_species,
    drop_unclassified,
    ...
  )

  data <- tealeaves_pack(
    data$copy(),
    normalization_function,
    plot_log,
    ...
    )

  data <- teabag_brew(
    data$copy(),
    dist_function,
    network,
    min_cluster_size,
    max_cluster_size,
    ...
    )

  data <- concotion_pour(
    data$copy(),
    inner_method,
    inter_method,
    inner_model,
    inter_model,
    probability_calculation,
    cluster_connection,
    ...
    )

  return(data)
}
