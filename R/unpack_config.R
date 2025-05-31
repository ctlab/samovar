#' @title Unpack SAMOVAR config
#' @description Unpack SAMOVAR config to a list of parameters
#' @param config_samovar SAMOVAR config
#' @return arguments for samovar_preprocess and samovar_boil
#' @import yaml
#' @example R/examples/config.R
#' @export
unpack_config <- function(config_samovar) {
    # Read SAMOVAR config file
    config_samovar_dict <- yaml::read_yaml(config_samovar, eval.expr=TRUE)
}
