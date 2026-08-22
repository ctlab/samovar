#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import methods
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr across arrange bind_rows coalesce count full_join group_by
#' @importFrom dplyr left_join mutate mutate_all n pull select summarise sym all_of
#' @importFrom tibble tibble as_tibble column_to_rownames rownames_to_column
#' @importFrom tidyr pivot_longer pivot_wider everything
#' @importFrom forcats fct_inorder
#' @importFrom stringr str_detect str_remove str_replace str_split
#' @importFrom plotly ggplotly hide_legend
#' @importFrom httr POST content timeout http_error set_config config status_code
#' @importFrom jsonlite fromJSON
#' @importFrom xml2 xml_text read_xml
#' @importFrom Matrix forceSymmetric
#' @importFrom distances distances
#' @importFrom htmlwidgets saveWidget
#' @importFrom progress progress_bar
#' @importFrom ggnewscale new_scale_color
#' @importFrom yaml read_yaml
#' @importFrom methods is new setClass setGeneric setMethod show slot slotNames
#' @importFrom stats as.dendrogram dist glm hclust lm na.omit prcomp predict
#' @importFrom stats runif uniroot
#' @importFrom utils read.csv write.csv head
#' @importFrom grDevices colorRampPalette
#' @import scclust
## usethis namespace: end
NULL

# dplyr / ggplot2 NSE symbols used in package code
utils::globalVariables(c(
  ".", "Freq", "N", "PC1", "PC2", "label_color", "n_bin", "name",
  "predicted taxID", "raw_value", "rect", "relative_abundance",
  "samples", "scientific_name", "sp", "true", "true taxID", "value",
  "x", "y"
))
