#' Misc functions
#'
#' @import tidyverse

minmaxscale <- function(x) {
  x <- (x - min(x))
  x / max(x)
}
