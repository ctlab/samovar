#' Progress bar
#'
#' @import progress

progress_function <- function (iters) {
  if (length(iters) != 1) iters = length(iters)

  progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                   total = iters,
                   complete = "=",
                   incomplete = "-",
                   current = ">",
                   clear = FALSE,
                   width = 100)
}
