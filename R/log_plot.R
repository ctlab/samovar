#' Print a log plot
#' @noRd

log_plot <- function(plot_log, postfix, gg, mode = "ggplot", write = F) {
  if(!write) {
    if(!isFALSE(plot_log)) {
      if(is.character(plot_log)) {
        if (mode == "ggplot") {
          ggsave(paste0(plot_log, postfix, ".png"), gg)
        } else if (mode == "plotly") {
          htmlwidgets::saveWidget(gg, paste0(plot_log, postfix, ".html"))
        }
      } else {
        print(gg)
      }
    }
  }
}
