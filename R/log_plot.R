#' Print a log plot
#'
#' @importFrom htmlwidgets saveWidget

log_plot <- function(plot_log, postfix, gg, mode = "ggplot") {
  if(!isFALSE(plot_log)) {
    if(is.character(plot_log)) {
      if (mode == "ggplot") {
        ggsave(paste0(plot_log, postfix, ".png"), gg)
      } else if (mode == "plotly") {
        htmlwidgets::saveWidget(as_widget(gg), paste0(plot_log, postfix, ".html"))
      }
    } else {
      print(gg)
    }
  }
}
