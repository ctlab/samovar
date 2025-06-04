#' Visualize composition
#'
#' @param data data.frame with dimensions of species * samples, or samovar objects: samovar_data, samovar_base, samovar_run or GMrepo_run (with data)
#' @param top integer, number of top-represented taxa to show, or FALSE to show all
#' @param reord_samples character, fpc, fpc_scaled, hcl, amount, tsne, or none reorder of samples on plot
#' @param reord_species character, same for reord_samples
#' @param interactive logical. ggplot or plotly object to return
#' @param ggplot_add functions to add to ggplot object, or FALSE.
#' @param bottom_legend vector length of samples to show on plot as a color legend, or FALSE
#' @param type character, column or tile for composition visualize, or donut (to implement) for mean composition visualization
#'
#' @example R/examples/data2viz.R
#' @export

viz_composition <- function(data,
                            reord_samples = "fpc",
                            reord_species = "amount",
                            type = "column",
                            top = 15, interactive = F,
                            ggplot_add = F, bottom_legend = F) {
  palette <- c("#FFFFFF","#FFFFFF", "#E7F09550", "#F5BFC7", "#90EE90", "#006400", "cadetblue", "darkblue")

  #configure
  data <- data %>% samovar2data

  #legend
  if(!isFALSE(bottom_legend)) {
    legend <- data.frame(samples = colnames(data), col = bottom_legend)
  }

  #reorder and subset
  if(!isFALSE(top)) {
    sumdata <- data %>% apply(1, sum) %>% order(decreasing = F)
    if (top < length(sumdata)) {
      data <- data[sumdata <= top,]
    }
  }
  data["other",] <- data %>% apply(2, function(x) 1 - sum(x))
  data_raw <- data

  if(type == "tile") {
    sp <- rownames(data)
    sm <- colnames(data)

    data <- data %>%
      apply(1, scale) %>%
      t %>%
      data.frame()

    rownames(data) <- sp
    colnames(data) <- sm
  }

  data <- data %>%
    reorder_df(reord_samples, dim = 2) %>%
    reorder_df(reord_species, dim = 1) %>%
    rownames_to_column("sp") %>%
    pivot_longer(names_to = "samples", values_to = "value", cols = -1) %>%
    mutate(sp = fct_inorder(sp),
           samples = fct_inorder(samples))

  data$value[is.na(data$value)] <- min(data$value, na.rm = T)

  data_raw <- data_raw %>%
    reorder_df(reord_samples, dim = 2) %>%
    reorder_df(reord_species, dim = 1) %>%
    rownames_to_column("sp") %>%
    pivot_longer(names_to = "samples", values_to = "raw_value", cols = -1) %>%
    mutate(sp = fct_inorder(sp),
           samples = fct_inorder(samples),
           raw_value = ifelse(
             raw_value < 0,
             0,
             raw_value
           )
    )

  data <- left_join(data, data_raw, by = c("sp", "samples"))
  data$value[data$raw_value == 0] <- min(data$value)

  if(!isFALSE(bottom_legend)) data <- left_join(data, legend, by = "samples")

  #!!!!! add adjusting colors to other taxa

  gg <- data %>%
    ggplot() +
    theme_minimal() +
    xlab("") + ylab("") +
    ggtitle("Composition, scaled values") +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      axis.text.y = element_text(),
      axis.text.x = element_text(angle = 90))

  if(type == "column") {
    gg <- gg +
      geom_col(aes(x = samples, y = value, fill = sp,
                   text = map(paste0(
                     "<b>Sample:</b> ", samples, "<br>",
                     "<b>Taxa:</b> ", "<i>", sp, "</i>", "<br>",
                     "<b>Value:</b> ", value
                   ),
                   ~HTML(.)
                   ))) +
      scale_fill_viridis_d("scaled values", direction = -1)

  } else if(type == "tile") {
    gg <- gg +
      geom_tile(aes(x = samples, fill = value, y = sp, text = raw_value)) +
      #coord_flip() +
      #scale_fill_viridis_c("scaled values", direction = -1)
      scale_fill_gradientn(colours = palette)
  }


  if(!isFALSE(bottom_legend)) {
    gg <- gg +
      geom_rug(aes(y = -0.01, col = col))
  }

  if(!isFALSE(ggplot_add)) {
    gg <- gg + ggplot_add
  }

  if(interactive) {
    gg <- gg +
      theme(legend.position = "none")

    plotly::ggplotly(gg)%>%
      plotly::layout(title="Composition",
                     hoverlabel = list(
                       bgcolor = "white",
                       font = list(size = 12),
                       align = "left"),
                     xaxis = list(title = '', showgrid = F, showticklabels = F),
                     yaxis = list(title = '', showgrid = F, tickfont = list(style = "italic")))

  } else {
    return(gg)
  }
}
