#' Scale species abundances
#'
#' @param samovar_data Samovar data object to rescale
#' @param normalization_function Function using for rescaling
#' @param plot_log Logical or path for log plots output
#' @example R/examples/preprocessing.R
# @importFrom plotly ggplotly
#' @export

tealeaves_pack = function(
    samovar_data,
    normalization_function = function(x) log10(x + 1),
    plot_log = T
) {

  data = samovar_data$copy()
  data$normalize(normalization_function)

  # prepare plots
  trun <- function(x, min_data) {
    x <- x %>% unlist
    x <- x[x > min_data]
    if (length(x) > 10^4 ) {
      warning("Only 10 000 random values represented on plot")
      x[runif(10^4, min = 1, max = length(x))]
    }
    x <- x %>%
      as.data.frame() %>%
      mutate(n_bin = cut(., breaks = 100, label = F, include.lowest = T))

    return(x)
  }

  ## assign data
  g1 <- samovar_data$data %>%
    trun(samovar_data$min_value)  %>%
    mutate(name = "before")

  g2 <- data$data %>%
    trun(data$min_value) %>%
    mutate(name = "after")


  g_count <- (g1) %>%
    rbind(g2) %>%
    count(ord = n_bin, name) %>%
    mutate(n = (n/sum(n)*2) )

  # check normalization function
  norm_df <- data.frame(N = seq(samovar_data$min_value, samovar_data$max_value, length.out = 100)) %>%
    mutate(before = N) %>%
    mutate(after = N %>% normalization_function %>% minmaxscale)

  if ((!all(norm_df$after == cummax(norm_df$after))) |
      (!all(abs(norm_df$after != Inf)))) {
    warning("Normalization function might be not monotone on [0,1] or produce Inf")
  }

  norm_df <- norm_df %>%
    pivot_longer(cols = -1) %>%
    mutate(ord = (order(N)+1) %/% 2) %>%
    left_join(g_count, by = c("ord", "name")) %>%
    mutate(name = fct_inorder(name))

  ggN <- ggplot(data = norm_df, aes(value, N)) +
    geom_point(aes(color = n), size = .5) +
    xlab("range") +
    ylab("value") +
    facet_grid(~name) +
    scale_color_viridis_c(name = "Amount", na.value = "transparent", direction = -1) +
    theme_minimal() +
    ggtitle("Filtration comparison")

  ## QQ-plot
  gg <- g1 %>% rbind(g2) %>%
    mutate(name = fct_inorder(name)) %>%
    ggplot() +
    facet_wrap(~name, scales = "free") +
    xlab("value") +
    theme_minimal()

  ggQQ <- gg +
    geom_qq(aes(sample = .)) +
    geom_qq_line(aes(sample = .))  +
    ggtitle("QQ-plot; only values > 0")

  ggD <- gg +
    geom_density(aes(.)) +
    ggtitle("Density plot; only values > 0")

  # render plots
  log_plot(plot_log, "norm", ggN)
  log_plot(plot_log, "dens", ggD)
  log_plot(plot_log, "QQ", ggQQ)

  return(data)
}
