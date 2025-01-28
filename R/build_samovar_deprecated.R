#' Build samovar object
#'
#' @description
#' Samovar network is a 2D-oriented graph with metadata and abundances of species per sample. Oriented graph could be used for network prediction, or for better generation some network could be used as initial (to be implemented)
#' For better understanding of building database and using it in generation, visit github source
#' @param samovar_data samovar data after preprocessing stages
#' @param network FALSE or graph that can be used for generation. To be implemented
#' @param distance_function function used for measuring distances between species based on samples
# @param k_means FALSE or number of k-means using for split. recommended for lazy-bootstrap
#' @param min_min_cluster_size FALSE or minimum number of species per cluster
#' @param max_min_cluster_size FALSE or minimum number of species per cluster
#' @param plot_log Logical or path for log plots output
#' @importFrom plotly ggplotly

build_samovar <- function(
    samovar_data,
    dist_function = function(x) dist(x),
    network = F,
    k_means = F,
    min_cluster_size = F,
    plot_log = T) {

  if (isFALSE(k_means)) {
    clustering <- function(df, clust_function) {
      df %>%
        distances() %>%
        clust_function() %>%
        as.numeric()
    }

    sc <- samovar_data$data %>%
      clustering(function(x) sc_clustering(x, size_constraint = min_cluster_size))

  } else {
    # deprecated in future
    ## functions
    recalculate_cluster <- function(a, min_cluster_size, distance) {
      # a - list of (clusters, cluster)
      # min_cluster_size - minimal amount of species in cluster
      clusters <- a[[1]]
      cluster <- a[[2]]

      cl <- clusters == cluster
      # if cluster is too small
      if (sum(cl) < min_cluster_size) {
        distance[cl, cl] <- Inf
        cluster <- clusters[distance[, cl] %>% which.min()]
        clusters[cl] <- cluster

        clusters <- recalculate_cluster(a = list(clusters, cluster), min_cluster_size)
      }
      return(clusters)
    }

    ## clustering
    distance <- samovar_data$data %>%
      dist_function()

    cluster <- distance %>%
      clust_function() %>%
      stats::cutree(k_means)

    ## recalculate clusters
    # NOT STABLE!!! NOT RECOMMENDED!!!
    for (k in 1:k_means) {
      if (k %in% cluster) {
        cluster <- recalculate_cluster(a = list(cluster, k), min_cluster_size, distance %>% as.matrix())
      }
    }
    sc <- cluster %>%
      factor() %>%
      levels() %>%
      as.numeric()
  }

  # re-calculate maximum clusters
  max_iter = 100
  max_cluster <- function(){
    tmp <- (tapply(sc, sc, length))
    return(c(max(tmp), which.max(tmp) %>% names %>% as.numeric()))
  }

  while ((max_cluster()[1] < max_cluster_size) | max_iter > 0){
    mc <-  max_cluster()[2]

    new_sc <- samovar_data$data %>%
      subset(sc == mc) %>%
      distances %>%
      hierarchical_clustering(min_cluster_size)

    if(sum(new_sc) == 0) {
      warning("Error in reducing maximum sample: to similar species!\nIf you want to use bootstrap for samovaR build, it might be better to re-filter species")
      break
    } else {
      sc[sc == mc] <- max(sc) + new_sc + 1
    }

    max_iter = max_iter - 1
  }

  # re-assign cluster names by abundance and viz
  sc_df <- as.data.frame(sc)
  sc_df$sp <- samovar_data$species

  sc_count <- sc_df %>%
    count(sc) %>%
    arrange(-n) %>%
    mutate(new_sc = 1:nrow(.))

  gg <- sc_count %>%
    ggplot(aes(new_sc, n, fill = new_sc)) +
    geom_col(show.legend = F) +
    theme_minimal() +
    scale_fill_viridis_c() +
    xlab("cluster") + ylab("amount") +
    ggtitle("Species per cluster")

  log_plot(plot_log, "cl", gg)

  sc_df <- sc_df %>%
    left_join(sc_count, by = "sc") %>%
    select(sp, new_sc)

  # PCA viz
  gg <- (samovar_data$data %>%
    apply(1, scale) %>%
    prcomp(center = T))[["rotation"]] %>%
    as.data.frame() %>%
    rownames_to_column("sp") %>%
    left_join(sc_df) %>%
    ggplot(aes(PC1, PC2, color = new_sc, text = sp)) +
    geom_point() +
    scale_color_viridis_c() +
    theme_minimal()

  log_plot(plot_log, "PCA", ggplotly(gg))


  # IMPLEMENT ON FUTURE
  if(F) {
  gg <- samovar_data$data %>%
    dist %>%
    hclust %>%
    as.dendrogram() %>%
    plotly::plot_dendro() %>%
    plotly::layout(split = ~sc_df$new_sc) %>%
    hide_legend()
  }

  samovar_data$cluster <- sc_df$new_sc %>% as.character()
  samovar_data$cluster_size <- sc_count$n
  names(samovar_data$cluster_size) <- sc_count$new_sc %>% as.character()

  return (samovar_data)
}
