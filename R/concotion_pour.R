#' Build samovar object
#'
#' @description
#' Samovar network is a 2D-oriented graph with metadata and abundances of species per sample. Oriented graph could be used for network prediction, or for better generation some network could be used as initial (to be implemented)
#' For better understanding of building database and using it in generation, visit github source
#' @param samovar_data samovar data after preprocessing stages
#' @param network FALSE or graph that can be used for generation. To be implemented
#' @param inner_method Character, glm, other to be implemented (bootstrap, bsPCA)
#' @param inter_method Character, glm, other to be implemented (bootstrap, bsPCA)
#' @param inner_model Character, model processed by glm(). For glm mode only. quasipoisson by default
#' @param inter_model Character, model processed by glm(). For glm mode only. quasipoisson by default
#' @param cooccurrence Character, co-occurrence calculation.
#' <br>
#' <br>If "simple", calculated as:
#' <br>P(A|B) = sum(A&B)/sum(A|B).
#' <br>
#' <br>If "oriented", calculated as
#' <br>P(A|B) = P(A&B|B)
#' <br>
#' <br>If "compositional", calculated as
#' <br>P(A|B) = P(A&B|B), and than sampled one of represented conditions of occurence
#' @param cluster_connection Character (mean, median), or function. The way of cluster connection. If function, way of summarize all samples of species cluster
#' @param ... Additional arguments, passed
#' @example R/examples/preprocessing.R
#' @export

concotion_pour <- function(
    samovar_data,
    inner_method = "glm",
    inter_method = "glm",
    inner_model = "gaussian",
    inter_model = "gaussian",
    probability_calculation = "oriented", # oriented or simple
    cluster_connection = "mean", # mean_oriented, mean, simple or PCA
    ...
) {
  # Additional functions re-import
  prob_calc <- function(...) prob_calc_general(..., probability_calculation = probability_calculation)
  summarise_cluster <- function(...) summarise_cluster_general(..., cluster_connection = cluster_connection)
  glm_calc <- function(...) glm_calc_general(..., summarise_cluster = summarise_cluster)
  rsq_calc <- function(...) rsq_calc_general(..., glm_calc = glm_calc)
  symmetric <- function(df) df %>% t %>% Matrix::forceSymmetric() %>% as.matrix

  # Parse data ----

  #df <- samovar_data$data
  min_value <- samovar_data$min_value
  clust_list_raw <- unique(samovar_data$cluster)

  # Remove clusters with 0 length
  clust_list <- c()
  for (clust in clust_list_raw){
    cl_len <- NULL
    suppressWarnings(
      try({
        df_cluster <- samovar_data$get_clean(clust)
        cl_len <- nrow(df_cluster)
      })
    )
    if(!is.null(cl_len)) clust_list <- c(clust_list, clust)
  }

  ####
  ## Build -----
  cat("---  Prepare tea leaves for different clusters  ---\n")

  #progress bar for all clusters calculation
  pb <- progress_function(clust_list)
  ### In different clusters ----

  if (inner_method == "glm") {
    #### glm ----

    Rs_cl <- list()
    Pr_cl <- list()

    for (cl in clust_list) {
      pb$tick()
      # prepare data
      df_cluster <- samovar_data$get_clean(cl)
      cl_len <- nrow(df_cluster)

      df_r2 <- prep_matrix(cl_len)
      df_pr <- prep_matrix(cl_len)

      # calculate
      for (i in 1:(cl_len-1)) {
        X <- df_cluster[i, ]
        for (j in (i+1):cl_len) {
          Y <- df_cluster[j, ]
          df_pr[j, i] <- prob_calc(X, Y, min_value)
          df_r2[j, i] <- rsq_calc(X, Y, inner_model)
        }
      }
      if(sum(df_pr) == 0) warning(paste0("May be problems with some species in cluster ", cl,"\nIt is better to re-filter or increase min_size of cluster"))

      Rs_cl[[cl]] <- df_r2 %>% symmetric
      Pr_cl[[cl]] <- df_pr %>% symmetric
    }
  } else {
    warning("Not implemented yet")
  }

  ### Between clusters ----
  cat("\n---  Prepare tea leaves for calculation between clusters  ---\n")

  #progress bar for connection between clusters calculation
  pb <- progress_function(length(clust_list)-1)

  if(samovar_data$cluster_n() > 1) {
    if(inter_method == "glm") {
      cl_len <- length(clust_list)
      Rs_il <- prep_matrix(cl_len)
      Pr_il <- prep_matrix(cl_len)

      for (i in 1:(cl_len-1)) {
        X <- samovar_data$get(clust_list[i])
        for (j in (i+1):cl_len) {
          Y <- samovar_data$get(clust_list[j])
          Pr_il[j,i] <- prob_calc(X, Y, min_value)
          Rs_il[j,i] <- rsq_calc(X, Y, inter_model)
        }
        pb$tick()
      }
    }
  } else {
    Rs_il <- matrix(1)
    Pr_il <- matrix(1)
  }


  # Return ----
  samovar_base_new <- new(
    "samovar_base",
    samovar_data = samovar_data,
    inner_cluster_graph_method = Rs_cl,
    inter_cluster_graph_method = Rs_il %>% symmetric,
    inner_cluster_graph_prob = Pr_cl,
    inter_cluster_graph_prob = Pr_il %>% symmetric,
    preferences = list(
      inner_method = inner_method,
      inter_method = inter_method,
      inner_model = inner_model,
      inter_model = inter_model,
      probability_calculation = probability_calculation,
      cluster_connection = cluster_connection
    )
  )
  cat("\n---  Samovar built  ---\n\n")
  return(samovar_base_new)
}
