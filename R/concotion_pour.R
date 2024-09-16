#' Build samovar object
#'
#' @description
#' Samovar network is a 2D-oriented graph with metadata and abundances of species per sample. Oriented graph could be used for network prediction, or for better generation some network could be used as initial (to be implemented)
#' For better understanding of building database and using it in generation, visit github source
#' @param samovar_data samovar data after preprocessing stages
#' @param network FALSE or graph that can be used for generation. To be implemented
#' @param inner_method Character, glm, other to be implemented (bootstrap, bsPCA)
#' @param inter_method Character, glm, other to be implemented (bootstrap, bsPCA)
#' @param inner_model Character, model processed by glm(). For glm mode only
#' @param inter_model Character, model processed by glm(). For glm mode only
#' @param cooccurrence Character, co-occurrence calculation.
#' <br>
#' <br>If "simple", calculated as:
#' <br>P(A|B) = sum(A&B)/sum(A|B).
#' <br>
#' <br>If "oriented", calculated as
#' <br>P(A|B) = P(A&B|B)
#' @param cluster_connection Character (mean, median), or function. The way of cluster connection. If function, way of summarize all samples of species cluster
#' @example R/examples/preprocessing.R
#' @export

concotion_pour <- function(
    samovar_data,
    inner_method = "glm",
    inter_method = "glm",
    inner_model = "gaussian",
    inter_model = "gaussian",
    minimal_cluster = 2,
    probability_calculation = "oriented", # oriented or simple
    cluster_connection = "mean" # mean_oriented, mean, simple or PCA
    ) {
  #__ ----
  # Additional functions ----

  ## misc ----
  prep_matrix = function(cl) {
    cl <- as.numeric(cl)
    matrix(nrow = cl, ncol = cl, data = 0)
  }

  ## Prob ----
  # calculate probabilities
  occurrence <- function(X, min_value = 0){

    #deal with matrix
    matrix_to_prob <- function(x){
      if(is.data.frame(x)) {
        x <- apply(x, 2, sum)
      }
      x %>% as.numeric()
    }

    X <- matrix_to_prob(X) > min_value

    return(X)
  }

  prob_calc <- function(X, Y, min_value = 0) {
    X <- occurrence(X)
    Y <- occurrence(Y)

    #calculate probs
    if (probability_calculation == "simple") {
      # PR(Y|X) = PR(Y)[PR(X)]
      Y1 <- Y[X] > min_value
      return(sum(Y1) / sum(X))
    } else if (probability_calculation == "oriented") {
      # get the prob = W(a&b)/W(a|b)
      return(sum(X & Y) / sum(X | Y))
    } else {
      warning("Not implemented")
      #break
    }
  }

  ## Cluster connection
  summarise_cluster <- function(x) {
    if (is.numeric(x)){
      return(x)
    } else if (nrow(x) == 1) {
      return(x %>% as.numeric)
    } else if (cluster_connection == "mean"){
      return(x %>% apply(2, mean))
    } else if (is.function(cluster_connection)) {
      return(x %>% apply(2, cluster_connection))
    } else {
      warning("Not implemented yet")
    }
  }

  ## GLM ----
  glm_calc <- function(X, Y, model) {
    prob = occurrence(X) & occurrence(Y)

    #deal with df
    X <- X %>% summarise_cluster()
    Y <- Y %>% summarise_cluster()

    # correct Rsq to number of cases
    if (sum(prob) == 0) {
      return(0)
    } else if (sum(prob) == 1) {
      return(1) # correct in future!!!!!!
    } else {
      return(glm(Y[prob] ~ X[prob], family = model))
    }
  }

  rsq_calc <- function(X, Y, model) {
    glmxy <- glm_calc(X, Y, model)
    if(is.numeric((glmxy))) return(glmxy)

    Rsq <- with(summary(glmxy), 1 - deviance / null.deviance)
    return(Rsq)
  }


  #__ ----
  # Parse data ----

  #df <- samovar_data$data
  min_value <- samovar_data$min_value
  clust_list <- unique(samovar_data$cluster)

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

      Rs_cl[[cl]] <- df_r2
      Pr_cl[[cl]] <- df_pr
    }
  } else {
    warning("Not implemented yet")
  }

  ### Between clusters ----
  cat("\n---  Prepare tea leaves for calculation between clusters  ---\n")

  #progress bar for connection between clusters calculation
  pb <- progress_function(length(clust_list)-1)

  if(inter_method == "glm") {
    cl_len <- length(clust_list)
    Rs_il <- prep_matrix(cl_len)
    Pr_il <- prep_matrix(cl_len)

    for (i in 1:(cl_len-1)) {
      X <- samovar_data$get(clust_list[i])
      for (j in (i:cl_len)) {
        Y <- samovar_data$get(clust_list[j])
        Pr_il[j,i] <- prob_calc(X, Y, min_value)
        Rs_il[j,i] <- rsq_calc(X, Y, inter_model)
      }
      pb$tick()
    }
  }

  # Return ----
  samovar_base_new <- new(
    "samovar_base",
    samovar_data = samovar_data,
    inner_cluster_graph_method = Rs_cl,
    inter_cluster_graph_method = Rs_il,
    inner_cluster_graph_prob = Pr_cl,
    inter_cluster_graph_prob = Pr_il,
    preferences = list(
      inner_method = inner_method,
      inter_method = inter_method,
      inner_model = inner_model,
      inter_model = inter_model,
      minimal_cluster = minimal_cluster,
      probability_calculation = probability_calculation,
      cluster_connection = cluster_connection
    )
  )
  cat("\n---  Samovar built  ---\n")
  return(samovar_base_new)
}
