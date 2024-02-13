build_samovar <- function(
    data_scaled,
    k_means = 1,
    inner_model = "gaussian",
    inter_model = "gaussian",
    minimal_cluster = 2,
    probability_calculation = "oriented", #oriented or simple
    cluster_connection = "mean" #mean_oriented, mean, simple or PCA
) {
  
  # Helpful functions ----
  
  recalculate_cluster <- function(a, minimal_cluster = 2) {
    #a - list of (clusters, cluster)
    #minimal_cluster - minimal amount of species in cluster
    clusters <- a[[1]] 
    cluster  <- a[[2]]
    
    cl = clusters == cluster
    #if cluster is too small
    if (sum (cl) < minimal_cluster) {
      distance[cl,cl] <- Inf
      cluster <- clusters[distance[,cl] %>% which.min]
      clusters[cl] <- cluster
      
      clusters <- recalculate_cluster(a = list(clusters, cluster), minimal_cluster)
    }
    return (clusters)
  }
  
  
  get_rsq <- function(X, Y, prob, model) {
    #correct Rsq to number of cases
    if (sum(prob) == 0) {
      return(0)
    } else  if (sum(prob) == 1) {
      return(1) #correct in future!!!!!!
    } else {
      glmij <- glm(Y[prob] ~ X[prob], family = model)
      Rsq <- with(summary(glmij), 1 - deviance/null.deviance)
      return(Rsq)
    }
  }
  
  ####
  # Parse data ----
  
  df <- data_scaled$data
  
  ####
  # Clusters calculation ----
  
  species = df %>% row.names
  distance = df %>% 
    dist(method = "euclidian") 
  
  clust <- distance %>% 
    hclust (method = "ward.D2") 
  
  distance <- distance %>% as.matrix
  
  clusters <- cutree(clust, k = k_means)
  
  ####
  # Recalculate clusters ----
  
  for (cluster in 1:k_means) {
    if (cluster %in% clusters) {
      clusters <- recalculate_cluster(a = list(clusters, cluster), minimal_cluster)
    }
  }
  clust_list <- clusters %>% factor %>% levels %>% as.numeric
  
  #rename clusters by order ----
  for (j in 1:length(clust_list)) {
    clusters[clusters == clust_list[j]] <- j
    clust_list[j] <- j
  }
  
  cat("---  Removed", k_means-j, "clusters, kept", j, " ---\n")
  k_means <- j
  
  ####
  # Build -----
  ## In different clusters ----
  cat("---  Prepare tea leaves for different clusters  ---\n\n")
  
  Rs_cl <- list()
  Pr_cl <- list()
  
  for (cluster in clust_list) {
    #prepare data
    cl = clusters == cluster
    df_cluster <- df[cl,]
    df_r2 <- matrix(nrow = sum(cl), ncol = sum(cl), data = 0)
    df_pr <- matrix(nrow = sum(cl), ncol = sum(cl), data = 0)
    
    #calculate
    for (i in 1:sum(cl)) {
      X = df[i,] %>% as.numeric
      eX = X > 0
      for (j in 1:sum(cl)) {
        if (i != j) {
          Y = df[j,] %>% as.numeric
          eY = Y > 0
          
          if (probability_calculation == "oriented") {
            #PR(Y|X) = PR(Y)[PR(X)]
            #calculate probability of case j if i = True
            eY1 = eY[eX] > 0
            df_pr[i,j] <- sum(eY1) / sum(eX)
            df_r2[i,j] <- get_rsq(X, Y, eX&eY, inner_model)

          } else {
            #get the prob = W(a&b)/W(a|b)
            df_pr[i,j] = sum (eX&eY) / sum (eX|eY)
            df_r2[i,j] <- get_rsq(X, Y, eX&eY, inner_model)
          }
        }
      }
    }
    
    Rs_cl[[cluster]] <- df_r2
    Pr_cl[[cluster]] <- df_pr
    
    cat("%", cluster, "%")
  }
  
  ## Between clusters ----
  cat("\n\n---  Prepare tea leaves for calculation between clusters  ---\n\n")
  
  Rs_il <- matrix(ncol = k_means, nrow = k_means)
  Pr_il <- matrix(ncol = k_means, nrow = k_means)
  
  for (i in clust_list) {
    for (j in clust_list) {
      if (i == j) {
        rs = 0
        prob = 0
      } else {
        X = df[clusters == i,]
        Y = df[clusters == j,]
        
        #check cluster connections
        if (cluster_connection == "mean") {
          
          X = X %>% apply(2, mean)
          Y = Y %>% apply(2, mean)
          eX = X > 0
          eY = Y > 0
          
          prob = sum (eX & eY) / sum (eX | eY)
          rs = get_rsq (X, Y, (eX & eY), inter_model)
          
        } else if (cluster_connection == "mean_oriented") {
          
          X = X %>% apply(2, mean)
          Y = Y %>% apply(2, mean)
          eX = X > 0
          eY = Y > 0
          eY1 = eY[eX] > 0
          
          prob = sum(eY1) / sum(eX)
          rs = get_rsq (X, Y, eX & eY, inter_model)
          
        } else if (cluster_connection == "PCA") {
          print("To be implemented")
          exit()
          
        } else {
          print("To be implemented")
          exit()
        }
      }
      Rs_il[i,j] <- rs
      Pr_il[i,j] <- prob
    }
    cat("%", i, "%")
  }
  
  cat("\n\n---  Samovar built  ---\n")
  
  return(
    list(
      "data" = df,
      "samovar" = list(
        "Rsq_inner" =  Rs_cl,
        "Prob_inner" = Pr_cl,
        "Rsq_inter" =  Rs_il,
        "Prob_inter" = Pr_il), 
      "parameters" = list(
        "trim" = data_scaled$parameters$trim,
        "normalize" = data_scaled$parameters$normalize,
        "samovar" = list(
          "inner_model" = inner_model, 
          "inter_model" = inter_model, 
          "k_means" = k_means,
          "clusters" = clusters,
          "probability_calculation" = probability_calculation,
          "cluster_connection" = cluster_connection
          )
        )
      )
    )
}



