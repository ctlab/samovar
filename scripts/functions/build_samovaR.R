build_samovar <- function(
    data_scaled,
    k_means = 1,
    inner_model = "gaussian",
    inter_model = "gaussian"
) {

  # calculate clusters ----
  
  species = rownames(data_scaled)
  clust = data_scaled %>% 
    dist(method = "euclidian") %>% 
    hclust (method = "ward.D2") 
  
  clusters <- cutree(clust, k = k_means)
  
  # BUILD -----
  ##### In different clusters ----
  cat("Calculating R square tables and GLMs for different clusters\n")
  
  Rsq_cl <- list()
  Rsq_pr <- list()
  
  # subsetting tables and make 
  for (k in 1:max(clusters)) {
    clust_num <- (clusters == k)
    if (sum(clust_num) > 1) {
      df_k <- data_scaled[clust_num,]
      ldfk <- length(df_k[,1])
      
      df_r2 <- matrix(ncol = ldfk, nrow = ldfk)
      df_pr <- matrix(ncol = ldfk, nrow = ldfk)
      
      for(i in 1:ldfk){
        j <- i
        while(j <= ldfk) {
          if(i == j) {
            df_r2[i,j] <- 0
          } else {
            
            dfkj <- df_k[j,]
            dfki <- df_k[i,]
            dfNA <- !((dfki>0) & (dfkj>0))
            dfBOTH <- ((dfki>0) | (dfkj>0)) %>% sum
            dfkj <- dfkj[!dfNA]
            dfki <- dfki[!dfNA]
            
            #get the prob = W(a&b)/W(a|b)
            probij <- (sum(!dfNA) / dfBOTH)
            if(sum(!dfNA) > 1) {
              glmij <- glm(dfki ~ dfkj, family = inner_model)
              df_r2[i,j] <- with(summary(glmij), 1 - deviance/null.deviance)
            } else {
              df_r2[i,j] <- 0
            }
            df_pr[i,j] <- probij
          }
          j <- j + 1
        }
      }
      df_r2 <- Matrix::forceSymmetric(df_r2)
      df_pr <- Matrix::forceSymmetric(df_pr)
    } else {
      df_r2 <- 0
      df_pr <- 0
    }
    
    Rsq_cl[[k]] <- df_r2
    Rsq_pr[[k]] <- df_pr
    
    cat("%", k, "%")
  }
  
  ##### Between clusters ----
  cat("\nCalculating R square tables and GLMs between clusters\n")
  
  Rsq_il <- matrix(ncol = k_means, nrow = k_means)
  
  for (i in 1:k_means){
    df_k1 <- data_scaled[clusters == i,]
    if(!is.null(ncol(df_k1))) df_k1 <- apply(df_k1, 2, mean) %>% unlist %>% as.numeric
    j <- i
    
    while (j <= k_means) {
      
      if(i == j) {
        Rsq_il[i,j] <- 0
        j <- j + 1
      } else {
        df_k2 <- data_scaled[clusters == j,]
        if(!is.null(ncol(df_k2))) df_k2 <- apply(df_k2, 2, mean) %>% unlist %>% as.numeric
        
        glmij <- glm(df_k1 ~ df_k2, family = inter_model)
        Rsq_il[i,j] <- with(summary(glmij), 1 - deviance/null.deviance)
        j <- j + 1
      }
    }
  }
  
  Rsq_il <- Matrix::forceSymmetric(Rsq_il) %>% as.matrix
  cat("Done")
  
  return(list(Rsq_cl,Rsq_pr,Rsq_il, list(inner_model, inter_model, k_means)))
}