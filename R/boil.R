boil <- function(samovar, 
                 init = F, 
                 init_level = F, 
                 n = 1, 
                 pref = F,
                 occurance_species = "simple",
                 occurance_clusters = "simple",
                 calculate_unclassified = F){
  # Samovar unpacking ----
  df <- samovar$data
  Rs_cl <- samovar$samovar$Rsq_inner
  Pr_cl <- samovar$samovar$Prob_inner
  Rs_il <- samovar$samovar$Rsq_inter
  Pr_il <- samovar$samovar$Prob_inter
  pn <- samovar$parameters$normalize
  clusters <- samovar$parameters$samovar$clusters
  k_means <- samovar$parameters$samovar$k_means
  inner_model <- samovar$parameters$samovar$inner_model
  inter_model <- samovar$parameters$samovar$inter_model
  probability_calculation <- samovar$parameters$samovar$probability_calculation
  cluster_connection <- samovar$parameters$samovar$cluster_connection
  
  # Helpful functions ----
  which.max.coord <- function(df, predicted) {
    df <- as.matrix(df)
    df[,predicted] <- 0
    df[-predicted,] <- 0
    
    if (max(df) == 0) return(c(NA,NA))
    
    dfl <- nrow(df)
    dfm <- which.max(df)[1]
    
    row <- dfm %% dfl
    if (row == 0) row <- dfl
    
    col <- dfm %/% dfl + 1
    if (row == dfl) col <- col - 1
    
    return(c(row,col))
  }
  
  # get accurate reverse function
  rf <- function(x) {
    if(is.na(x)) return (0)
    if(x < 0) return (0)
    if(x > 1) return (1)
    
    x = x * (pn$max_data - pn$min_data) +  pn$min_data
    x = pn$reverse_function(x)
    return(x)
  }
  
  #get glm
  get_newval <- function(X, Y, coex = "both", model, yval){
    X <- as.numeric(X)
    Y <- as.numeric(Y)
    yval <- as.numeric(yval)
    
    if (coex == "both") {
      eX <- X > 0
      eY <- Y > 0
      XY <- eX&eY
    } else {
      XY <- TRUE
    }
    
    if(sum(XY) == 0) return(0)
    if(sum(XY) == 1) Y[XY]/X[XY] * yval
    
    glm_df <- tibble("X" = X[XY], "Y" = Y[XY])
    glmij <- glm(formula = X~Y, data = glm_df, model)
    xval <- predict(glmij, data.frame(Y = yval, X = NA))
    return (as.numeric(xval))
  }
  
  # Get initial ----
  ## If initials not specified, mean of most occurred species is used
  if (isFALSE(init)) {
    init <- (df > 0) %>% apply(1, sum) %>% which.max %>% names
  }
  if (isFALSE(init_level)) {
    init_level <- df[rownames(df) == init, ] %>% 
      as.numeric %>% 
      mean %>% 
      rf
  }
  
  #CHECK!!!!!
  init_level <- (samovar$parameters$normalize$normalisation_function(init_level) -
                 samovar$parameters$normalize$min_data) /
                (samovar$parameters$normalize$max_data - 
                 samovar$parameters$normalize$min_data)
  
  inits <- c(init, init_level)
  
  # Generate ----
  ## Prepare ----
  cat("---  Generation started  ---\n\n")
  results <- data.frame("sp" = rownames(df))
  
  ## Main loop ----
  for (counter in 1:n) {
    ### initialize ----
    init <- inits[1]
    init_level <- as.numeric(inits[2])
    res <- data.frame()
    cluster <- clusters[which(init == rownames(df))]
    cluster_len <- sum(clusters == cluster)
    cl_done <- c()
    
    ### calculate cluster occurrence in stream ----
    if (occurance_clusters == "simple") {
      xPr_il <- Pr_il[cluster,]
      xPr_il[cluster] <- 1
      cl_0 <- which(xPr_il < runif(k_means))
      xRs_il <- Rs_il
      xRs_il[cl_0,] <- 0
      xRs_il[,cl_0] <- 0
      cl_done <- c(cl_done, cl_0)
    } else {
      cat ("To be implemented")
      exit()
    }
    
    ### calculations species amounts ----
    while(length(cl_done) < k_means) {
      #### prediction of the new initial level ----
      # for any cluster except of the first initial
      if(is.na(init)){
        #### get best cluster to predict ----
        pr <- which.max.coord(xRs_il, cl_done)
        #### check ----
        if(is.na(pr[1])){
          cat("---  Warning: probably incomplete data on generation", 
              cluster, "-", length(cl_done), "of", k_means, "clusters were generated  ---")
          cl_done <- 1:k_means
          #all probabilities turn to 0
          cluster <- 0
        } else {
          cluster <- pr[2]
        }
        #### predict ----
          if (cluster_connection %in% c("mean", "mean_oriented")) {
            X <- df[clusters == pr[1],] %>% apply(2, mean)
            Y <- df[clusters == pr[2],] %>% apply(2, mean)
            newmean <- get_newval(X,Y, 
                                  model = inter_model,  
                                  yval = mean(res[res$cl == pr[1],2]))
            
            #get random species from new cluster and its mean
            init <- sample(rownames(df)[clusters == cluster], 1)
            init_mean <- mean(df[rownames(df) == init, (X>0)&(Y>0)] %>% as.numeric)
            clus_mean <- mean(Y[(X>0)&(Y>0)])
            init_level <- init_mean/clus_mean * newmean
            
          } else {
            cat ("To be implemented")
            exit()
            }
        }
      
      #### generate cluster ----
      #if not error in prediction graph
      if (cluster != 0) {
        ##### unpack ----
        df_cl  <- df[clusters == cluster,]
        xRs_cl <- Rs_cl[[cluster]]
        xPr_cl <- Pr_cl[[cluster]]
        
        res_cl <- data.frame("sp" = rownames(df_cl), "l" = NA) #results
        res_cl$cl <- cluster
        res_pr <- data.frame("sp" = rownames(df_cl), "l" = 0)  #occurrence, logical
        
        ##### initialize ----
        sp_done <- which(res_cl$sp == init)
        sp_0 <- c(0)
        res_cl[sp_done,2] <- init_level
        res_pr[sp_done,2] <- 1
        
        ##### generate occurrence table ----
        if (occurance_species == "simple") {
          X <- xPr_cl[sp_done,]
          for (i in (1:cluster_len)[-sp_done]) {
            # КОСТЫЛЬ!!!----
            if (i > length(X)) i <- length(X)
            if (X[i] > runif(1)) {
              res_pr[i,2] <- 1
            } else {
              sp_0 <- c(sp_0, i)
            }
          }
        } else {
          cat ("To be implemented")
          exit()
        }
        sp_done <- c (sp_done, sp_0[sp_0 != 0])
        
        ##### generate prediction graph ----
        pred <- data.frame(from = which(res_cl$sp == init), to = NA)
        pr <- 0
        xRs_cl[,sp_0] <- 0
        xRs_cl[sp_0,] <- 0
        
        while((length(sp_done) < cluster_len) & !is.na(pr[1])) {
          #probs prediction graph
          pr <- which.max.coord(xRs_cl, pred$from)
          pred <- rbind(pred, pr)
          sp_done <- c(sp_done, pr[2])
        }
        pred <- pred[!is.na(pred$to),]
        
        ##### generate data ----
        # might be improved if factorisation of pred$from is used
        if(nrow(pred)>0){
          for (i in 1:nrow(pred)) {
            X = df_cl[pred[i,1],] %>% as.numeric
            Y = df_cl[pred[i,2],] %>% as.numeric
            if (probability_calculation %in% c("oriented", "simple")) {
              res_cl[pred[i,2],2] <- get_newval(X, Y, 
                                                model = inner_model,  
                                                yval = res_cl[pred[i,1],2])
            }
          }
        }
        ##### merge predictions ----
        res_cl[is.na(res_cl)] <- 0
        cl_done <- c(cl_done, cluster)
        res <- rbind(res, res_cl)
        init <- NA
        init_level <- NA
      }
    }
    ### merge generations ----
    
    if(!(isFALSE(pref))) { 
      colnames(res)[2] <- paste0(pref,counter)
    } else {
      colnames(res)[2] <- counter
    }
    
    cat("%", counter, "%")
    results <- full_join(results, res[,-3], by = "sp")
  }
  ## Return and correct data ----
  results[,-1] <- results[,-1] %>% sapply(rf)
  sumres <- apply(as.data.frame(results[,-1]), 2, sum)
  
  for (i in which(sumres > 1) + 1) {
    results[,i] <- results[,i] / sumres[i-1] }
  
  if (calculate_unclassified) {
    sumres <- apply(results, 2, sum)
    results[nrow(results) + 1,-1] <- 1 - sumres
    results[nrow(results), 1] <- "Unclassified"
  }

  cat("\n\n---  Tea done!  ---\n")
  return(results)
}