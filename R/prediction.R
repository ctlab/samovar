prediction <- function (init, init_level,
                        Rsq_cl, Rsq_il, Rsq_pr,
                        kmeans,
                        res_sp,
                        res_scale,
                        res_sp_scale,
                        trsh,
                        species,
                        min_res,
                        inner_model = "gaussian",
                        inter_model = "gaussian") {
  res <- data.frame()
  cl_todo <- 1:max(kmeans)
  cl_done <- c()
  
  #check is there some cluster else to do and start with init cluster
  cl <- kmeans[which(init == species)]
  
  #rescale init
  init_level <- res_scale(init_level) - min_res_sp_scale
  
  which.max.coord <- function(df, pred) {
    df <- as.data.frame(df)
    df[,-pred] <- 0
    df[pred,] <- 0
    
    df_row <- apply(df, 1, max)
    df_row <- which.max(df_row)[1]
    
    df_col <- apply(df, 2, max)
    df_col <- which.max(df_col)[1]
    
    return(c(df_col, df_row))
  }
  
  prob_cl <- function(a, b){
    and <- ((a > 0) & (b > 0)) %>% sum
    or  <- ((a > 0) | (b > 0)) %>% sum
    return(and/or)
  }
  
  while (length(cl_todo)>0) {
    #make predictions for relevant cluster
    #for now: easiest prediction just based on one group. 
    #in future, find the best path must be based on ML
    
    #first: make cluster tables for R2 and means 
    sp_cl <- species[cl == kmeans]
    df_k <- res_sp_scale[cl == kmeans,]
    df_r_cl <- Rsq_cl[[cl]] %>% as.matrix
    df_r_pr <- Rsq_pr[[cl]] %>% as.matrix
    df_r_pr[is.na(df_r_pr)] <- 0
    
    res_cl <- data.frame("sp" = sp_cl,
                         "res_scale" = NA)
    
    #prediction of initial level 
    if(!is.na(init)){
      #for represented cluster
      pred <- data.frame(from = which(init == sp_cl), to = which(init == sp_cl))
      res_cl[which(sp_cl == init),2] <- init_level
      
      #second: start predictions based on biggest R2
      while(sum(is.na(res_cl[,2]))) {
        
        #prediction graph
        pr <- which.max.coord(df_r_pr, pred$to)
        pred <- rbind(pred, pr)
        
        #КОСТЫЛЬ!!!
        if(sum(pr == c(1,1)) == 2) {
          res_cl[is.na(res_cl)] <- 0
        }
        
        #prediction of presence
        if(runif(1) <= df_r_pr[pr[1], pr[2]]) {
          
          df <- df_r_cl
          df[-pred$from,] <- 0
          best <- which.max(df[,pr[2]])
          
          #fit GLM
          dfi <- df_k[best,]
          dfj <- df_k[pr[2],]
          
          dfij <- ((dfi > 0) & (dfj > 0)) %>% which
          
          if(length(dfij) < 2){
            #edit in future; for now if there are less than 2 connected observations for 2 groups level of 2 group <- 0
            res_cl[pr[2],2] <- 0
            
          } else {
            dfij <- data.frame("i" = as.numeric(dfi[dfij]), 
                               "j" = as.numeric(dfj[dfij]))
            glmij <-glm(j~i, data = dfij, family = inner_model)
            res_cl[pr[2],2] <- predict.glm(glmij, newdata = data.frame(i = res_cl[pr[1],2]))
          }
          
        } else {
          #the best way is bayesian correction if prob <- 0, so we should re-celculate:
          # prob <- w(!a&b)/w(!a|b)
          
          #but for now, lets just remove that probabilities
          df_r_pr[pr[2],] <- 0
          df_r_pr[,pr[2]] <- 0
          
          res_cl[pr[2],2] <- 0
        }
        
      }
      
    } else {
      #for not represented cluster
      res_cl$res_scale <- 0
    }
    
    #add results
    res <- rbind(res, res_cl)
    cl_todo <- cl_todo[-cl]
    cl_done <- c(cl_done, cl)
    
    cat("%",cl,"cl %")
    
    #select new cluster for generation and predict initial levels
    #it is also better to use bayesian correction
    if(length(cl_todo) > 0) {
      cl_new <- which.max.coord(Rsq_il, cl_done)
      
      #костыль
      if((cl_new[1] == 1) & (cl_new[2] == 1)) break
      
      df_k1 <- df_k
      if(!is.null(ncol(df_k1))) df_k1 <- apply(df_k1, 2, mean)
      
      df_k <- res_sp_scale[cl_new[2] == kmeans,]
      df_k2 <- df_k
      if(!is.null(ncol(df_k2))) df_k2 <- apply(df_k2, 2, mean)
      
      df_i12 <- (df_k1 > 0) & (df_k2 > 0)
      df_o12 <- (df_k1 > 0) | (df_k2 > 0)
      
      #random generation of cluster presence based on Bayesian model
      if(runif(1) < (sum(df_i12)/sum(df_o12))) {
        
        #mean of cluster level
        dfk12 <- data.frame("i" = as.numeric(df_k2[df_i12]),
                            "j" = as.numeric(df_k1[df_i12]))
        glmij <- glm(i~j, data = dfk12, family = inter_model)
        cl_level <- predict.glm(glmij, newdata = data.frame(j = mean(res_cl$res_scale)), type = "response")
        
        #select init
        sp_cl <- species[cl_new[2] == kmeans]
        
        if(length(sp_cl) == 1) {
          res <- rbind(res, c(sp_cl, cl_level))
        } else {
          probs <- apply(df_k, 1, prob_cl, df_k1)
          probs <- probs / sum(probs)
          init <- sample(sp_cl, size = 1, prob = probs)
          df_ki <- df_k[which(init == sp_cl),] %>% as.numeric
          df_k2i <- data.frame("i" = as.numeric(df_ki),
                               "j" = as.numeric(df_k2))
          glmij <- glm(i~j, data = dfk12, family = inter_model)
          init_level <- predict.glm(glmij, newdata = data.frame(j = mean(res_cl$res_scale)), type = "response")
        }
        
      } else {
        sp_cl <- species[cl_new[2] == kmeans]
        init <- NA
      }
      
      cl <- cl_new[2]
    }
  }
  
  #КОСТЫЛЬ!!!
  res <- res[!duplicated(res$sp),]
  
  #rescale
  res[is.na(res)] <- 0
  res$res_scale <- as.numeric(res$res_scale)
  
  cat("\nrescaling...")
  
  #if(sum(res$res_scale) > -min_res_sp_scale) res$res_scale <- -min_res_sp_scale * res$res_scale / sum(res$res_scale)
  res$res_scale <- sapply(res$res_scale, res_unscale)
  
  res$res_scale <- res$res_scale %>% as.numeric
  
  res <- rbind(res, c("unclassified", 1 - sum(res$res_scale, na.rm = T)))
  res_uncl <- res[res$sp %in% c("Unknown", "unclassified"),]
  res <- res[!(res$sp %in% c("Unknown", "unclassified")),]
  res <- res[order(res$res_scale),]
  res <- rbind(res, c("unclassified", sum(as.numeric(res_uncl$res_scale))))
  
  res[is.na(res)] <- 0
  res$res_scale <- res$res_scale %>% as.numeric
  
  cat("\ndone...")
  
  return(res)
}
