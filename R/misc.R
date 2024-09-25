#' Misc functions
#'
#' @import tidyverse

## Min-max scaling ----
minmaxscale <- function(x) {
  x <- (x - min(x))
  x / max(x)
}

## Best transform ----
best_transform <- function(normalization_function) {
  if(is.function(normalization_function)) {
    return(normalization_function)
  } else {
    warning("Not implemented auto-detection of normalization function")
  }
}

## Reorder data.frame for viz ----
reorder_df <- function(data, dim = 1, reord = "fpc") {
  if(dim == 1) {
    tmp <- data
  } else {
    tmp <- data %>% t
  }

  if (reord == "hcl") {
    ord <- (tmp %>% dist %>% hclust)$order
  } else if (reord == "tsne") {
    if (requireNamespace("tsne")) {
      ord <- tmp %>% apply(2, scale) %>%  tsne::tsne(k = 1, max_iter = 500) %>% order
    }
  } else if (reord == "fpc_scaled") {
    ord <- (tmp %>% apply(2, scale) %>% prcomp())[["x"]][1,]
  } else if (reord == "fpc") {
    ord <- (tmp %>% prcomp())[["x"]][,1] %>% order
  } else if (reord == "amount") {
    ord <- tmp %>% apply(1, sum) %>% order(decreasing = T)
  } else if (reord %in% c("none", F)) {
    ord <- 1:nrow(tmp)
  }

  if(dim == 1) {
    return(data[ord,])
  } else {
    return(data[,ord])
  }
}

# Transform to data ----
samovar2data <- function(data) {
  data <- switch(class(data),
                 "data.frame" = data,
                 "tibble" = data,
                 "matrix" = data,
                 "samovar_run" = data$data,
                 "samovar_data" = data$data,  #%>% data$reverse_normalize_df,
                 "samovar_base" = data$samovar_data$data) %>%  #%>% data$samovar_data$reverse_normalize_df) %>%
    data.frame

  data[is.na(data)] <- 0

  data[,apply(data, 2, sum) > 1] <- data[,apply(data, 2, sum) > 1] /
    apply(data[,apply(data, 2, sum) > 1], 2, sum)

  data["unclassified",] <- apply(data, 2, function(x) 1 - sum(x))
  return(data)
}

#___ ----
# concotion_pour ----

## Misc ----
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

prob_calc_general <- function(X, Y, min_value = 0, probability_calculation) {
  X <- occurrence(X)
  Y <- occurrence(Y)

  #calculate probs
  if (probability_calculation == "oriented") {
    # PR(Y|X) = PR(Y)[PR(X)]
    Y1 <- Y[X] > min_value
    return(sum(Y1) / sum(X))
  } else if (probability_calculation == "simple") {
    # get the prob = W(a&b)/W(a|b)
    return(sum(X & Y) / sum(X | Y))
  } else if (probability_calculation == "compositional"){
    return(sum(X & Y) / sum(X | Y))
  } else {
    warning("Not implemented")
    #break
  }
}

## Cluster connection
summarise_cluster_general <- function(x, cluster_connection) {
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

## GLM calculation ----
glm_calc_general <- function(X, Y,
                             model, summarise_cluster,
                             mode = "glm") {
  prob = occurrence(X) & occurrence(Y)

  #deal with df
  X <- X %>% summarise_cluster()
  Y <- Y %>% summarise_cluster()

  # correct Rsq to number of cases
  if (sum(prob) == 0) {
    return(0)
  } else if (sum(prob) == 1) {
    return(X[prob] / Y[prob])
  } else {
    if(mode == "glm"){
      xy <- data.frame(x = X[prob], y = Y[prob])
      return(glm(y~x, data = xy, family = model))
    } else if(mode == "lm") {
      xy <- data.frame(x = X[prob], y = Y[prob])
      return(lm(y~x, data = xy))
    }
  }
}

# Rsq calculation based on glm ----
rsq_calc_general <- function(X, Y, model, glm_calc) {
  glmxy <- glm_calc(X, Y, model)
  if(is.numeric((glmxy))) return(0.00001) # correct in future!!!!!!

  Rsq <- with(summary(glmxy), 1 - deviance / null.deviance)
  return(Rsq)
}

#__ ----

# samovar_boil() ----
## maximum coordinate ----
which.max.coord <- function(df, predicted) {
  df[,predicted] <- 0
  df[-predicted,] <- 0

  if (max(df) == 0) return(c(NA,NA))

  which(df == max(df), arr.ind=T)[1,] %>% as.numeric()
}

## Get reverse function ----
rf <- function(x) {
  if(is.na(x)) return (0)
  if(x < 0) return (0)
  if(x > 1) return (1)

  x = x * (pn$max_data - pn$min_data) +  pn$min_data
  x = pn$reverse_function(x)
  return(x)
}


## Get new value ----
get_newval_general <- function(X, Y,
                               xval,
                               model,
                               mode = "glm",
                               summarise_cluster,
                               glm_calc){

  if (mode == "glm") {
    glmxy <- glm_calc(X, Y, model)
    xval <- summarise_cluster(xval)

    if(is.numeric(glmxy)) return(glmxy*xval)

    yval <- try(predict(glmxy, newdata = data.frame(x = xval)),
                silent = T) #!!!!!!
    if(is.character(yval)) yval <- 0
  } else {
    warning("Not implemented")
  }

  return (as.numeric(yval))
}

# Presence calculation ----
presence_calculation_general <- function(xPr_cl, sp_done,
                                         df_cl, init_sp,
                                         probability_calculation,
                                         force_simple = F,
                                         from_distribution = T) {
  if(force_simple) if(probability_calculation %in% "compositional") probability_calculation <- "simple"

  if(from_distribution) {
    edf <- xPr_cl %>% unlist
    runif_cluster <- sample(edf, nrow(xPr_cl)-1)
  } else {
    runif_cluster <- runif(nrow(xPr_cl)-1)
  }

  switch(probability_calculation,
         "simple" = {
           xPr_cl[sp_done,-sp_done] > runif_cluster
         },

         "oriented" = {
           xPr_cl[sp_done,-sp_done] > runif_cluster
         },

         "compositional" = {
           sm <- (df_cl > 0)[-sp_done,df_cl[init_sp,] > 0]
           sm[,sample(1:ncol(sm), 1)]
         },

         "bayesian" = warning("Not implemented yet...")
  )
}
