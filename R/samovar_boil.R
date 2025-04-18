#' Generate artificial data
#'
#' @description
#' Use pre-built samovar_data with its parameters
#'
#' @param samovar_base samovar data after preprocessing and building stages
#' @param N number of artificial samples to generate
#' @param init_sp species vector for initializing data generation, or FALSE for usage most common taxa, auto for choosing random taxa
#' @param init_ab species amount vector (values from 0 to 1) for initializing data generation, or FALSE for mean initial taxa assignment, or auto for usage from known edf for each species from init_sp
#' @param avoid_zero_generations logical, avoid zero-based generations or not. FALSE might results in under-distributed communities, while TRUE in over-represented with species from different clusters possibly come from different samples groups
#' @param seed initial seed for the seeds generation
#' @example R/examples/processing.R
#' @export

samovar_boil <- function(samovar_base, N = 1,
                         init_sp = F, init_ab = F,
                         avoid_zero_generations = T,
                         seed = 42) {
  #__ ----
  # Additional functions ----
  prob_calc <- function(...) prob_calc_general(..., probability_calculation = samovar_base$preferences$probability_calculation)
  summarise_cluster <- function(...) summarise_cluster_general(..., cluster_connection = samovar_base$preferences$cluster_connection)
  glm_calc <- function(...) glm_calc_general(..., summarise_cluster = summarise_cluster)
  get_newval <- function(...) get_newval_general(..., summarise_cluster = summarise_cluster, glm_calc = glm_calc, mode = "glm")
  presence_calculation <- function(...) presence_calculation_general(..., probability_calculation = samovar_base$preferences$probability_calculation)
  #__ ----

  # Main ----
  ## Initializing all ----
  if(N < max(length(init_sp), length(init_ab))) N <- max(length(init_sp), length(init_ab))
  pb <- progress_function(N)

  ## If initials not specified, mean of most occurred species is used
  if (isFALSE(init_sp)) {
    init_sp <- samovar_base$samovar_data$species %>% sample(N, replace = T)
  } else if(init_sp == "auto"){
    init_sp <- (samovar_base$samovar_data$data > 0) %>%
      apply(1, sum) %>%
      which.max %>% names
  }

  if (isFALSE(init_ab)) {
    for (sp in seq_along(init_sp)) {
      tmp <- samovar_base$samovar_data$data %>%
        subset(rownames(.) %in% init_sp[sp]) %>%
        as.numeric
      init_ab[sp] <- tmp[tmp > samovar_base$samovar_data$min_value] %>%
        sample(1)
    }
  } else if (init_ab == "auto") {
    for (sp in seq_along(init_sp)) {
      init_ab[sp] <- samovar_base$samovar_data$data %>%
        subset(rownames(.) %in% init_sp[sp]) %>%
        as.numeric %>%
        mean(na.rm = T)
    }
    init_ab <- init_ab %>% as.numeric
  }

  ## vectorizing inits
  init_sp <- rep(init_sp, length.out = N)
  init_ab <- rep(init_ab %>% samovar_base$samovar_data$normalization_function(.), length.out = N)

  #generate seeds
  seed_list <- runif(n = N)*100000

  cat("---  Generation started  ---\n\n")
  results <- new("samovar_run")

  ## make iterations ----
  for(iter in 1:N) {
    set.seed(seed_list[iter])

    iter_error <- try(silent = T, {

      init <- init_sp[iter]
      init_level <- init_ab[iter]

      ## Initializing iter ----
      cluster <- samovar_base$get_cluster(init)
      cluster_len <- samovar_base$get_cluster_len(cluster)
      cluster_todo <- rep(0, samovar_base$samovar_data$cluster_n())
      cluster_todo[cluster] <- 1

      ## Calculate cluster representation ----
      xPr_il <- samovar_base$inter_cluster_graph_prob
      xRs_il <- samovar_base$inter_cluster_graph_method
      cluster_todo[-cluster] <- presence_calculation(xPr_il, cluster,
                                                     force_simple = T)

      # avoid zero-based generations
      cl0 <- which(cluster_todo == 0)

      if(avoid_zero_generations) {
        xRs_il[cl0,] <- 0
        xRs_il[,cl0] <- 0
      }

      cl_done <- cl0

      ## Generation loop ----
      while(length(cl_done) < samovar_base$samovar_data$cluster_n()) {

        ### Re-initializing cluster ----
        # for any cluster except of the first initial
        if(is.na(init)) {

          # get new cluster for generation
          from_to <- which.max.coord(xRs_il, cl_done[
            !(cl_done %in% cl0) #!!!!!
          ]) #samovar_base$get_max_inter(cl_done)
          cluster <- from_to[2]
          cluster_len <- samovar_base$get_cluster_len(cluster)

          if ((from_to[1] == from_to[2])|(is.na(cluster))) {
            warning(paste0("Exit: possible incomplete generation on cluster: ",
                           cl_done[length(cl_done)], " at iteration ", iter))
            break
          }

          # predict new value of the cluster
          Y <- samovar_base$samovar_data$get(from_to[2])
          X <- samovar_base$samovar_data$get(from_to[1])
          xval <- suppressWarnings(results$get_cluster(from_to[1]) %>% mean(na.rm = T))

          cluster_init_level <- get_newval(X = X, Y = Y, xval = xval,
                                           model = samovar_base$preferences$inter_model)
          #!!!
          tmp1 <- samovar_base$samovar_data$get_clean(from_to[1]) %>% colnames
          tmp <- samovar_base$samovar_data$get(from_to[2])

          tmp2 <- tmp[,tmp1]
          tmp2 <- tmp2[,apply(tmp2, 2, function(x) sum(x) > 0)]
          sp_cl <- rownames(tmp2)
          tmp2 <- tmp2[,sample(colnames(tmp2), 1)]

          if(sum(tmp2) != 0) {
            init <- sp_cl[which(tmp2 > 0)] %>% sample(.,1)
            means <- suppressWarnings(tmp %>% apply(1, mean, na.rm = T))
            init_level <- means[init] / mean(means) * cluster_init_level
          } else {
            init <- NA
            results$new_sp(cluster, sp_cl, 0, iter)
          }
        }

        ### Generating cluster ----
        #### Parse data ----
        if(!is.na(init)) {
          df_cl <- samovar_base$samovar_data$get(cluster)
          xRs_cl <- samovar_base$inner_cluster_graph_method[[as.character(cluster)]]
          xPr_cl <- samovar_base$inner_cluster_graph_prob[[as.character(cluster)]]

          species <- rownames(df_cl)
          species_abundance = rep(0, cluster_len)
          species_presence = rep(0, cluster_len)

          sp_done <- which(init == species)
          species_abundance[sp_done] <- init_level
          species_presence[sp_done] <- 1

          if(cluster_len > 1) {
            #### Generate occurrence data
            sp_done <- which(init == rownames(df_cl))
            species_presence[-sp_done] <- presence_calculation(xPr_cl, sp_done,
                                                               df_cl, init)

            #### Make generation graph ----
            prediction_graph <- tibble(from = sp_done, to = sp_done)

            # assign to 0 all removed species
            w0 <- which(species_presence == 0)
            sp_done <- c(sp_done, w0)

            if(avoid_zero_generations){
              xRs_cl[w0,] <- 0
              xRs_cl[,w0] <- 0
            }

            while(length(sp_done) < cluster_len) {
              pr <- which.max.coord(xRs_cl, prediction_graph$to)
              if(!is.na(pr[2])) prediction_graph <- prediction_graph %>% rbind(pr)
              sp_done <- c(sp_done, pr[2])
            }
            prediction_graph <- prediction_graph[-1,] %>% na.omit
            #circlize::chordDiagram(prediction_graph)

            if (nrow(prediction_graph) > 0)
              for (i in 1:nrow(prediction_graph)) {
                Y <- df_cl[prediction_graph$to[i],]
                X <- df_cl[prediction_graph$from[i],]
                xval <- species_abundance[prediction_graph$from[i]]
                yval <- get_newval(X = X, Y = Y, xval = xval,
                             model = samovar_base$preferences$inner_model)
                #cat(yval, "\n")
                species_abundance[prediction_graph$to[i]] <- yval
              }
          }
          #### Combine data ----
          results$new_sp(as.character(cluster), species, species_abundance, to_run = iter)
        }

        cl_done <- c(cl_done, cluster)
        init <- NA
      }
    })
    if(is.character(iter_error)) {
      warning(iter_error)
      warning(paste("Possible incompete generation on iter", iter))
      results$new_sp("error", "error", samovar_base$samovar_data$max_value, to_run = iter)
    }

    pb$tick()
  }
  set.seed(seed)


  #### Postprocessing ----
  #!!!!!
  results$data[results$data < samovar_base$samovar_data$min_value] <- samovar_base$samovar_data$min_value
  results$data[results$data > samovar_base$samovar_data$max_value] <- samovar_base$samovar_data$max_value
  results$data <- results$data %>%
    data.frame() %>%
    apply(c(1,2), samovar_base$samovar_data$reverse_normalization_function) %>%
    apply(2, function(line) {
      if((sum(line)) > 1) {line/sum(line)} else{line}
    }) %>%
    data.frame()

  results$data["unclassified",] <- apply(results$data, 2, function(x) 1 - sum(x) )

  colnames(results$data) <- paste0("gen", 1:ncol(results$data))
  ###
  cat("---  Generation done  ---\n\n")
  #results$data <- results$data %>%  samovar$samovar_data$reverse_normalize_df
  return(results$copy())
}
