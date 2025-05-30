# GMrepo runs ----

#' GMrepo run data class
#'
#' @slot metadata metadata DataFrame
#' @slot run character
#' @slot data data.frame
#' @name GMrepo_run
#' @rdname GMrepo_run
#' @docType class
#' @export
setClass(
  "GMrepo_run",
  slots = list(
    metadata = "data.frame",
    run = "character",
    data = "data.frame"
  )
)

#' @export
setMethod("show", "GMrepo_run", function(object) {
  cat("GMrepo run object:")
  cat(object@run, sep = "\n")
})

#' @export
setMethod("bind", "GMrepo_run", function(object, GMrepo_run2) {
  object@metadata <- rbind(object@metadata, GMrepo_run2@metadata)
  object@run <- c(object@run, GMrepo_run2@run)
  return(object)
})

#' @export
setMethod("filter", "GMrepo_run", function(object, field, value) {
  tmp <- which(object@metadata[, field] == value)
  object@metadata <- object@metadata[tmp, ]
  object@run <- object@run[tmp]
  return(object)
})

# External runs data ----

#' Samovar run data class
#'
#' @slot metadata metadata DataFrame
#' @slot data data
#' @slot run character, sample IDs
#' @slot cluster character
#' @slot species character
#' @name samovar_run
#' @rdname samovar_run
#' @docType class
#' @export
setClass(
  "samovar_run",
  slots = list(
    metadata = "data.frame",
    data = "data.frame",
    run = "character",
    cluster = "character",
    species = "character"
  )
)

#' @export
setMethod("show", "samovar_run", function(object) {
  cat("Samovar run object, data:\n")
  print(object@data %>% rownames_to_column("sp") %>% as_tibble)
})

#' @export
setMethod("bind", "samovar_run", function(object, samovar_data2) {
  object@metadata <- rbind(object@metadata, samovar_data2@metadata)
  object@run <- c(object@run, samovar_data2@run)

  new_species <- which(!(samovar_data2@species %in% object@species))
  object@species <- c(object@species, samovar_data2@species[new_species])
  object@cluster <- c(object@cluster, samovar_data2@cluster[new_species])

  object@data[samovar_data2@species[new_species],] <- 0
  object@data <- cbind(object@data, samovar_data2@data)
  return(object)
})

#' @export
setMethod("add_species", "samovar_run", function(object, cluster_sp, sp, abundance, to_run = 1) {
  new_species <- sp[!(sp %in% object@species)]
  object@data[new_species,] <- 0
  object@data[sp,to_run] <- abundance
  object@cluster <- c(object@cluster, rep(cluster_sp, length(sp)))
  object@species <- c(object@species, new_species)
  object@data[is.na(object@data)] <- 0
  return(object)
})

#' @export
setMethod("filter", "samovar_run", function(object, field, value) {
  tmp <- which(object@metadata[, field] == value)
  object@metadata <- object@metadata[tmp, ]
  object@data <- object@data[tmp,]
  object@run <- object@run[tmp]
  return(object)
})

#' @export
setMethod("get_cluster", "samovar_run", function(object, cl) {
  return(object@data[object@cluster == cl,])
})

#' @export
setMethod("export_data", "samovar_run", function(object) {
  return(object@data)
})

# Samovar data ----

#' samovar data class
#'
#' @slot description metadata DataFrame
#' @slot data DataFrame with species abundances. No NA pass
#' @slot run character, runs
#' @slot species character, runs
#' @slot normalization_function normalization function for samples
#' @slot reverse_normalization_function reverse normalization function
#' @slot min_value minimal value after scaling
#' @slot max_value maximal value after scaling
#' @slot cluster character vector, enumerated clusters for each species
#' @slot cluster_size named numeric, cluster sizes per cluster
#' @docType class
#' @name samovar_data
#' @rdname samovar_data
#' @export
setClass(
  "samovar_data",
  slots = list(
    metadata = "data.frame",
    data = "data.frame",
    run = "character",
    species = "character",
    normalization_function = "function",
    reverse_normalization_function = "function",
    min_value = 'numeric',
    max_value = 'numeric',
    cluster = 'character',
    cluster_size = 'numeric'
  )
)

#' @export
setMethod("show", "samovar_data", function(object) {
  cat("samovar_data object:\n")
  cat(length(object@run), "samples\n")
  cat(length(object@species), "species\n")
  print(tibble(object@data))
})

#' @export
setMethod("bind", "samovar_data", function(object, bind_with) {
  if (class(bind_with) == "GMrepo_run") {
    object@metadata <- rbind(object@metadata, bind_with@metadata)
    object@run <- c(object@run, bind_with@run)
    object@species <- unique(object@species, bind_with@species)
    object@data <- object@data %>% cbind(bind_with@data)
  } else if(class(bind_with) == "samovar_run") {
    object@metadata <- rbind(object@metadata, bind_with@metadata)
    object@run <- c(object@run, bind_with@run)
    object@data <- object@data %>% cbind(bind_with@data)

    new_species <- which(!(bind_with@species %in% object@species))
    object@species <- c(object@species, bind_with@species[new_species])
    object@cluster <- c(object@cluster, bind_with@cluster[new_species])
    object@cluster_size <- table(object@cluster)
  }
  return(object)
})

#' @export
setMethod("filter", "samovar_data", function(object, field, value) {
  tmp <- which(object@metadata[, field] %in% value)
  object@metadata <- object@metadata[tmp, ]
  object@data <- object@data[tmp, ] %>%
    subset(apply(1, sum, na.rm = T) > 0)
  object@run <- object@run[tmp]

  if (is.character(object@cluster)) {
    object@cluster <- object@cluster[tmp]
    warning("Clusters after filtration is not reassigned; may cause breaks in future steps")
  }
  return(object)
})

#' @export
setMethod("rescale", "samovar_data", function(object, scale_function = function(x) {
  x / sum(x)
}) {
  object@data[is.na(object@data)] <- 0
  object@data <- apply(object@data, 2, scale_function) %>% as.data.frame()
  object@min_value <- 0
  object@max_value <- 1
  return(object)
})

#' @export
setMethod("normalize", "samovar_data", function(object, normalization_function = object@normalization_function) {
  object@data <- apply(object@data, 2, normalization_function) %>% as.data.frame()
  object@data[is.na(object@data)] <- 0

  # assign normalization and inverse functions
  object@normalization_function <- normalization_function

  inverse <- function(f, lower, upper) {
    function(y) {
      uniroot((function(x) f(x) - y),
              lower = lower,
              upper = upper)[1] %>%
        unlist() %>%
        as.numeric()
    }
  }

  object@reverse_normalization_function <- inverse(function(x) normalization_function(x),
                                               lower = 0, upper = 1)

  object@min_value <- normalization_function(object@min_value)
  object@max_value <- normalization_function(object@max_value)
  return(object)
})

#' @export
setMethod("reverse_normalize", "samovar_data", function(object) {
  object@data <- apply(object@data, 2, object@reverse_normalization_function)
  return(object)
})

#' @export
setMethod("reverse_normalize_df", "samovar_data", function(object, x) {
  if(!is.null(object@reverse_normalization_function)) {
    x[x < object@min_value] <- object@min_value
    x[x > object@max_value] <- object@max_value

    x <- x %>%
      data.frame() %>%
      apply(c(1,2), object@reverse_normalization_function) %>%
      apply(2, function(line) {
        if((sum(line)) > 1) {line/sum(line)} else{line}
      }) %>%
      data.frame()

    x["unclassified",] <- apply(x, 2, function(x) 1 - sum(x) )
  }
  return(x)
})

#' @export
setMethod("rebuild", "samovar_data", function(object, min_sp = 1, min_samp = 1) {
  cat("Initial:", nrow(object@data), "x", ncol(object@data), "\n")
  # filter
  above0 <- function(x) sum(x > object@min_value)
  object@data <- object@data[apply(object@data, 1, above0) >= min_samp, ]
  object@data <- object@data[, apply(object@data, 2, above0) >= min_sp]

  # rebuild
  tmp <- which(object@species %in% rownames(object@data))
  object@run <- colnames(object@data)
  object@species <- rownames(object@data)
  object@metadata <- object@metadata[rownames(object@metadata) %in% object@run, ]
  cat("After filtering:", nrow(object@data), "x", ncol(object@data), "\n")

  if (length(object@cluster)>0) {
    object@cluster <- object@cluster[tmp]
    warning("Clusters after filtration is dropped; may cause breaks in future steps")
  }
  return(object)
})

#' @export
setMethod("get_clean", "samovar_data", function(object, cl) {
  df <- object@data[object@cluster == cl,]
  #filter
  df <- df[,apply(df,2, sum) > object@min_value]
  return(df)
})

#' @export
setMethod("get_data", "samovar_data", function(object, cl) {
  df <- object@data[object@cluster == cl,]
  return(df)
})

#' @export
setMethod("cluster_count", "samovar_data", function(object) {
  return(length(unique(object@cluster)))
})

# Samovar base ----

#' samovar base class
#'
#' @slot samovar_data samovar_data object
#' @slot inner_cluster_graph_method list of graphs in matrix form of inner cluster connections
#' @slot inter_cluster_graph_method list of graphs in matrix form of inter cluster connections
#' @slot inner_cluster_graph_prob list of co-occurrence probabilities in matrix form of inner cluster members
#' @slot inter_cluster_graph_prob list of co-occurrence probabilities in matrix form between clusters
#' @slot preferences concotion_pour() properties
#' @docType class
#' @name samovar_base
#' @rdname samovar_base
#' @export
setClass(
  "samovar_base",
  slots = list(
    samovar_data = 'samovar_data',
    inner_cluster_graph_method = 'list',
    inter_cluster_graph_method = 'matrix',
    inner_cluster_graph_prob = 'matrix',
    inter_cluster_graph_prob = 'matrix',
    preferences = 'list'
  )
)

#' @export
setMethod("show", "samovar_base", function(object) {
  cat("Prepared samovar:")
  cat("\nClusters:", object@samovar_data@cluster)
  print(object@samovar_data)
  print(as.data.frame(object@preferences) %>% t)
})

#' @export
setMethod("get_cluster_info", "samovar_base", function(object, sp) {
  return(as.numeric(object@samovar_data@cluster[object@samovar_data@species == sp]))
})

#' @export
setMethod("get_cluster_length", "samovar_base", function(object, cl) {
  return(object@samovar_data@cluster_size[cl])
})

#' @export
setMethod("get_max_inter_cluster", "samovar_base", function(object, predicted, type = "prob") {
  if(type != "prob") {
    df <- object@inter_cluster_graph_method
  } else {
    df <- object@inter_cluster_graph_prob
  }
  return(which.max.coord(df, predicted))
})

#' @export
setMethod("get_max_inner_cluster", "samovar_base", function(object, cl, predicted, type = "prob") {
  if(type != "prob") {
    df <- object@inter_cluster_graph_method[[cl]]
  } else {
    df <- object@inner_cluster_graph_prob[[cl]]
  }
  return(which.max.coord(df, predicted))
})
