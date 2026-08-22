# S4 class definitions (from ctlab/samovar@5d37936), with prototypes,
# dplyr-safe generic names, and $ / $<- slot accessors for compatibility.

.identity <- function(x) x

.slot_get <- function(x, name) {
  if (identical(name, "copy")) {
    return(function() x)
  }
  if (name %in% methods::slotNames(x)) {
    return(methods::slot(x, name))
  }
  stop(sprintf("'%s' has no slot '%s'", class(x)[1], name), call. = FALSE)
}

.slot_set <- function(x, name, value) {
  if (!(name %in% methods::slotNames(x))) {
    stop(sprintf("'%s' has no slot '%s'", class(x)[1], name), call. = FALSE)
  }
  methods::slot(x, name) <- value
  x
}

# Generics ----

#' S4 generics for samovaR objects
#'
#' @param object A samovaR S4 object
#' @param bind_with Object to bind to \code{object}
#' @param field Metadata column name
#' @param value Filter value
#' @param cluster_sp Cluster label for new species
#' @param sp Species names
#' @param abundance Numeric abundances
#' @param to_run Generated-sample column index
#' @param cl Cluster identifier
#' @param scale_function Column rescaling function
#' @param normalization_function Abundance transform
#' @param x Matrix or data frame to reverse-normalize
#' @param min_sp Minimum species per sample
#' @param min_samp Minimum samples per species
#' @param predicted Indices already generated
#' @param type Graph type, usually \code{"prob"}
#' @name bind_samovar
#' @rdname samovar-s4
#' @export
setGeneric("bind_samovar", function(object, bind_with) standardGeneric("bind_samovar"))

#' @rdname samovar-s4
#' @export
setGeneric("filter_samovar", function(object, field, value) standardGeneric("filter_samovar"))

#' @rdname samovar-s4
#' @export
setGeneric("add_species", function(object, cluster_sp, sp, abundance, to_run = 1) {
  standardGeneric("add_species")
})

#' @rdname samovar-s4
#' @export
setGeneric("export_data", function(object) standardGeneric("export_data"))

#' @rdname samovar-s4
#' @export
setGeneric("get_cluster_data", function(object, cl) standardGeneric("get_cluster_data"))

#' @rdname samovar-s4
#' @export
setGeneric("rescale_samovar", function(object, scale_function = function(x) x / sum(x)) {
  standardGeneric("rescale_samovar")
})

#' @rdname samovar-s4
#' @export
setGeneric("normalize_samovar", function(object, normalization_function = object@normalization_function) {
  standardGeneric("normalize_samovar")
})

#' @rdname samovar-s4
#' @export
setGeneric("reverse_normalize_samovar", function(object) standardGeneric("reverse_normalize_samovar"))

#' @rdname samovar-s4
#' @export
setGeneric("reverse_normalize_df_samovar", function(object, x) standardGeneric("reverse_normalize_df_samovar"))

#' @rdname samovar-s4
#' @export
setGeneric("rebuild_samovar", function(object, min_sp = 1, min_samp = 1) {
  standardGeneric("rebuild_samovar")
})

#' @rdname samovar-s4
#' @export
setGeneric("get_clean_samovar", function(object, cl) standardGeneric("get_clean_samovar"))

#' @rdname samovar-s4
#' @export
setGeneric("get_data_samovar", function(object, cl) standardGeneric("get_data_samovar"))

#' @rdname samovar-s4
#' @export
setGeneric("cluster_count", function(object) standardGeneric("cluster_count"))

#' @rdname samovar-s4
#' @export
setGeneric("get_cluster_info", function(object, sp) standardGeneric("get_cluster_info"))

#' @rdname samovar-s4
#' @export
setGeneric("get_cluster_length", function(object, cl) standardGeneric("get_cluster_length"))

#' @rdname samovar-s4
#' @export
setGeneric("get_max_inter_cluster", function(object, predicted, type = "prob") {
  standardGeneric("get_max_inter_cluster")
})

#' @rdname samovar-s4
#' @export
setGeneric("get_max_inner_cluster", function(object, cl, predicted, type = "prob") {
  standardGeneric("get_max_inner_cluster")
})

# GMrepo runs ----

#' GMrepo run data class
#'
#' @slot metadata metadata DataFrame
#' @slot run character
#' @slot data data.frame
#' @name GMrepo_run
#' @aliases GMrepo_run-class
#' @rdname GMrepo_run
#' @docType class
#' @export
setClass(
  "GMrepo_run",
  slots = list(
    metadata = "data.frame",
    run = "character",
    data = "data.frame"
  ),
  prototype = list(
    metadata = data.frame(),
    run = character(),
    data = data.frame()
  )
)

#' @export
setMethod("$", "GMrepo_run", function(x, name) .slot_get(x, name))

#' @export
setMethod("$<-", "GMrepo_run", function(x, name, value) .slot_set(x, name, value))

#' @export
setMethod("show", "GMrepo_run", function(object) {
  cat("GMrepo run object:")
  cat(object@run, sep = "\n")
})

#' @export
setMethod("bind_samovar", "GMrepo_run", function(object, bind_with) {
  object@metadata <- rbind(object@metadata, bind_with@metadata)
  object@run <- c(object@run, bind_with@run)
  object
})

#' @export
setMethod("filter_samovar", "GMrepo_run", function(object, field, value) {
  tmp <- which(object@metadata[, field] == value)
  object@metadata <- object@metadata[tmp, , drop = FALSE]
  object@run <- object@run[tmp]
  object
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
#' @aliases samovar_run-class
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
  ),
  prototype = list(
    metadata = data.frame(),
    data = data.frame(),
    run = character(),
    cluster = character(),
    species = character()
  )
)

#' @export
setMethod("$", "samovar_run", function(x, name) .slot_get(x, name))

#' @export
setMethod("$<-", "samovar_run", function(x, name, value) .slot_set(x, name, value))

#' @export
setMethod("show", "samovar_run", function(object) {
  cat("Samovar run object, data:\n")
  print(object@data %>% rownames_to_column("sp") %>% as_tibble)
})

#' @export
setMethod("bind_samovar", "samovar_run", function(object, bind_with) {
  object@metadata <- rbind(object@metadata, bind_with@metadata)
  object@run <- c(object@run, bind_with@run)

  new_species <- which(!(bind_with@species %in% object@species))
  object@species <- c(object@species, bind_with@species[new_species])
  object@cluster <- c(object@cluster, bind_with@cluster[new_species])

  object@data[bind_with@species[new_species], ] <- 0
  object@data <- cbind(object@data, bind_with@data)
  object
})

#' @export
setMethod("add_species", "samovar_run", function(object, cluster_sp, sp, abundance, to_run = 1) {
  sp <- as.character(sp)
  abundance <- as.numeric(abundance)
  if (length(abundance) == 1 && length(sp) > 1) {
    abundance <- rep(abundance, length(sp))
  }

  if (nrow(object@data) == 0 && ncol(object@data) == 0) {
    object@data <- as.data.frame(
      matrix(0, nrow = length(sp), ncol = to_run, dimnames = list(sp, NULL)),
      stringsAsFactors = FALSE
    )
    object@species <- sp
    object@cluster <- rep(as.character(cluster_sp), length(sp))
    object@data[sp, to_run] <- abundance
    object@data[is.na(object@data)] <- 0
    return(object)
  }

  if (ncol(object@data) < to_run) {
    extra <- to_run - ncol(object@data)
    for (i in seq_len(extra)) {
      object@data[[ncol(object@data) + 1]] <- 0
    }
  }

  new_species <- sp[!(sp %in% object@species)]
  if (length(new_species) > 0) {
    object@data[new_species, ] <- 0
    object@cluster <- c(object@cluster, rep(as.character(cluster_sp), length(new_species)))
    object@species <- c(object@species, new_species)
  }
  object@data[sp, to_run] <- abundance
  object@data[is.na(object@data)] <- 0
  object
})

#' @export
setMethod("filter_samovar", "samovar_run", function(object, field, value) {
  tmp <- which(object@metadata[, field] == value)
  object@metadata <- object@metadata[tmp, , drop = FALSE]
  object@data <- object@data[tmp, , drop = FALSE]
  object@run <- object@run[tmp]
  object
})

#' @export
setMethod("get_cluster_data", "samovar_run", function(object, cl) {
  object@data[object@cluster == cl, , drop = FALSE]
})

#' @export
setMethod("export_data", "samovar_run", function(object) {
  object@data
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
#' @aliases samovar_data-class
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
    min_value = "numeric",
    max_value = "numeric",
    cluster = "character",
    cluster_size = "numeric"
  ),
  prototype = list(
    metadata = data.frame(),
    data = data.frame(),
    run = character(),
    species = character(),
    normalization_function = .identity,
    reverse_normalization_function = .identity,
    min_value = 0,
    max_value = 1,
    cluster = character(),
    cluster_size = numeric()
  )
)

#' @export
setMethod("$", "samovar_data", function(x, name) .slot_get(x, name))

#' @export
setMethod("$<-", "samovar_data", function(x, name, value) .slot_set(x, name, value))

#' @export
setMethod("show", "samovar_data", function(object) {
  cat("samovar_data object:\n")
  cat(length(object@run), "samples\n")
  cat(length(object@species), "species\n")
  print(tibble(object@data))
})

#' @export
setMethod("bind_samovar", "samovar_data", function(object, bind_with) {
  if (is(bind_with, "GMrepo_run")) {
    object@metadata <- rbind(object@metadata, bind_with@metadata)
    object@run <- c(object@run, bind_with@run)
    object@species <- unique(c(object@species, bind_with@species))
    object@data <- object@data %>% cbind(bind_with@data)
  } else if (is(bind_with, "samovar_run")) {
    object@metadata <- rbind(object@metadata, bind_with@metadata)
    object@run <- c(object@run, bind_with@run)
    object@data <- object@data %>% cbind(bind_with@data)

    new_species <- which(!(bind_with@species %in% object@species))
    object@species <- c(object@species, bind_with@species[new_species])
    object@cluster <- c(object@cluster, bind_with@cluster[new_species])
    object@cluster_size <- as.numeric(table(object@cluster))
  }
  object
})

#' @export
setMethod("filter_samovar", "samovar_data", function(object, field, value) {
  tmp <- which(object@metadata[, field] %in% value)
  object@metadata <- object@metadata[tmp, , drop = FALSE]
  object@data <- object@data[tmp, , drop = FALSE] %>%
    subset(apply(., 1, sum, na.rm = TRUE) > 0)
  object@run <- object@run[tmp]

  if (is.character(object@cluster) && length(object@cluster) > 0) {
    object@cluster <- object@cluster[tmp]
    warning("Clusters after filtration is not reassigned; may cause breaks in future steps")
  }
  object
})

#' @export
setMethod("rescale_samovar", "samovar_data", function(object, scale_function = function(x) {
  x / sum(x)
}) {
  object@data[is.na(object@data)] <- 0
  object@data <- apply(object@data, 2, scale_function) %>% as.data.frame()
  object@min_value <- 0
  object@max_value <- 1
  object
})

#' @export
setMethod("normalize_samovar", "samovar_data", function(object, normalization_function = object@normalization_function) {
  object@data <- apply(object@data, 2, normalization_function) %>% as.data.frame()
  object@data[is.na(object@data)] <- 0

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

  object@reverse_normalization_function <- inverse(
    function(x) normalization_function(x),
    lower = 0, upper = 1
  )

  object@min_value <- normalization_function(object@min_value)
  object@max_value <- normalization_function(object@max_value)
  object
})

#' @export
setMethod("reverse_normalize_samovar", "samovar_data", function(object) {
  object@data <- apply(object@data, 2, object@reverse_normalization_function)
  object
})

#' @export
setMethod("reverse_normalize_df_samovar", "samovar_data", function(object, x) {
  if (!is.null(object@reverse_normalization_function)) {
    x[x < object@min_value] <- object@min_value
    x[x > object@max_value] <- object@max_value

    x <- x %>%
      data.frame() %>%
      apply(c(1, 2), object@reverse_normalization_function) %>%
      apply(2, function(line) {
        if ((sum(line)) > 1) {
          line / sum(line)
        } else {
          line
        }
      }) %>%
      data.frame()

    x["unclassified", ] <- apply(x, 2, function(x) 1 - sum(x))
  }
  x
})

#' @export
setMethod("rebuild_samovar", "samovar_data", function(object, min_sp = 1, min_samp = 1) {
  cat("Initial:", nrow(object@data), "x", ncol(object@data), "\n")
  above0 <- function(x) sum(x > object@min_value)
  object@data <- object@data[apply(object@data, 1, above0) >= min_samp, , drop = FALSE]
  object@data <- object@data[, apply(object@data, 2, above0) >= min_sp, drop = FALSE]

  tmp <- which(object@species %in% rownames(object@data))
  object@run <- colnames(object@data)
  object@species <- rownames(object@data)
  if (nrow(object@metadata) > 0) {
    object@metadata <- object@metadata[rownames(object@metadata) %in% object@run, , drop = FALSE]
  }
  cat("After filtering:", nrow(object@data), "x", ncol(object@data), "\n")

  if (length(object@cluster) > 0) {
    object@cluster <- object@cluster[tmp]
    warning("Clusters after filtration is dropped; may cause breaks in future steps")
  }
  object
})

#' @export
setMethod("get_clean_samovar", "samovar_data", function(object, cl) {
  df <- object@data[object@cluster == cl, , drop = FALSE]
  df <- df[, apply(df, 2, sum) > object@min_value, drop = FALSE]
  df
})

#' @export
setMethod("get_data_samovar", "samovar_data", function(object, cl) {
  object@data[object@cluster == cl, , drop = FALSE]
})

#' @export
setMethod("cluster_count", "samovar_data", function(object) {
  length(unique(object@cluster))
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
#' @aliases samovar_base-class
#' @rdname samovar_base
#' @export
setClass(
  "samovar_base",
  slots = list(
    samovar_data = "samovar_data",
    inner_cluster_graph_method = "list",
    inter_cluster_graph_method = "matrix",
    inner_cluster_graph_prob = "list",
    inter_cluster_graph_prob = "matrix",
    preferences = "list"
  ),
  prototype = list(
    samovar_data = new("samovar_data"),
    inner_cluster_graph_method = list(),
    inter_cluster_graph_method = matrix(1),
    inner_cluster_graph_prob = list(),
    inter_cluster_graph_prob = matrix(1),
    preferences = list()
  )
)

#' @export
setMethod("$", "samovar_base", function(x, name) .slot_get(x, name))

#' @export
setMethod("$<-", "samovar_base", function(x, name, value) .slot_set(x, name, value))

#' @export
setMethod("show", "samovar_base", function(object) {
  cat("Prepared samovar:")
  cat("\nClusters:", object@samovar_data@cluster)
  print(object@samovar_data)
  print(as.data.frame(object@preferences) %>% t)
})

#' @export
setMethod("get_cluster_info", "samovar_base", function(object, sp) {
  as.numeric(object@samovar_data@cluster[object@samovar_data@species == sp])
})

#' @export
setMethod("get_cluster_length", "samovar_base", function(object, cl) {
  object@samovar_data@cluster_size[cl]
})

#' @export
setMethod("get_max_inter_cluster", "samovar_base", function(object, predicted, type = "prob") {
  if (type != "prob") {
    df <- object@inter_cluster_graph_method
  } else {
    df <- object@inter_cluster_graph_prob
  }
  which.max.coord(df, predicted)
})

#' @export
setMethod("get_max_inner_cluster", "samovar_base", function(object, cl, predicted, type = "prob") {
  if (type != "prob") {
    df <- object@inter_cluster_graph_method[[cl]]
  } else {
    df <- object@inner_cluster_graph_prob[[cl]]
  }
  which.max.coord(df, predicted)
})
