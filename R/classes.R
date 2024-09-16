# GMrepo runs ----

#' GMrepo run data class
#'
#' @slot description metadata DataFrame
#' @slot run character
#' @method bind GMrepo_run
#' @method filter GMrepo_run
#' @name GMrepo_run
#' @rdname GMrepo_run
#' @export

setRefClass("GMrepo_run",
  fields = list(
    metadata = "data.frame",
    run = "character"
  ),
  methods = list(
    show = function() {
      cat("GMrepo run object:")
      cat(run, sep = "\n")
    },
    bind = function(GMrepo_run2) {
      metadata <<- rbind(metadata, GMrepo_run2$metadata)
      run <<- c(run, GMrepo_run2$run)
    },
    filter = function(field, value) {
      tmp <- which(metadata[, field] == value)
      metadata <<- metadata[tmp, ]
      run <<- run[tmp]
    }
  )
)

# Samovar data ----

#' samovar data class
#'
#' @slot description metadata DataFrame
#' @slot data DataFrame with species abundances. No NA pass
#' @slot run character, runs
#' @slot species character, runs
#' @method bind samovar_data
#' @method filter samovar_data
#' @method rescale samovar_data
#' @method rebuild samovar_data
#' @method normalize samovar_data
#' @method reverse_normalize samovar_data
#' @method get samovar_data
#' @name samovar_data
#' @rdname samovar_data
#' @export


setRefClass("samovar_data",
  fields = list(
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
  ),
  methods = list(
    show = function() {
      cat("samovar_data object:\n")
      cat(length(run), "samples\n")
      cat(length(species), "species\n")
      print(tibble(data))
    },
    bind = function(GMrepo_run2) {
      metadata <<- rbind(metadata, GMrepo_run2$metadata)
      run <<- c(run, GMrepo_run2$run)
      species <<- unique(species, GMrepo_run2$species)
      data <<- data %>% cbind(GMrepo_run2$data)
    },
    filter = function(field, value) {
      tmp <- which(metadata[, field] %in% value)
      metadata <<- metadata[tmp, ]
      data <<- data[tmp, ] %>%
        subset(apply(1, sum, na.rm = T) > 0)
      run <<- run[tmp]

      if (is.character(cluster)) {
        cluster <<- cluster[tmp]
        warning("Clusters after filtration is not reassigned; may cause breaks in future steps")
      }
    },
    rescale = function(scale_function = function(x) {
                         x / sum(x)
                       }) {
      data[is.na(data)] <<- 0
      data <<- apply(data, 2, scale_function) %>% as.data.frame()
      min_value <<- 0
      max_value <<- 1
    },
    normalize = function(normalization_function = normalization_function) {
      data <<- apply(data, 2, normalization_function) %>% as.data.frame()
      data[is.na(data)] <<- 0

      # assign normalization and inverse functions
      normalization_function <<- normalization_function

      inverse <- function(f, lower, upper) {
        function(y) {
          uniroot((function(x) f(x) - y),
                  lower = lower,
                  upper = upper)[1] %>%
            unlist() %>%
            as.numeric()
        }
      }

      reverse_normalization_function <<- inverse(function(x) normalisation_function(x),
        lower = 0, upper = 1
      )

      min_value <<- normalization_function(min_value)
      max_value <<- normalization_function(max_value)
    },
    reverse_normalize = function() {
      data <<- apply(data, 2, reverse_normalization_function)
    },
    rebuild = function(min_sp = 1, min_samp = 1) {
      cat("Initial:", nrow(data), "x", ncol(data), "\n")
      # filter
      above0 <- function(x) sum(x > min_value)
      data <<- data[apply(data, 1, above0) >= min_samp, ]
      data <<- data[, apply(data, 2, above0) >= min_sp]

      # rebuild
      tmp <- which(species %in% rownames(data))
      run <<- colnames(data)
      species <<- rownames(data)
      metadata <<- metadata[rownames(metadata) %in% run, ]
      cat("After filtering:", nrow(data), "x", ncol(data), "\n")

      if (length(cluster)>0) {
        cluster <<- cluster[tmp]
        warning("Clusters after filtration is dropped; may cause breaks in future steps")
      }
    },
    get_clean = function(cl) {
      df <- data[cluster == cl,]

      #filter
      df <- df[,apply(df,2, sum) > min_value]
      return(df)
    },
    get = function(cl) {
      df <- data[cluster == cl,]
      return(df)
    },
    cluster_n = function() return(length(unique(cluster)))
  )
)


# Samovar base----

#' samovar base class
#'
#' @slot samovar_dataa samovar_data object
#' @slot method method to obtain samovar_data
#' @slot inner_cluster_graph_method list of graphs in matrix form of inner cluster connections
#' @slot inter_cluster_graph_method list of graphs in matrix form of inter cluster connections
#' @slot inner_cluster_graph_prob list of co-occurrence probabilities in matrix form of inner cluster members
#' @slot inter_cluster_graph_prob list of co-occurrence probabilities in matrix form between clusters
#' @method bind samovar_data
#' @method filter samovar_data
#' @method rescale samovar_data
#' @method rebuild samovar_data
#' @method normalize samovar_data
#' @method reverse_normalize samovar_data
#' @name samovar_data
#' @rdname samovar_data
#' @export

setRefClass("samovar_base",
  fields = list(
    samovar_data = 'samovar_data',
    inner_cluster_graph_method = 'list',
    inter_cluster_graph_method = 'matrix',
    inner_cluster_graph_prob = 'list',
    inter_cluster_graph_prob = 'matrix',
    preferences = 'list'
  ),
  methods = list(
    show = function() {
      cat("Prepared samovar:")
      cat("\nClusters:", samovar_data$cluster)
      print(samovar_data)
      print(str(preferences))
    }
  )
)
