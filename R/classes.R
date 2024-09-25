# GMrepo runs ----

#' GMrepo run data class
#'
#' @slot metadata metadata DataFrame
#' @slot run character
#' @method bind GMrepo_run
#' @method filter GMrepo_run
#' @name GMrepo_run
#' @rdname GMrepo_run
#' @docType class
#' @export

setRefClass("GMrepo_run",
  fields = list(
    metadata = "data.frame",
    run = "character",
    data = "data.frame"
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

# External runs data ----

#' Samovar run data class
#'
#' @slot metadata metadata DataFrame
#' @slot data data
#' @slot run character, samle IDs
#' @method bind samovar_run
#' @method add_species samovar_run
#' @method filter samovar_run
#' @name samovar_run
#' @rdname samovar_run
#' @docType class
#' @export

setRefClass("samovar_run",
            fields = list(
              metadata = "data.frame",
              data = "data.frame",
              run = "character",
              cluster = "character",
              species = "character"
            ),
            methods = list(
              show = function() {
                cat("Samovar run object, data:\n")
                print(data %>% rownames_to_column("sp") %>% as_tibble)
              },
              bind = function(samovar_data2) {
                metadata <<- rbind(metadata, samovar_data2$metadata)
                run <<- c(run, samovar_data2$run)

                new_species <- which( !(samovar_data2$species %in% species) )
                species <<- c(species, samovar_data2$species[new_species])
                cluster <<- c(cluster, samovar_data2$cluster[new_species])

                data[samovar_data2$species[new_species],] <<- 0
                data <<- cbind(data, samovar_data2$data)
              },
              new_sp = function(cluster_sp, sp, abundance, to_run = 1) {
                new_species <- sp[!(sp %in% species)]
                data[new_species,] <<- 0
                data[sp,to_run] <<- abundance
                cluster <<- c(cluster, rep(cluster_sp, length(sp)))
                species <<- c(species, new_species)
                data[is.na(data)] <<- 0
              },
              filter = function(field, value) {
                tmp <- which(metadata[, field] == value)
                metadata <<- metadata[tmp, ]
                data <<- data[tmp,]
                run <<- run[tmp]
              },
              get_cluster = function(cl) {
                data[cluster == cl,]
              },
              export = function() return(data)
            )
)

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
#' @method bind samovar_data
#' @method filter samovar_data
#' @method rescale samovar_data
#' @method rebuild samovar_data
#' @method normalize samovar_data
#' @method reverse_normalize samovar_data
#' @method get samovar_data
#' @method get_clean samovar_data
#' @method cluster_n samovar_data
#' @docType class
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
    bind = function(bind_with) {
      if (class(bind_with) == "GMrepo_run") {
        metadata <<- rbind(metadata, bind_with$metadata)
        run <<- c(run, bind_with$run)
        species <<- unique(species, bind_with$species)
        data <<- data %>% cbind(bind_with$data)
      } else if(class(bind_with) == "samovar_run") {
        metadata <<- rbind(metadata, bind_with$metadata)
        run <<- c(run, bind_with$run)
        data <<- data %>% cbind(bind_with$data)

        new_species <- which( !(bind_with$species %in% species) )
        species <<- c(species, bind_with$species[new_species])
        cluster <<- c(cluster, bind_with$cluster[new_species])
        cluster_size <<- table(cluster)
      }

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

      reverse_normalization_function <<- inverse(function(x) normalization_function(x),
        lower = 0, upper = 1)

      min_value <<- normalization_function(min_value)
      max_value <<- normalization_function(max_value)
    },
    reverse_normalize = function() {
      data <<- apply(data, 2, reverse_normalization_function)
    },
    reverse_normalize_df = function(x) {
      if(!is.null(reverse_normalization_function)) {
        x[x < min_value] <- min_value
        x[x > max_value] <- max_value

        x <- x %>%
          data.frame() %>%
          apply(c(1,2), reverse_normalization_function) %>%
          apply(2, function(line) {
            if((sum(line)) > 1) {line/sum(line)} else{line}
          }) %>%
          data.frame()

        x["unclassified",] <- apply(x, 2, function(x) 1 - sum(x) )
      }
      return(x)
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
#' @slot samovar_base samovar_data object
#' @slot method method to obtain samovar_base
#' @slot inner_cluster_graph_method list of graphs in matrix form of inner cluster connections
#' @slot inter_cluster_graph_method list of graphs in matrix form of inter cluster connections
#' @slot inner_cluster_graph_prob list of co-occurrence probabilities in matrix form of inner cluster members
#' @slot inter_cluster_graph_prob list of co-occurrence probabilities in matrix form between clusters
#' @slot properties concotion_pour() properties
#' @method get_cluster samovar_base
#' @method get_cluster_len samovar_base
#' @method get_max_inner samovar_base
#' @method get_max_inter samovar_base
#' @docType class
#' @name samovar_base
#' @rdname samovar_base
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
      print(as.data.frame(preferences) %>% t)
    },
    get_cluster = function(sp) {
      samovar_data$cluster[samovar_data$species == sp] %>%
        as.numeric()
    },
    get_cluster_len = function(cl) {
      samovar_data$cluster_size[cl]
    },
    get_max_inter = function(predicted,
                               type = "prob") {
      if(type != "prob") {
        df <- samovar$inter_cluster_graph_method
      } else {
        df <- inter_cluster_graph_prob
      }
      return(which.max.coord(df, predicted))
    },
    get_max_inner = function(cl, predicted,
                              type = "prob") {
      if(type != "prob") {
        df <- inter_cluster_graph_method[[cl]]
      } else {
        df <- inner_cluster_graph_prob[[cl]]
      }
      return(which.max.coord(df, predicted))
    }

  )
)
