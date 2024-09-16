#' Filter species and samples from samovar_data object
#'
#' @param samovar_data Samovar data object to filter
#' @param metadata_filter False, character or data.frame with 2 columns: first contain metadata names for filtering, and second values per column
#' @param treshhold_amount Minimum value to conclude as not the noise.
#' @param treshhold_samples Minimum number of representing samples to keep species.
#' @param treshhold_species Minimum number of representing species to keep samples.
#' @param drop_unclassified Drop unknown and unclassified ranks. True by default
#' @example R/examples/preprocessing.R
#' @export

teatree_trim <- function(samovar_data,
                     metadata_filter = F,
                     treshhold_amount = 10^(-5),
                     treshhold_samples = 1,
                     treshhold_species = 1,
                     drop_species = F,
                     drop_unclassified = T) {

  # init
  data <- samovar_data$copy()

  # unclassified filtering
  if (drop_unclassified) {
    data$data <- data$data %>%
      subset(rownames(.) %>% str_detect("([Uu]nclassified)|([Uu]nknown)", negate = T))
  }

  # species filtering
  if(!isFALSE(drop_species)){
    data$data <- data$data[!(rownames(data$data) %in% drop_species),]
  }

  # metadata filtering
  if (!isFALSE(metadata_filter)) {
    if (is.character(metadata_filter)) {
      data$filter(metadata_filter[1], metadata_filter[2])
    } else {
      for (i in 1:nrow(metadata_filter)) {
        data$filter(metadata_filter[i,1], metadata_filter[i,2])
      }
    }
  }

  # treshhold filtering
  data$data[data$data <= treshhold_amount] <- 0
  data$rebuild(min_samp = treshhold_samples, min_sp = treshhold_species)

  return(data)
}
