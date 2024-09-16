#' Build samovar data object from file or environment
#'
#' @param data Data.frame or path to abundance file. Row names is using as species list, column names as sample list. Unique names required
#' @param metadata Data.frame or path to metadata file
#' @param ... Parameters processed by read.table()
#' @export
# Need tests!

read.abundance <- function(data, metadata) {
  if (is.character(data)) data <- read.table(data, header = T, row.names = 1, ...)
  if (is.character(metadata)) metadata <- read.table(metadata, header = T, row.names = 1, ...)
  data[is.na(data)] <- 0
  data_samovar = new("samovar_data",
                     data = data,
                     metadata = metadata,
                     run = colnames(data),
                     species = rownames(data))
  data_samovar$rebuild

  return(data_samovar)
}

