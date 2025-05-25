#' Build samovar data object from file or environment
#'
#' @param data Data.frame or path to abundance file. Row names is using as species list, column names as sample list. Unique names required
#' @param metadata Data.frame or path to metadata file
#' @param ... Parameters processed by read.csv()
#' @export
# Need tests!

read_annotation_table <- function(data, metadata = F, ...) {
  if (is.character(data)) data <- read.csv(data, header = T, row.names = 1, ...)
  if (is.character(metadata)) metadata <- read.csv(metadata, header = T, row.names = 1, ...)
  data[is.na(data)] <- 0
  data_samovar = new("samovar_data",
                     data = data,
                     metadata = metadata,
                     run = colnames(data),
                     species = rownames(data))
  data_samovar$rebuild

  return(data_samovar)
}

#' Build samovar annotation object from directory
#'
#' @param data_dir Path to abundance table. Row names: sequence IDs, column names: annotators; true_annotation
#' @param metadata Data.frame or path to metadata file
#' @param ... Parameters processed by read.csv()
#' @export
# Need tests!

read_annotation_dir <- function(data_dir, metadata = F, sample_name = 0, ...) {
  for (data in dir(data_dir()))
  if (is.character(data)) data <- read.csv(data, header = T, row.names = 1, ...)
  if (!isFALSE(metadata)) metadata <- read.csv(metadata, header = T, row.names = 1, ...)
  data[is.na(data)] <- 0
  data_samovar = new("samovar_data",
                     data = data,
                     metadata = metadata,
                     run = colnames(data),
                     species = rownames(data))
  data_samovar$rebuild

  return(data_samovar)
}



