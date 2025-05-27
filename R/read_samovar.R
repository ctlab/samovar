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

#' Read annotations produced with samovar pipeline from the directory
#'
#' @param data_dir Path to abundance table. Row names: sequence IDs, column names: annotators; true for true annotation
#' @param sample_name_position Position to split the file basename by "." to extract sample names. 0 by default
#' @param ... Parameters processed by read.csv()
#' @importFrom dplyr bind_rows mutate
#' @importFrom stringr str_split str_remove
#' @example R/examples/check_samovar.R
#' @export

read_annotation_dir <- function(data_dir,  sample_name_position = 0, ...) {
  results <- tibble()
  for (data_path in dir(data_dir, pattern = ".csv$", full.names = T)) {
    sample_name <- (basename(data_path) %>%
      stringr::str_split("\\."))[[1]][1:(sample_name_position+1)]

    tmp <- read.csv(data_path, ...) %>%
      mutate(sample = sample_name)

    # fix colnames
    colnames(tmp) <- colnames(tmp) %>%
      stringr::str_remove("_[0-9]*$") %>%
      stringr::str_replace("taxid", "taxID")

    results <- dplyr::bind_rows(
      results,
      tmp
    )
  }

  return(results)
}


#' Process annotation data.frame to SamovaR
#'
#' @param data Processed abundance table. Row names: sequence IDs, column names: annotators; true for true annotation
#' @example R/examples/check_samovar.R
#' @return list of samovar data objects
#' @export


annotation2samovar <- function(data) {
 return()
}
