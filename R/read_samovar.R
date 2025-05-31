#' Build samovar object from the abundance matrix
#'
#' @param data abundance matrix
##' Row names: organisms/OTUs/ASVs
##' Column names: samples
#' @param metadata metadata data.frame in format: ADD, or FALSE
#' @param min_sp data_samovar$rebuild() option, minimal number of species to filter. 0 by default
#' @param min_samp data_samovar$rebuild() option, minimal number of species to filter. 0 by default
#' @example R/examples/load_samovar.R
#' @return samovar object
#' @export
# Need tests!

table2samovar <- function(data, metadata = F, min_sp = 0, min_samp = 0){
  if(isFALSE(metadata)) metadata <- data.frame()
  data <- as.data.frame(data)
  data[is.na(data)] <- 0
  data_samovar = new("samovar_data",
                     data = data,
                     metadata = metadata,
                     run = colnames(data),
                     species = rownames(data))
  data_samovar$rebuild(min_sp, min_samp)
  data_samovar$rescale()

  return(data_samovar)
}


#' Build samovar data object from file or environment
#'
#' @param data path to abundance matrix
##' Row names: organisms/OTUs/ASVs
##' Column names: samples
#' @param metadata Data.frame or path to metadata file
#' @param ... Parameters processed by read.csv()
#' @example R/examples/load_samovar.R
#' @export
#' @return samovar object
# Need tests!

read_samovar <- function(data, metadata = F, ...) {
  if (is.character(data)) data <- read.csv(data, header = T, row.names = 1, ...)
  if (is.character(metadata)) metadata <- read.csv(metadata, header = T, row.names = 1, ...)

  data_samovar <- table2samovar(data, metadata, min_samp = 0, min_sp = 0)

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
#' @return samovar object
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
#' @param data Processed abundance table.
##' Row names: sequence IDs,
##' Column names:
##' - annotators: (starting with `taxID_`);
##' - `true`: for true annotation
##' - `length`: length of sequence
##' - sample
##'
#' @example R/examples/check_samovar.R
#' @return list of samovar data objects
#' @export


annotation2samovar <- function(data) {
  # fix colnames
  colnames(data) <- colnames(data) %>%
    stringr::str_remove("_[0-9]*$") %>%
    stringr::str_replace("taxid", "taxID")

  selected_columns <- (str_detect(colnames(data) , "^taxID_|^N_") &
                         str_detect(colnames(data) , "confidence", negate = T)) %>%
    which

  # build samovar
  res <- list()
  for (colname in colnames(data)[selected_columns]) {
    data_tmp <- data[,c(colname, "sample")] %>%
      summarise(value = n(), .by = c(!!sym(colname), sample)) %>%
      pivot_wider(values_from = value, names_from = !!sym(colname), id_cols = sample, values_fill = 0) %>%
      column_to_rownames("sample")

    data_samovar <- data_tmp %>% t %>%
      table2samovar(min_samp = 0, min_sp = 0)

    res[[colname]] <- data_samovar
  }

 return(res)
}
