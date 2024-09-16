#' Get data from GMrepo_run object
#'
#' @param runs GMrepo_run object got by GMrepo_type2run or created by user with `new('GMrepo_run', metadata = data.frame(), run = run_list)`
#' @param number_to_out False by default, maximum number of obtained data
#' @param at_level "species" by default. level to obtain classification from GMrepo
# @param keep_metadata To be implemented. Keep metadata from query
#' @param QC_filter QCStatus by default. Perform auto QC filtering based on metadata column, or False for no checking.
#' @import tidyverse
#' @importFrom xml2 xml_text
#' @importFrom jsonlite fromJSON
#' @importFrom httr content
#' @example R/examples/GMrepo.R
#' @export

GMrepo_run2data <- function(run,
                            number_to_out = F,
                            at_level = "species",
                            QC_filter = "QCStatus") {

  # set config to ignore certificate
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  options(RCurlOptions = list(ssl_verifypeer = FALSE))
  options(rsconnect.check.certificate = FALSE)

  # initialize
  if (!isFALSE(QC_filter)) {
    run$filter(QC_filter, 1) # get only runs with annotated data
  }
  pb <- progress_function(length(run$run))

  # init resulted tables
  res_sp <- tibble(taxa = character())

  for (RUN in run$run) {
    pb$tick()

    # get run info from GMrepo
    query <- httr::POST("https://gmrepo.humangut.info/api/getFullTaxonomicProfileByRunID",
      body = list("run_id" = RUN), encode = "json"
    ) %>%
      httr::content() %>%
      xml2::xml_text() %>%
      jsonlite::fromJSON()

    if (!is.null(query[[at_level]])) {
      df_sp <- query[[at_level]][, -2] %>%
        summarise(N = sum(relative_abundance), .by = scientific_name)
      colnames(df_sp) <- c("taxa", RUN)
      res_sp <- full_join(res_sp, df_sp, by = "taxa")
    }
  }

  res_sp <- res_sp %>% column_to_rownames("taxa")

  GMrepo <- new("samovar_data",
    run = run$run,
    metadata = run$metadata,
    data = res_sp,
    species = rownames(res_sp)
  )

  GMrepo$rescale()

  return(GMrepo)
}
