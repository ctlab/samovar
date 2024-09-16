#' Get runs from GMrepo by meshID
#'
#' @param mesh_ids Character. All types of meshID to use. List of relations between meshID and phenotype could be obtained using `GMrepo_meshID()`
#' @param number_to_process False by default, or maximum number of runs per meshID
#' @importFrom xml2 xml_text
#' @importFrom jsonlite fromJSON
#' @import httr
#' @example R/examples/GMrepo.R
#' @export

GMrepo_type2run <- function(mesh_ids = c("D006262"),
                            number_to_process = F) {

  # set config to ignore certificate
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  options(RCurlOptions = list(ssl_verifypeer = FALSE))
  options(rsconnect.check.certificate = FALSE)

  run <- new("GMrepo_run")
  #iterate throw mesh_ids
  for (iter in seq_along(mesh_ids)) {

    #get variables for iteration
    mesh_id = mesh_ids[iter]
    if (!isFALSE(number_to_process)) {
      if (length(iter) <= length(number_to_process)) {
        N = number_to_process[iter]
      }
    } else {
        N = 10^10
    }

    # get runIDs metadata table
    GMrepo_run_data <- httr::POST("https://gmrepo.humangut.info/api/getAssociatedRunsByPhenotypeMeshIDLimit",
      body = list("mesh_id" = mesh_id, "skip" = 0, "limit" = N),
      encode = "json"
    ) %>%
      httr::content() %>%
      xml2::xml_text() %>%
      jsonlite::fromJSON() %>%
      column_to_rownames("run_id") %>%
      new("GMrepo_run", metadata = ., run = rownames(.))

    run$bind(GMrepo_run_data)
  }

  return(run)
}
