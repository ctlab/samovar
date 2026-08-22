#' Get runs from GMrepo by meshID
#'
#' @param mesh_ids Character. All types of meshID to use. List of relations between meshID and phenotype could be obtained using `GMrepo_meshID()`
#' @param number_to_process False by default, or maximum number of runs per meshID
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
    GMrepo_run_data <- gmrepo_post(
      "/api/getAssociatedRunsByPhenotypeMeshIDLimit",
      body = list("mesh_id" = mesh_id, "skip" = 0, "limit" = N)
    ) %>%
      gmrepo_parse() %>%
      column_to_rownames("run_id") %>%
      new("GMrepo_run", metadata = ., run = rownames(.))

    run <- bind_samovar(run, GMrepo_run_data)
  }

  return(run)
}
