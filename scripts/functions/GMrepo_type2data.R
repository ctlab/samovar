GMrepo_type2data <- function(
  mesh_ids = c("D006262"), #character of all types to use
  number_to_process = F, #numeric vector for number of runs per meshID
  experiment_type = F, #if logical F, use all exp.types; else - only chars
  test = F
){
  if (test) {
    #get current dir
    getSrcDirectory(function(x) {x})
    cur_dir <- here(getSrcDirectory(GMrepo_type2data))
    #read test file
    t1 <- read.table(paste0(cur_dir, "/../test/test_data_main.txt"))
    t2 <- read.table(paste0(cur_dir, "/../test/test_data_meta.txt"))
    t3 <- read_json(paste0(cur_dir, "/../test/test_data_meta.json"))
    return(list(t1,t2,t3))
  } else {
    GMrepo_type2run(mesh_ids, number_to_process, experiment_type) %>% GMrepo_run2data()
  }
}
