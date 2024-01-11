GMrepo_type2run <- function(mesh_ids = c("D006262"), #character of all types to use
                            number_to_process = F, #numeric vector for number of runs per meshID
                            experiment_type = "Metagenomics", #if logical F, use all exp.types; else - only chars
                            test = F
                            ) {
  
  if(test){
    #get current dir
    getSrcDirectory(function(x) {x})
    cur_dir <- here(getSrcDirectory(GMrepo_type2run))
    #read test file
    read.table(paste0(cur_dir, "/../test/health_runIDs.txt")) %>% unlist
  } else {
  
  runIDs <- c()
    for (mesh_id in mesh_ids) {
      #number of associated runs
      pheno_07 <- POST("https://gmrepo.humangut.info/api/countAssociatedRunsByPhenotypeMeshID", body = list( "mesh_id" = "D006262"), encode = "json");
      pheno_07_cont <- content( pheno_07 );
      phenotyp_nr_assoc_runs <- fromJSON( xml_text( pheno_07_cont ));
      
      #get runIDs
      
      if((number_to_process[which(mesh_id==mesh_ids)] == F)|
         (number_to_process[which(mesh_id==mesh_ids)] > phenotyp_nr_assoc_runs)){
        number_to_process[which(mesh_id==mesh_ids)] <- phenotyp_nr_assoc_runs
      }
      
      params <- list( "mesh_id" = mesh_id, "skip" = 0, "limit" = number_to_process)
      
      pheno_08 <- POST("https://gmrepo.humangut.info/api/getAssociatedRunsByPhenotypeMeshIDLimit", 
                       body = params, encode = "json")
      
      pheno_08_cont <- content( pheno_08 )
      phenotyp_a_page_of_assoc_runs <- fromJSON( xml_text( pheno_08_cont ))
      
      #remove irrelevant
      if (experiment_type != F) {
        phenotyp_a_page_of_assoc_runs <- 
          phenotyp_a_page_of_assoc_runs[phenotyp_a_page_of_assoc_runs$experiment_type == experiment_type,]
      }
      #QC
      phenotyp_a_page_of_assoc_runs <- 
        phenotyp_a_page_of_assoc_runs[!is.na(phenotyp_a_page_of_assoc_runs$QCStatus),]
      phenotyp_a_page_of_assoc_runs <- 
        phenotyp_a_page_of_assoc_runs[phenotyp_a_page_of_assoc_runs$QCStatus != 1,]
      
      runIDs <- c(runIDs, phenotyp_a_page_of_assoc_runs$run_id)
    }
    return(runIDs %>% unlists)
  }
}



