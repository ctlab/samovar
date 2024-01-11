GMrepo_run2data <- function(runs, #run names from GMrepo
                            number_to_process, #amount of runs to be processed
                            test = F #test or not
) {
  
  if(length(runs) > number_to_process) runs <- runs[1:number_to_process]
  
  rlenprint <- runs %>% length %/% 10
  rlenprint <- rlenprint * 1:10
  
  #init resulted tables
  res_sp <- data.frame()
  res_ph <- data.frame()
  res_mt <- list()
  
  for (RUN in runs) {
    
    #get run info from GMrepo
    query <- POST("https://gmrepo.humangut.info/api/getFullTaxonomicProfileByRunID", 
                  body = list( "run_id" = RUN ), encode = "json");
    retrieved_contents <- content( query );
    df <- fromJSON( xml_text( retrieved_contents ));
    
    if(!is.null(df[["species"]])) {
      dfsp <- df[["species"]][,-2]
      colnames(dfsp)[2] <- RUN
      
      #init res_sp
      if (is.null(res_sp[1,1])) {
        res_sp <- dfsp[1,-2]
        res_sp <- res_sp[-1,]
      }
      
      #find duplicated taxa
      df_dup <- (dfsp$ncbi_taxon_id %>% duplicated) %>% which
      
      if(!is.null(df_dup)) {
        for (i in df_dup) {
          df_nd <- which(dfsp$ncbi_taxon_id == dfsp[df_dup,1])
          df_str <- c(dfsp[df_dup, 1], 
                      sum(dfsp[df_nd, 2]),
                      dfsp[df_dup, 3])
          dfsp <- dfsp[-df_nd,]
          dfsp <- dfsp %>% rbind (df_str)
          dfsp[,1:2] <- sapply(dfsp[,1:2], as.numeric)
        }
      }
      
      
      #process data to table
      res_sp <- full_join(res_sp, dfsp, by = c("ncbi_taxon_id", "scientific_name"))
      res_ph <- rbind(res_ph, c(RUN, df[["phenotypes"]]$term))
      res_mt[[RUN]] <- df[["run"]]
      
    }
    
    if (which(RUN == runs) %in% rlenprint){
      cat("%", which(which(RUN == runs) == rlenprint) * 10, "%")
    }
  }
  
  return(list(res_sp, res_ph, res_mt))
}