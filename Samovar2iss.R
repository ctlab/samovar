dict <- c(
  "root" = "Phix",
  "unclassified" = "Scer"
)

catp <- function(...) {
  paste0(..., "\n", collapse = "") %>% 
    cat(file = "iss_list.bash", append = T)
}

cnames <- function(table_samovar, 
                   dict, output) {
  for (i in names(dict)) {
    rownames(table_samovar)[rownames(table_samovar) == i] <- dict[i]
  }
  return(table_samovar)
}


table2iss <- function(..., 
                      dict = dict,
                      output = "benchmarking/",
                      n_reads = 10000,
                      path_to_genomes = "/mnt/tank/scratch/dsmutin/samovar/test_genomes/",
                      ext = ".fna"
                        ) {
  df <- samovar(...) %>% 
    #cnames(dict) %>% 
    rownames_to_column("sp") %>% 
    mutate(sp = paste0(sp, ext)) %>% 
    subset(
      sp %in% dir(path_to_genomes, pattern = paste0(ext, "$"))
    ) %>% 
    "rownames<-"(c(NULL)) %>% 
    as.data.frame %>% 
    column_to_rownames("sp") %>% 
    apply(2, function(x) x / sum(x)) %>% 
    apply(2, function(x) round(x * n_reads)) %>% 
    as.data.frame() %>% 
    rownames_to_column("sp")
  
  #generate using iss
  catp('mkdir ', output)
  catp('cd ', output)
  
  catp('
    iss_generate () {
      iss generate --genomes $1 --model hiseq --output ${2}_cont --n_reads $3
    } 
  ')
  
  acs3 <- '
      cat *${i}*R1* > ${i}_full_R1.fq
      cat *${i}*R2* > ${i}_full_R2.fq
      
      rm *.fastq
      rm *abundance*
      rm *vcf
      
      echo "done " ${i}
  '
  
  for (i in 2:ncol(df)){
    for (j in 1:nrow(df)) {
      acs2 <- paste0 (
        "iss_generate ",
        path_to_genomes, "/", df$sp[j], " ",
        df$sp[j], ".", i-1, " ",
        df[j,i] 
        )
      
      catp(acs2)
    }
    catp(
         acs3 %>% 
           str_replace_all("\\$\\{i\\}", as.character(i-1)))
  }
  
  #temp
  catp("cd ../../")
  cat ("\nCommand written...\n")
}