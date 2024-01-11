#library installation ----

liblist <- c("tidyverse",
             "corrplot",
             "viridis",
             "plotly",
             "ggdendro",
             "httr",
             "jsonlite",
             "xml2",
             "tsne",
             "cluster",
             "Matrix",
             "shiny",
             "here")

for (i in liblist) {
  if (!require(i, quietly = T, character.only = T))
    install.packages(i)
  library(i, character.only = T)
  print(paste0(i, " loaded"))
}

#get current dir
getSrcDirectory(function(x) {x})
unnamed <- function(){}
cur_dir <- here(getSrcDirectory(unnamed))

#source functions

functions <- c(str_c(cur_dir, "scripts/functions/", dir(str_c(cur_dir, "scripts/functions"))), 
               str_c(cur_dir, "scripts/plots/", dir(str_c(cur_dir, "scripts/plots"))))

for (i in functions) source(i)

rm(functions, liblist, i, cur_dir, unnamed)
