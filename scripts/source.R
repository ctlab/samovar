#library installation ----

liblist <- c("tidyverse",
             "corrplot",
             "viridis",
             "plotly",
             "ggdendro",
             "tsne",
             "jsonlite",
             "here")

for (i in liblist) {
  if (!require(i, quietly = T, character.only = T))
    install.packages(i)
  library(i, character.only = T)
  cat("--- ", i, " loaded  ---\n")
}

#get current dir
getSrcDirectory(function(x) {x})
unnamed <- function(){}
cur_dir <- here(getSrcDirectory(unnamed))

#source functions

functions <- c(str_c(cur_dir, "/functions/", dir(str_c(cur_dir, "/functions"))), 
               str_c(cur_dir, "/plots/", dir(str_c(cur_dir, "/plots"))))

for (i in functions) source(i)

rm(functions, liblist, i, cur_dir, unnamed)

cat("---------\n\n")
