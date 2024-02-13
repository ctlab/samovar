# get and split arguments ----
args <- commandArgs(TRUE)

if (sum(args == "--help") != 0) {
  cat("\n")
  wd <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[4])
  wd <- gsub(".R", ".read", wd)
  help_info <- readLines(wd)
  for (i in help_info) cat(i, "\n")
  cat("\n")
  stop()
}

wd <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[4])
wd <- gsub("SH_samova.R", "", wd)

#source functions
liblist <- c("tidyverse")
for (i in liblist) {
  if (!require(i, quietly = T, character.only = T))
    install.packages(i)
  library(i, character.only = T)
  cat("--- ", i, " loaded  ---\n")
}
mkp <- function(wd, path) {
  ps <- paste0(wd, path)
  paths <- paste0(ps, dir(ps))
  return(paths)
}
for (i in mkp(wd, "scripts/functions/")) source(i)


if (sum(args == "--test") != 0) {
  teatree <- mkp(wd, "scripts/test/test_data_main.txt") %>% read.table
  cat("\n---  Data loaded  ---\n")
} else {
  #load
  cat("now only test data is supported")
  cat("To be implemented")
  stop()
}

wn <- function(args, a) {
  pos1 <- args[which(args == a) + 1]
  if (length(pos1) == 0) return()
  
  pos2 <- args[(which(args == a) + 1):length(args)] 
  pos <- pos2 %>% str_detect("-")
  if (sum(pos) == 0){
    res <- paste(pos2, collapse = " ")
  } else {
    res <- pos2[1:(which(pos)[1]-1)] %>% paste(collapse = " ")
  }
  
  return(res)
}

#i
trA <- args[which(args == "-trA") + 1] %>% as.numeric
trP <- args[which(args == "-trA") + 1] %>% as.numeric
km <- args[which(args == "-k") + 1] %>% as.numeric
mc <- args[which(args == "-mc") + 1] %>% as.numeric
nf <- args[which(args == "-nf") + 1]
o <- args[which(args == "-o") + 1]
pref <- args[which(args == "-pref") + 1]
n <- args[which(args == "-n") + 1] %>% as.numeric
inn <- args[which(args == "-inn") + 1]
int <- args[which(args == "-int") + 1]
is <- wn(args, "-is")
il <- args[which(args == "-il") + 1] %>% as.numeric

if(length(trA) == 0) trA <- 10^(-5)
if(length(trP) == 0) trP <- 2
if(length(km) == 0) km <- 10
if(length(mc) == 0) mc <- 2
if(length(nf) == 0) nf <- function(x) log10(x+1) else as.call(nf)
if(length(o) == 0) o <- "generated_data.csv"
if(length(pref) == 0) pref <- F
if(length(n) == 0) n <- 1
if(length(inn) == 0) inn <- "gaussian"
if(length(int) == 0) int <- "gaussian"
if(length(il) == 0) il <- F
if(length(is) == 0) is <- F

add <- ifelse(sum(args == "--add") != 0, T, F)
du <- ifelse(sum(args == "--not_drop_unclassified") != 0, F, T)
cu <- ifelse(sum(args == "--calc_unclassified") != 0, T,F)
  
cat("---  Arguments parsed  ---\n")

tea <- teatree %>% 
  res_normalize(treshhold_amount = trA,
                treshhold_presence = trP,
                normalisation_function = nf, 
                drop_unclassified = du) %>%
  build_samovar(k_means = km,
                inner_model = inn,
                inter_model = int,
                minimal_cluster = mc) %>%
  boil(init = il,
       init_level = is, 
       n = n, 
       pref = pref, 
       calculate_unclassified = cu)

write.csv(tea, o)