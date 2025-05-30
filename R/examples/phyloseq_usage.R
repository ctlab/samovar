library(tidyverse)

# load phyloseq
library(phyloseq)
data(GlobalPatterns)
GP <- GlobalPatterns
ps <- phyloseq(otu_table(GP))

# transformations
library(samovaR)

# transform to phyloseq
samovar_data <- phyloseq2samovar(ps)
# and backwards
phyloseq_data <- samovar2phyloseq(samovar_data)

# regenerate data
samovar <- samovar_data %>%
  teatree_trim(
    treshhold_amount = 10^(-3),
    treshhold_samples = .1,
    treshhold_species = .5) %>%
  tealeaves_pack() %>%
  teabag_brew() %>%
  concotion_pour()

samovar_new <- samovar %>%
  samovar_boil(N = 100)

# add new data to phyloseq
phyloseq_new <- samovar2phyloseq(samovar_new)
