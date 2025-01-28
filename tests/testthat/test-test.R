library(tidyverse)

test_that("installed", {
  library(samovaR)
})

test_that("GMrepo API", {
  GMrepo_type2data(number_to_process = 1)
})

test_that("processing stages", {
  teatree <- GMrepo_type2data(number_to_process = 100)
  tealeaves <- teatree %>%
    teatree_trim(treshhold_species = 3, treshhold_samples = 3, treshhold_amount = 10^(-3))
  teabag <- tealeaves %>%
    tealeaves_pack(normalization_function = function(x) log10(x+1))
  concotion <- teabag %>%
    teabag_brew(min_cluster_size = 4, max_cluster_size = 6)
})

test_that("full pipeline", {
  samovar <- GMrepo_type2data(number_to_process = 100) %>%
    teatree_trim(treshhold_species = 3, treshhold_samples = 3, treshhold_amount = 10^(-3)) %>%
    tealeaves_pack(normalization_function = function(x) log10(x+1)) %>%
    teabag_brew(min_cluster_size = 4, max_cluster_size = 6) %>%
    concotion_pour() %>%
    samovar_boil(N = 1)
})

test_that("viz", {
  GMrepo_type2data(number_to_process = 100) %>%
    viz_composition()
})
