library(tidyverse)
library(samovaR)

make_test_samovar_fixtures <- function() {
  teatree <<- example_samovar_data(n_species = 4, n_samples = 4, seed = 1)
  tealeaves <<- teatree %>%
    teatree_trim(
      treshhold_species = 1,
      treshhold_samples = 1,
      treshhold_amount = 10^(-3)
    )
  teabag <<- tealeaves %>%
    tealeaves_pack(normalization_function = function(x) log10(x + 1))
  concotion <<- teabag %>%
    teabag_brew(min_cluster_size = 2, max_cluster_size = 4)
  data.samovar <<- concotion %>%
    concotion_pour(probability_calculation = "simple")
}

make_test_samovar_fixtures()
