library(tidyverse)


test_that("Samova.R: preprocessing", {
  expect_no_error({
    tealeaves_test <- teatree %>%
      teatree_trim(
        treshhold_species = 3,
        treshhold_samples = 3,
        treshhold_amount = 10^(-3)
      )
  })

  expect_no_error({
    teabag_test <- tealeaves %>%
      tealeaves_pack(
        normalization_function = function(x)
          log10(x + 1)
      )
  })

  expect_no_error({
    concotion_test <- teabag %>%
      teabag_brew(min_cluster_size = 4,
                  max_cluster_size = 6)
  })

})
