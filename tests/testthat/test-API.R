library(tidyverse)

test_that("Samova.R: GMrepo API", {
  expect_no_error(GMrepo_type2data(number_to_process = 10))
})
