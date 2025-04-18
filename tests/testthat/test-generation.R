library(tidyverse)

test_that("Samova.R: build", {
  expect_no_error({
    concotion %>%
      concotion_pour(probability_calculation = "simple")
  })
})

test_that("Samova.R: generation", {
  expect_no_error({
    data.samovar %>%
      samovar_boil(N = 1, avoid_zero_generations = T)
  })
})
