test_that("Samova.R: build", {
  expect_no_error({
    concotion %>%
      concotion_pour(probability_calculation = "simple")
  })
})

test_that("Samova.R: generation", {
  expect_no_error({
    generated <- data.samovar %>%
      samovar_boil(N = 2, avoid_zero_generations = T)
    expect_equal(ncol(generated$data), 2)
    expect_false("error" %in% rownames(generated$data))
  })
})
