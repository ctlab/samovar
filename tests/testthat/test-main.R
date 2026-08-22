test_that("Samova.R: installation", {
  expect_no_error(library(samovaR))
})

test_that("Samova.R object", {
  expect_equal(class(data.samovar)[1], "samovar_base")
  expect_equal(class(data.samovar$samovar_data)[1], "samovar_data")
  expect_true(isS4(data.samovar))
  expect_true(isS4(data.samovar$samovar_data))
  expect_equal(typeof(data.samovar$samovar_data$data), "list")
  expect_gte(data.samovar$samovar_data$max_value, 0)
  expect_lte(data.samovar$samovar_data$min_value, data.samovar$samovar_data$max_value)
  expect_equal(typeof(data.samovar$samovar_data$metadata), "list")
  expect_equal(typeof(data.samovar$samovar_data$normalization_function), "closure")
  expect_equal(typeof(data.samovar$samovar_data$reverse_normalization_function),
               "closure")
  expect_equal(typeof(data.samovar$samovar_data$species), "character")
  expect_true(is.function(get_clean_samovar))
})
