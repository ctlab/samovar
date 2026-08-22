test_that("Samova.R: GMrepo API", {
  skip_on_cran()
  skip_if_offline()
  result <- tryCatch(
    GMrepo_type2data(number_to_process = 10),
    error = function(e) NULL
  )
  skip_if(is.null(result), "GMrepo API unavailable or returned invalid response")
  expect_s3_class(result, "samovar_data")
})
