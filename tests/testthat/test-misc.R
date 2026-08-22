test_that("Samova.R: vizualization", {
  gg <- concotion %>%
    viz_composition()
  expect_s3_class(gg, "ggplot")
})
