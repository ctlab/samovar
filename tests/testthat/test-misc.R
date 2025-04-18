library(tidyverse)

test_that("Samova.R: vizualization", {
  gg <- concotion %>%
    viz_composition()
  expect_equal(class(gg), c("gg", "ggplot"))
})
