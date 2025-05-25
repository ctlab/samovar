# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(samovaR)

# Create a new environment for test data
test_data_env <- new.env()

# Function to load test data
load_test_data <- function() {
  test_data_path <- system.file("testdata/data_test.RData", package = "samovaR")
  if (!file.exists(test_data_path) || test_data_path == "") {
    stop("Test data file not found at: ", test_data_path)
  }
  load(test_data_path, envir = test_data_env)
  attach(test_data_env, name = "test_data_env", pos = 2)
}

# Load test data when the package is loaded
load_test_data()

# Run tests
test_check("samovaR")
