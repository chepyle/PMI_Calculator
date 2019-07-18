library(testthat)
library(shinytest)
# Many thanks: https://github.com/rstudio/shinytest-ci-example, 
test_that("Application works", {
  expect_pass(testApp(".", compareImages = FALSE))
})