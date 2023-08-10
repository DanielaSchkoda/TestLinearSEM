source("./R/constraints_to_test.R")
library(testthat)

test_that("A_T is constructed correctly", {
  p <- 4
  A_T <- construct_A_T(p)
  expect_equal(A_T[[1]][[1]], list(mom = c(1,1,2), coeff=1))
  expect_equal(A_T[[1]][[2]], list(mom = c(1,1,1), coeff=-1))
  expect_equal(A_T[[3]][[2]], list(mom = numeric(0), coeff=1))
})