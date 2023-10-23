# imports
library(testthat)
library(wrapr)
library(MASS)
library(Rcpp)
source("./R/CRT_statistic.R")
sourceCpp("./src/CRT_statistic.cpp")
source("./R/constraints_to_test.R")

test_that("Matrix is estimated correctly", {
    n <- 10000000
    X <- mvrnorm(n = n,
                 mu = c(0,0,0), 
                 Sigma = diag(x = c(1,2,0.5), nrow = 3, ncol = 3))
    # Pi with true rank 2
    Pi <- matrix(c("m1", "m2", "m3", "m1_1", "m1_2", "m1_3", "m2_2", "m2_3", "m3_3"), nrow=3)
    occuring_moments <- list_of_moments(1:2, 3)

    Pi_hat <- CRT(X, Pi, 2, moment_expressions = occuring_moments)$estimated_Pi
    Pi_expected <- matrix(c(0, 0, 0, 1, 0, 0, 2, 0, 0.5), nrow = 3)
    tolerance <- 1e-2
    expect_true(all(abs(Pi_hat - Pi_expected) < tolerance))
})

test_that("Matrix with coefficients is estimated correctly", {
  n <- 10000000
  X <- mvrnorm(n = n,
               mu = c(0,0,0), 
               Sigma = diag(x = c(1,2,0.5), nrow = 3, ncol = 3))
  # Pi with true rank 2
  Pi <- matrix(c("m1", "m2", "m3", "-m1_1", "m1_2", "m1_3", "m2_2", "m2_3", "m3_3"), nrow=3)

  occuring_moments <- list_of_moments(1:2, 3)
  Pi_hat <- CRT(X, Pi, 2, moment_expressions = occuring_moments)$estimated_Pi
  Pi_expected <- matrix(c(0, 0, 0, -1, 0, 0, 2, 0, 0.5), nrow = 3)
  tolerance <- 1e-2
  expect_true(all(abs(Pi_hat - Pi_expected) < tolerance))
})

test_that("Rank is estimated correctly by RS test", {
  n_simulations <- 100
  n <- 1000
  # Pi with true rank 2
  Pi <- matrix(c("m1", "m2", "m3", "m1_1", "m1_2", "m1_3", "m2_2", "m2_3", "m3_3"), nrow=3)
  occuring_moments <- list_of_moments(1:2, 3)
  # Test if test accepts most of the time when testing for rank = 2
  p_values <- rep(0, 100)
  for (sim in 1:100) {
    X <- mvrnorm(n = n,
                 mu = c(0,0,0), 
                 Sigma = diag(x = c(1,2,0.5), nrow = 3, ncol = 3))
    p_values[sim] <- CRT(X, Pi, 2, moment_expressions = occuring_moments)$PVAL
  }
  
  nominal_size <- 0.05
  rejection_rate <- length(which(p_values <= nominal_size))/100
  tolerance <- 0.05
  expect_true(abs(rejection_rate - nominal_size) < tolerance)
  
  # Test if test rejects most of the time when testing for rank = 1
  p_values <- rep(0, 100)
  for (sim in 1:100) {
    X <- mvrnorm(n = n,
                 mu = c(0,0,0), 
                 Sigma = diag(x = c(1,2,0.5), nrow = 3, ncol = 3))
    p_values[sim] <- CRT(X, Pi, 1, moment_expressions = occuring_moments)$PVAL
  }
  
  nominal_size <- 0.05
  rejection_rate <- length(which(p_values <= nominal_size))/100
  expect_true(rejection_rate > 0.5)
})