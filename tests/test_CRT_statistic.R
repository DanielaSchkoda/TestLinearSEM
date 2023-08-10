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
    Pi <- matrix(list(list(mom=c(1), coeff=1), list(mom=c(2), coeff=1), list(mom=c(3), coeff=1),
                      list(mom=c(1,1), coeff=1), list(mom=c(1,2), coeff=1), list(mom=c(1,3), coeff=1),
                      list(mom=c(2,2), coeff=1), list(mom=c(2,3), coeff=1), list(mom=c(3,3), coeff=1)), 
                 nrow=3)
    Pi <- as.list(as.data.frame(Pi))
    
    Pi_hat <- CRT(X, Pi, 2)$estimated_Pi
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
  Pi <- matrix(list(list(mom=c(1), coeff=1), list(mom=c(2), coeff=1), list(mom=c(3), coeff=1),
                    list(mom=c(1,1), coeff=-1), list(mom=c(1,2), coeff=1), list(mom=c(1,3), coeff=1),
                    list(mom=c(2,2), coeff=1), list(mom=c(2,3), coeff=1), list(mom=c(3,3), coeff=1)), 
               nrow=3)
  Pi <- as.list(as.data.frame(Pi))
  
  Pi_hat <- CRT(X, Pi, 2)$estimated_Pi
  Pi_expected <- matrix(c(0, 0, 0, -1, 0, 0, 2, 0, 0.5), nrow = 3)
  tolerance <- 1e-2
  expect_true(all(abs(Pi_hat - Pi_expected) < tolerance))
})

test_that("Rank is estimated correctly by RS test", {
  n_simulations <- 100
  n <- 1000
  # Pi with true rank 2
  Pi <- matrix(list(list(mom=c(1), coeff=1), list(mom=c(2), coeff=1), list(mom=c(3), coeff=1),
                    list(mom=c(1,1), coeff=-1), list(mom=c(1,2), coeff=1), list(mom=c(1,3), coeff=1),
                    list(mom=c(2,2), coeff=1), list(mom=c(2,3), coeff=1), list(mom=c(3,3), coeff=1)), 
               nrow=3)
  Pi <- as.list(as.data.frame(Pi))
  
  # Test if test accepts most of the time when testing for rank = 2
  p_values <- rep(0, 100)
  for (sim in 1:100) {
    X <- mvrnorm(n = n,
                 mu = c(0,0,0), 
                 Sigma = diag(x = c(1,2,0.5), nrow = 3, ncol = 3))
    p_values[sim] <- CRT(X, Pi, 2)$PVAL
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
    p_values[sim] <- CRT(X, Pi, 1)$PVAL
  }
  
  nominal_size <- 0.05
  rejection_rate <- length(which(p_values <= nominal_size))/100
  expect_true(rejection_rate > 0.5)
})