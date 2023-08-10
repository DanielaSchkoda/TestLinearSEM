library(Rcpp)
source("./R/methods_polynomial_constraints.R")
sourceCpp("./src/methods_polynomial_constraints.cpp")
library(testthat)

test_that("U stat estimates correctly", {
  n = 100
  n_sim = 100
  
  X = cbind(rnorm(n, 0, 1), rnorm(n, 0, 1))
  
  ineq_satisfied = list(
    list(coef=-1, moms=matrix(c(1,1,0, 2,2,0), ncol=3, byrow=TRUE))
  )
  ineq_not_satisfied = list(
    list(coef=1, moms=matrix(c(1,1,0, 2,2,0), ncol=3, byrow=TRUE))
  )
  eq_satisfied = list(
    list(coef=1, moms=matrix(c(1,2,0, 1,1,0), ncol=3, byrow=TRUE))
  )
  eq_not_satisfied = list(
    list(coef=-1, moms=matrix(c(1,1,0, 2,2,0), ncol=3, byrow=TRUE))
  )
  
  # Test should accept for ineq_satisfied
  pvals = sapply(1:n_sim, function(sim) incomplete_U_stat(X, ineq_constraints = list(ineq_satisfied))$PVAL)
  ratio_rejected = length(which(pvals <= 0.05))/n_sim
  expect_true(ratio_rejected <= 0.1)
  
  # Test should accept for eq_satisfied
  pvals = sapply(1:n_sim, function(sim) incomplete_U_stat(X, equality_constraints = list(eq_satisfied))$PVAL)
  ratio_rejected = length(which(pvals <= 0.05))/n_sim
  expect_true(ratio_rejected <= 0.1)
  
  # Test should accept for ineq_satisfied and eq_satisfied together
  pvals = sapply(1:n_sim, function(sim) incomplete_U_stat(X, equality_constraints = list(eq_satisfied), ineq_constraints = list(ineq_satisfied))$PVAL)
  ratio_rejected = length(which(pvals <= 0.05))/n_sim
  expect_true(ratio_rejected <= 0.1)
  
  # Test should reject for ineq_not_satisfied 
  pvals = sapply(1:n_sim, function(sim) incomplete_U_stat(X, ineq_constraints = list(ineq_not_satisfied))$PVAL)
  ratio_rejected = length(which(pvals <= 0.05))/n_sim
  expect_true(ratio_rejected > 0.5)
  
  # Test should reject for eq_not_satisfied 
  pvals = sapply(1:n_sim, function(sim) incomplete_U_stat(X, equality_constraints = list(eq_not_satisfied))$PVAL)
  ratio_rejected = length(which(pvals <= 0.05))/n_sim
  expect_true(ratio_rejected > 0.5)
  
  # Test should reject for ineq_not_satisfied and eq_not_satisfied 
  pvals = sapply(1:n_sim, function(sim) incomplete_U_stat(X, equality_constraints = list(eq_not_satisfied), ineq_constraints = list(ineq_not_satisfied))$PVAL)
  ratio_rejected = length(which(pvals <= 0.05))/n_sim
  expect_true(ratio_rejected > 0.5)
})
