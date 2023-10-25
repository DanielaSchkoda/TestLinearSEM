library(Rcpp)
#source("./R/Incomplete_U_statistic.R")
#sourceCpp("./src/Incomplete_U_statistic.cpp")
library(testthat)

test_that("U stat estimates correctly", {
  n = 100
  n_sim = 100
  
  # Sample from two independent standard normal variables
  X = cbind(rnorm(n, 0, 1), rnorm(n, 0, 1))
  
  # Inequality -m11*m22 <= 0 is satisfied. (m11 = m22 = 1 for infinite sample size)
  ineq_satisfied = list(
    list(coef=-1, moms=list(c(1,1), c(2,2)))
  )
  ineq_not_satisfied = list(
    list(coef=1, moms=list(c(1,1), c(2,2)))
  )
  eq_satisfied = list(
    list(coef=1, moms=list(c(1,2), c(1,1)))
  )
  eq_not_satisfied = list(
    list(coef=-1, moms=list(c(1,1), c(2,2)))
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
