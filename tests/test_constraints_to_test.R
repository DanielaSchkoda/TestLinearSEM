#source("./R/constraints_to_test.R")
#source("./R/CRT_statistic.R")
library(testthat)
library(Matrix)

test_that("Y_T is constructed correctly", {
  p <- 4
  k <- 3
  Y_T <- construct_Y_T(p, k)
  expect_equal(Y_T[1,1], "1 * m1_1_2")
  expect_equal(Y_T[5,1], "-1 * m1_1_1")
  expect_equal(Y_T[9,2], "0")
})

evaluate_poly <- function(t, poly) {
  result = 0
  for (summand in poly) {
    summand_evaluated <- summand$coef * prod(sapply(summand$moms, function(ind) {
      t[matrix(ind, 1)]
    }))
    result = result + summand_evaluated
  }
  return(result)
}

test_that("Aronhold invariant is correct", {
  p <- 3
  indices <- gtools::permutations(n=p, r=3, v=1:p, repeats.allowed=TRUE)
  indices <- apply(indices, 1, unlist, simplify = FALSE)
  
  # Tensor with lower rank: Aronhold should always be zero up to small numerical error
  for (sim in 1:100) {
    C <- array(runif(p^2, -2,2), c(p,p)) 
    t <- array(NA, c(p,p,p)) 
    for (ind in indices) {
      t[matrix(ind, 1)] <- sum(sapply(1:p, function(i) prod(C[i, ind])))
    }
    # Evaluate Aronhold invariant at t
    expect_true(evaluate_poly(t, Aronhold_invariant) < 1e-10)
  }
  

  # General sym. tensor: Aronhold should almost always be non-zero
  results <- sapply(1:100, function(sim){
    t <- array(runif(p^3, -10, 10), c(p,p,p)) 
    for (ind in indices) {
      ind_sorted <- sort(ind)
      t[matrix(ind, 1)] <- t[matrix(ind_sorted, 1)]
    }
    evaluate_poly(t, Aronhold_invariant)
  })
  expect_true(sum(abs(results) > 1e-10)>=95)
})

test_that("Cumulants are caluclated correctly", {
  cum_expr <- c(1,2,2,4,5,6)
  expect_false(grepl("3", get_cumulant_formula(cum_expr), fixed = TRUE))
  expect_true(grepl("1 * m1_2_2_4_5_6", get_cumulant_formula(cum_expr), fixed = TRUE))
  
  class(get_cumulant_formula(cum_expr))
  
  cum_expr <- c(1,1,2)
  expect_equal(get_cumulant_formula(cum_expr), "m1_1_2")
})
