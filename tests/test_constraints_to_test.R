source("./R/constraints_to_test.R")
library(testthat)

test_that("A_T is constructed correctly", {
  p <- 4
  A_T <- construct_A_T(p)
  expect_equal(A_T[1,1], "1 * m112")
  expect_equal(A_T[2,1], "-1 * m111")
  expect_equal(A_T[2,3], "0")
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
  for (i in 1:100) {
    C <- array(runif(p^2, -2,2), c(p,p)) 
    t <- array(NA, c(p,p,p)) 
    for (ind in indices) {
      t[matrix(ind, 1)] <- sum(sapply(1:p, function(i) C[i,ind[1]]*C[i,ind[2]]*C[i, ind[3]]))
    }
    # Evaluate Aronhold invariant at t
    expect_true(evaluate_poly(t, Aronhold_invariant) < 1e-10)
  }
  

  # General sym. tensor: Aronhold should almost always be non-zero
  results <- sapply(1:100, function(i){
    t <- array(runif(p^3, -10, 10), c(p,p,p)) 
    for (ind in indices) {
      ind_sorted <- sort(ind)
      t[matrix(ind, 1)] <- t[matrix(ind_sorted, 1)]
    }
    evaluate_poly(t, Aronhold_invariant)
  })
  expect_true(sum(abs(results) > 1e-10)>=95)
})

test_that("Polynomials n=2 l=1 are correct", {
  p <- 2
  r <- 3
  k <- 4
  indices <- gtools::permutations(n=p, r=k, v=1:p, repeats.allowed=TRUE)
  indices <- apply(indices, 1, unlist, simplify = FALSE)
  
  # Tensor with lower rank: One of the polys should be <= zero up to small numerical error
  for (i in 1:100) {
    C <- array(runif(p*r, -2,2), c(r,p)) 
    t <- array(NA, c(p,p,p,p)) 
    for (ind in indices) {
      t[matrix(ind, 1)] <- sum(sapply(1:r, function(i) C[i,ind[1]]*C[i,ind[2]]*C[i, ind[3]]*C[i, ind[4]])) 

    }
    # Evaluate polys at t
    results <- c(evaluate_poly(t, disc), evaluate_poly(t, D), evaluate_poly(t, P))
    expect_true(min(results) < 1e-10)
  }
  
  
  # General sym. tensor: Sometimes all polys should be > 0
  results <- sapply(1:1000, function(i){
    t <- array(runif(p^4, -10, 10), c(p,p,p,p)) 
    for (ind in indices) {
      ind_sorted <- sort(ind)
      t[matrix(ind, 1)] <- t[matrix(ind_sorted, 1)]
    }
    min(c(evaluate_poly(t, disc), evaluate_poly(t, D), evaluate_poly(t, P)))
  })
  expect_true(sum(results > 1e-10)>=150)
})
