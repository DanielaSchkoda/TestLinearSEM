#' CRT statistic
#' 
#' @param X Matrix of size n x p containing the observed data.
#' @param Pi Matrix for that we test whether rank(Pi) <= 4r.
#' @param r Integer, rank bound.
#' @param J Matrix, Jacobian of Pi needed for calculating the asymtotic covariance of Pi using the delta method. Optional.
#' @param moment_expressions List of moment expressions appearing in Pi. Optional.
#' 
#' @return List containing the p-value, the test statistic, the estimated Pi and its estimated covariance matrix.
#' 
#' @export
CRT <- function(X, Pi, r, J=NULL, moment_expressions=NULL) {
  n <- nrow(X)
  p <- ncol(X) 

  if (is.null(moment_expressions)) {
    highest_order <- max(c(apply(Pi, c(1,2), function(entry) {
      # Extract all m1_2_3 and count number of indices 
      if (grepl("m", entry)) {
        occuring_moments <- stringr::str_extract_all(entry, "m[0-9]+(_[0-9]+)*")[[1]]
        max(sapply(occuring_moments, function(mom) length(strsplit(mom, "_")[[1]])))
      }
      else 0
    }
    )))
    moment_expressions <- list_of_moments(2:highest_order, p)
  }
  
  if (is.null(J)) {
    J <- calculus::jacobian(c(Pi), var = names(moment_expressions))
  }
  
  estimate_Pi_and_W(Pi, moment_expressions, J, X) %.>% to(
    Pi_hat <- Pi_hat,
    W <- W
  )
  
  svd(Pi_hat, nu = nrow(Pi_hat), nv = ncol(Pi_hat)) %.>% to(
    Sigma <- d,
    U <- u,
    V <- v
  )

  CRT <- n* sum(Sigma[(r + 1):min(nrow(Pi_hat), ncol(Pi_hat))]^2)
  
  U_2 <- U[,(r+1):nrow(Pi_hat)]
  V_2 <- V[,(r+1):ncol(Pi_hat)]
  Omega <- (t(V_2) %x% t(U_2)) %*% W %*% (V_2 %x% U_2)
  nr_summands <- (nrow(Pi_hat) - r) * (ncol(Pi_hat) - r)
  weights <- rev(sort(eigen(Omega)$values))[1:nr_summands]
  # In case there's a small complex part in one of the eigenvalues due to numerical errors, disregard it:
  weights <- Re(weights) 
  
  # Sample weighted chi2 distribution to calculate p-value
  E <- 10000
  sample_limiting_distr <- rowSums(sapply(1:nr_summands, function(i) weights[i]*rchisq(E, df=1)))
  
  pval <- length(which(c(sample_limiting_distr) >= CRT))/E
  
  return(list("PVAL" = pval, "TSTAT" = CRT, "estimated_Pi" = Pi_hat, "estimated_W" = W))
}

estimate_Pi <- function(Pi, moment_expressions, X) {
  estimated_moments <- c(estimate_Moments(moment_expressions, X))
  names(estimated_moments) <- names(moment_expressions)
  calculus::evaluate(Pi, estimated_moments)
}

estimate_Pi_and_W <- function(Pi, moment_expressions, J, X) {
  estimated_moments <- c(estimate_Moments(moment_expressions, X))
  names(estimated_moments) <- names(moment_expressions)
  Pi_hat <- calculus::evaluate(Pi, estimated_moments)
  
  Sigma_moments <- estimate_Covs(moment_expressions, X)
  Jacobian <- calculus::evaluate(J, estimated_moments) 
  W <- Jacobian %*% Sigma_moments %*% t(Jacobian)
  return(list(Pi_hat=Pi_hat, W=W))
}

list_of_moments <- function(orders, p) {
  moms <- list()
  for (k in setdiff(orders, c(0))) { 
    order_k_moms <- gtools::combinations(p, k, 1:p, repeats.allowed = TRUE)
    order_k_moms <- lapply(split(order_k_moms, seq(nrow(order_k_moms))), function(ind) unlist(ind, use.names = FALSE))
    names <- lapply(order_k_moms, function(ind) paste0("m", paste0(ind, collapse = "_")))
    order_k_moms <- lapply(order_k_moms, function(ind) {
      list(coeff=1, mom=ind)
    })
    names(order_k_moms) <- names
    moms <- append(moms, order_k_moms)
  }
  return(moms)
}