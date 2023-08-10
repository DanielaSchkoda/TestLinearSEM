#' @export
CRT <- function(X, Pi, r, upper_bound_s=Inf) {
  n <- nrow(X)
  Pi_hat <- estimate_Pi(Pi, X)
  svd(Pi_hat, nu = nrow(Pi_hat), nv = ncol(Pi_hat)) %.>% to(
    Sigma <- d,
    U <- u,
    V <- v
  )
  CRT <- n* sum(Sigma[(r + 1):min(nrow(Pi_hat), ncol(Pi_hat))]^2)
  
  Pi_vectorized <- unlist(Pi, recursive = FALSE)
  W <-  calculateW(Pi_vectorized, X)

  U_2 <- U[,(r+1):nrow(Pi_hat)]
  V_2 <- V[,(r+1):ncol(Pi_hat)]
  Omega <- (t(V_2) %x% t(U_2)) %*% W %*% (V_2 %x% U_2)
  nr_summands <- min(upper_bound_s, (nrow(Pi_hat) - r) * (ncol(Pi_hat) - r))
  weights <- rev(sort(eigen(Omega)$values))[1:nr_summands]
  # In case there's a complex part in one of the eigenvalues due to numerical errors, disregard it:
  weights <- Re(weights) 
  
  # Sample weighted chi2 distribution to calculate p-value
  E <- 10000
  sample_limiting_distr <- rowSums(sapply(1:nr_summands, function(i) weights[i]*rchisq(E, df=1)))
  
  pval <- length(which(c(sample_limiting_distr) >= CRT))/E
  
  return(list("PVAL" = pval, "TSTAT" = CRT, "estimated_Pi" = Pi_hat, "estimated_Omega" = Omega))
}