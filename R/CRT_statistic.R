#' @export
CRT <- function(X, Pi, r, upper_bound_s=Inf, use_delta_method=FALSE, J=NULL, moment_expressions=NULL) {
  n <- nrow(X)
  
  if (use_delta_method) {
    delta_method(Pi, X, J, moment_expressions) %.>% to(
      Pi_hat <- M_hat,
      W <- W
    )
  } else {
    Pi_hat <- estimate_Pi(Pi, X)
    Pi_vectorized <- unlist(Pi, recursive = FALSE)
    W <- estimate_W(Pi_vectorized, X)
  }
  
  svd(Pi_hat, nu = nrow(Pi_hat), nv = ncol(Pi_hat)) %.>% to(
    Sigma <- d,
    U <- u,
    V <- v
  )
  CRT <- n* sum(Sigma[(r + 1):min(nrow(Pi_hat), ncol(Pi_hat))]^2)
  
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
  
  return(list("PVAL" = pval, "TSTAT" = CRT, "estimated_Pi" = Pi_hat, "estimated_W" = W))
}

#' @export
delta_method <- function(M_in_terms_of_moments, X, J, moment_expressions) {
  # Estimate all moments up to order four and their covariances.
  p <- ncol(X)
  # First we specify the formulas for the cumulants in terms of the moments
  #M_in_terms_of_moments <- apply(M, c(1,2), function(entry) if (grepl("c", entry)) get_cumulant_formula(entry) else as.character(entry))
  
  # Estimate M 
  # occuring_orders <- unique(c(apply(M_in_terms_of_moments, c(1,2), function(entry)
  #   if (grepl("m", entry)) nchar(sub('.*m', '', entry))
  #   else 0
  # )))
  # moment_expressions <- list_of_moments(occuring_orders, p) 
  estimated_moments <- c(estimate_Pi(list(moment_expressions), X))
  names(estimated_moments) <- names(moment_expressions)
  M_hat <- calculus::evaluate(M_in_terms_of_moments, estimated_moments)
  
  # Estimate W
  Sigma_moments <- estimate_W(moment_expressions, X)
  #J <- jacobian(c(M_in_terms_of_moments), var = names(moment_expressions))
  Jacobian <- calculus::evaluate(J, estimated_moments)
  W <- Jacobian %*% Sigma_moments %*% t(Jacobian)
  return(list(M_hat=M_hat, W=W))
}

list_of_moments <- function(orders, p) {
  moms <- list()
  for (k in setdiff(orders, c(0))) { 
    order_k_moms <- gtools::combinations(p, k, 1:p, repeats.allowed = TRUE)
    order_k_moms <- lapply(split(order_k_moms, seq(nrow(order_k_moms))), function(ind) unlist(ind, use.names = FALSE))
    names <- lapply(order_k_moms, function(ind) paste(c("m", ind), collapse = '') )
    order_k_moms <- lapply(order_k_moms, function(ind) {
      list(coeff=1, mom=ind)
    })
    names(order_k_moms) <- names
    moms <- append(moms, order_k_moms)
  }
  return(moms)
}