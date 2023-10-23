# Random signed permutation matrix
rand_signed_perm <- function(p) {
  signs <- sample(c(-1,1), p, replace=TRUE)
  mat <- matrix(0,p,p)
  mat[cbind(seq_len(p),sample(p))] <- 1
  return(mat %*% diag(signs))
}

ICA <- function(X) {
  # Use dcovICA for p=2 since the function is faster than steadyICA
  if (ncol(X) == 2) {
    # X %*% whitener approx. Z
    steadyICA::whitener(X) %.>% to(
      Z <- Z,
      O <- whitener
    )
    # S_hat %*% t(W_hat) approx. Z
    steadyICA::dcovICA(Z, theta.0 = 0) %.>% to(
      W_hat <- W,
      S_hat <- S
    )
    # X approx. epsilon_hat %*% B_hat 
    return(list("B_hat"=t(W_hat) %*% solve(O), "epsilon_hat"=S_hat))
  } else {
    steadyICA::steadyICA(X, whiten = TRUE, method = 'Cpp') %.>% to(
      B_hat <- M,
      epsilon_hat <- S
    )
    return(list("B_hat"=B_hat, "epsilon_hat"=epsilon_hat))
  }
}

independence_stat <- function(epsilon) {
  if (ncol(epsilon) == 2) {
    n * steadyICA::dcovustat(1/n*rank(epsilon[,1]), 1/n*rank(epsilon[,2]))
  } else {
    n * steadyICA::multidcov(1/n*apply(epsilon, 2, rank), symmetric = FALSE)
  }
}

#' dCovICA
#' 
#' @param X Matrix of size n x p containing the observed data.
#' 
#' @return PVAL, the p-value. 
#' 
#' @export
dCovICA <- function(X) {
  # Initial ICA
  ICA(X) %.>% to(
    B_hat <- B_hat,
    epsilon_hat <- epsilon_hat
  )
  # Bootstrapping
  R = 1000
  n = nrow(X)
  bootstrap_samples <- parallel::mclapply(1:R, function(sim) {
    epsilon_star <- apply(epsilon_hat, 2, function(col) col[sample(n)])
    X_star <- epsilon_star %*% B_hat
    ICA(X_star) %.>% to(
      B_star_hat <- B_hat,
      epsilon_star_hat <- epsilon_hat
    )
    P <- rand_signed_perm(p)
    independence_stat(epsilon_star_hat %*% P)
  })
  pval <- (1+length(which(c(bootstrap_samples) >= independence_stat(epsilon_hat))))/(1+R)
  return(list("PVAL"=pval))
}