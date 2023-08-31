# The code is adapted from Nils Sturma, https://github.com/NilsSturma/TestGGM

###########
# helpers #
###########

permutations <- function(n){
  return(do.call(rbind,combinat::permn(seq(n))))
}

findn <- function(N,D){
  rem <- N%%D
  if (rem==0){
    return(N)
  } else {
    return(N-rem)
  }
}

random_combs <- function(n,k,nr){
  v <- seq(n)
  res <- matrix(0, nr, k)
  for (i in 1:nr){
    res[i,] <- sample(v,k)
  }
  # replacement is possible but very unlikely
  return(res)
}

degree <- function(poly) {
    degrees_monomials <- sapply(poly, function(mon) length(mon[["moms"]]))
    return(max(degrees_monomials))
}

###########
## indep ##
###########

#' @export
indep_stat <- function(X, equality_constraints=list(), ineq_constraints=list(), E=1000){
    all_constraints <- append(equality_constraints, ineq_constraints)
    nr_eqs = length(equality_constraints)
    nr_constraints <- length(all_constraints)
    order_kernel <- max(sapply(all_constraints, degree))
    N <- findn(nrow(X),order_kernel)
    indices <- matrix(1:N, ncol=order_kernel)

    H <- H_c(X, indices, all_constraints)
    n = dim(H)[1]

    # Mean and centering
    H_mean <- Rfast::colmeans(H)
    H_centered <- Rfast::transpose(Rfast::transpose(H) - H_mean) # Centering: H_i <- (H_i - H_mean)
    
    # Diagonal of the sample covariance of H
    H_cov <- Rfast::colsums(H_centered**2) / n
    
    # Vector for studentizing
    studentizer <- H_cov**(-1/2)
    
    # Test statistic
    marginal_stats <- studentizer * sqrt(n) * H_mean
    marginal_stats[1:nr_eqs] <- abs(marginal_stats[1:nr_eqs])
    test_stat <-  max(marginal_stats)
    
    # Bootstrapping 
    W <- matrix(0, E, nr_constraints)
    for (i in 1:E){
        epsilons <- rnorm(n,0,1)
        for (j in 1:nr_constraints){
            W[i,j] <- 1/sqrt(n) * sum(H_centered[,j] * epsilons)
        }
    }
    W[,1:nr_eqs] <- abs(W[,1:nr_eqs])
    W_studentized <- Rfast::transpose(Rfast::transpose(W) * studentizer)
    results <- Rfast::rowMaxs(W_studentized, value = TRUE)
    
    # pval
    pval <- (1 + sum(results >= test_stat)) / (1+E)
    
    return(list("PVAL"=pval, "TSTAT"=test_stat, "Results"=results, "H"=H))
}

############
## U-stat ##
############

#' @export
incomplete_U_stat <- function(X, equality_constraints=list(), ineq_constraints=list(), combine=c("and", "or"), E=1000, n1=min(nrow(X),500), N=2*nrow(X)){
  n <- nrow(X)
  all_constraints <- append(equality_constraints, ineq_constraints)
  nr_constraints <- length(all_constraints)
  nr_eqs <- length(equality_constraints)
  order_kernel <- max(sapply(all_constraints, degree))
  N <- min(0.7*choose(n,order_kernel), N)

  # Determine N_hat by Bernoulli sampling
  N_hat <- stats::rbinom(1, choose(n,order_kernel), (N / choose(n,order_kernel)))

  # Choose randomly N_hat unique subsets with cardinality r of {1,...,n}
  indices <- matrix(unlist(random_combs_c(n,order_kernel,N_hat)[[1]]), ncol = order_kernel, byrow = TRUE)

  # Compute matrix H
  H <- H_c(X, indices, all_constraints)
  H_mean <- colMeans(H)
  H_centered <- t(t(H) - H_mean)

  # Compute matrix G on subset of samples
  X_G <- X[sample(n, size=min(n1,n), replace=FALSE),]
  G <- G_c(X_G, order_kernel, all_constraints)
  G_mean <- colMeans(G)
  G_centered <- t(t(G) - G_mean)

  # Diagonal of the approximate variance of H
  cov_H_diag <- colSums(H_centered**2) / N_hat
  cov_G_diag <- colSums(G_centered**2) / n1
  cov_diag <- order_kernel**2 * cov_G_diag + (n/N) * cov_H_diag

  # Vector for studentizing
  studentizer <- cov_diag**(-1/2)

  # Test statistic
  marginal_stats <- sqrt(n) * H_mean
  # Take absolute values for equality constraints
  if (nr_eqs > 0) marginal_stats[1:nr_eqs] <- abs(marginal_stats[1:nr_eqs])
  # studentizing
  test_stat <-  max(studentizer * marginal_stats)

  # Bootstrap
  U_B <- matrix(0, E, nr_constraints)
  for (i in 1:E){
    epsilons <- rnorm(N_hat,0,1)
    for (j in 1:nr_constraints){
        U_B[i,j] <- 1/sqrt(N_hat) * sum(H_centered[,j] * epsilons)
    }
  }
  U_A <- matrix(0, E, nr_constraints)
  for (i in 1:E){
    epsilons <- rnorm(n1,0,1)
    for (j in 1:nr_constraints){
        U_A[i,j] <- 1/sqrt(n1) * sum(G_centered[,j] * epsilons)
    }
  }
  U <- order_kernel * U_A + sqrt(n/N) * U_B
  # Take absolute values for equality constraints
  if (nr_eqs > 0) U[,1:nr_eqs] <- abs(U[,1:nr_eqs]) 
  # studentizing
  U_studentized <- t(t(U) * studentizer)
  results <- Rfast::rowMaxs(U_studentized, value = TRUE)

  # pval
  pval <- (1 + sum(results >= test_stat)) / (1+E)

  return(list("PVAL"=pval, "TSTAT"=test_stat))
}
