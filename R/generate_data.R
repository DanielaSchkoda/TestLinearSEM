#' Generates data under the null hypothetis and different alternatives
#' 
#' @param n Integer, sample size.
#' @param p Integer, number of nodes.
#' @param noise_distr String, distribution of noise terms.
#' @param scenario String, scenario to sample from.
#' @param delta Real Number, delta for local alternatives.
#' 
#' @return Matrix with \code{n} rows and \code{p} columns. Each row corresponds 
#' to an independent sample of the distribution of X.
#' 
#' @export
generate_data = function(n, p, noise_distr, scenario = c("H_0", "alt_cos", "alt_quadr", "alt_confounding"), delta = NULL) {
  scenario <- match.arg(scenario)
  if ((scenario %in% c("alt_cos", "alt_quadr", "alt_confounding")) & is.null(delta)) {
    stop("for alt_cos, alt_quadr, and alt_confounding delta has to be specified")
  }
  
  eps <- calculate_eps(n, p, noise_distr)
  B <- matrix(runif(p^2, -1, 1), nrow = p)
  diag(B) <- 1

  data <- if (scenario == "H_0") {
      eps %*% solve(B)
    } else if (scenario == "alt_cos") {
      apply(eps %*% solve(B), c(1,2), function(x) (1-delta)*x+delta*cos(x))
    } else if (scenario == "alt_quadr") {
      data <- matrix(, n, 2)
      data[,1] <- eps[,1]
      data[,2] <- sapply(eps[,1], function(x) delta * x^2 + x) + eps[,2]
      data
    } else if (scenario == "alt_confounding") {
      Z = calculate_eps(n, 1, noise_distr)
      Z_copied = do.call(cbind, replicate(p, Z, simplify=FALSE))
      eps %*% solve(B) + delta * Z_copied
    }
  return(data)
}

calculate_eps = function(n, p, distr_eps = c("gamma", "beta", "overlapping_gaussian", "gaussian")) {
  distr_eps <- match.arg(distr_eps)
  
  if (distr_eps == "gamma") {
    shapes <- runif(p, 2, 3)
    rates <- runif(p, 1, 5)
    eps <- mapply(function(n, shape, rate) rgamma(n, shape, rate) - shape/rate, n, shapes, rates)
  } else if (distr_eps == "beta") {
    alphas <- runif(p, 1.5, 2)
    betas <- runif(p, 2, 10)
    eps <- mapply(function(n, beta, alpha) rbeta(n, alpha, beta) - alpha/(alpha+beta), n, alphas,  betas)
  } else if (distr_eps == "overlapping_gaussian") {
    sds_1 <- runif(p, 0.5, 1)
    sds_2 <- runif(p, 1, 3)
    mean_diffs <- runif(p, 3, 4)
    eps <- mapply(sample_overlapping_gaussian, n, sds_1, sds_2, mean_diffs)
  } else if (distr_eps == "gaussian") {
    sds <- runif(p, 0.5, 2)
    eps <- mapply(function(n, sd) rnorm(n, 0, sd), n, sds)
  }
  return(eps)
}

sample_overlapping_gaussian = function(n, sd_1, sd_2, mean_diff) {
  # Sample such that order is random
  return(sample(c(rnorm(n/2, 0, sd_1), rnorm(n/2, mean_diff, sd_2)) - 1/2 * mean_diff))
}