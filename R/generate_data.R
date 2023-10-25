#' Generates data under the null hypothesis and different alternatives
#' 
#' @param n Integer, sample size.
#' @param p Integer, number of nodes.
#' @param l Integer, number of confounders
#' @param noise_distr String, distribution of noise terms.
#' @param scenario String, scenario to sample from.
#' @param delta Real Number, delta for local alternatives.
#' 
#' @return Matrix with \code{n} rows and \code{p} columns. Each row corresponds 
#' to an independent sample of the distribution of X.
#' 
#' @export
generate_data = function(n, p, l, noise_distr, scenario = c("H_0", "alt_cos", "alt_quadr", "alt_conf"), delta = NULL) {
  scenario <- match.arg(scenario)
  if ((scenario %in% c("alt_cos", "alt_quadr", "alt_conf")) & is.null(delta)) {
    stop("For alt_cos, alt_quadr, and alt_conf, delta has to be specified")
  }
  
  data <- if (scenario == "H_0") {
      generate_data_H_0(n, p, l, noise_distr)
    } else if (scenario == "alt_cos") {
      apply(generate_data_H_0(n, p, l, noise_distr), c(1,2), function(x) (1-delta)*x+delta*cos(x))
    } else if (scenario == "alt_conf") {
      eps <- calculate_eps(n, p, noise_distr)
      B <- matrix(runif(p^2, -1, 1), nrow = p)
      diag(B) <- 1
      Gamma <- matrix(runif(p*(l+1), -1, 1), nrow = l+1)
      # Multiply last col, which corresponds to additional confounder, by delta to obtain local alternatives
      L = calculate_eps(n, l+1, noise_distr) 
      L[, l+1] <- delta * L[, l+1]
      (eps + L %*% Gamma) %*% solve(B)
    } 
  return(data)
}

generate_data_H_0 <- function(n, p, l, noise_distr) {
  eps <- calculate_eps(n, p, noise_distr)
  B <- matrix(runif(p^2, -1, 1), nrow = p)
  diag(B) <- 1
  if (l > 0) {
    Gamma <- matrix(runif(p*l, -1, 1), nrow = l)
    L <-  calculate_eps(n, l, noise_distr)
    (eps + L %*% Gamma) %*% solve(B)
  } else {
    eps %*% solve(B)
  }
}

calculate_eps = function(n, p, distr_eps = c("gamma", "beta", "overlapping_gaussian", "gaussian", "uniform")) {
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
  } else if (distr_eps == "uniform") {
    as <- runif(p, 0.5, 2)
    eps <- mapply(function(n, a) runif(n, -a, a), n, as)
  }
  return(eps)
}

sample_overlapping_gaussian = function(n, sd_1, sd_2, mean_diff) {
  # Sample such that order is random
  return(sample(c(rnorm(n/2, 0, sd_1), rnorm(n/2, mean_diff, sd_2)) - 1/2 * mean_diff))
}