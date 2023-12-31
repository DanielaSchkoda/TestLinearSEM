# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

estimateMoment <- function(moment, X) {
    .Call(`_TestLinearSEM_estimateMoment`, moment, X)
}

estimateCumulant <- function(cumulant, X) {
    .Call(`_TestLinearSEM_estimateCumulant`, cumulant, X)
}

estimateCovBtwMoments <- function(moment_1, moment_2, X) {
    .Call(`_TestLinearSEM_estimateCovBtwMoments`, moment_1, moment_2, X)
}

estimate_Moments <- function(Pi, X) {
    .Call(`_TestLinearSEM_estimate_Moments`, Pi, X)
}

estimate_Covs <- function(Pi_vectorized, X) {
    .Call(`_TestLinearSEM_estimate_Covs`, Pi_vectorized, X)
}

random_combs_c <- function(n, k, nr) {
    .Call(`_TestLinearSEM_random_combs_c`, n, k, nr)
}

permutations <- function(n) {
    .Call(`_TestLinearSEM_permutations`, n)
}

estimate_polynomial <- function(poly, L) {
    .Call(`_TestLinearSEM_estimate_polynomial`, poly, L)
}

H_c <- function(X, indices, polynomials) {
    .Call(`_TestLinearSEM_H_c`, X, indices, polynomials)
}

G_c <- function(X, order_kernel, polynomials) {
    .Call(`_TestLinearSEM_G_c`, X, order_kernel, polynomials)
}

