#' Calulate basis function normalizing constant
#'
#' Calculate the basis function normalizing constant so that all basis functions have the same marginal variance
#'
#' @param W The sparse basis function matrix
#' @param Q_inv The inverse of the sparse precision matrix
#'
#' @return A set of normalizing weights for the basis functions
#' @export
#'
calc_omega_norm <- function(W, Q_inv) {
    omega <- sqrt(apply(W * (W %*% Q_inv), 1, sum))
    return(omega)
}

#' Calculate normalized basis function
#'
#' Calculate the normalized basis function matrix so that all basis functions have the same marginal variance
#'
#' @param W The sparse basis function matrix
#' @param omega The normalizing weights for the basis functions which is the output of \code{\link{calc_omega_norm()}}
#' @param tau2 The basis function resolution precision
#'
#' @return A normalized basis function matrix
#' @export
#'
calc_W_star <- function(W, omega, tau2) {
    W_star <- spam_diag(sqrt(tau2) / omega) %*% W
    return(W_star)
}
