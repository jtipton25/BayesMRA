
#' Transform the beta parameters to restore the parameters to the original scale.
#' The function `mcmc_mra()` fits the regression coefficients using centered and
#' scaled response and covariate values. This function undoes that transformation
#' to scale the parameter estimates to their usual (uncentered and unscaled) scale
#' estimates
#'
#' @param out The output of `mcmc_mra()`
#'
#' @return A matrix of the size of parameters in out$beta of transformed
#' regression parameters equivalent to fitting the model on the original scale
#' @export
#'
unscale_beta <- function(out) {
    K <- out$model$params$n_mcmc / out$model$params$n_thin
    p <- ncol(out$data$X)
    mu_y <- out$data$mu_y
    sd_y <- out$data$sd_y
    mu_X <- out$data$mu_X
    sd_X <- out$data$sd_X
    beta <- matrix(0, K, p)
    for (k in 1:K) {
        # intercept
        beta[k, 1] <- out$beta[k, 1] * sd_y + mu_y + sd_y * out$beta[k, -1] * (mu_X / sd_X)
        if (p > 1) {
            beta[k, -1] <- sd_y * out$beta[k, -1] / sd_X
        }
    }
    return(beta)
}


#' Transform the sigma2 parameter to restore the parameter to the original scale.
#' The function `mcmc_mra()` fits the regression variance using a centered and
#' scaled response. This function undoes that transformation to scale the
#' parameter estimate to its usual (uncentered and unscaled) scale
#'
#' @param out The output of `mcmc_mra()`
#'
#' @return A vector of the size of out$sigma2 of transformed regression
#' variance equivalent to fitting the model on the original scale
#' @export
#'
unscale_sigma2 <- function(out) {
    K <- out$model$params$n_mcmc / out$model$params$n_thin
    p <- ncol(out$data$X)
    sd_y <- out$data$sd_y
    sigma2 <- out$sigma2 * sd_y^2
    return(sigma2)
}

