#' dmvn_smw_mra
#'
#' @param y is a \eqn{n}{n} vector of Gaussian data.
#' @param X is a \eqn{n \times p}{n x p} matrix of fixed effects (like latitude, elevation, etc)
#' @param beta is the regression parameter beta
#' @param tW is the transpose of the sparse Wendland basis matrix of class spam
#' @param tWW is the transpose of the sparse Wendland basis matrix multiplied by itself of class spam
#' @param Q_alpha_tau2 is the prior precision matrix for random effects alpha
#' @param sigma2 is the residual error
#' @param Rstruct is the Cholesky prior precision matrix for random effects alpha
#' @param logd is a logical value of whether to calculate the log density (`logd = TRUE`) or the density on the data scale (`logd = FALSE`)
#'
#' @return The (log) density of a normal distribution with mean \eqn{\mathbf{X} \boldsymbol{\beta}} and covariance matrix \eqn{\sigma^2 \mathbf{I} + \mathbf{W} \mathbf{Q}_{\alpha_{\tau^2}} \mathbf{W}'}
#' @import spam
#' @import spam64
#' @export
dmvn_smw <- function(y, X, beta, tW, tWW, Q_alpha_tau2, sigma2, Rstruct = NULL, logd = TRUE) {
    n <- length(y)

    devs <- y - X %*% beta
    ## calculate the Sherman-Morrison-Woodbury inverse

    G <- Q_alpha_tau2 + tWW / sigma2
    G_chol <- NULL
    if (is.null(Rstruct)) {
        G_chol <- chol(G)
    } else {
        G_chol <- update.spam.chol.NgPeyton(Rstruct, G)
    }

    ## calculate the Sherman-Morrison-Woodbury determinant
    Sigma_det <- 2.0 * sum(log(diag(G_chol))) -
        determinant(Q_alpha_tau2, logartithm = TRUE)$modulus[1] +
        n * log(sigma2)

    ## calculate the Sherman-Morrison-Woodbury log likelihood
    ll <- - 0.5 * n * log(2 * base::pi) -
                 0.5 * Sigma_det -
                 0.5 / sigma2 * (sum(devs^2) - sum(forwardsolve(G_chol, transpose = TRUE, tW %*% devs, upper.tri = TRUE)^2) / sigma2)

    if (logd) {
        return(ll)
    } else {
        return(exp(ll))
    }
}
