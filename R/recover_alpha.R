
#' Recover the latent spatial parameters alpha
#'
#' @param out The output from `mcmc_mra_integrated()`
#' @param n_message A positive integer number of frequency of iterations
#' @param constraint What constraint should be applied to the spatial process? Options include no constraint (`constraint = "unconstrained"`), a constraint so the entire MRA process sums to 0 (`constraint = "overall"`), a constraint so that each of the M levels of the MRA process sum to 0 (`constraint = "resolution"`), or whether the predicte process must sum to 0 (`constraint = "predicted"`). Note: `constraint = "predicted"` is NOT currently supported.
#'
#' @return
#' @export
#'
#' @import spam
#' @import spam64
#'

recover_alpha <- function(out, n_message = 100, constraint = "unconstrained") {

    if (class(out) != "mcmc_mra_integrated") {
        stop('out must be of class "mcmc_mra_integrated" which is the output of the function mcmc_mra_integrated()')
    }
    ## check the constraints on alpha
    if (!(constraint %in% c("unconstrained", "overall", "resolution", "predicted"))) {
        stop('constraint must be either "unconstrained", "overall", "resolution", or "predicted"')
    }
    if (constraint == "predicted") {
        stop('constraint = "predicted" is not currently supported -- developer note: add W_pred to function call to enable this in future results')
    }

    # initialize the process for compositional sampling of alpha
    tW           <- t(out$MRA$W)
    tWW          <- tW %*% out$MRA$W
    Q_alpha      <- make_Q_alpha_2d(sqrt(out$MRA$n_dims), rep(1, out$MRA$M))
    Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, out$tau2[1, ])
    A_alpha      <- 1 / out$sigma2[1] * tWW + Q_alpha_tau2
    ## precalculate the sparse cholesky structure for faster Gibbs updates
    Rstruct      <- chol(A_alpha)
    alpha        <- matrix(0, nrow(out$beta), sum(out$MRA$n_dims))
    constraints  <- make_constraint(out$MRA, constraint = constraint)
    A_constraint <- constraints$A_constraint
    a_constraint <- constraints$a_constraint

    for (j in 1:nrow(out$beta)) {
        # can parallelize this if desired
        if (j %% n_message == 0) {
            message("Compsition sampling recovery of alpha for Iteration ", j, " out of ", nrow(out$beta))
        }
        Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, out$tau2[j, ])
        A_alpha      <- 1 / out$sigma2[j] * tWW + Q_alpha_tau2
        b_alpha      <- 1 / out$sigma2[j] * tW %*% ((out$data$y - out$data$mu_y) / out$data$sd_y - out$data$X %*% out$beta[j, ])
        # alpha   <- as.vector(rmvnorm.canonical(1, b_alpha, A_alpha, Rstruct = Rstruct))
        ## sample constrained to sum to 0
        if (constraint == "unconstrained") {
            alpha[j, ] <- as.vector(rmvnorm.canonical(1, b_alpha, A_alpha, Rstruct = Rstruct))
        } else if (constraint %in% c("resolution", "overall", "predicted")) {
            alpha[j, ] <- as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha, Rstruct = Rstruct,
                                                            A = A_constraint, a = a_constraint))
        }
    }
    return(alpha)
}
