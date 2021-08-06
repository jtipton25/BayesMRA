# TODO check the sum-to-zero constraint for each mini-grid
# TODO check the normalizing of the basis functions


#' Bayesian Multi-resolution Spatial Regression
#'
#' this function runs Markov Chain Monte Carlo to estimate the Bayesian multi-resolution spatial regression model.
#' @param y is a \eqn{n}{n} vector of Gaussian data.
#' @param X is a \eqn{n \times p}{n x p} matrix of fixed effects (like latitude, elevation, etc)
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param params is a list of parameter settings. The list
#' \code{params} must contain the following values:
#' * \code{n_adapt}: A positive integer number of adaptive MCMC iterations.
#' * \code{n_mcmc}: A positive integer number of total MCMC iterations
#' post adaptation.
#' * \code{n_thin}: A positive integer number of MCMC iterations per saved
#' sample.
#' * \code{n_message}: A positive integer number of frequency of iterations
#'  to output a progress message. For example, \code{n_message = 50}
#'  outputs progress messages every 50 iterations.
#' @param priors is the list of prior settings.
#' @param M The number of resolutions.
#' @param n_neighbors The expected number of neighbors for each interior basis function. This determines the basis radius parameter.
#' @param n_coarse_grid The number of basis functions in one direction (e.g. \code{n_coarse_grid = 10} results in a \eqn{10 \times 10}{10x10} course grid which is further extended by the number of additional padding basis functions given by \code{n_padding}.
#' @param n_padding The number of additional boundary points to add on each boundary. For example, n_padding = 5 will add 5 boundary knots to the both the left  and right side of the grid).
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the initial values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param verbose Should verbose output be printed? Typically this is only useful for troubleshooting.
#' @param constraint What constraint should be applied to the spatial process? Options include no constraint (`constraint = "unconstrained"`), a constraint so the entire MRA process sums to 0 (`constraint = "overall"`), a constraint so that each of the M levels of the MRA process sum to 0 (`constraint = "resolution"`), or whether the predicted process must sum to 0 (`constraint = "predicted"`). Note: `constraint = "predicted"` is NOT currently supported.
#' @param n_chain is the MCMC chain id. The default is 1.
#'
#' @examples
#' set.seed(111)
#' ## genereate the locations
#' locs <- matrix(runif(20), 10, 2)
#' ## generate some covariates and regression coefficients
#' X <- cbind(1, matrix(rnorm(30), 10, 3))
#' beta <- rnorm(ncol(X))
#'
#' ## simulate the MRA process
#' M <- 2
#' MRA <- mra_wendland_2d(locs, M = M, n_coarse_grid = 4)
#' W <- MRA$W
#' n_dims <- MRA$n_dims
#' dims_idx <- MRA$dims_idx
#'
#' ## set up the process precision matrices
#' Q_alpha <- make_Q_alpha_2d(sqrt(n_dims), c(0.9, 0.8))
#' Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2 = c(2, 4))
#'
#' ## add in constraints so each resolution has random effects that sum to 0
#' constraints <- make_constraint(MRA, constraint = "resolution", joint = TRUE)
#' A_constraint <- constraints$A_constraint
#' a_constraint <- constraints$a_constraint
#' alpha <- as.vector(spam::rmvnorm.prec.const(
#'     n = 1,
#'     mu = rep(0, nrow(W)),
#'     Q = Q_alpha_tau2,
#'     A = A_constraint,
#'     a = a_constraint))
#' ## define the data
#' y <-  as.vector(X %*% beta + W %*% alpha + rnorm(10))
#'
#' ## define the params for MCMC fitting
#' params <- list(
#'     n_mcmc    = 5,
#'     n_adapt   = 5,
#'     n_thin    = 1,
#'     n_message = 5)
#'
#' ## define the model priors
#' priors <- list(
#'     alpha_tau2   = 1,
#'     beta_tau2    = 1,
#'     alpha_sigma2 = 1,
#'     beta_sigma2  = 1,
#'     mu_beta      = rep(0, ncol(X)),
#'     Sigma_beta   = 5 * diag(ncol(X)))
#'
#' ## Fit the MRA model using MCMC
#' out     <- mcmc_mra_vecchia(
#'     y             = y,
#'     X             = X,
#'     locs          = locs,
#'     params        = params,
#'     priors        = priors,
#'     M             = 2,
#'     n_coarse_grid = 4
#' )
#'
#' @export
#'
#' @importFrom mvnfast rmvn
#' @importFrom fields fields.rdist.near
#' @importFrom Matrix Cholesky
#' @importFrom stats lm rgamma sd
#' @import spam
#' @import spam64

mcmc_mra_vecchia <- function(
    y,
    # X,
    locs,
    params,
    M             = 3,
    n_neighbors   = 68,
    n_coarse_grid = 10,
    n_padding     = 5L,
    inits         = NULL,
    normalize     = TRUE,
    verbose       = FALSE,
    config        = NULL,
    constraint    = "overall",
    n_chain       = 1
) {


    ##
    ## Run error checks
    ##

    if (!is_numeric_vector(y, length(y)))
        stop("y must be a numeric vector of length N.")
    # if (length(y) != nrow(X))
    #     stop("X must have the same number of rows as the length of y.")
    # if (!is_numeric_matrix(X, length(y), ncol(X)))
    #     stop("X must be a numeric matrix with N rows.")
    if (!is_numeric_matrix(locs, length(y), 2))
        stop("locs must be a numeric matrix with N rows and 2 columns.")

    ## check the params list
    if (!is_positive_integer(params$n_adapt, 1))
        stop("params must contain a positive integer n_adapt.")
    if (!is_positive_integer(params$n_mcmc, 1))
        stop("params must contain a positive integer n_mcmc.")
    if (!is_positive_integer(params$n_thin, 1))
        stop("params must contain a positive integer n_thin.")
    if (!is_positive_integer(params$n_message, 1))
        stop("params must contain a positive integer n_message.")

    params$n_adapt   <- as.integer(params$n_adapt)
    params$n_mcmc    <- as.integer(params$n_mcmc)
    params$n_thin    <- as.integer(params$n_thin)
    params$n_message <- as.integer(params$n_message)

    ## check the constraints on alpha
    if (!(constraint %in% c("unconstrained", "overall", "resolution", "predicted"))) {
        stop('constraint must be either "unconstrained", "overall", "resolution", or "predicted"')
    }
    if (constraint == "predicted") {
        stop('constraint = "predicted" is not currently supported -- developer note: add W_pred to function call to enable this in future results')
    }

    # check MCMC input
    if (!is.logical(normalize) || length(normalize) != 1 || is.na(normalize)) {
        stop("normalize must be either TRUE or FALSE.")
    }

    if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
        stop("verbose must be either TRUE or FALSE.")
    }

    if (!is_positive_integer(n_chain, 1)) {
        stop("n_chain must be a positive integer")
    }

    n_chain <- as.integer(n_chain)

    ##
    ## setup config
    ##

    ## do we sample the functional relationship parameters? This is primarily
    ## used to troubleshoot model fitting using simulated data

    sample_beta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_beta']])) {
            sample_beta <- config[['sample_beta']]
            if (!is.logical(sample_beta) | is.na(sample_beta))
                stop('If specified, sample_beta must be TRUE or FALSE')
        }
    }

    ## do we sample the MRA variance parameter? This is primarily
    ## used to troubleshoot model fitting using simulated data
    sample_tau2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_tau2']])) {
            sample_tau2 <- config[['sample_tau2']]
            if (!is.logical(sample_tau2) | is.na(sample_tau2))
                stop('If specified, sample_tau2 must be TRUE or FALSE')
        }
    }

    ## do we sample the nugget variance parameter
    sample_sigma2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma2']])) {
            sample_sigma2 <- config[['sample_sigma2']]
            if (!is.logical(sample_sigma2) | is.na(sample_sigma2))
                stop('If specified, sample_sigma2 must be TRUE or FALSE')
        }
    }

    ## do we sample the latent intensity parameter alpha
    sample_alpha <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_alpha']])) {
            sample_alpha <- config[['sample_alpha']]
            if (!is.logical(sample_alpha) | is.na(sample_alpha))
                stop('If specified, sample_alpha must be TRUE or FALSE')
        }
    }

    ##
    ## setup the MRA process -- turn this into a function
    ##

    W_list <- vector(mode = "list", length = M)
    tau2 <- 100 / (2^(1:M - 1))
    Walpha <- matrix(0, N, M)
    alpha <- vector(mode = "list", length = M)
    sigma2 <- 5

    for (m in 1:M) {
        W_list[[m]] <- vector(mode = "list", length = 4^(m-1))

        if (m == 1) {
            grid_idx <- partition_observations(locs, m)
            MRA <- mra_wendland_2d(locs, M = 1, n_padding = n_padding)
            constraints <- make_constraint(MRA, constraint = constraint, joint = TRUE)
            A_constraint <- constraints$A_constraint
            a_constraint <- constraints$a_constraint
            Q_alpha <- list(make_Q_alpha_2d(sqrt(MRA$n_dims), phi = 0.999))
            class(Q_alpha) <- "Q_alpha"
            class(Q_alpha[[1]]) <- "spam"
            Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2[m])
            W  <- MRA$W
            tW <- t(MRA$W)
            tWW <- t(MRA$W) %*% MRA$W
            G <- tWW + Q_alpha_tau2
            Rstruct <- chol(G)

            Q_alpha_inv <- NULL
            omega       <- NULL
            W_star      <- NULL
            tW_star     <- NULL
            tWW_star    <- NULL
            if (normalize) {
                Q_alpha_inv <- solve.spam(Q_alpha[[1]])
                omega <- calc_omega_norm(W, Q_alpha_inv)
                W_star <- calc_W_star(W, omega, tau2[m])
                tW_star <- t(W_star)
                tWW_star <- tW_star %*% W_star
            }
            if(normalize) {
                W_list[[m]] <- list(MRA = MRA, W = W, tW = tW, tWW = tWW,
                                    W_star = W_star, tW_star = tW_star, tWW_star = tWW_star,
                                    Rstruct = Rstruct, Q_alpha_inv = Q_alpha_inv, omega = omega,
                                    Q_alpha = Q_alpha, Q_alpha_tau2 = Q_alpha_tau2, grid_idx = grid_idx,
                                    A_constraint = A_constraint, a_constraint = a_constraint)
            } else {
                W_list[[m]] <- list(MRA = MRA, W = W, tW = tW, tWW = tWW, Rstruct = Rstruct,
                                    Q_alpha = Q_alpha, Q_alpha_tau2 = Q_alpha_tau2, grid_idx = grid_idx,
                                    A_constraint = A_constraint, a_constraint = a_constraint)
            }
        } else {
            grid_idx <- partition_observations(locs, m)
            for (j in 1:(4^(m-1))) {
                MRA <- mra_wendland_2d(locs[grid_idx == j, ], M = 1, n_padding = n_padding)
                constraints <- make_constraint(MRA, constraint = constraint, joint = TRUE)
                A_constraint <- constraints$A_constraint
                a_constraint <- constraints$a_constraint
                Q_alpha <- list(make_Q_alpha_2d(sqrt(MRA$n_dims), phi = 0.999))
                class(Q_alpha) <- "Q_alpha"
                class(Q_alpha[[1]]) <- "spam"
                Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2[m])
                W  <- MRA$W
                tW <- t(MRA$W)
                tWW <- t(MRA$W) %*% MRA$W
                G <- tWW + Q_alpha_tau2
                Rstruct = chol(G)


                Q_alpha_inv <- NULL
                omega       <- NULL
                W_star      <- NULL
                tW_star     <- NULL
                tWW_star    <- NULL
                if (normalize) {
                    Q_alpha_inv <- solve.spam(Q_alpha[[1]])
                    omega <- calc_omega_norm(W, Q_alpha_inv)
                    W_star <- calc_W_star(W, omega, tau2[m])
                    tW_star <- t(W_star)
                    tWW_star <- tW_star %*% W_star
                }

                if (normalize) {
                    W_list[[m]][[j]] <- list(MRA = MRA, W = W, tW = tW, tWW = tWW,
                                             W_star = W_star, tW_star = tW_star, tWW_star = tWW_star,
                                             Rstruct = Rstruct, Q_alpha_inv = Q_alpha_inv, omega = omega,
                                             Q_alpha = Q_alpha, Q_alpha_tau2 = Q_alpha_tau2, grid_idx = grid_idx,
                                             A_constraint = A_constraint, a_constraint = a_constraint)
                } else {
                    W_list[[m]][[j]] <- list(MRA = MRA, W = W, tW = tW, tWW = tWW, Rstruct = Rstruct,
                                             Q_alpha = Q_alpha, Q_alpha_tau2 = Q_alpha_tau2, grid_idx = grid_idx,
                                             A_constraint = A_constraint, a_constraint = a_constraint)
                }
            }
        }
    }

    ##
    ## initialize the random variables ----
    ##

    # can I think of a more effective way to store the alpha parameters?
    for (m in 1:M) {
        if (m == 1) {
            alpha[[m]] <- drop(rmvnorm.prec.const(1, rep(0, W_list[[m]]$MRA$n_dims), W_list[[m]]$Q_alpha_tau2 + 1e-4 * diag(W_list[[m]]$MRA$n_dims),
                                                  A = W_list[[m]]$A_constraint, a = W_list[[m]]$a_constraint))
            if (normalize) {
                Walpha[, m] <- W_list[[m]]$W_star %*% alpha[[m]]
            } else {
                Walpha[, m] <- W_list[[m]]$MRA$W %*% alpha[[m]]
            }
        } else {
            for (j in 1:(4^(m-1))) {
                alpha[[m]][[j]] <- drop(rmvnorm.prec(1, rep(0, W_list[[m]][[j]]$MRA$n_dims), W_list[[m]][[j]]$Q_alpha_tau2 + 1e-4 * diag(W_list[[m]][[j]]$MRA$n_dims),
                                                     A = W_list[[m]][[j]]$A_constraint, a = W_list[[m]][[j]]$a_constraint))
                if (normalize) {
                    Walpha[W_list[[m]][[j]]$grid_idx == j, m] <- W_list[[m]][[j]]$W_star %*% alpha[[m]][[j]]
                } else {
                    Walpha[W_list[[m]][[j]]$grid_idx == j, m] <- W_list[[m]][[j]]$MRA$W %*% alpha[[m]][[j]]
                }
            }
        }
    }

    ##
    ## check for initial values ----
    ##

    ## initial values for beta
    # if (!is.null(inits[['beta']])) {
    #     if(!is_numeric_vector(inits[['beta']], p))
    #         stop("initial value for beta must be a numeric vector of length p")
    #     if (all(!is.na(inits[['beta']]))) {
    #         beta <- inits[['beta']]
    #     }
    # }
    # X_beta <- X %*% beta

    ## intial values for sigma2
    if (!is.null(inits[['sigma2']])) {
        if(!is_positive_numeric(inits[['sigma2']], 1))
            stop("initial value for sigma2 must be a positive numeric value")
        if (all(!is.na(inits[['sigma2']]))) {
            sigma2 <- inits[['sigma2']]
        }
    }

    ## initial values for alpha
    if (!is.null(inits[['alpha']])) {
        # add in appropriate size checks for alpha
        # if(!is_numeric_vector(inits[['alpha']], length(alpha)))
            # stop("initial value for alpha must be positive numeric vector of length equal to the number of all grid points")
        if (all(!is.na(inits[['alpha']]))) {
            alpha <- inits[['alpha']]
        }
    }

    for (m in 1:M) {
        if (m == 1) {
            if (normalize) {
                Walpha[, m] <- W_list[[m]]$W_star %*% alpha[[m]]
            } else {
                Walpha[, m] <- W_list[[m]]$MRA$W %*% alpha[[m]]
            }
        } else {
            for (j in 1:(4^(m-1))) {
                if (normalize) {
                    Walpha[W_list[[m]][[j]]$grid_idx == j, m] <- W_list[[m]][[j]]$W_star %*% alpha[[m]][[j]]
                } else {
                    Walpha[W_list[[m]][[j]]$grid_idx == j, m] <- W_list[[m]][[j]]$MRA$W %*% alpha[[m]][[j]]
                }
            }
        }
    }

    ## initial values for tau2
    if (!is.null(inits[['tau2']])) {
        if(!is_positive_numeric(inits[['tau2']], M) | !is.vector(inits[['tau2']]))
            stop("initial value for tau2 must be a positive numeric vector of length M")
        if (all(!is.na(inits[['tau2']]))) {
            tau2 <- inits[['tau2']]
        }
    }
    for (m in 1:M) {
        if (m == 1) {
            W_list[[m]]$Q_alpha_tau2 <- make_Q_alpha_tau2(W_list[[m]]$Q_alpha, tau2[m])
            if (normalize) {
                W_list[[m]]$W_star <- calc_W_star(W_list[[m]]$W, W_list[[m]]$omega, tau2[m])
                W_list[[m]]$tW_star <- t(W_list[[m]]$W_star)
                W_list[[m]]$tWW_star <- W_list[[m]]$tW_star %*% W_list[[m]]$W_star
            }
        } else {
            for (j in 1:(4^(m-1))) {
                W_list[[m]][[j]]$Q_alpha_tau2 <- make_Q_alpha_tau2(W_list[[m]][[j]]$Q_alpha, tau2[m])
                if (normalize) {
                    W_list[[m]][[j]]$W_star <- calc_W_star(W_list[[m]][[j]]$W, W_list[[m]][[j]]$omega, tau2[m])
                    W_list[[m]][[j]]$tW_star <- t(W_list[[m]][[j]]$W_star)
                    W_list[[m]][[j]]$tWW_star <- W_list[[m]][[j]]$tW_star %*% W_list[[m]][[j]]$W_star
                }
            }
        }
    }


    ##
    ## set up the save variables ----
    ##

    n_save <- params$n_mcmc / params$n_thin

    alpha_save <- vector(mode = "list", length = M)
    for (m in 1:M) {
        if (m == 1) {
            alpha_save[[m]] <- matrix(0, n_save, length(alpha[[m]]))
        } else {
            alpha_save[[m]] <- array(0, dim = c(n_save, length(alpha[[m]][[1]]), 4^(m-1)))
        }
    }
    Walpha_save <- array(0, dim = c(n_save, dim(Walpha)))
    tau2_save <- matrix(0, n_save, M)
    mu_save <- matrix(0, n_save, N)
    sigma2_save <- rep(0, n_save)


    message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations")
    message("Normalization = ", normalize)

    for (k in 1:(params$n_adapt + params$n_mcmc)) {
        if (k == params$n_adapt + 1) {
            message("Starting MCMC fitting for chain ", n_chain, ", running for ", params$n_mcmc, " iterations")
        }
        if (k %% params$n_message == 0) {
            if (k <= params$n_adapt) {
                message("MCMC adaptation iteration ", k, " for chain ", n_chain)
            } else {
                message("MCMC fitting iteration ", k - params$n_adapt, " for chain ", n_chain)
            }
        }

        ##
        ## sample alpha ----
        ##

        if (sample_alpha) {
            if (verbose)
                message("sampling alpha")
            for (m in 1:M) {
                if (m == 1) {
                    if (normalize) {
                        A <- 1 / sigma2 * W_list[[m]]$tWW_star + W_list[[m]]$Q_alpha_tau2
                        b <- 1 / sigma2 * W_list[[m]]$tW_star %*% (y - apply(Walpha[, -m, drop = FALSE], 1, sum))
                        alpha[[m]] <- drop(rmvnorm.canonical.const(1, b, A,
                                                                   A = W_list[[m]]$A_constraint,
                                                                   a = W_list[[m]]$a_constraint,
                                                                   Rstruct = W_list[[m]]$Rstruct))
                        Walpha[, m] <- W_list[[m]]$W_star %*% alpha[[m]]
                    } else {
                        A <- 1 / sigma2 * W_list[[m]]$tWW + W_list[[m]]$Q_alpha_tau2
                        b <- 1 / sigma2 * W_list[[m]]$tW %*% (y - apply(Walpha[, -m, drop = FALSE], 1, sum))
                        alpha[[m]] <- drop(rmvnorm.canonical.const(1, b, A,
                                                                   A = W_list[[m]]$A_constraint,
                                                                   a = W_list[[m]]$a_constraint,
                                                                   Rstruct = W_list[[m]]$Rstruct))
                        Walpha[, m] <- W_list[[m]]$W %*% alpha[[m]]
                    }
                } else {
                    for (j in 1:(4^(m-1))) {
                        if (normalize) {
                            A <- 1 / sigma2 * W_list[[m]][[j]]$tWW_star + W_list[[m]][[j]]$Q_alpha_tau2
                            b <- 1 / sigma2 * W_list[[m]][[j]]$tW_star %*% (y[W_list[[m]][[j]]$grid_idx == j] - apply(Walpha[W_list[[m]][[j]]$grid_idx == j, -m, drop = FALSE], 1, sum))
                            alpha[[m]][[j]] <- drop(rmvnorm.canonical.const(1, b, A,
                                                                            A = W_list[[m]][[j]]$A_constraint,
                                                                            a = W_list[[m]][[j]]$a_constraint,
                                                                            Rstruct = W_list[[m]]$Rstruct))
                            Walpha[W_list[[m]][[j]]$grid_idx == j, m] <- W_list[[m]][[j]]$W_star %*% alpha[[m]][[j]]
                        } else {
                            A <- 1 / sigma2 * W_list[[m]][[j]]$tWW + W_list[[m]][[j]]$Q_alpha_tau2
                            b <- 1 / sigma2 * W_list[[m]][[j]]$tW %*% (y[W_list[[m]][[j]]$grid_idx == j] - apply(Walpha[W_list[[m]][[j]]$grid_idx == j, -m, drop = FALSE], 1, sum))
                            alpha[[m]][[j]] <- drop(rmvnorm.canonical.const(1, b, A,
                                                                            A = W_list[[m]][[j]]$A_constraint,
                                                                            a = W_list[[m]][[j]]$a_constraint,
                                                                            Rstruct = W_list[[m]]$Rstruct))
                            Walpha[W_list[[m]][[j]]$grid_idx == j, m] <- W_list[[m]][[j]]$W %*% alpha[[m]][[j]]
                        }
                    }
                }
            }
        }

        ##
        ## sample sigma2 ----
        ##

        if (sample_sigma2) {
            if (verbose)
                message("sampling sigma2")
            sigma2 <- 1 / rgamma(1, 0.1 + 0.5 * N, 0.1 + 0.5 * sum((y - apply(Walpha, 1, sum))^2))
        }

        ##
        ## sample tau2 ----
        ##

        if (sample_tau2) {
            if (verbose)
                message("sampling tau2")
            for (m in 1:M) {
                if (m == 1) {
                    SS <- t(alpha[[m]]) %*% W_list[[m]]$Q_alpha[[1]] %*% alpha[[m]]
                    tau2[m] <- rgamma(1, 0.1 + 0.5 * length(alpha[[m]]), 0.1 + 0.5 * SS)
                    W_list[[m]]$Q_alpha_tau2 <- make_Q_alpha_tau2(W_list[[m]]$Q_alpha, tau2[m])
                    if (normalize) {
                        W_list[[m]]$W_star <- calc_W_star(W_list[[m]]$W, W_list[[m]]$omega, tau2[m])
                        W_list[[m]]$tW_star <- t(W_list[[m]]$W_star)
                        W_list[[m]]$tWW_star <- W_list[[m]]$tW_star %*% W_list[[m]]$W_star
                    }
                } else {
                    SS <- 0
                    for (j in 1:(4^(m-1))) {
                        SS <- SS + t(alpha[[m]][[j]]) %*% W_list[[m]][[j]]$Q_alpha[[1]] %*% alpha[[m]][[j]]
                    }
                    tau2[m] <- rgamma(1, 0.1 + 0.5 * length(alpha[[m]][[j]]) * (4^(m-1)), 0.1 + 0.5 * SS)
                    for (j in 1:(4^(m-1))) {
                        W_list[[m]][[j]]$Q_alpha_tau2 <- make_Q_alpha_tau2(W_list[[m]][[j]]$Q_alpha, tau2[m])
                        if (normalize) {
                            W_list[[m]][[j]]$W_star <- calc_W_star(W_list[[m]][[j]]$W, W_list[[m]][[j]]$omega, tau2[m])
                            W_list[[m]][[j]]$tW_star <- t(W_list[[m]][[j]]$W_star)
                            W_list[[m]][[j]]$tWW_star <- W_list[[m]][[j]]$tW_star %*% W_list[[m]][[j]]$W_star
                        }
                    }
                }
            }
        }

        ##
        ## save the variables ----
        ##

        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx <- (k - params$n_adapt) / params$n_thin
                for (m in 1:M) {
                    if (m == 1) {
                        alpha_save[[m]][save_idx, ] <- alpha[[m]]
                    } else {
                        for (j in 1:(4^(m-1))) {
                            alpha_save[[m]][save_idx, , j] <- alpha[[m]][[j]]
                        }
                    }
                }
                Walpha_save[save_idx, , ] <- Walpha
                mu_save[save_idx, ]      <- apply(Walpha, 1, sum)
                sigma2_save[save_idx]    <- sigma2
                tau2_save[save_idx, ]    <- tau2
            }
        }

        # End of MCMC loop ----
    }

    ##
    ## save MCMC output ----
    ##

    out <- NULL
    if (normalize) {
        out <- list(alpha  = alpha_save,
                    Walpha = Walpha_save,
                    mu     = mu_save,
                    sigma2 = sigma2_save,
                    tau2   = tau2_save,
                    omega  = omega)
    } else {
        out <- list(alpha  = alpha_save,
                    Walpha = Walpha_save,
                    mu     = mu_save,
                    sigma2 = sigma2_save,
                    tau2   = tau2_save)
    }

    class(out) <- "mcmc_mra_vecchia"
    return(out)
}
