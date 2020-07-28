#' Bayesian Polya-gamma regression
#'
#' this function runs the Bayesian multinomial regression using Polya-gamma data augmentation
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
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the initial values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param n_chain is the MCMC chain id. The default is 1.
#'
#' @export
#'
#' @importFrom mvnfast rmvn
#' @importFrom fields rdist
#' @importFrom Matrix Cholesky
#' @importFrom stats lm rgamma
#' @import spam

## polya-gamma spatial linear regression model
mcmc_mra_integrated <- function(
    y,
    X,
    locs,
    params,
    priors,
    M = 4,
    n_coarse_grid = 10,
    # n_max_fine_grid = 2^6,
    n_cores = 1L,
    inits   = NULL,
    config  = NULL,
    verbose = FALSE,
    use_spam = TRUE, ## use spam or Matrix for sparse matrix operations
    n_chain       = 1
    # pool_s2_tau2  = true,
    # file_name     = "DM-fit",
    # corr_function = "exponential"
) {

    ##
    ## Run error checks
    ##

    # check_input_spatial(Y, X, locs)
    # check_params(param)
    # check_inits_pgLM(params, inits)
    # check_config(params, config)

    ##
    ## setup config
    ##

    ## do we sample the functional relationship parameters? This is primarily
    ## used to troubleshoot model fitting using simulated data

    sample_beta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_beta']])) {
            sample_beta <- config[['sample_beta']]
        }
    }


    # sample_rho <- TRUE
    # if (!is.null(config)) {
    #     if (!is.null(config[['sample_rho']])) {
    #         sample_rho <- config[['sample_rho']]
    #     }
    # }

    ## do we sample the MRA variance parameter? This is primarily
    ## used to troubleshoot model fitting using simulated data
    sample_tau2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_tau2']])) {
            sample_tau2 <- config[['sample_tau2']]
        }
    }

    sample_lambda <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_lambda']])) {
            sample_lambda <- config[['sample_lambda']]
        }
    }

    ## do we sample the nugger variance parameter
    sample_sigma2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma2']])) {
            sample_sigma2 <- config[['sample_sigma2']]
        }
    }

    ## add in a counter for the number of regularized Cholesky
    num_chol_failures <- 0

    N      <- length(y)
    # n_time <- dim(Y)[3]
    p      <- ncol(X)

    ##
    ## setup MRA spatial basis
    ##

    MRA      <- mra_wendland_2d(locs, M, n_coarse_grid = n_coarse_grid, use_spam = use_spam)
    # MRA      <- mra_wendland_2d(locs, M, n_max_fine_grid = n_max_fine_grid, use_spam = use_spam)
    W_list   <- MRA$W
    tW_list  <- vector(mode = 'list', length = M)
    tWW_list <- vector(mode = 'list', length = M)
    for (m in 1:M) {
        if (use_spam) {
            tW_list[[m]] <- t(W_list[[m]])
        } else {
            tW_list[[m]] <- Matrix::t(W_list[[m]])
        }
        tWW_list[[m]] <- tW_list[[m]] %*% W_list[[m]]
    }
    n_dims <- rep(NA, length(W_list))
    dims_idx <- c()
    for (i in 1:M) {
        n_dims[i] <- ncol(W_list[[i]])
        dims_idx <- c(dims_idx, rep(i, n_dims[i]))
    }
    W <- do.call(cbind, W_list)

    tW <- NULL
    if (use_spam) {
        tW <- t(W)
    } else {
        stop ('Only support use_spam = TRUE')
        tW <- Matrix::t(W)
    }

    tWW <- tW %*% W


    ##
    ## initial values
    ##

    tX  <- t(X)
    tXX <- tX %*% X

    ##
    ## priors for beta
    ##

    mu_beta        <- rep(0, p)
    Sigma_beta     <- 10 * diag(p)

    ## check if priors for mu_beta are specified
    if (!is.null(priors[['mu_beta']])) {
        if (all(!is.na(priors[['mu_beta']]))) {
            mu_beta <- priors[['mu_beta']]
        }
    }

    ## check if priors for Sigma_beta are specified
    if (!is.null(priors[['Sigma_beta']])) {
        if (all(!is.na(priors[['Sigma_beta']]))) {
            Sigma_beta <- priors[['Sigma_beta']]
        }
    }
    Sigma_beta_chol <- tryCatch(
        chol(Sigma_beta),
        error = function(e) {
            if (verbose)
                message("The Cholesky decomposition of the prior covariance Sigma_beta was ill-conditioned and mildy regularized.")
            chol(Sigma_beta + 1e-8 * diag(p))
        }
    )
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)

    ##
    ## initialize beta
    ##

    # beta   <- as.vector(rmvn(1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    beta   <- as.vector(lm (y ~ X - 1)$coeff)
    X_beta <- X %*% beta

    ##
    ## intialize a proper CAR structure to initialize the parameter alpha
    ##

    Q_alpha <- make_Q_alpha_2d(sqrt(n_dims), rep(0.999, length(n_dims)), use_spam = use_spam)

    ## prior for lambda is scale_lambda
    ## fix this so it can be set by the user
    scale_lambda <- 0.5

    lambda <- rgamma(M, 0.5, scale_lambda)
    tau2   <- rep(1, M)
    #100 * pmin(pmax(1 / rgamma(M, 0.5, lambda), 1), 100)

    ##
    ## initialize rho
    ##

    # rho      <- runif(1, -1, 1)

    ##
    ## initialize sigma2
    ##

    sigma2  <- pmax(pmin(1 / rgamma(1, priors$alpha_sigma2, priors$beta_sigma2), 5), 0.1)
    sigma   <- sqrt(sigma2)

    ## intialize an ICAR structure for fitting alpha
    Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(1, length(n_dims)), use_spam = use_spam)
    Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)


    ##
    ## sampler config options -- to be added later
    ##

    ##
    ## check for initial values
    ##

    ## initial values for beta
    if (!is.null(inits[['beta']])) {
        if (all(!is.na(inits[['beta']]))) {
            beta <- inits[['beta']]
            if(!is_numeric_vector(beta, p))
                stop("initial value for beta must be a numeric vector of length p")
        }
    }

    ## intial values for sigma2
    if (!is.null(inits[['sigma2']])) {
        if (all(!is.na(inits[['sigma2']]))) {
            sigma2 <- inits[['sigma2']]
            if(!is_positive_numeric(sigma2, 1))
                stop("initial value for sigma2 must be a positive numeric value")
        }
    }

    ## initial values for alpha
    if (!is.null(inits[['alpha']])) {
        if (all(!is.na(inits[['alpha']]))) {
            alpha <- inits[['alpha']]
        }
    }

    ## initial values for tau2
    if (!is.null(inits[['tau2']])) {
        if (all(!is.na(inits[['tau2']]))) {
            tau2 <- inits[['tau2']]
            if(!is_positive_numeric(tau2, M))
                stop("initial value for tau2 must be a positive numeric vector of length M")
        }

    }

    ## initial values for rho
    # if (!is.null(inits[['rho']])) {
    #     if (all(!is.na(inits[['rho']]))) {
    #         rho <- inits[['rho']]
    #     }
    # }

    ##
    ## setup save variables
    ##

    n_save       <- params$n_mcmc / params$n_thin
    beta_save    <- matrix(0, n_save, p)
    # rho_save     <- rep(0, n_save)
    tau2_save    <- matrix(0, n_save, M)
    sigma2_save  <- rep(0, n_save)
    lambda_save  <- matrix(0, n_save, M)

    ##
    ## initialize tuning
    ##

    ##
    ## tuning variables for adaptive MCMC
    ##

    # ## tuning for rho
    # rho_accept       <- 0
    # rho_accept_batch <- 0
    # rho_tune         <- 0.025

    ##
    ## Starting MCMC chain
    ##

    message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations \n")

    for (k in 1:(params$n_adapt + params$n_mcmc)) {
        if (k == params$n_adapt + 1) {
            message("Starting MCMC fitting for chain ", n_chain, ", running for ", params$n_mcmc, " iterations \n")
        }
        if (k %% params$n_message == 0) {
            if (k <= params$n_adapt) {
                message("MCMC adaptation iteration ", k, " for chain ", n_chain)
            } else {
                message("MCMC fitting iteration ", k - params$n_adapt, " for chain ", n_chain)
            }
        }

        ##
        ## sample sigma2
        ##

        if (sample_sigma2) {
            if (verbose)
                message("sample sigma2")

            devs   <- y - X_beta - W_alpha
            SS     <- sum(devs^2)
            sigma2 <- 1 / rgamma(1, N / 2 + priors$alpha_sigma2, SS / 2 + priors$beta_sigma2)
            sigma  <- sqrt(sigma2)
        }

        ##
        ## sample beta -- double check these values
        ##

        if (sample_beta) {
            if (verbose)
                message("sample beta")

            ## only use the modern climate state to update beta
            A      <- 1 / sigma2 * tXX + Sigma_beta_inv
            b      <- 1 / sigma2 * tX %*% (y - W_alpha) + Sigma_beta_inv %*% mu_beta
            ## guarantee a symmetric matrix
            A      <- (A + t(A)) / 2
            beta   <- rmvn_arma(A, b)
            X_beta <- X %*% beta
        }

        ##
        ## sample alpha
        ##

        # if (sample_alpha) {
        #     if (verbose)
        #         message("sample alpha")
        #     A_alpha <- 1 / sigma2 * tWW + Q_alpha_tau2
        #     b_alpha <- 1 / sigma2 * tW %*% (y - X_beta)
        #     alpha   <- as.vector(rmvnorm.canonical(1, b_alpha, A_alpha, Rstruct = Rstruct))
        # }
        #
        # ## update W_alpha
        # W_alpha <- W %*% alpha

        ##
        ## sample rho
        ##

        # if (sample_rho) {
        #     if (verbose)
        #         message("sample rho")
        #
        #     rho_vals <- rowSums(
        #         sapply(2:n_time, function(tt) {
        #             t_alpha_Q <- t(alpha[, tt-1]) %*% Q_alpha_tau2
        #             c(
        #                 t_alpha_Q %*% alpha[, tt-1],
        #                 t_alpha_Q %*% alpha[, tt]
        #             )
        #         })
        #     )
        #
        #     a_rho <- rho_vals[1]
        #     b_rho <- rho_vals[2]
        #     rho   <- rtrunc(1, "norm", a = -1, b = 1, mean = b_rho / a_rho, sd = sqrt(1 / a_rho))
        # }

        ##
        ## sample tau2
        ##

        if (sample_tau2) {
            if (verbose)
                message("sample tau2")
            for (m in 1:M) {
                devs    <- alpha[dims_idx == m]
                SS      <- as.numeric(devs %*% (Q_alpha[[m]] %*% devs))
                tau2[m] <- 1 / rgamma(1, 0.5 + n_dims[m] / 2, lambda[m] + SS / 2)
                # tau2_inv[m] <- 1 / rgamma(1, 0.5 + n_dims[m] / 2, lambda[m] + SS / 2)
                # tau2[m]     <- 1 / tau2_inv[m]
            }
        }

        Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)

        ##
        ## sample lambda
        ##

        if (sample_lambda) {
            if (verbose)
                message("sample lambda")
            for (m in 1:M) {
                # lambda[m]     <- rgamma(1, 1, scale_lambda + 1 / tau2[m])
                lambda[m]     <- rgamma(1, 1, scale_lambda + tau2[m])
            }
        }

        ##
        ## save variables
        ##

        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, ]   <- beta
                # rho_save[save_idx]       <- rho
                tau2_save[save_idx, ]   <- tau2
                sigma2_save[save_idx]   <- sigma2
                alpha_save[save_idx, ]  <- alpha
                lambda_save[save_idx, ] <- lambda
            }
        }

        ##
        ## End of MCMC loop
        ##
    }

    ## print out acceptance rates -- no tuning in this model


    if (num_chol_failures > 0)
        warning("The Cholesky decomposition for the update of alpha was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore.")

    ##
    ## return the MCMC output -- think about a better way to make this a class
    ##

    out <- list(
        beta     = beta_save,
        # rho      = rho_save,
        tau2     = tau2_save,
        sigma2   = sigma2_save,
        alpha    = alpha_save,
        lambda   = lambda_save,
        MRA      = MRA
    )

    class(out) <- "mcmc_mra"

    return(out)
}
