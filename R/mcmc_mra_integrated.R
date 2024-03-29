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
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the initial values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param verbose Should verbose output be printed? Typically this is only useful for troubleshooting.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE} (see spam and Matrix packages for details).
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
#' out     <- mcmc_mra(
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
#' @importFrom mvnfast rmvn dmvn
#' @importFrom fields rdist
#' @importFrom Matrix Cholesky
#' @importFrom stats lm rgamma dgamma sd rnorm runif
#' @import spam

mcmc_mra_integrated <- function(
    y,
    X,
    locs,
    params,
    priors        = NULL,
    M             = 4,
    n_neighbors   = 68,
    n_coarse_grid = 10,
    n_padding     = 5L,
    n_cores       = 1L,
    inits         = NULL,
    config        = NULL,
    verbose       = FALSE,
    use_spam      = TRUE, ## use spam or Matrix for sparse matrix operations
    n_chain       = 1
) {

    ##
    ## Run error checks
    ##

    if (!is_numeric_vector(y, length(y)))
        stop("y must be a numeric vector of length N.")
    if (length(y) != nrow(X))
        stop("X must have the same number of rows as the length of y.")
    if (!is_numeric_matrix(X, length(y), ncol(X)))
        stop("X must be a numeric matrix with N rows.")
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

    ## check the priors list


    ## check mcmc input
    if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
        stop("verbose must be either TRUE or FALSE.")
    }

    if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
        stop("use_spam must be either TRUE or FALSE.")
    }

    if (!is_positive_integer(n_cores, 1)) {
        stop("n_cores must be a positive integer")
    }

    n_cores <- as.integer(n_cores)

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

    ## do we sample the nugger variance parameter
    sample_sigma <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma']])) {
            sample_sigma <- config[['sample_sigma']]
            if (!is.logical(sample_sigma) | is.na(sample_sigma))
                stop('If specified, sample_sigma must be TRUE or FALSE')
        }
    }

    ## add in a counter for the number of regularized Cholesky
    num_chol_failures <- 0

    N      <- length(y)
    # n_time <- dim(Y)[3]
    p      <- ncol(X)

    ##
    ## center and scale the input and covariates
    ##

    sd_y <- sd(y)
    mu_y <- mean(y)
    y <- (y - mu_y) / sd_y

    if (ncol(X) >= 2) {
        mu_X <- apply(X[, -1, drop = FALSE], 2, mean)
        sd_X <- apply(X[, -1, drop = FALSE], 2, sd)
        for(i in 2:ncol(X)) {
            X[, i] <- (X[, i] - mu_X[i-1]) / sd_X[i-1]
        }
    }

    ##
    ## setup MRA spatial basis
    ##

    MRA <- mra_wendland_2d(
        locs          = locs,
        M             = M,
        n_coarse_grid = n_coarse_grid,
        n_padding     = n_padding,
        n_neighbors   = n_neighbors,
        use_spam      = use_spam
    )

    W        <- MRA$W
    n_dims   <- MRA$n_dims
    dims_idx <- MRA$dims_idx

    tW <- NULL
    if (use_spam) {
        tW <- t(W)
    } else {
        stop ('Only support use_spam = TRUE')
        # tW <- Matrix::t(W)
    }

    tWW <- tW %*% W


    ##
    ## initial values
    ##

    tX  <- t(X)
    tXX <- tX %*% X

    ##
    ## initialize sigma2
    ##

    alpha_sigma2 <- 0.01
    beta_sigma2  <- 0.01
    ## check if priors for alpha_sigma2 are specified
    if (!is.null(priors[['alpha_sigma2']])) {
        if (!is_positive_numeric(priors[['alpha_sigma2']], 1))
            stop("If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
        if (all(!is.na(priors[['alpha_sigma2']]))) {
            alpha_sigma2 <- priors[['alpha_sigma2']]
        }
    }
    ## check if priors for beta_sigma2 are specified
    if (!is.null(priors[['beta_sigma2']])) {
        if (!is_positive_numeric(priors[['beta_sigma2']], 1))
            stop("If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
        if (all(!is.na(priors[['beta_sigma2']]))) {
            beta_sigma2 <- priors[['beta_sigma2']]
        }
    }

    sigma  <- pmax(pmin(1 / rgamma(1, alpha_sigma2, beta_sigma2), 5), 0.1)
    sigma2   <- sigma^2

    ##
    ## priors for beta
    ##

    mu_beta        <- rep(0, p)
    Sigma_beta     <- 100 * diag(p)

    ## check if priors for mu_beta are specified
    if (!is.null(priors[['mu_beta']])) {
        if(!is_numeric_vector(priors[['mu_beta']], p))
            stop("If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
        if (all(!is.na(priors[['mu_beta']]))) {
            mu_beta <- priors[['mu_beta']]
        }
    }

    ## check if priors for Sigma_beta are specified
    if (!is.null(priors[['Sigma_beta']])) {
        if(!is_sympd_matrix(priors[['Sigma_beta']], p))
            stop("If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
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

    tau2   <- 10*2^(1:M - 1)

    alpha_tau2 <- 0.01
    beta_tau2  <- 0.01
    ## check if priors for alpha_tau2 are specified
    if (!is.null(priors[['alpha_tau2']])) {
        if (!is_positive_numeric(priors[['alpha_tau2']], 1))
            stop("If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
        if (all(!is.na(priors[['alpha_tau2']]))) {
            alpha_tau2 <- priors[['alpha_tau2']]
        }
    }
    ## check if priors for beta_tau2 are specified
    if (!is.null(priors[['beta_tau2']])) {
        if (!is_positive_numeric(priors[['beta_tau2']], 1))
            stop("If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
        if (all(!is.na(priors[['beta_tau2']]))) {
            beta_tau2 <- priors[['beta_tau2']]
        }
    }

    Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)

    ## initialize the cholesky structure
    G       <- Q_alpha_tau2 + tWW / sigma2
    Rstruct <- chol(G)

    ## intialize an ICAR structure for fitting alpha
    Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(0.9, length(n_dims)), use_spam = use_spam)
    # Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(1, length(n_dims)), use_spam = use_spam)
    Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)

    ##
    ## sampler config options -- to be added later
    ##

    ##
    ## check for initial values
    ##

    ## initial values for beta
    if (!is.null(inits[['beta']])) {
        if(!is_numeric_vector(inits[['beta']], p))
            stop("initial value for beta must be a numeric vector of length p")
        if (all(!is.na(inits[['beta']]))) {
            beta <- inits[['beta']]
        }
    }

    ## intial values for sigma2
    if (!is.null(inits[['sigma2']])) {
        if(!is_positive_numeric(inits[['sigma2']], 1))
            stop("initial value for sigma2 must be a positive numeric value")
        if (all(!is.na(inits[['sigma2']]))) {
            sigma2 <- inits[['sigma2']]
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

    ##
    ## initialize the log likelihood
    ##

    ll_current <- dmvn_smw(y, X, beta, tW, tWW, Q_alpha_tau2, sigma2, Rstruct = Rstruct)

    ##
    ## setup save variables
    ##

    n_save       <- params$n_mcmc / params$n_thin
    beta_save    <- matrix(0, n_save, p)
    tau2_save    <- matrix(0, n_save, M)
    sigma2_save  <- rep(0, n_save)
    ll_save      <- rep(0, n_save)

    ##
    ## initialize tuning
    ##

    # tuning for sigma2
    sigma2_accept       <- 0
    sigma2_accept_batch <- 0
    sigma2_tune         <- 0.25

    # tuning for tau2
    tau2_accept       <- 0
    tau2_accept_batch <- 0
    tau2_batch        <- matrix(0, 50, M)
    lambda_tau2_tune     <- 1.0 / 3.0^0.8
    Sigma_tau2_tune      <- diag(M)
    Sigma_tau2_tune_chol <- chol(Sigma_tau2_tune)

    # tuning for beta
    beta_accept          <- 0.0
    beta_accept_batch    <- 0.0
    beta_batch           <- matrix(0, 50, p)
    lambda_beta_tune     <- 1.0 / 3.0^0.8
    Sigma_beta_tune      <- diag(p)
    Sigma_beta_tune_chol <- chol(Sigma_beta_tune)

    ##
    ## tuning variables for adaptive MCMC
    ##

    ##
    ## Starting MCMC chain
    ##

    message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations")

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
        ## sample sigma2
        ##

        if (sample_sigma) {
            if (verbose)
                message("sample sigma2")

            # devs <- y - X_beta - W_alpha
            # SS        <- sum(devs^2)
            # sigma2 <- 1 / rgamma(1, N / 2 + alpha_sigma2, SS / 2 + beta_sigma2)
            # sigma     <- sqrt(sigma2)

            ##
            ## regular scale
            ##

            # sigma2_star <- rnorm(1, sigma2, sigma_tune)
            # if (sigma2_star > 0) {
            #     ll_proposal <- dmvn_smw(y, X, beta, tW, tWW, Q_alpha_tau2, sigma2_star, Rstruct = Rstruct)
            #     mh1 <- ll_proposal +
            #         dgamma(1 / sigma2_star, alpha_sigma2, beta_sigma2, log = TRUE)
            #     mh2 <- ll_current +
            #         dgamma(1 / sigma2, alpha_sigma2, beta_sigma2, log = TRUE)
            #     mh <- exp(mh1 - mh2)
            #     if (mh > runif(1, 0.0, 1.0)) {
            #         sigma2     <- sigma2_star
            #         ll_current <- ll_proposal
            #         if (k <= params$n_adapt) {
            #             sigma_accept_batch <- sigma_accept_batch + 1.0 / 50.0
            #         } else {
            #             sigma_accept <- sigma_accept + 1.0 / params$n_mcmc
            #         }
            #     }
            # }

            ##
            ## log scale proposal for sigma2
            ##

            sigma2_star  <- exp(rnorm(1, log(sigma2), sigma2_tune))
            sigma_star <- sqrt(sigma2_star)
            # if (sigma2_star > 0) {
            ll_proposal <- dmvn_smw(y, X, beta, tW, tWW, Q_alpha_tau2, sigma2_star, Rstruct = Rstruct)
            mh1 <- ll_proposal +
                dgamma(1 / sigma2_star, alpha_sigma2, beta_sigma2, log = TRUE) +
                log(sigma_star)
            mh2 <- ll_current +
                dgamma(1 / sigma2, alpha_sigma2, beta_sigma2, log = TRUE) +
                log(sigma)
            mh <- exp(mh1 - mh2)
            if (mh > runif(1, 0.0, 1.0)) {
                sigma     <- sigma_star
                sigma2    <- sigma2_star
                ll_current <- ll_proposal
                if (k <= params$n_adapt) {
                    sigma2_accept_batch <- sigma2_accept_batch + 1.0 / 50.0
                } else {
                    sigma2_accept <- sigma2_accept + 1.0 / params$n_mcmc
                }
                # }
            }

            ## update tuning
            if (k <= params$n_adapt) {
                if (k %% 50 == 0){
                    out_tuning <- update_tuning(
                        k,
                        sigma2_accept_batch,
                        sigma2_tune
                    )
                    sigma2_tune         <- out_tuning$tune
                    sigma2_accept_batch <- out_tuning$accept
                }
            }
        }

        ##
        ## sample beta using block Metropolis Hastings
        ##

        if (sample_beta) {
            if (verbose)
                message("sample beta")

            beta_star <- c(
                rmvn(
                    n      = 1,
                    mu     = beta,
                    sigma  = lambda_beta_tune * Sigma_beta_tune_chol,
                    isChol = TRUE
                )
            )

            ll_proposal <- dmvn_smw(y, X, beta_star, tW, tWW, Q_alpha_tau2, sigma2, Rstruct = Rstruct)
            mh1 <- ll_proposal +
                dmvn(beta_star, mu_beta, Sigma_beta_chol, isChol = TRUE, log = TRUE)
            mh2 <- ll_current +
                dmvn(beta, mu_beta, Sigma_beta_chol, isChol = TRUE, log = TRUE)
            mh <- exp(mh1 - mh2)
            if (mh > runif(1, 0, 1)) {
                beta       <- beta_star
                ll_current <- ll_proposal
                if (k <= params$n_adapt) {
                    beta_accept_batch <- beta_accept_batch + 1 / 50
                } else {
                    beta_accept <- beta_accept + 1 / params$n_mcmc
                }
            }

            ##
            ## Update tuning for beta
            ##

            if (k <= params$n_adapt) {
                ## update tuning
                beta_batch[k %% 50, ] <- t(beta)
                if (k %% 50 == 0) {
                    out_tuning <- update_tuning_mv(
                        k               = k,
                        accept          = beta_accept_batch,
                        lambda          = lambda_beta_tune,
                        batch_samples   = beta_batch,
                        Sigma_tune      = Sigma_beta_tune,
                        Sigma_tune_chol = Sigma_beta_tune_chol
                    )
                    beta_accept_batch    <- out_tuning$accept
                    lambda_beta_tune     <- out_tuning$lambda
                    beta_batch           <- out_tuning$batch_samples
                    Sigma_beta_tune      <- out_tuning$Sigma_tune
                    Sigma_beta_tune_chol <- out_tuning$Sigma_tune_chol
                }
            }
        }

        ##
        ## sample tau2
        ##

        if (sample_tau2) {
            if (verbose)
                message("sample tau2")

            tau2_star <- c(
                rmvn(
                    n      = 1,
                    mu     = tau2,
                    sigma  = lambda_tau2_tune * Sigma_tau2_tune_chol,
                    isChol = TRUE
                )
            )

            if (all(tau2_star > 0)) {
                Q_alpha_tau2_star <- make_Q_alpha_tau2(Q_alpha, tau2_star, use_spam = use_spam)

                ll_proposal <- dmvn_smw(y, X, beta, tW, tWW, Q_alpha_tau2_star, sigma2, Rstruct = Rstruct)

                mh1 <- ll_proposal +
                    sum(dgamma(tau2_star, alpha_tau2, beta_tau2, log = TRUE))
                mh2 <- ll_current +
                    sum(dgamma(tau2, alpha_tau2, beta_tau2, log = TRUE))
                mh <- exp(mh1 - mh2)
                if (mh > runif(1, 0.0, 1.0)) {
                    tau2         <- tau2_star
                    Q_alpha_tau2 <- Q_alpha_tau2_star
                    ll_current   <- ll_proposal
                    if (k <= params$n_adapt) {
                        tau2_accept_batch <- tau2_accept_batch + 1.0 / 50.0
                    } else {
                        tau2_accept <- tau2_accept + 1.0 / params$n_mcmc
                    }
                }
            }

            ##
            ## Update tuning for tau2
            ##

            if (k <= params$n_adapt) {
                ## update tuning
                tau2_batch[k %% 50, ] <- t(tau2)
                if (k %% 50 == 0) {
                    out_tuning <- update_tuning_mv(
                        k               = k,
                        accept          = tau2_accept_batch,
                        lambda          = lambda_tau2_tune,
                        batch_samples   = tau2_batch,
                        Sigma_tune      = Sigma_tau2_tune,
                        Sigma_tune_chol = Sigma_tau2_tune_chol
                    )
                    tau2_accept_batch    <- out_tuning$accept
                    lambda_tau2_tune     <- out_tuning$lambda
                    tau2_batch           <- out_tuning$batch_samples
                    Sigma_tau2_tune      <- out_tuning$Sigma_tune
                    Sigma_tau2_tune_chol <- out_tuning$Sigma_tune_chol
                }
            }
        }

        # Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)

        ##
        ## save variables
        ##

        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, ]   <- beta
                tau2_save[save_idx, ]   <- tau2
                sigma2_save[save_idx]   <- sigma2
                ll_save[save_idx]       <- ll_current
            }
        }

        ##
        ## End of MCMC loop
        ##
    }

    ## print out acceptance rates -- no tuning in this model
    message("Average acceptance rate for beta is ", mean(beta_accept))
    message("Average acceptance rate for sigma2 is ", mean(sigma2_accept))
    message("Average acceptance rate for tau2 is ", mean(tau2_accept))


    if (num_chol_failures > 0)
        warning("The Cholesky decomposition for the update of alpha was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore.")

    ##
    ## return the MCMC output -- think about a better way to make this a class
    ##

    out <- list(
        beta     = beta_save,
        tau2     = tau2_save,
        sigma2   = sigma2_save,
        ll       = ll_save,
        MRA      = MRA,
        data   = list(
            y    = y * sd_y + mu_y,
            mu_y = mu_y,
            sd_y = sd_y,
            X    = X,
            mu_X = mu_X,
            sd_X = sd_X,
            locs = locs),
        model  = list(
            params     = params,
            priors     = priors,
            inits      = inits,
            config     = config,
            verbose    = verbose,
            use_spam   = use_spam,
            n_chain    = n_chain)
    )

    class(out) <- "mcmc_mra_integrated"

    return(out)
}
