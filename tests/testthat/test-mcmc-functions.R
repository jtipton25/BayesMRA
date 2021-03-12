context("mcmc functions")

# mcmc_mra ---------------------------------------------------------------------

test_that("mcmc_mra", {
    y <- rnorm(10)
    X <- matrix(1:20, 10, 2)
    params <- list(
        n_mcmc    = 5,
        n_adapt   = 5,
        n_thin    = 1,
        n_message = 5
    )
    locs <- matrix(1:20, 10, 2)

    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )
    inits <- list(
        beta = rep(0, ncol(X)),
        sigma2 = 1,
        alpha = rep(0, 520), ## use M = 2 and n_coarse_grid = 4
        tau2  = 1:2
    )

    out <- mcmc_mra(y, X, locs, params, priors, inits = inits, M = 2, n_coarse_grid = 4)
    expect_s3_class(out, "mcmc_mra")
    expect_s3_class(out$MRA, "mra_wendland_2d")
    expect_error(mcmc_mra(y, X, locs, params, M = -3), "the number of resolutions M must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_coarse_grid = -3), "n_coarse_grid must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_neighbors = -3), "n_neighbors must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_padding = -3), "n_padding must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_cores = -3), "n_cores must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_cores = 1:4), "n_cores must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_cores = NA), "n_cores must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_cores = NULL), "n_cores must be a positive integer")

    expect_error(mcmc_mra(y, X, locs, params, n_chain = 1:4), "n_chain must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_chain = NA), "n_chain must be a positive integer")
    expect_error(mcmc_mra(y, X, locs, params, n_chain = NULL), "n_chain must be a positive integer")

    expect_error(mcmc_mra(y, X, locs, params, use_spam = FALSE), "The Matrix package is not currently supported")

    expect_error(mcmc_mra(y, X), 'argument "locs" is missing, with no default')
    expect_error(mcmc_mra(y, X, locs), 'argument "params" is missing, with no default')
    y <- matrix(1:10, 5, 2)
    expect_error(mcmc_mra(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep(NA, 10)
    expect_error(mcmc_mra(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep("aaa", 10)
    expect_error(mcmc_mra(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep(TRUE, 10)
    expect_error(mcmc_mra(y, X, locs, params), "y must be a numeric vector of length N.")

    y <- rnorm(10)
    X <- array(0, dim=c(2, 2, 2))
    expect_error(mcmc_mra(y, X, locs, params), "X must have the same number of rows as the length of y.")
    X <- matrix(0, 5, 3)
    expect_error(mcmc_mra(y, X, locs, params), "X must have the same number of rows as the length of y.")
    X <- array(0, dim=c(10, 2, 2))
    expect_error(mcmc_mra(y, X, locs, params), "X must be a numeric matrix with N rows.")
    X <- matrix(NA, 10, 3)
    expect_error(mcmc_mra(y, X, locs, params), "X must be a numeric matrix with N rows.")

    X <- matrix(1:20, 10, 2)
    locs <- matrix(0, 10, 3)
    expect_error(mcmc_mra(y, X, locs, params), "locs must be a numeric matrix with N rows and 2 columns.")
    locs <- matrix(NA, 10, 2)
    expect_error(mcmc_mra(y, X, locs, params), "locs must be a numeric matrix with N rows and 2 columns.")
    locs <- matrix(1:20, 10, 2)

    params$n_mcmc <- -10
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- NA
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- NULL
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- "aaa"
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- 50.5
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_mcmc.")

    params$n_mcmc <- 500
    params$n_adapt <- -10
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- NA
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- NULL
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- "aaa"
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- 50.5
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_adapt.")

    params$n_adapt <- 500
    params$n_thin <- -10
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- NA
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- NULL
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- "aaa"
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- 50.5
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_thin.")


    params$n_thin <- 1
    params$n_message <- -10
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- NA
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- NULL
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- "aaa"
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- 50.5
    expect_error(mcmc_mra(y, X, locs, params), "params must contain a positive integer n_message.")

    params$n_message = 100
    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )

    priors$mu_beta <- rep(0, ncol(X) + 1)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep(NA, ncol(X))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep("aa", ncol(X))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- matrix(rep(0, ncol(X)))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep(0, ncol(X))

    priors$Sigma_beta <- diag(ncol(X) + 1)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- matrix(1, ncol(X), ncol(X))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigam_beta <- diag(NA, ncol(X))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- matrix(rep(0, ncol(X), ncol(X)))
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- diag(ncol(X))

    priors$alpha_sigma2 <- rep(1, 2)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- -1
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- "aa"
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- 0
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- NA
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- 1

    priors$beta_sigma2 <- rep(1, 2)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- -1
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- "aa"
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- 0
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- NA
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- 1

    priors$alpha_tau2 <- rep(1, 2)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- -1
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- "aa"
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- 0
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- NA
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- 1

    priors$beta_tau2 <- rep(1, 2)
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- -1
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- "aa"
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- 0
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- NA
    expect_error(mcmc_mra(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- 1

    priors <- NULL

    expect_error(mcmc_mra(y, X, locs, params, use_spam = "aaa"), "use_spam must be either TRUE or FALSE.")
    expect_error(mcmc_mra(y, X, locs, params, use_spam = 3), "use_spam must be either TRUE or FALSE.")
    expect_error(mcmc_mra(y, X, locs, params, use_spam = NA), "use_spam must be either TRUE or FALSE.")

    expect_error(mcmc_mra(y, X, locs, params, verbose = "aaa"), "verbose must be either TRUE or FALSE.")
    expect_error(mcmc_mra(y, X, locs, params, verbose = 3), "verbose must be either TRUE or FALSE.")
    expect_error(mcmc_mra(y, X, locs, params, verbose = NA), "verbose must be either TRUE or FALSE.")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_beta = "TRUE")), "If specified, sample_beta must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_beta = NA)), "If specified, sample_beta must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_beta = 3)), "If specified, sample_beta must be TRUE or FALSE")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_tau2 = "TRUE")), "If specified, sample_tau2 must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_tau2 = NA)), "If specified, sample_tau2 must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_tau2 = 3)), "If specified, sample_tau2 must be TRUE or FALSE")

    # expect_error(mcmc_mra(y, X, locs, params, config = list(sample_lambda = "TRUE")), "If specified, sample_lambda must be TRUE or FALSE")
    # expect_error(mcmc_mra(y, X, locs, params, config = list(sample_lambda = NA)), "If specified, sample_lambda must be TRUE or FALSE")
    # expect_error(mcmc_mra(y, X, locs, params, config = list(sample_lambda = 3)), "If specified, sample_lambda must be TRUE or FALSE")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_sigma2 = "TRUE")), "If specified, sample_sigma2 must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_sigma2 = NA)), "If specified, sample_sigma2 must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_sigma2 = 3)), "If specified, sample_sigma2 must be TRUE or FALSE")

    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_alpha = "TRUE")), "If specified, sample_alpha must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_alpha = NA)), "If specified, sample_alpha must be TRUE or FALSE")
    expect_error(mcmc_mra(y, X, locs, params, config = list(sample_alpha = 3)), "If specified, sample_alpha must be TRUE or FALSE")


    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = 3)), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = rep("aaa", ncol(X)))), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = rep(NA, ncol(X)))), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = matrix(1:ncol(X)))), "initial value for beta must be a numeric vector of length p")

    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = 1:3)), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = matrix(1:ncol(X)))), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = NA)), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = -3)), "initial value for sigma2 must be a positive numeric value")

    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(alpha = 2:ncol(out$alpha))), "initial value for alpha must be positive numeric vector of length equal to the number of all grid points")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(alpha = matrix(1:ncol(out$alpha)))), "initial value for alpha must be positive numeric vector of length equal to the number of all grid points")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(alpha = c(2:ncol(out$alpha), NA))), "initial value for alpha must be positive numeric vector of length equal to the number of all grid points")

    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = 1:3)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = 1)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = matrix(1:ncol(X)))), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = NA)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = -3)), "initial value for tau2 must be a positive numeric vector of length M")

    # mcmc_mra(y, X, locs, params, priors,  M = 2, n_coarse_grid = 4)

    # expect_error(mcmc_mra(y, X, params), 'argument "priors" is missing, with no default')
    # priors <- default_priors_pg_lm(y, X)
    # out <- pg_lm(y, X, params, priors)
    # expect_true(class(out) == "pg_lm")

})


# mcmc_mra_integrated ----------------------------------------------------------

test_that("mcmc_mra_integrated", {
    y <- rnorm(10)
    X <- matrix(1:20, 10, 2)
    params <- list(
        n_mcmc    = 5,
        n_adapt   = 5,
        n_thin    = 1,
        n_message = 5
    )
    locs <- matrix(1:20, 10, 2)

    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )
    inits <- list(
        beta = rep(0, ncol(X)),
        sigma2 = 1,
        alpha = rep(0, 520), ## use M = 2 and n_coarse_grid = 4
        tau2  = 1:2
    )

    out <- mcmc_mra_integrated(y, X, locs, params, priors, inits = inits, M = 2, n_coarse_grid = 4)
    expect_s3_class(out, "mcmc_mra_integrated")
    expect_s3_class(out$MRA, "mra_wendland_2d")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = -3), "the number of resolutions M must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_coarse_grid = -3), "n_coarse_grid must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_neighbors = -3), "n_neighbors must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_padding = -3), "n_padding must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_cores = -3), "n_cores must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_cores = 1:4), "n_cores must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_cores = NA), "n_cores must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_cores = NULL), "n_cores must be a positive integer")

    expect_error(mcmc_mra_integrated(y, X, locs, params, n_chain = 1:4), "n_chain must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_chain = NA), "n_chain must be a positive integer")
    expect_error(mcmc_mra_integrated(y, X, locs, params, n_chain = NULL), "n_chain must be a positive integer")

    expect_error(mcmc_mra_integrated(y, X, locs, params, use_spam = FALSE), "The Matrix package is not currently supported")

    expect_error(mcmc_mra_integrated(y, X), 'argument "locs" is missing, with no default')
    expect_error(mcmc_mra_integrated(y, X, locs), 'argument "params" is missing, with no default')
    y <- matrix(1:10, 5, 2)
    expect_error(mcmc_mra_integrated(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep(NA, 10)
    expect_error(mcmc_mra_integrated(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep("aaa", 10)
    expect_error(mcmc_mra_integrated(y, X, locs, params), "y must be a numeric vector of length N.")
    y <- rep(TRUE, 10)
    expect_error(mcmc_mra_integrated(y, X, locs, params), "y must be a numeric vector of length N.")

    y <- rnorm(10)
    X <- array(0, dim=c(2, 2, 2))
    expect_error(mcmc_mra_integrated(y, X, locs, params), "X must have the same number of rows as the length of y.")
    X <- matrix(0, 5, 3)
    expect_error(mcmc_mra_integrated(y, X, locs, params), "X must have the same number of rows as the length of y.")
    X <- array(0, dim=c(10, 2, 2))
    expect_error(mcmc_mra_integrated(y, X, locs, params), "X must be a numeric matrix with N rows.")
    X <- matrix(NA, 10, 3)
    expect_error(mcmc_mra_integrated(y, X, locs, params), "X must be a numeric matrix with N rows.")

    X <- matrix(1:20, 10, 2)
    locs <- matrix(0, 10, 3)
    expect_error(mcmc_mra_integrated(y, X, locs, params), "locs must be a numeric matrix with N rows and 2 columns.")
    locs <- matrix(NA, 10, 2)
    expect_error(mcmc_mra_integrated(y, X, locs, params), "locs must be a numeric matrix with N rows and 2 columns.")
    locs <- matrix(1:20, 10, 2)

    params$n_mcmc <- -10
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- NULL
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- "aaa"
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_mcmc.")
    params$n_mcmc <- 50.5
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_mcmc.")

    params$n_mcmc <- 500
    params$n_adapt <- -10
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- NULL
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- "aaa"
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_adapt.")
    params$n_adapt <- 50.5
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_adapt.")

    params$n_adapt <- 500
    params$n_thin <- -10
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- NULL
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- "aaa"
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_thin.")
    params$n_thin <- 50.5
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_thin.")


    params$n_thin <- 1
    params$n_message <- -10
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- NULL
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- "aaa"
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_message.")
    params$n_message <- 50.5
    expect_error(mcmc_mra_integrated(y, X, locs, params), "params must contain a positive integer n_message.")

    params$n_message = 100
    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )

    priors$mu_beta <- rep(0, ncol(X) + 1)
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep(NA, ncol(X))
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep("aa", ncol(X))
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- matrix(rep(0, ncol(X)))
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter mu_beta in priors must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- rep(0, ncol(X))

    priors$Sigma_beta <- diag(ncol(X) + 1)
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- matrix(1, ncol(X), ncol(X))
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigam_beta <- diag(NA, ncol(X))
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- matrix(rep(0, ncol(X), ncol(X)))
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter Sigma_beta in priors must be a symmetric positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- diag(ncol(X))

    priors$alpha_sigma2 <- rep(1, 2)
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- -1
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- "aa"
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- 0
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma2 in priors must be a positive numeric value.")
    priors$alpha_sigma2 <- 1

    priors$beta_sigma2 <- rep(1, 2)
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- -1
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- "aa"
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- 0
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma2 in priors must be a positive numeric value.")
    priors$beta_sigma2 <- 1

    priors$alpha_tau2 <- rep(1, 2)
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- -1
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- "aa"
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- 0
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_tau2 in priors must be a positive numeric value.")
    priors$alpha_tau2 <- 1

    priors$beta_tau2 <- rep(1, 2)
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- -1
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- "aa"
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- 0
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_tau2 in priors must be a positive numeric value.")
    priors$beta_tau2 <- 1

    priors <- NULL

    expect_error(mcmc_mra_integrated(y, X, locs, params, use_spam = "aaa"), "use_spam must be either TRUE or FALSE.")
    expect_error(mcmc_mra_integrated(y, X, locs, params, use_spam = 3), "use_spam must be either TRUE or FALSE.")
    expect_error(mcmc_mra_integrated(y, X, locs, params, use_spam = NA), "use_spam must be either TRUE or FALSE.")

    expect_error(mcmc_mra_integrated(y, X, locs, params, verbose = "aaa"), "verbose must be either TRUE or FALSE.")
    expect_error(mcmc_mra_integrated(y, X, locs, params, verbose = 3), "verbose must be either TRUE or FALSE.")
    expect_error(mcmc_mra_integrated(y, X, locs, params, verbose = NA), "verbose must be either TRUE or FALSE.")

    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_beta = "TRUE")), "If specified, sample_beta must be TRUE or FALSE")
    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_beta = NA)), "If specified, sample_beta must be TRUE or FALSE")
    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_beta = 3)), "If specified, sample_beta must be TRUE or FALSE")

    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_tau2 = "TRUE")), "If specified, sample_tau2 must be TRUE or FALSE")
    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_tau2 = NA)), "If specified, sample_tau2 must be TRUE or FALSE")
    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_tau2 = 3)), "If specified, sample_tau2 must be TRUE or FALSE")

    # expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_lambda = "TRUE")), "If specified, sample_lambda must be TRUE or FALSE")
    # expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_lambda = NA)), "If specified, sample_lambda must be TRUE or FALSE")
    # expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_lambda = 3)), "If specified, sample_lambda must be TRUE or FALSE")

    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_sigma2 = "TRUE")), "If specified, sample_sigma2 must be TRUE or FALSE")
    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_sigma2 = NA)), "If specified, sample_sigma2 must be TRUE or FALSE")
    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_sigma2 = 3)), "If specified, sample_sigma2 must be TRUE or FALSE")


    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = 3)), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = rep("aaa", ncol(X)))), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = rep(NA, ncol(X)))), "initial value for beta must be a numeric vector of length p")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(beta = matrix(1:ncol(X)))), "initial value for beta must be a numeric vector of length p")

    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = 1:3)), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = matrix(1:ncol(X)))), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = NA)), "initial value for sigma2 must be a positive numeric value")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(sigma2 = -3)), "initial value for sigma2 must be a positive numeric value")

    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = 1:3)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = 1)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = matrix(1:ncol(X)))), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = NA)), "initial value for tau2 must be a positive numeric vector of length M")
    expect_error(mcmc_mra_integrated(y, X, locs, params, M = 2, n_coarse = 4, inits = list(tau2 = -3)), "initial value for tau2 must be a positive numeric vector of length M")

    # mcmc_mra_integrated(y, X, locs, params, priors,  M = 2, n_coarse_grid = 4)

    # expect_error(mcmc_mra_integrated(y, X, params), 'argument "priors" is missing, with no default')
    # priors <- default_priors_pg_lm(y, X)
    # out <- pg_lm(y, X, params, priors)
    # expect_true(class(out) == "pg_lm")

})

# update_tuning ----------------------------------------------------------------

test_that("update_tuning", {
    # update_tuning
    expect_equal(
        update_tuning(50, 0.6, 0.5),
        list(accept = 0, tune = 0.575954955084454))
    expect_equal(
        update_tuning(250, 0.6, 0.5),
        list(accept = 0, tune = 0.532644196047304))
    expect_equal(
        update_tuning(50, 0.3, 0.5),
        list(accept = 0, tune = 0.434061722697292))
    expect_equal(
        update_tuning(50, 0.3, 0.1),
        list(accept = 0, tune = 0.0868123445394585))

    # update_tuning_mv
    set.seed(111)
    k <- 50
    accept <- 0.6
    lambda <- 0.2
    batch_samples <- matrix(rnorm(150), 50, 3)
    Sigma_tune <- stats::rWishart(1, 5, diag(3))[, , 1]
    Sigma_tune_chol <- chol(Sigma_tune)
    expect_equal(
        update_tuning_mv(k, accept, lambda, batch_samples, Sigma_tune, Sigma_tune_chol),
        list(accept = 0, lambda = 0.668921761531105,
             batch_samples = structure(c(0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0), .Dim = c(50L, 3L)),
             Sigma_tune = structure(c(2.99149062979637,
                                      0.401601806709832, 0.424753903084598, 0.401601806709832, 3.88366262056398,
                                      0.594329214360549, 0.424753903084598, 0.594329214360549, 4.71125545103594
             ), .Dim = c(3L, 3L)),
             Sigma_tune_chol = structure(c(1.72959261960624,
                                           0, 0, 0.232194449812848, 1.95697428650457, 0, 0.245580316584201,
                                           0.274559983527489, 2.13905646830272), .Dim = c(3L, 3L))))
    accept <- 0.1
    expect_equal(
        update_tuning_mv(k, accept, lambda, batch_samples, Sigma_tune, Sigma_tune_chol),
        list(accept = 0, lambda = 0.12854540862658,
             batch_samples = structure(c(0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0), .Dim = c(50L, 3L)),
             Sigma_tune = structure(c(2.99149062979637,
                                      0.401601806709832, 0.424753903084598, 0.401601806709832, 3.88366262056398,
                                      0.594329214360549, 0.424753903084598, 0.594329214360549, 4.71125545103594
             ), .Dim = c(3L, 3L)),
             Sigma_tune_chol = structure(c(1.72959261960624,
                                           0, 0, 0.232194449812848, 1.95697428650457, 0, 0.245580316584201,
                                           0.274559983527489, 2.13905646830272), .Dim = c(3L, 3L))))

    # update_tuning_vec
    expect_equal(
        update_tuning_vec(50, c(0.1, 0.6), c(0.5, 0.5)),
        list(accept = c(0, 0), tune = c(0.434061722697292, 0.575954955084454)))
    expect_equal(
        update_tuning_vec(250, c(0.1, 0.6), c(0.25, 0.75)),
        list(accept = c(0, 0), tune = c(0.234678235354129, 0.798966294070957)))
    expect_equal(
        update_tuning_vec(50, c(0.6, 0.1), c(0.5, 0.5)),
        list(accept = c(0, 0), tune = c(0.575954955084454, 0.434061722697292)))
    expect_equal(
        update_tuning_vec(250, c(0.6, 0.1), c(0.25, 0.75)),
        list(accept = c(0, 0), tune = c(0.266322098023652, 0.704034706062386)))

    # update_tuning_mat
    expect_equal(
        update_tuning_mat(50, matrix(c(0.1, 0.6, 0.6, 0.1), 2, 2), matrix(c(0.5, 0.25, 1, 0.5), 2, 2)),
        list(accept = structure(c(0, 0, 0, 0), .Dim = c(2L, 2L)), tune = structure(c(0.434061722697292,
                                                                                     0.287977477542227, 1.15190991016891, 0.434061722697292), .Dim = c(2L,
                                                                                                                                                   2L))))
    # update_tuning_mv_mat
    set.seed(111)
    k <- 50
    accept <- c(0.1, 0.6, 0.6)
    lambda <- c(0.1, 0.2, 0.5)
    batch_samples <- array(rnorm(50*3*4), dim = c(50, 3, 4))
    Sigma_tune <- stats::rWishart(4, 5, diag(3))
    Sigma_tune_chol <- Sigma_tune
    for (j in 1:4) {
        Sigma_tune_chol[, , j] <- chol(Sigma_tune[, , j])
    }


    expect_equal(
        update_tuning_mv_mat(k, accept, lambda, batch_samples, Sigma_tune, Sigma_tune_chol),
        list(accept = c(0, 0, 0), lambda = c(0.0483971455170536, 0.50369599627459, 1.25923999068647),
             batch_samples = structure(c(0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0, 0), .Dim = c(50L, 3L, 4L)),
             Sigma_tune = structure(c(0.861843460847109,
                                      0.0350297118298129, -0.205904229686601, 0.0350297118298129, 5.23938194426575,
                                      0.113512965023776, -0.205904229686601, 0.113512965023776, 1.11610391509234,
                                      12.819834315426, 0.0307410097016759, -3.91159618091204, 0.0307410097016759,
                                      0.606063359369128, 0.0753128883635706, -3.91159618091204, 0.0753128883635706,
                                      3.74553987550437, 3.58825485621149, 1.11700532441496, 0.294659706303476,
                                      1.11700532441496, 6.08365678995394, 2.20167448094668, 0.294659706303476,
                                      2.20167448094668, 6.46870007854512, 5.0065339819057, 2.07428171851721,
                                      -3.28335619632131, 2.07428171851721, 2.42602278480844, -0.181070805667884,
                                      -3.28335619632131, -0.181070805667884, 4.16265859287012), .Dim = c(3L, 3L, 4L)),
             Sigma_tune_chol = structure(c(0.928355244961275, 0,
                                           0, 0.0377330898058039, 2.28865859363066, 0, -0.221794653290498,
                                           0.0532547593310908, 1.03154009978899, 3.5804796208645, 0, 0,
                                           0.00858572396908477, 0.778453367076702, 0, -1.09247826970388,
                                           0.108795990628847, 1.59379877592884, 1.89426895033717, 0, 0,
                                           0.5896762042244, 2.39498199661823, 0, 0.155553257762594, 0.880987176243409,
                                           2.38083280761989, 2.23752854326055, 0, 0, 0.927041455969336,
                                           1.25164568617588, 0, -1.46740304440397, 0.94217769624022, 1.05909777020423
             ), .Dim = c(3L, 3L, 4L))))
})


# predict_mra ------------------------------------------------------------------

test_that("predict_mra", {
    set.seed(111)
    y <- rnorm(10)
    X <- matrix(1:20, 10, 2)
    params <- list(
        n_mcmc    = 5,
        n_adapt   = 5,
        n_thin    = 1,
        n_message = 5
    )
    locs <- matrix(1:20, 10, 2)

    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )
    inits <- list(
        beta = rep(0, ncol(X)),
        sigma2 = 1,
        alpha = rep(0, 520), ## use M = 2 and n_coarse_grid = 4
        tau2  = 1:2
    )

    out <- mcmc_mra(y, X, locs, params, priors, inits = inits, M = 2, n_coarse_grid = 4)
    new_data <- list(locs_pred = locs, X_pred = X)
    # test equality here
    # expect_snapshot(predict_mra(out, new_data))
    # expect_equal(predict_mra(out, new_data),
    new_data <- list(locs_pred = locs, X_pred = NULL)
    expect_error(predict_mra(out, new_data), "new_data must be a list that contains a n_pred by p matrix of covariates in the variable X_pred which are used to make predictions")
    new_data <- list(locs_pred = NULL, X_pred = X)
    expect_error(predict_mra(out, new_data), "new_data must be a list that contains a n_pred by 2 matrix of locations in the variable locs_pred at which to make predictions")

    new_data <- list(locs_pred = locs, X_pred = cbind(X, 1))
    expect_error(predict_mra(out, new_data), "The columns of X_pred in the list new_data must be the same number as X used to fit the model")
    new_data <- list(locs_pred = cbind(locs, 1), X_pred = X)
    expect_error(predict_mra(out, new_data), "locs_pred must be a numeric matrix with n_pred rows and 2 columns")

    new_data <- list(locs_pred = locs, X_pred = rbind(X, 1))
    expect_error(predict_mra(out, new_data), "The number of rows in locs_pred in the list new_data must be equal to the number of rows in X_pred in the list new_data")
    new_data <- list(locs_pred = rbind(locs, 1), X_pred = X)
    expect_error(predict_mra(out, new_data), "The number of rows in locs_pred in the list new_data must be equal to the number of rows in X_pred in the list new_data")

    new_data <- list(locs_pred = locs, X_pred = X)
    class(out) <- "adsfads"
    expect_error(predict_mra(out, new_data), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")


})

# recover_alpha ----------------------------------------------------------------
test_that("mcmc_mra_integrated", {
    set.seed(111)
    y <- rnorm(10)
    X <- matrix(1:20, 10, 2)
    params <- list(
        n_mcmc    = 5,
        n_adapt   = 5,
        n_thin    = 1,
        n_message = 5
    )
    locs <- matrix(1:20, 10, 2)

    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )
    inits <- list(
        beta = rep(0, ncol(X)),
        sigma2 = 1,
        alpha = rep(0, 520), ## use M = 2 and n_coarse_grid = 4
        tau2  = 1:2
    )

    out <- mcmc_mra_integrated(y, X, locs, params, priors, inits = inits, M = 2, n_coarse_grid = 4)
    expect_equal(
        recover_alpha(out),
        structure(c(-64.7333866463099, -130.839970750028, -81.7630451607594,
                    20.0650235474602, -379.294611826064, -64.2990510786463, -132.505404176397,
                    -82.0737465128427, 21.2962810968424, -379.298300087779, -63.9213331419306,
                    -132.55783679027, -82.3390894186505, 19.9795739955749, -379.968919682045,
                    -64.6828490684561, -132.568570152654, -82.6835164334958, 20.3553658356635,
                    -379.67823724831, -64.7807972874185, -132.602382774529, -82.445784152543,
                    20.1570128074586, -380.390929076791, -64.8533195496013, -132.65806059458,
                    -82.6023096615292, 20.4647788931254, -380.733266296382, -65.4251412251807,
                    -133.03853273358, -81.9548225985245, 19.675779571268, -380.442634045589,
                    -66.32932836786, -132.644328946591, -82.021524244391, 19.2959511825487,
                    -380.560571940578, -66.3013426209743, -133.054606139242, -82.6593432452581,
                    19.7774166021312, -381.171941871622, -66.5935395630529, -133.611556235203,
                    -83.2873107598587, 20.3342038807348, -381.548583799912, -65.7838014176614,
                    -132.177909702956, -81.8625442076821, 19.5348856650733, -381.184887576088,
                    -65.7848031834382, -131.22610895694, -81.1951199512631, 19.0319687228632,
                    -380.693599554924, -66.4041016312409, -132.26666538173, -82.1769251374973,
                    18.3079632852364, -381.119733477059, -66.6541374623996, -131.470488019429,
                    -81.2405325170505, 18.2171991390801, -381.191616158872, -65.5065822472886,
                    -132.024660419568, -81.3294347122307, 20.2001359505615, -378.977505163057,
                    -64.8094008690011, -131.686726748763, -81.8090010657046, 20.9235629364437,
                    -379.498315081341, -65.3856052047179, -132.446738885159, -81.4492251794448,
                    19.7200792421361, -379.73110790472, -65.4963933669284, -132.778928637775,
                    -82.2878607846389, 20.2556772956924, -379.512531633246, -65.0462608121672,
                    -132.555251020851, -82.6941005288126, 20.3936624170144, -380.340326831371,
                    -65.6505221966375, -132.417378918041, -82.2128607018779, 20.1490440125241,
                    -380.032573256183, -66.2377649024772, -132.736673230609, -82.5065648079132,
                    19.60166454804, -380.477356195764, -66.6847354304493, -133.084881793987,
                    -82.9474380981397, 18.920809514848, -380.785352181432, -66.2004481944652,
                    -132.686144959527, -82.5131197476966, 19.9344280038371, -379.973550533895,
                    -65.8976614634331, -132.457422215091, -82.6447942934494, 20.2066397097815,
                    -381.236719713631, -66.1981903858843, -132.553756538668, -81.9832785830469,
                    19.2268947199064, -381.222604327969, -66.691476019507, -132.069762251667,
                    -82.3286732482982, 18.7884474310926, -380.446144275303, -66.9802399613722,
                    -132.76743780836, -81.6851245027929, 18.8135914863086, -380.96427699799,
                    -66.77147188876, -131.997603332844, -82.1129965675474, 18.6219884634683,
                    -380.610732705762, -65.713390479068, -131.948320886952, -82.3653365433017,
                    19.0167045919134, -379.660159879674, -65.2733765245733, -132.300365129082,
                    -81.4174391803004, 19.9523927732042, -379.471858325386, -65.9292030247524,
                    -132.01480822595, -81.894697024825, 19.9048230710124, -379.573042965516,
                    -65.8230684797929, -132.136333402033, -82.0246730700588, 20.2344142666813,
                    -379.506265092402, -65.7269473347927, -132.424556771408, -82.2342470105524,
                    19.9739257004719, -379.952188946703, -65.9133431243844, -132.678485237909,
                    -82.310239904739, 19.4801699653784, -380.45213757111, -66.5178601401541,
                    -133.300868713239, -82.3642666294538, 19.7493458210508, -380.047437522508,
                    -66.7478888141312, -132.492815321465, -82.9552406277819, 20.1685155686455,
                    -380.435833819445, -66.2528990476359, -132.190264094076, -82.2379644844429,
                    20.2895876428491, -380.000986755193, -66.5653774147342, -132.187938023323,
                    -82.4704181584185, 19.7058441231079, -380.689638976827, -67.1224661887427,
                    -131.954349870809, -82.743960617802, 18.863384304064, -380.44874050863,
                    -66.6791621501202, -132.200407258653, -82.4420698326459, 18.8143472903174,
                    -380.491077045769, -67.0780566924263, -132.149892369048, -82.1866980512954,
                    17.8366864162686, -380.55210022875, -67.6934555023236, -132.952881041791,
                    -82.7166206409655, 17.5659363100807, -380.110590449647, -65.880903224764,
                    -132.116495131276, -82.7253638427203, 20.5483946160073, -379.412214726035,
                    -65.5199075362484, -132.456727411969, -82.1928506827471, 19.8727973657607,
                    -378.80249695608, -65.1912292205822, -132.607105606092, -81.7593172931353,
                    20.22997113013, -380.447532650584, -65.9975661385125, -132.668645660318,
                    -80.8766404448296, 19.5089412801238, -379.987707857334, -66.4420641773468,
                    -132.651863570371, -82.0875942565958, 19.0531403772499, -379.70653037801,
                    -65.3692163466998, -132.941872037034, -82.63004517416, 19.2658983356222,
                    -380.075014800951, -65.8245205529665, -132.380237141115, -82.5567108971806,
                    19.2149155952069, -380.034374900876, -66.3980191923077, -132.537262369728,
                    -82.6513617198326, 19.8314831641682, -379.919889085478, -65.6998939531443,
                    -131.613468856041, -82.8389923679438, 20.3805776107045, -380.329284578916,
                    -66.4441335439505, -132.117976594902, -83.5278576581782, 19.7423018432639,
                    -380.780291410166, -67.124298900969, -131.539653704668, -83.0528246298891,
                    18.614269163689, -380.722857478884, -67.3131274237298, -132.287985569612,
                    -82.4046981083764, 18.3804472840114, -380.820890820141, -67.1480688572288,
                    -131.708374863723, -82.2196925381818, 18.3939823931431, -381.279849991764,
                    -67.6708917802159, -131.968446525506, -81.7060593435994, 18.521240830075,
                    -380.89988149194, -65.836097005632, -132.907359912127, -81.7358479979875,
                    21.407339240689, -379.457720263451, -65.7713094010804, -132.41878084056,
                    -81.4863814930496, 19.949223547624, -379.105297852774, -66.2839049131087,
                    -132.30847743888, -82.1423898155921, 19.8015208865529, -379.526880873468,
                    -66.1199801342712, -132.748571989045, -81.3659685155472, 19.2682119152725,
                    -380.127167639954, -65.859971659463, -133.020342138633, -81.1547094899581,
                    19.156569022675, -379.903190679781, -65.1890661121249, -132.072623705795,
                    -82.0227438550416, 18.9932992182367, -379.784600460923, -65.3120259335474,
                    -132.449793216616, -82.2922622852684, 18.7467962715778, -380.241018275154,
                    -65.9489528797688, -133.321583433862, -82.7345983998696, 19.5816127210072,
                    -379.979925999884, -65.5253893260436, -131.856503798804, -81.8350185400294,
                    19.653695074964, -379.620282675166, -65.9474839787633, -131.756529101015,
                    -81.959079430592, 19.657129877153, -379.322763423516, -66.5626756359745,
                    -132.603280041532, -82.8399948192813, 19.3310669023219, -380.121014623415,
                    -66.4278731778762, -132.697480269586, -82.3429796603879, 19.5649941531548,
                    -380.922757313944, -67.3973318206266, -132.045857006918, -82.2461066410157,
                    19.0971034225778, -380.7795061377, -67.1421847629561, -132.959058149658,
                    -82.9993205764984, 18.9150052722769, -381.769547675662, -66.4059952161801,
                    -133.47080547673, -81.940628696685, 20.1752540911826, -379.294682114424,
                    -66.4299373530627, -131.941827411045, -82.1497971159092, 19.9638693339821,
                    -379.363220199417, -65.803868671551, -132.599575931566, -81.7049234862595,
                    19.5397257966571, -379.996177805486, -65.5994442942466, -132.257210104199,
                    -81.6248226808921, 19.3515611508026, -380.245395788874, -65.3980638994467,
                    -132.171095858392, -81.324461074561, 19.424512624031, -379.113139637284,
                    -65.1947750179536, -133.219863459553, -81.5124604056215, 19.6280633792385,
                    -379.481637932685, -66.5041441737846, -133.421028964756, -82.0058260850508,
                    19.5721378962794, -380.930208892531, -65.9222995263907, -133.177174687117,
                    -82.2362937432109, 19.598129162873, -380.39155545032, -65.7539660307476,
                    -132.741703092581, -81.7692523254293, 19.6876216780357, -379.523976842391,
                    -66.7384385081534, -132.361013509983, -81.5441192728544, 19.8632532265858,
                    -380.00185426125, -66.0486018531511, -132.349102065764, -81.9494731445754,
                    19.3334647657841, -380.816164151337, -66.4434389413797, -132.910199692857,
                    -82.1304340314163, 19.4178074785828, -381.267811331318, -66.6156873491401,
                    -132.941814939897, -82.2651577612214, 19.1742112046019, -381.000614268321,
                    -66.1985795118594, -133.536963943496, -82.2669512752132, 20.0147682174226,
                    -381.813312258644, -64.8278743264624, -132.613487242594, -82.7416393090175,
                    19.473419077087, -379.303825158261, -65.6750324456138, -132.285030650005,
                    -81.813400978595, 20.3643198437234, -379.248009418564, -66.6453129233769,
                    -132.242459545685, -81.8721843342773, 19.5468477192287, -380.243761188476,
                    -65.9587064829364, -133.100718811766, -82.4113445428063, 19.8529832381438,
                    -379.987021997859, -65.8276594161259, -133.060574384091, -81.5449074699627,
                    19.0843729341814, -379.949770894203, -66.3412704669797, -133.310175197865,
                    -81.8787044110034, 19.1441800807098, -379.744016652084, -67.228485035035,
                    -133.440023015148, -81.5027890651238, 19.6297574070519, -380.600127327589,
                    -66.4770358008617, -133.641836445141, -81.2080480201034, 19.4332306070616,
                    -381.077096776414, -66.2375283013442, -132.30196004782, -81.7769912888836,
                    19.8878331769561, -380.537250975646, -66.5745334440057, -132.652538450129,
                    -82.177956441671, 19.3330934918156, -381.084680113761, -66.3698327812703,
                    -132.100849464109, -82.5896334258823, 19.4823656034557, -380.780988801122,
                    -66.2609757584939, -133.703074859928, -82.1530142704922, 18.8504579408502,
                    -380.605938204604, -66.6912973531879, -133.415227461074, -82.6032463917606,
                    19.8681085383855, -381.531625521799, -66.1814773903082, -133.662750188768,
                    -82.0339572272199, 19.8653832307748, -381.661432396412, -65.6490784641345,
                    -132.288178838581, -82.4599541433933, 20.0296379133268, -379.776128303751,
                    -65.6177198026924, -132.640204527223, -82.2054845291714, 19.9333644707704,
                    -379.281597940816, -65.4902155901642, -131.877535769718, -81.8874655153821,
                    19.8650649820188, -379.574468934244, -66.0008096926726, -132.701404970888,
                    -81.905843470735, 19.4053025899134, -379.935128180115, -65.7260277610004,
                    -132.089899029748, -81.3899871786933, 18.4760343138254, -380.283216109907,
                    -65.5890336745042, -132.201050879425, -81.6874668521098, 18.4438926030932,
                    -380.587867642879, -66.1579521640004, -132.501057538354, -81.9422800071342,
                    19.3350446517658, -380.237524359188, -66.8723269837005, -132.902846782961,
                    -81.3428541518959, 19.729113746533, -380.974740536346, -66.4816738357528,
                    -132.650571970728, -81.8180155758109, 19.5774070693933, -380.947036917177,
                    -66.341959425319, -133.014571648656, -82.0007061120978, 19.0542721493426,
                    -381.390477756426, -67.1994314547267, -132.56344979443, -81.8739181530027,
                    18.4837810766601, -381.637255519443, -66.4589621625191, -133.048323338255,
                    -82.3470228154387, 19.3126312111408, -380.736756681681, -66.5265646830769,
                    -133.722157540374, -82.609876045602, 18.7551604429939, -381.514319562434,
                    -66.1741616437411, -133.982237912153, -82.1650154609846, 19.6427877094475,
                    -381.963993357546, -66.3500873815926, -132.007706008563, -82.2827549957247,
                    19.8678101208787, -380.032162388621, -65.6871494735241, -133.069886807366,
                    -82.0329556978025, 19.7768397607313, -379.80310060373, -65.8391345894385,
                    -132.736237681211, -82.066368094126, 19.6948559820524, -380.183925381702,
                    -66.136705592665, -132.504107520988, -82.5836246863942, 20.1191047510299,
                    -380.639437963188, -65.4715780237787, -132.365597067214, -82.0818312170591,
                    18.8392942851166, -380.443985830389, -65.6571480926068, -132.067689560286,
                    -82.1783897331859, 18.7821651008084, -380.755552287473, -65.9192961396058,
                    -131.860473806809, -81.8605783897813, 19.3425194832386, -380.602464249134,
                    -66.1222133797698, -132.153102737533, -80.8771950546791, 19.7895615635421,
                    -380.818931950391, -66.0213905355554, -132.807756539424, -81.8466319369892,
                    19.4603098554749, -381.233388627435, -65.7835767048016, -133.342195464273,
                    -81.8532577066913, 19.0779845324929, -381.307365182347, -65.8644067855247,
                    -133.33718923199, -81.6332704494257, 18.793238055203, -381.226070853815,
                    -66.7748912742162, -132.881893810011, -81.992030599389, 18.863368912739,
                    -380.889193634131, -66.1638744370529, -133.42406638794, -82.3412655010897,
                    18.8367055214151, -381.299016566655, -66.1850189084969, -133.775284466302,
                    -82.8913951375134, 19.3692036729103, -381.397394011055, -66.2264338196865,
                    -132.511776794302, -82.0212758779192, 20.5262825325252, -381.313095305772,
                    -66.0676691577529, -132.486564685738, -82.1026933751608, 20.0759098230413,
                    -380.554051699289, -66.5654582260959, -132.447189231652, -82.2905698337624,
                    19.3276744688769, -380.938171846227, -66.0579319339412, -133.271482990971,
                    -82.251067767922, 19.7024887098266, -380.603910401184, -65.6670286813248,
                    -132.058230282238, -82.0485328728934, 19.6305122438859, -380.637896045432,
                    -66.0218377041978, -132.345870501019, -82.1400425447046, 19.3913146125617,
                    -381.047065815775, -65.95736438701, -132.072461680456, -82.5163902390617,
                    19.85488399728, -380.951338350335, -66.5872290987351, -133.084474832483,
                    -81.6472613844665, 19.8289330092485, -381.142524698465, -66.7864792331684,
                    -132.973760825565, -81.7051906259743, 19.4671613092978, -381.227517600973,
                    -66.3997010716032, -133.516705700999, -81.3771461813747, 19.5574633560373,
                    -381.645379331646, -66.6432518678396, -133.317948202828, -81.8631092729427,
                    18.4795873404545, -380.967224471299, -66.4346114993049, -132.923710517894,
                    -82.4344901112947, 19.0230406529467, -381.049087330541, -66.9576614988353,
                    -133.683576134099, -81.9504652200102, 19.0070638892433, -381.527333470107,
                    -66.4361170589338, -132.806335430213, -82.6313806758864, 19.6750873146379,
                    -382.051191233197, -66.0122860180657, -132.662370761959, -81.9480403559937,
                    20.2510281278273, -380.621337261195, -66.1896157146921, -132.519164655192,
                    -82.7098679619878, 20.3319424170074, -380.96716190463, -66.8251722143138,
                    -132.631840369975, -82.8488307549176, 19.4466889348711, -380.582759071729,
                    -65.8479562530063, -133.081078719865, -82.4202892977358, 19.2351311445767,
                    -381.328536798029, -65.528729964641, -132.497037482726, -83.0536491029278,
                    19.4966161395307, -381.210626806817, -65.799515219772, -132.816249448779,
                    -82.87335556896, 19.2562021186805, -381.117525794756, -66.4256930469674,
                    -132.664443715387, -82.8759328685872, 18.7394344611145, -380.569314892906,
                    -66.335475839367, -132.672135582481, -82.1286657487665, 18.9936308615013,
                    -380.467456008129, -66.4690025752553, -132.666940634947, -81.9499000865003,
                    18.526306212037, -381.3841318255, -66.3611757354151, -133.08049269184,
                    -81.890493875942, 18.8903936910787, -381.61601279833, -66.5249604968849,
                    -132.832332616368, -82.0208490765502, 19.0650278017944, -381.740555836251,
                    -66.619315693022, -132.96429332575, -82.600317697886, 19.1758323052918,
                    -381.753869393166, -66.5186062303762, -132.926667410465, -82.5841734393235,
                    19.0304582780512, -381.827156349371, -66.2825071320324, -132.805203383011,
                    -82.5932255337141, 19.7849985994987, -381.71412172744, -66.3099897279833,
                    -132.927875603134, -82.7987060586476, 18.8550657992683, -380.766939340768,
                    -66.1920255914018, -132.607996789233, -82.9983271761923, 18.9386913774976,
                    -380.468336793914, -65.6103507824756, -132.164695363328, -82.3455728587155,
                    19.5946119608305, -381.087744999535, -66.5958213737661, -132.623377400292,
                    -83.028907050001, 19.4697879578797, -380.288183846497, -66.3581926290083,
                    -132.752420701739, -82.8133276102771, 19.5651683917856, -381.309535023446,
                    -65.8194500194001, -133.041971099143, -83.0287011835432, 19.5007300473075,
                    -381.022090136707, -66.0497464597091, -133.154438020395, -82.5876159827057,
                    19.5284088308081, -381.411125249213, -66.2334959597712, -132.603746650735,
                    -82.2588959520652, 18.4534614793384, -381.250656629135, -65.7267231159311,
                    -132.922233181206, -82.4403790327815, 18.9203296072425, -381.467590082991,
                    -66.0759783779837, -132.893176546462, -82.332864635005, 19.1011867653148,
                    -381.863489235828, -66.5735885465328, -132.530819735661, -83.107143768134,
                    19.6642265865804, -380.71818699739, -65.8832889535286, -132.323214426807,
                    -82.4413961055603, 18.9190479194201, -381.449024610422, -66.2988733143792,
                    -132.887085050178, -82.01060658455, 19.2008730470359, -381.925361334289,
                    -65.50215282087, -132.548299895827, -82.0718688131589, 19.2071612432254,
                    -381.044216246218, -65.7546057070637, -132.044157713314, -82.985393129818,
                    18.8004167801283, -379.956108115869, -65.1885980187607, -131.90457459863,
                    -82.9771392318849, 19.7241941428136, -380.582118552101, -65.3967998103729,
                    -132.470621969339, -82.4031374629877, 19.7364164527255, -380.792579808565,
                    -65.4410960795612, -132.461641433587, -82.0056515303486, 19.3556194411952,
                    -380.332317436955, -65.96421759819, -132.548266509328, -82.983149577819,
                    18.7602161843045, -381.080841545258, -65.8182344100582, -133.096589439061,
                    -83.1854697971797, 18.8358792983501, -381.478899049619, -65.8024155768475,
                    -132.425207543394, -82.7359874654743, 18.8716059192002, -380.872969957056,
                    -66.1622803376563, -132.780258546847, -82.9572709678734, 18.8522089768715,
                    -381.999651221335, -66.1950283854871, -132.860917824542, -81.8736036904068,
                    18.5761440104347, -382.075187243925, -66.6912708860711, -133.10999697673,
                    -81.6406692155662, 18.9742209679265, -381.301137979265, -66.56604111797,
                    -132.231206423437, -82.8846752562088, 20.2094270988705, -381.488973783857,
                    -66.3271460498179, -132.551397093904, -82.2211552151397, 19.7449133578117,
                    -381.615489131704, -66.1382189309666, -132.377490015497, -81.9650385008984,
                    19.2806166218043, -381.412162987297, -66.2603502106959, -132.669700067453,
                    -81.9859144777386, 19.3849074708351, -381.127691111336, -65.9870561735533,
                    -131.950149745379, -82.7693886730415, 18.5623705421995, -379.810102981466,
                    -65.5226824622527, -131.767971869722, -82.7238927278092, 19.1647713247451,
                    -380.172250622394, -65.530484015311, -131.633295021058, -82.7502381289401,
                    19.985096296218, -381.537471144802, -65.1977682298028, -131.67330048182,
                    -82.5703819872817, 19.1743742335044, -380.885091720591, -65.6462640800884,
                    -132.421338430912, -82.8730628491436, 18.1012916927332, -381.060007406444,
                    -64.629307204558, -131.900186777381, -84.1256577435145, 18.4780704538172,
                    -381.138791765205, -65.5961488727808, -132.204139051117, -83.0681846016583,
                    19.3029593604472, -380.930240066225, -65.9924260380845, -132.579443950369,
                    -82.8760866305988, 18.8028517807194, -381.374610607246, -65.9649658755645,
                    -133.13240471499, -81.8660830820966, 18.2156752860896, -382.124082362735,
                    -65.8390765588748, -133.191547896599, -82.2942885177675, 18.8108402747453,
                    -381.635893074002, -65.8872171485627, -132.252305640683, -81.9469334923309,
                    19.9324244157831, -382.686837504776, -65.5350634003839, -132.854804812008,
                    -81.8209209657005, 19.9383637641604, -381.33199996752, -65.5366364387338,
                    -132.829261806379, -82.3971869409698, 20.0622998429979, -380.031374561301,
                    -65.8864308088384, -132.767424946856, -81.0792516021539, 18.4168160114924,
                    -380.378760919719, 49.3028008735474, 98.2831569924198, 60.7136377754325,
                    -15.2986444522901, 279.067654263811, 48.4088096679486, 97.684431147785,
                    60.5880588072864, -15.3642275649942, 278.922008685268, 48.1822045277791,
                    98.6421542821381, 60.8783940447515, -15.1358157867625, 278.695115011595,
                    48.256564945364, 97.9845053017217, 61.1412152064364, -15.6828674606716,
                    278.651600412158, 48.5539250242178, 97.7693698283139, 61.1199460424027,
                    -14.3594701026283, 278.420677306826, 48.2716004407096, 98.2076497277083,
                    60.7662491853016, -13.7052419242282, 279.751775973467, 48.1139069279855,
                    98.2508687431069, 59.9567302480981, -14.1710107216006, 279.626933836959,
                    48.5596197716669, 97.8474337728921, 60.1687943898807, -13.9258024630201,
                    278.949406804217, 48.1625170734771, 97.8735973571606, 59.935073816915,
                    -13.7402327557122, 280.529748778334, 48.6915728142271, 97.7039676567692,
                    59.9725047447554, -14.3532724434975, 280.961958308393, 48.9893889391641,
                    97.6849230065188, 59.4684634097451, -14.5409616629179, 279.229125249728,
                    48.9476639410266, 98.6761115614305, 59.1392252713454, -14.8817448308068,
                    278.904674431402, 48.3403702427108, 97.9140836564932, 58.5537452276433,
                    -15.0537571022589, 279.257609233479, 48.7667120852905, 98.5671401496776,
                    59.0412051800476, -15.104289223228, 278.22762140563, 49.5692460072504,
                    97.8227609276182, 59.0547222871638, -15.1879042565551, 277.834518614067,
                    49.8516009143782, 98.0335059398787, 58.3934138364672, -16.4732835615162,
                    278.723167709046, 49.1398665129579, 97.7611093615146, 58.9177023297015,
                    -16.3348007794866, 277.252737203567, 48.1570322152172, 97.564183332559,
                    58.7398947477139, -15.935814103852, 278.294117013405, 48.4489729811681,
                    98.8860344907204, 60.239990894151, -15.2360965203358, 278.48934311936,
                    47.8450347696696, 98.1699287161337, 61.0121194552868, -15.4246856350875,
                    278.919607669456, 47.8068093212106, 98.4142497406312, 61.6422459019233,
                    -15.6323778718682, 278.60659728966, 48.3056599759211, 97.7805162470627,
                    60.8961100017111, -15.1082029176727, 278.346940209615, 49.1008937752129,
                    97.5755028888485, 60.6401625302417, -14.1447924102401, 278.496627421852,
                    48.6769194784168, 97.8838823786397, 60.3996916017933, -13.3087094231011,
                    279.461137190393, 48.1479129101346, 97.7801774293946, 59.4065207773641,
                    -13.9100329179259, 279.363644258185, 48.4308464996356, 98.4639846033722,
                    60.4124020906968, -13.1808613313174, 279.500169999552, 48.4235123825108,
                    97.7903788831933, 60.1659295679604, -13.9151702878544, 280.324315396664,
                    48.1956864841488, 97.5893148078219, 59.8178779357455, -14.4927505621841,
                    280.745040366454, 49.2151401808909, 97.4962195944655, 59.2851254299583,
                    -14.2224715168544, 279.735931659633, 49.6899441340138, 98.3726560382711,
                    59.0118766695063, -14.753912119789, 279.33285424761, 49.3499499778824,
                    98.5198231620917, 59.1349384975347, -14.9055201346394, 278.528812313986,
                    49.248713590253, 98.6142816806982, 59.7677983090355, -14.9755451386662,
                    279.296429291781, 49.1762225699947, 98.6735605357129, 59.7066965531679,
                    -15.5110479564517, 278.495865001843, 49.9854039462944, 98.4299276327906,
                    57.9190670582748, -16.2019029801602, 277.891006403058, 49.5237245247224,
                    97.5146209880072, 57.8466777562984, -15.190372355331, 277.568864116663,
                    48.7083522193392, 96.9437893826555, 58.0655122906115, -15.2381518617733,
                    278.265701659628, 48.6232793047363, 98.6512613751814, 60.3381763949894,
                    -14.9153076677568, 278.825056010902, 48.1637445047237, 98.5944775923875,
                    61.4818459781196, -15.6745060004487, 279.138306242321, 48.1985633128069,
                    98.3112998561894, 61.5807528183939, -15.4757337467208, 278.299309103168,
                    48.1688719567705, 98.7990876676039, 61.4234108858098, -15.3977391396203,
                    279.217730847577, 49.271919310867, 97.5600954072673, 60.8697119532662,
                    -14.4372898769797, 280.043709432395, 48.6434052269884, 97.2166991776499,
                    60.8571054085654, -14.2608829904764, 280.323570998795, 48.7246141271962,
                    97.4947514297053, 60.6562908819888, -14.1147615756466, 279.795058522921,
                    47.9040183095278, 97.9307168509498, 59.952814563243, -14.2005227074787,
                    280.242142414046, 48.3025576161348, 97.46453701493, 60.3381395759146,
                    -14.5106067811261, 280.377648386728, 48.4644491851326, 97.5468006027509,
                    59.8976927675351, -14.5621902907204, 281.120551001332, 48.2887581776524,
                    97.2334404138962, 59.5277602195243, -14.8288932648063, 279.674154089115,
                    48.6463962999045, 97.8377543120298, 59.5972504162489, -14.3101543275938,
                    279.443882485115, 48.4103421058502, 98.4311240425387, 59.5819272010349,
                    -14.2568467682125, 279.441463988412, 48.9697860164047, 97.9059344949847,
                    59.3759845651419, -13.7515551414069, 279.814222617059, 48.7291612696071,
                    97.8487209858474, 59.9807711762262, -15.2202351661952, 278.427455338055,
                    48.6587509794865, 98.2900057558964, 58.7442766170572, -15.4386709444963,
                    278.254420609128, 48.9427462397121, 97.5217131234522, 58.1724052673128,
                    -14.4126972667183, 278.020194111517, 48.4101353089549, 97.4386784134001,
                    58.3365028288976, -14.0450965929556, 276.811282545209, 47.9216188534202,
                    97.9213669045246, 61.0540194606595, -15.5591873527011, 278.389631626757,
                    47.9654788336618, 98.4057321626429, 60.808506206456, -15.1264626523874,
                    278.375189876736, 48.09570425388, 97.8182626025007, 61.3372046202423,
                    -15.5287220071531, 279.026713991755, 48.4284703715427, 97.9297593673671,
                    61.080540184268, -15.569867171715, 278.597335990416, 48.8679933388718,
                    97.1483081667327, 60.7722205108568, -15.3235323190126, 278.509368474647,
                    48.7465760222957, 96.9305581309927, 60.7941694048863, -14.8049919388036,
                    279.170376199335, 48.2584845173995, 96.5964787043145, 61.1707587224243,
                    -14.3378943622611, 280.147749817254, 47.5244305616132, 97.5263234063039,
                    59.8913087889177, -14.1641923943627, 280.900423516486, 48.255263006235,
                    97.4580723828364, 60.0334944418271, -14.4685364191985, 280.415523316491,
                    48.4053432946502, 97.4862029353527, 59.6658085183036, -14.8639981167832,
                    280.786760514119, 48.6292951019125, 96.9942476793213, 59.9379142754495,
                    -14.8210292635893, 279.569205935335, 48.490558061476, 97.6855993530924,
                    59.4450006118029, -14.8328638789143, 279.96240404507, 48.1227150223047,
                    98.0046807189551, 59.6012099076582, -14.3742201558673, 278.983976338782,
                    49.1611078110976, 97.1344215331171, 60.2610762632205, -14.4305454298676,
                    278.855015940463, 49.2117842451992, 97.5539857922151, 59.6284464532399,
                    -14.9125495638476, 278.281314696263, 48.0203271783461, 98.0243325619898,
                    58.7450194815865, -14.3799626974192, 278.608291741533, 48.6395496217159,
                    98.1545763674514, 57.4212166296682, -14.1228066438825, 278.269220995529,
                    48.313064704311, 98.5579019496846, 57.7416186531504, -14.1602344182264,
                    277.608161101123, 48.0926747267398, 97.4606986306699, 59.9851822634933,
                    -15.5574613608747, 278.808279768173, 48.4104969667232, 97.5545776020882,
                    59.8624298482433, -14.9411462524032, 278.054407427876, 48.0569934874961,
                    97.9459536151915, 60.6801821502089, -15.6626054464874, 278.914292419631,
                    48.565827420883, 97.8805336983869, 60.5153443885988, -15.683409372791,
                    279.076905470548, 48.8655487631577, 97.695800592966, 60.559560251452,
                    -15.446508447223, 279.173784642811, 48.6456199657554, 97.7886861721463,
                    61.0757419211678, -15.1643548594808, 279.17232233408, 48.4828673865016,
                    96.7743093083303, 60.7636760949494, -14.2090825356259, 279.62378424214,
                    48.1411723823686, 96.9228781856068, 60.109336224617, -14.2070743466851,
                    279.850342401705, 48.0227318088698, 96.9832564370281, 60.1809789658444,
                    -14.4059663327105, 280.371340185303, 47.933773663183, 96.7094250283695,
                    59.7363793212387, -15.1809694954031, 280.336679054766, 48.5585498638513,
                    97.1705409840714, 59.7682612180812, -13.6021660705767, 279.129665673207,
                    47.9596635229339, 97.0150871460811, 59.5089859540878, -13.9860053273124,
                    279.831362647917, 48.5996603198735, 97.0815021629039, 59.2828201270886,
                    -14.5131210977778, 279.827790390181, 49.410017453331, 96.9577377018142,
                    59.8014796756074, -14.5641279137692, 279.395440393066, 49.0801516841827,
                    96.9709520533435, 59.8960360312753, -14.0396951226752, 278.707935264295,
                    48.3148651048234, 97.0102692997759, 59.849765921975, -14.4519856523061,
                    278.504616603112, 49.0189070530201, 97.2949794404534, 59.2672862019686,
                    -14.1709440436942, 277.896976850198, 48.8220093398706, 98.0362824743411,
                    59.1181043148762, -14.1769811963869, 276.714751250799, 47.9776955377991,
                    97.846751373194, 59.6607460151567, -13.9100154466252, 277.844988495041,
                    47.9756439626157, 97.6734791665651, 60.1750219217112, -14.5478658134923,
                    277.740853895223, 47.9437303620599, 97.671962403655, 60.0341228528923,
                    -14.8421525332255, 278.694239806202, 48.4500982245158, 97.5357731128485,
                    60.0486978763852, -14.8506980231448, 278.500406672235, 49.0448055595228,
                    98.2357909449728, 60.1031728794927, -14.7975709623774, 279.634810237248,
                    49.1028257021275, 98.4766441165896, 60.2412352259548, -13.8314994221989,
                    278.741200644336, 48.5814570655667, 97.7210554710506, 60.1509041144823,
                    -14.2154483892319, 279.221365959348, 48.3082516043465, 98.0643050179157,
                    59.5792102263618, -14.3122383834911, 280.222816729696, 48.4366634927194,
                    97.5657053037181, 59.6718784216831, -14.6291022006214, 279.81006080904,
                    48.7861144956883, 96.9026132251528, 59.7856166596317, -14.2316323185227,
                    279.158148334173, 49.0838348671631, 97.17061720831, 60.4042675554738,
                    -14.627346787134, 279.658836171266, 48.828476197397, 96.6888102528528,
                    60.5261515788747, -15.0913618557352, 279.681679791186, 49.045830019137,
                    97.007206003599, 59.7337293673847, -15.073794000572, 279.805678345021,
                    49.3395028393367, 97.3333296401596, 60.1337269481269, -14.8183287479182,
                    279.284164530527, 48.8708328648053, 97.5535939889652, 59.7629703119637,
                    -14.3110852866884, 279.594716792846, 49.0642194996671, 97.7445745559083,
                    59.461432236034, -15.1112740028622, 278.601895587802, 48.6345208958563,
                    96.5115813309105, 59.6399277991965, -14.3710285962932, 277.621932356268,
                    48.7869420108598, 96.5502110988607, 59.6914977816206, -14.2753068980698,
                    277.634941576865, 47.8511536475648, 97.7207000101357, 60.1263424284041,
                    -15.3322078749861, 277.987929435182, 48.2765675214067, 97.9119061653952,
                    60.0468673741302, -14.3852181014412, 278.621252360381, 47.8927017332889,
                    98.0752460062335, 60.0547548732477, -14.5761880526453, 278.262823262437,
                    48.1412527349139, 98.3458735184761, 59.9785749200618, -14.2765625475226,
                    279.012694239586, 48.2377111896222, 97.8686571576334, 60.5331906335877,
                    -13.9276302944821, 279.338675578553, 49.1052327115149, 97.9076187990807,
                    60.3213757023287, -13.2435818231935, 279.017171443466, 48.6259280018017,
                    97.4612588776874, 60.0572447060036, -13.9471924283908, 279.004418295797,
                    48.6614428483529, 97.2506909106654, 60.237261516644, -14.6208509564116,
                    280.314117285137, 48.602248339548, 98.0829713902568, 59.8188978291468,
                    -14.156029461388, 280.22712620967, 48.9528766984942, 97.3035596099359,
                    59.2786060421949, -13.3598067186039, 279.875132191874, 48.4595986413267,
                    97.3315469491606, 60.1778217882289, -14.1904591467292, 279.811734228108,
                    49.1290810571486, 97.4042931170373, 60.0232008422757, -13.99091287386,
                    279.642301244579, 48.9545112568068, 97.4503787861159, 59.8374151822767,
                    -14.1945609455406, 280.003116851042, 49.0459167216543, 97.4147173431096,
                    59.9920259941873, -15.0539261785835, 280.106162265282, 49.1196103638113,
                    97.5629620950074, 59.7281710769744, -15.1880710444247, 279.26961518243,
                    48.9277435089751, 97.0399050760508, 59.4474950319381, -14.5434539510404,
                    278.742015942859, 48.5019674854533, 97.1882406615935, 59.6331953959122,
                    -16.0135375674004, 279.115669461088, 48.5042701321072, 96.5108259690773,
                    59.0298523817462, -14.5701503254545, 279.480127389151, 48.5867077490783,
                    97.9751470684763, 59.556339257617, -14.6416295766927, 279.320616091134,
                    49.1655870573359, 97.8656652712667, 59.3290518785331, -13.8179952299152,
                    278.956883242828, 48.8446538548758, 97.9385862894036, 59.7258996717183,
                    -13.5247513558659, 278.818232809417, 49.5493026508715, 98.2791835200918,
                    59.7866669093042, -14.0225387721478, 278.813328688955, 48.0770354942038,
                    97.7937068957575, 60.1011954434296, -13.6642993017091, 278.772591099256,
                    49.0113282140847, 98.1254155603362, 59.6653001432416, -14.0888954848188,
                    279.289052120892, 49.103658124896, 97.9024592934548, 59.830962257612,
                    -14.1041053337412, 279.634781489199, 47.8235358613305, 97.5654210127966,
                    59.8115155936653, -13.5757232371707, 278.663882883029, 48.4716335556639,
                    97.7751956043502, 60.4172538913412, -13.7908255716516, 279.943001406174,
                    49.1493311021328, 98.1277994145833, 59.8666183457312, -15.3744405889241,
                    280.082964922352, 48.5213527568592, 98.0559013929092, 60.2812119962879,
                    -14.7247900980992, 279.666492258127, 49.1644524602334, 98.1843130917019,
                    60.1262735580337, -14.0735561161115, 279.526829717878, 49.0245056124365,
                    97.6399136929191, 60.1739797280763, -14.8664893394275, 279.852709058128,
                    49.7692681199617, 97.4710860999943, 60.1578571099032, -15.4727999213681,
                    279.524920470935, 49.2127085676517, 97.8498789310205, 59.9284175402599,
                    -15.482958334686, 279.178463700761, 48.8509605275286, 97.3132829691604,
                    59.8621167188605, -15.4356888563022, 279.49995083115, 48.4445826625353,
                    96.9937685705276, 59.2792961924308, -15.2759232265824, 278.974621394222,
                    48.6846676245178, 97.1173777606937, 59.8272071890921, -15.1264198406348,
                    278.012476533665, 48.8855769571129, 98.4715914708794, 58.7241508550642,
                    -14.2534806613908, 279.54384389434, 48.9821001059548, 98.8616492319778,
                    59.9015017824024, -13.6254041975233, 279.485712425088, 49.1987493824879,
                    97.4072925112223, 59.7479071261097, -13.5424821550897, 278.141996257446,
                    48.7875556215625, 97.8270454938237, 59.4356519700523, -13.4649521235864,
                    278.567225177256, 49.2534608747114, 97.9712414978681, 59.5572493784504,
                    -14.0100067431848, 278.529887635131, 48.4827275482262, 97.0917385777873,
                    59.6882748442452, -14.2006918642001, 279.231118334172, 48.0673982668962,
                    97.7393714762762, 59.7296743501589, -13.7929471893809, 280.00261831538,
                    48.0336421131373, 97.5715090465627, 59.6273580369146, -14.2324593104863,
                    279.755987545544, 48.3245352729202, 98.1809982740895, 60.0997960647958,
                    -14.0901558971227, 279.845571964633, 49.1052295404216, 98.1537397772888,
                    60.1425448730927, -13.5690638294622, 279.792317843521, 48.4329236432458,
                    97.4295895650463, 59.76541382611, -14.7013843992053, 279.101347194496,
                    48.5757675486614, 98.4980360192451, 59.3793941015085, -14.641498246589,
                    280.650993674405, 49.0758601610437, 97.8719782167563, 59.7954332249531,
                    -14.5369285599682, 279.776109688759, 49.2293833852416, 97.4699193119291,
                    60.1586495172839, -15.969672295411, 280.143840946011, 48.5306364873722,
                    96.7140772274842, 60.2575581234072, -15.836727940833, 279.72273584618,
                    47.6891224744451, 96.8290701279327, 59.4784684421358, -14.5193979892544,
                    279.201425892261, 48.2582456496533, 97.1145003575432, 59.4328945184944,
                    -15.2615978803923, 278.551073643634, 48.2453603536073, 96.840889083604,
                    59.7144826897883, -14.8404954140905, 278.525386333614, 48.046589694204,
                    97.1958626389487, 58.960347832169, -14.0843600786081, 278.852178151292,
                    48.8077648055057, 97.4892412012363, 59.7514134530993, -14.3458822150119,
                    278.423486178195, 49.6964277916147, 97.8687624116452, 60.0118376584365,
                    -13.2641220608804, 277.889826793666, 48.5747503843614, 97.6953406677784,
                    59.7662978934864, -13.7720210290439, 278.6376340712, 48.0943010596662,
                    97.5209562767469, 59.0567829349593, -13.4849055830123, 278.354829116614,
                    48.2296776030847, 96.9973895212768, 59.36831173379, -13.8740899847199,
                    278.434995832946, 47.8069080213637, 96.8703297769241, 59.4167256173779,
                    -14.3566041953283, 279.995252444599, 48.5126741220586, 97.4693053275464,
                    59.8256234405549, -15.5266662713007, 280.08022132938, 49.152586913742,
                    97.6215740805715, 59.4383800587825, -14.1058959526654, 279.280769781798,
                    48.5519563477806, 97.7250627014984, 59.6972545015412, -13.502241290595,
                    280.530502288871, 48.6513828892331, 97.2937677340927, 60.2607737488322,
                    -13.8633497233733, 279.897562248308, 49.1197733774639, 97.7592830010918,
                    59.2065034847155, -14.350692698943, 279.790871467915, 48.8574595870132,
                    97.5612338750014, 59.7757865111078, -14.3524558308179, 279.780078668908,
                    48.965640123362, 96.5346025311563, 59.6600481822799, -14.4719013145298,
                    280.133290743178, 48.1481878361008, 96.6910688229572, 60.089745699849,
                    -15.2399113833263, 279.662786913129, 47.1234724861387, 96.6372689646273,
                    59.719631812805, -14.8882337038914, 279.164812309041, 47.9190006158293,
                    97.1567894452263, 59.6973105408096, -15.2718571873526, 278.690286445078,
                    48.3381740655053, 96.949586346049, 59.395896131019, -14.305432869354,
                    279.832345947411, 48.5436349348749, 97.4202033388967, 59.4052802285069,
                    -14.6232923268946, 279.007184536379, 48.5582902574623, 97.416254976866,
                    60.2350972114622, -14.1260048320956, 278.507190464641, 48.4948631949238,
                    97.6414630870764, 60.0991538353101, -13.742231145687, 278.369643798732,
                    48.2063867691782, 97.7320265852312, 60.440407803598, -13.9673741451282,
                    278.326561448323, 47.8699319932666, 97.7006883906716, 59.1417705053674,
                    -14.1486605688568, 278.586411713708, 48.0940475811022, 97.3804575603496,
                    59.8058165833039, -14.7385990380862, 278.347856190575, 48.2395365435353,
                    96.4991698798531, 59.5090982136575, -14.6503197312303, 280.001149963005,
                    49.1871623535983, 97.4972301602228, 59.7230097900366, -15.503379905484,
                    280.104075465508, 49.0184283723649, 97.5038343037346, 58.8245889541643,
                    -14.1771273770923, 280.19057290168, 49.12060199119, 98.133914241634,
                    59.8427651658107, -13.6362982798759, 280.107507351522, 48.3043036175069,
                    97.487085521172, 59.6724735927962, -13.6073472298006, 279.417896728455,
                    48.803179692952, 97.2347617883145, 59.3825069887675, -14.1907210910545,
                    279.522070495147, 48.8068611451819, 97.9441183932129, 59.713862206163,
                    -14.3753457923441, 279.041672177792, 48.6903834704043, 97.6631840961595,
                    59.832082796904, -15.1427550837608, 279.818456668846, 48.9162517429495,
                    97.8901738367445, 59.8735506780372, -14.9602278688963, 279.743622888784,
                    47.7280993249542, 96.8687551739816, 59.7370438097091, -15.0107267091134,
                    279.737308958817, 48.2320207556445, 97.2415589806331, 58.8321420277244,
                    -14.9421049551474, 279.53981220288, 48.7505704781962, 97.3130173507157,
                    59.8805518388089, -13.9433954425485, 280.161414864922, 47.5523114064447,
                    97.1705437582971, 60.1003876319794, -13.7690278471365, 279.562124523005,
                    48.8167451626367, 97.2032876586476, 59.5996970323527, -13.8291177985402,
                    278.535858204412, 48.4103203292248, 97.3062694011965, 59.9204937946395,
                    -14.885966870417, 278.634269590318, 48.4837330755607, 97.4328151027867,
                    60.3618868290766, -14.0707262632024, 278.640050101032, 48.3133062257185,
                    97.0465992801544, 59.4024992891478, -14.8682925600219, 278.847855189406,
                    48.4074643226129, 96.6688934252667, 59.6784578083765, -15.1421043232289,
                    278.002922648523, 48.6121045362569, 97.0124107173384, 59.9260508433478,
                    -15.1894044285544, 278.829063882812, 48.6743359628293, 97.8889986501278,
                    59.9631690731773, -14.3805039896819, 278.489495705417, 49.0687444104009,
                    97.8960179357801, 59.5506743389235, -14.2231150237392, 279.810972033628,
                    48.4679718085972, 98.3191510141501, 59.806520625785, -13.7861980614024,
                    280.635711226507, 48.2452894814097, 97.658697959148, 60.1844873899197,
                    -14.1162785483844, 280.246996583016, 48.5030607014267, 97.5816997268662,
                    59.7629461931932, -13.742997794949, 280.019852562373, 47.8107420252643,
                    97.3368435201636, 60.1573931665956, -14.3283276173152, 279.296873589884,
                    48.3472212454411, 97.5210724086108, 59.9467382064826, -14.4100387915736,
                    279.709351137992, 49.3979597073869, 96.837394866266, 59.2997625927818,
                    -14.5875915513913, 279.472874705632, 49.3968417729186, 96.9403299643282,
                    59.7665132443134, -14.3911837883224, 279.311136355912, 48.7401760817514,
                    97.4944534785759, 59.5872171889245, -13.774944160236, 280.778349190634,
                    49.3028577374834, 96.557169920354, 59.8953104315157, -14.0737865285237,
                    280.154177434912, 47.9640406213591, 97.608418647434, 59.529970846659,
                    -13.7482255298683, 278.424986851695, 49.0619001642088, 97.5813689079045,
                    59.6668938873598, -14.1128515270485, 278.895709435945, 48.16167062025,
                    97.2606822954838, 60.3029137067826, -14.870027113075, 278.509190044085,
                    48.1621754544847, 96.8838771814341, 60.3001321506294, -15.095894217802,
                    277.896529938782, 48.2633987768064, 96.825133828992, 59.4375207483803,
                    -15.1113201578075, 278.624500871625, 48.5820087628209, 97.282262401719,
                    59.9215624499024, -14.6816132439795, 279.014288468581, 48.484661396108,
                    97.351899626969, 59.9684937212877, -14.8290286516003, 278.556312751995,
                    49.5083754764168, 98.1260519392207, 59.7543010943062, -14.6926208192525,
                    278.667580663654, 48.7781102529604, 98.0984659836996, 59.4030590588453,
                    -13.8897788930045, 279.458868014587, 49.0309129434192, 97.7449571977136,
                    59.9876930294756, -13.4627908376604, 279.360333754005, 48.4723539815964,
                    97.900710245738, 60.308171861163, -14.3284986462856, 279.155409934611,
                    48.8949845500324, 97.5774101577909, 60.0454321254852, -14.9579087654596,
                    280.209993430342, 49.0019886593323, 97.2234082422252, 60.5542613704129,
                    -14.8163318256184, 280.091201636446, 48.4396051932132, 97.2631727680428,
                    59.8962799450637, -14.2842750110796, 280.332065653405, 48.9082310353263,
                    96.6111165397301, 60.2425247064802, -13.3753529293978, 279.862470590553,
                    48.7595825030567, 97.0007971199989, 60.7189246925772, -13.2587280656824,
                    279.882156057352, 49.2152393414327, 97.0495513592932, 60.112777800639,
                    -13.6056885937196, 279.924304649299, 48.9285789847335, 96.677186261381,
                    59.7881101776016, -13.4297471888597, 279.101543664522, 49.0255756322736,
                    96.9543375111851, 59.890819056085, -15.4220433850501, 278.409172074697,
                    48.9499395095535, 97.2030350878234, 58.8009825413863, -15.1957432776165,
                    278.444828807207, 48.6807658619257, 96.9727352028413, 59.2299269608772,
                    -15.7548614507631, 278.314772499013, 48.3190931595406, 96.69342605095,
                    59.9787510998721, -14.7106610887439, 277.85741389624, 48.7133655761082,
                    97.0027778988681, 59.8488034803994, -14.2254959378706, 278.567559524997,
                    48.1264332892853, 97.4690314186817, 59.5484752946036, -14.8835965921661,
                    279.119635319568, 47.9830003987974, 97.3016519499355, 59.279733005908,
                    -14.861147793335, 278.637814945225, 49.1464781270611, 96.5386201525334,
                    59.9224642765663, -15.0357141679377, 279.128920692825, 49.9056425101201,
                    97.2734606559821, 59.9888480891183, -14.0746092751251, 279.844706970805,
                    48.8376075226499, 96.9835027451299, 59.9993075889473, -14.0017373516679,
                    279.942227035671, 48.7386920488409, 96.5665926764391, 60.5626970272253,
                    -14.237156133873, 279.019233249266, 49.1856210183298, 97.1782282866581,
                    59.7654145680149, -14.3874785033833, 280.502587410258, 48.6645013626141,
                    96.8486676275348, 60.1443545419114, -14.6718465916035, 279.649285836033,
                    48.6715563238932, 97.416265974181, 60.3547974638254, -13.6931950307669,
                    279.607547014321, 49.22776386279, 97.5807821693952, 61.0422014094497,
                    -13.4323012631888, 280.479579757574, 48.7551172003622, 96.667648939858,
                    60.8197905134455, -13.2707222997222, 279.480850648299, 48.6193704515822,
                    96.7424888911933, 60.374314071677, -14.0799133818897, 280.35114291828,
                    47.9868548542796, 97.5317908270252, 59.8557532540409, -13.4010192687618,
                    279.868298793018, 49.4178982119053, 96.3658843065426, 59.666848676814,
                    -15.0534442140427, 277.654002553628, 49.2823799373461, 96.1763436034672,
                    59.5975301514112, -14.5904911328056, 278.016186878779, 49.004027919348,
                    96.6367046042564, 59.0196516139009, -15.1052983754891, 277.586817682999,
                    48.7819171658221, 96.1338484172483, 59.5689506960717, -14.8323809669417,
                    277.223128461708, 48.5586819719752, 96.4960554804752, 59.0679270342836,
                    -15.047723498979, 278.369593225389, 48.5562448749988, 97.4111738794239,
                    59.3716858638047, -15.1937039753178, 278.347589245316, 49.0780097929567,
                    97.4354378896155, 59.4558462985573, -14.1407047924192, 278.249653505127,
                    48.930547980398, 97.3596448089542, 59.9703061094095, -15.2134900494677,
                    278.768620022705, 49.8032865604279, 97.3968032560939, 59.9544847476936,
                    -14.3229947952463, 278.92219182074, 49.6126007302502, 97.1106955028247,
                    60.0739508431389, -14.1722283428636, 278.560769969158, 49.197777036373,
                    97.422424620717, 59.7996468942215, -14.4962361765267, 279.700626184725,
                    49.4305060469248, 98.0224796210579, 59.5908540282876, -14.3364890390941,
                    280.664187711263, 49.3118603203349, 96.7853960744986, 60.4733478168296,
                    -13.6927198944052, 279.648936301041, 48.7196037643002, 97.5229047779254,
                    59.7735953183164, -14.3327419934457, 281.346431608506, 49.1879790224308,
                    98.4768659150732, 60.3717806689484, -15.3912451303533, 280.064288797803,
                    48.5980478456115, 96.5086334610583, 60.3264673982328, -13.4467587453207,
                    279.709577280154, 48.6234133850637, 96.3944703099594, 59.6478985790904,
                    -13.4573341252972, 280.170904037853, 47.8895623944643, 96.8371020139791,
                    58.9584495505765, -12.975017724177, 279.956630052618, 48.8613430132395,
                    96.2606144673752, 59.6470442308376, -14.2182044426494, 277.816039972417,
                    49.5753379388893, 96.5619789115284, 59.810856349152, -14.3876445368705,
                    278.195867459504, 48.6042603996025, 97.2244393473607, 59.5628326352742,
                    -15.1284836805655, 277.422453396082, 48.930555081529, 97.0536705327795,
                    59.4578111129358, -15.7664321672432, 277.124097826873, 48.7908028030992,
                    97.3402425348105, 59.5387706822792, -14.9287245967733, 278.038059633559,
                    48.7966969214223, 97.6011483464514, 59.5425658765957, -15.1016062418095,
                    279.086257389469, 48.8666310241664, 97.5037097721768, 59.8999550162942,
                    -15.0958065503186, 278.140465619617, 49.2721473116522, 97.2509837461926,
                    60.4520479834141, -14.7336801090022, 278.996398076382, 49.25693589893,
                    98.0921138332367, 60.3097960653703, -14.9936805414867, 278.790277169177,
                    49.2906266654099, 97.05897029274, 61.0010917531496, -14.3043473097685,
                    278.909102961761, 48.9254420917753, 97.6914686645694, 59.9634672275646,
                    -13.9997785415992, 279.451343077727, 49.3941389074086, 97.3000613375621,
                    60.2047748449979, -14.5718931315842, 279.045129363443, 48.8573328083539,
                    96.889858471022, 60.3933864362151, -13.6551418220951, 279.920905065007,
                    49.4454382169863, 96.6788831341715, 60.1972268027725, -13.8695135885362,
                    280.782663308365, 48.9843813412258, 97.3649642447697, 60.6494980471045,
                    -14.3280442079325, 280.646508251991, 48.1498554426189, 97.093196054041,
                    59.9954863347766, -13.6292603681989, 280.274396887673, 47.5378839251762,
                    96.8679758787299, 59.9034630976028, -13.0631534480374, 280.647593967785,
                    48.4903986933697, 97.2752521857302, 59.540806217279, -12.7975594871259,
                    281.150550466538, 48.7959824845697, 96.4740189184726, 59.2435976350377,
                    -15.0543351111662, 278.361713610698, 49.1417040924754, 96.0421597987044,
                    59.6678301739561, -14.4071148360617, 278.702796009137, 48.7758101534446,
                    97.2858508887275, 59.1933034283883, -14.984629566377, 278.749393710077,
                    48.6740558821183, 97.4851238789092, 59.7642101345605, -15.6479402433612,
                    278.064464097274, 48.8554385349751, 98.0502280531398, 59.5704849325438,
                    -13.883274482572, 277.220968515208, 49.3955325427019, 97.984635140748,
                    60.182616081631, -14.3170049565402, 277.640890159785, 48.8849961731103,
                    98.353435156065, 59.9715082998309, -14.3087738959495, 278.188406823557,
                    48.891443469175, 98.1012844005214, 60.621940208721, -14.1842602303765,
                    279.0195954502, 48.790910660089, 97.8556987583295, 60.4390383749602,
                    -14.0966696406302, 279.513241800978, 48.6100714871294, 97.8748611413199,
                    60.1796332431344, -13.6044336316626, 280.003985204375, 49.4727973105652,
                    97.1274195010033, 59.753951069525, -12.9655845680632, 280.312136537874,
                    49.3037452040012, 97.4562913890612, 60.2514651357266, -12.8983901953326,
                    279.56280379995, 49.3551001622448, 97.2825310522394, 60.3067644103981,
                    -12.7087003081665, 280.114536366436, 48.8479318797756, 97.5858378477522,
                    59.9499082144628, -12.6062316553059, 281.095437146205, 48.3745886085886,
                    98.6293598345858, 60.6740394282528, -12.2985659896985, 281.465925986836,
                    48.0989518897909, 97.6070827631157, 60.3313383536298, -11.8790873238958,
                    282.023999417286, 47.5112142709431, 96.849749245932, 60.3624948491637,
                    -12.6833602460927, 281.342999063852, 48.2493907175739, 97.5494322534648,
                    60.4286928152387, -12.9165537646267, 282.014652127013, 48.3843379818875,
                    95.6023363929294, 58.2199002208906, -14.1804892486205, 278.444270498807,
                    48.576901464409, 96.5439522233737, 58.9538591688646, -14.5275238460738,
                    279.68631319789, 48.1293163217126, 97.5200147258201, 58.4568660967323,
                    -14.9222314470374, 278.113660347507, 48.1681957688855, 97.4974746456782,
                    58.5781702208922, -14.9114103331779, 279.245874664855, 48.2156337791874,
                    97.7490474971579, 59.9764963800755, -14.2057209225895, 278.451473178658,
                    49.2093462007454, 98.5161282432589, 60.3892243366928, -14.580726599361,
                    277.938214733748, 48.8373489988828, 98.8247777626408, 59.8959070198674,
                    -14.6187932488796, 278.484493145365, 49.1302963273878, 98.0265774943306,
                    60.6780731207781, -13.8637270116427, 277.428073869741, 48.7348225067038,
                    97.901208709225, 61.1999782585916, -13.5874675919385, 279.43108587129,
                    49.1563306106806, 98.4092586843123, 59.7073797815027, -12.1824831503215,
                    279.841987261719, 49.7070961011754, 98.8190005925303, 60.5834592372583,
                    -12.9103183643926, 280.772764646982, 49.8563925229714, 98.6849420307137,
                    59.9669376179267, -12.6165207815599, 280.536598583505, 49.4637261810402,
                    98.6289095111027, 60.2770032136712, -13.1148290836223, 280.531279904608,
                    49.1778082139945, 97.9082943696372, 60.3480360286735, -12.9391484032632,
                    280.448484068951, 48.6605889114446, 97.1674551069144, 60.6429095407463,
                    -12.3242773454389, 279.987845388293, 48.9799867509356, 96.9457629500706,
                    60.4679057878933, -11.4766853888861, 281.401428719133, 48.1572892981079,
                    97.4948830789382, 60.0570994198688, -10.9537389699499, 280.49290904478,
                    47.9837453435929, 97.4349109421368, 61.0525058601533, -13.0643824826933,
                    280.73789254596), .Dim = c(5L, 520L)))

    expect_equal(
        recover_alpha(out, n_message = 1, constraint = "resolution"),
        structure(c(-0.0599404486331423, 0.916412114376179, 0.224072472768398,
                    -0.0513652844057617, -0.127764758834701, 0.87129464121719, 0.589541849431702,
                    -0.878647960895478, -0.897009108230634, -0.164726346221073, 0.647176044702803,
                    -0.0576206480996242, 0.150847710176947, -1.03028019734748, -0.458264868445895,
                    0.516919602349361, 0.355250435020487, 0.0450862479774514, -0.742702138519206,
                    -1.10551713707781, 0.227159635120444, 0.517586566887453, 0.012209444392397,
                    -1.06130547445707, -0.965383685938235, -0.224892773185928, 0.991897440544474,
                    0.460668555648283, -1.53468331935721, -0.728342445049634, -0.0771355472511175,
                    -0.0743067828006025, -0.466480845139699, -1.16625966589473, -0.646010061815161,
                    -0.724533280068471, 0.359256822852152, -0.459650078113953, -0.81181843597858,
                    -0.499152681731545, 0.225677905837614, 0.249860821931178, 0.171372651125182,
                    -1.103072840094, 0.0715066555670774, -0.749328687135574, 0.121935469901547,
                    0.585341560056747, -1.56542536308376, 1.68181089684984, 0.208887603770563,
                    -0.447527049641719, 0.743964180178182, -1.30426871292062, 2.34136920742191,
                    1.33503036927792, -0.554746506863408, 0.349734490767446, -0.872095756224155,
                    2.55875067771288, 0.795042895328749, -0.312453192928399, 0.708466670537462,
                    -0.310989400767625, 2.16934703820147, 0.997429495224139, -0.443513239541016,
                    0.400912373739914, 0.472609780458555, 2.03736254098388, -0.0196719226400717,
                    0.928868760758013, 0.285513588607387, -0.460152854102034, -0.257112003067959,
                    0.697999366069013, 1.03559565468763, -0.395075853252365, -0.726057395667709,
                    -0.0381328678781472, 0.445689983564165, 0.127971259619147, 0.30295052716383,
                    -0.988545557928802, -0.389815399640241, 0.14096488940487, 0.12728864920044,
                    -0.260105268242341, -0.682387727574337, -0.953920795969566, -0.174934941901171,
                    0.452654843190277, -0.718004202536861, -0.778262590971224, -0.829634287365565,
                    0.0838992262904412, 0.803750279516712, 0.366942549773853, -1.17733438589045,
                    -0.498520717061496, 0.383211748838448, 0.330525159777187, -0.024376671566813,
                    -1.05766078936952, -0.171596443678453, -0.2893429166038, 0.444064607069816,
                    0.104467567078018, -1.01656001860263, -0.761589241490611, -0.192425463260577,
                    0.95080704705109, -0.340385003879987, -0.791881145861112, -0.0178915586013204,
                    -1.01181892460079, -0.0891475657817162, 0.352103368120027, -1.16616844073113,
                    1.43587161859773, 0.551860149959026, -0.405748993184194, 0.590125624238851,
                    -1.3807950650623, 1.88497933533415, 1.28956121384246, -0.56211843358112,
                    0.0432840647505373, -0.817302778243544, 2.08160687111399, 0.834654093660362,
                    -0.442807825471419, 0.733544161177065, -0.16034658386539, 2.11938763641547,
                    1.1249058587496, -0.998530669422205, 0.581432553170885, -0.0226615929998104,
                    2.24685244150808, 0.838047543294977, 0.0996067106603875, -0.120587274484251,
                    -0.942301861980951, 0.873871389378252, 0.43872480197345, 0.370781675804273,
                    -1.10908563144068, -0.811420690585209, 0.583754051013415, -0.340631477648657,
                    -0.742915953251553, -1.06522653460451, -0.778796186699566, -0.649477129773686,
                    0.0266985133829678, 0.0529808428932421, -0.116665394391703, -1.21603158969225,
                    -0.980964474262777, 0.125479766261662, 0.13818715803645, 0.120599810157159,
                    -0.103543759727273, -0.719258039948414, 0.61595581105405, -0.373179195752553,
                    -0.0474139219629706, -0.937689378384647, 0.0740302098353709,
                    1.54151514558731, -0.333637614491997, -0.136838833451122, -0.508111350143196,
                    0.432588964071897, 1.29650753043498, 0.227996103073565, -0.735188504623181,
                    -1.18738724528538, 0.334385777279039, 0.763981528148804, 1.06265820295719,
                    -0.202163141988365, -1.19935576293665, 0.960844827062914, -0.500810431469034,
                    1.16875177759385, 0.0970421297354847, -1.90059445086828, 1.04238978317022,
                    0.348991901658763, -0.113462597513171, 0.396650527831866, -1.82201706186621,
                    1.49044350060169, 0.115190241698912, -0.580856876643498, 0.121145704052537,
                    -1.23233525795106, 1.85235386856868, -0.0640744235441275, -0.95537820567327,
                    0.988944204874144, -0.498097170368283, 2.13465321889669, 1.19940522870772,
                    -1.06004124447726, 1.02103222988563, -0.921433004712128, 2.16549432986318,
                    0.188981974611977, -0.0956021575113732, 0.194800918186502, -0.104899798401135,
                    -0.124787594766048, 0.278944713643227, -0.152745545055431, -0.581190565383565,
                    -0.571542630697849, 0.514568336548706, 0.741238815604135, -0.834197918464126,
                    -0.198934004763032, -0.302744216580464, -0.232163479980983, 0.127276739896871,
                    -0.746924056681741, 0.0777525872273372, 0.102552875958651, 0.0608082790667765,
                    -0.0252657562349441, 0.0723871569002954, -0.517820227539801,
                    0.598022330808561, -0.30890376802455, 0.204591887552567, -0.290334287739814,
                    -0.0599634077416908, 0.14717747207203, -0.154797624330826, 1.33877170607431,
                    0.0763113398805766, 0.0366481666481491, -0.187279835810415, -0.2383464390808,
                    0.912361467718085, 0.240417866124778, -0.20345594701763, -0.0991119557896241,
                    0.159189051069035, 1.09355091601387, 0.541824462981265, 0.0565822483988541,
                    0.173537476297355, 0.31474936840911, 0.70844986743063, 0.0490780103152,
                    0.336006227265329, -0.750164807981236, 0.804463423588942, 0.387184959051041,
                    0.351187689794144, 0.656504962437964, -1.0400173252429, 0.380453905394859,
                    0.297263383206115, -0.236532841654376, 0.332726179640758, -0.349336958665333,
                    1.36346990185578, -0.126469011079934, -0.42225432899329, 0.406438410707949,
                    -0.782304854326014, 1.57020485109837, 0.185199325217269, 0.154892862177853,
                    0.301046720510385, -0.414181273704713, 1.77984670299162, 0.158396886839881,
                    0.61214297203807, -0.23443816405802, -0.0557446709494229, -0.394041871005214,
                    -0.0798830330124076, -0.812293582437238, -0.0117341092825427,
                    0.222999279260478, 0.835557210518203, 0.808975789928155, -0.514462700555538,
                    -0.389734701180998, 0.0201033753017441, -0.480381987780731, 0.859863354022167,
                    0.123266381094069, 0.535580120365012, 0.908323549274286, -0.144588011252551,
                    0.294219700765694, -0.0714042781357054, 0.146107755095727, 0.863154394998503,
                    -0.228737398085883, -0.829706861928969, 0.551867963013114, 0.148376727982566,
                    -0.176723500944206, 0.115013659477569, -0.358979704016434, 0.00798488482850068,
                    -0.110189002575737, 0.185333198340118, -0.251296403723501, -0.0208034173713258,
                    0.471456804216267, -0.174031753009331, -0.357906457336242, 0.0890399520026222,
                    0.349603588343626, 0.291999320266542, 0.247919831845763, -0.147082730396676,
                    -1.06811435101673, 0.399190781784, -0.180470916694723, -0.254410818485269,
                    -0.72572670226981, 0.184718645144045, 0.466801630224325, -0.205777684343545,
                    0.346296364233581, 0.0928458169684632, 0.133746792488949, 0.358833457128213,
                    0.490301722232326, 0.539874609676815, -0.112733553986232, 1.82796583767703,
                    0.497394955039681, 0.235065938344293, -0.0115237505998493, 0.172081285643685,
                    1.20181681396343, -0.023298792653577, 0.803671606187893, 0.237447823796913,
                    -0.283288160043014, 1.36639872318506, 0.548354262027885, 0.291399638972155,
                    0.297219153655874, -0.561105400472002, 0.485159887668857, 0.647144964019191,
                    -0.216309348490541, 0.331568767518661, 0.200972534624043, 0.723483701680948,
                    1.01061384676208, 0.388146384004685, -0.292875267479019, 0.214014714128325,
                    0.688085978591204, 1.00534191545302, -0.139536868909295, 0.347764693888394,
                    0.195331860276269, 0.261667313290332, 0.563281624901606, -0.883603175054304,
                    0.302089504876903, -0.36024288268294, -0.512522497702662, -0.597917113416941,
                    -0.518448038592078, 0.0337515063412752, 0.759662174059372, -0.972210608789027,
                    -0.0402209808918599, 0.349727153271249, 0.122899987847802, 0.161732696357205,
                    -0.289589837237173, 0.443712304913163, 0.687192184652665, 0.694769313110328,
                    -0.203269785003386, -0.296720708396293, -0.127341105537568, 0.24018386864052,
                    -0.038501563439354, 0.165226885856669, 0.00269338693311738, 0.261475165186567,
                    0.456466193752978, 0.74647481976524, 0.46382302682963, 0.0598514731706246,
                    0.39884325203775, 0.473600205445416, 0.145382268966671, 0.456770054951903,
                    0.945131223354736, 0.676716056863199, 0.0182856808703544, 1.02269630483892,
                    0.461277619287159, 1.27629643982377, -0.191290875935977, 0.677620717379848,
                    0.754662301806945, -0.255123788708289, 1.2831926706582, 0.25024574448625,
                    0.608145784227617, 0.130210624836309, 0.37697828143672, 1.32175918460445,
                    0.7046538464962, 0.15034433332894, -0.484756345724747, 0.489893626034416,
                    -0.312920959542737, 0.756105246654215, -0.917806737692874, -0.247964525854513,
                    0.474187443565086, -0.149077241642516, 0.251967517866206, -0.445571768515812,
                    0.068525895706216, 0.0402329025291692, 0.701470576823056, -0.505426707351489,
                    0.523677467081088, 0.278983043394078, 0.047051537963597, 0.691651067506825,
                    -0.693491599382639, -0.438184401291543, 0.555682642242175, -0.168531588402999,
                    0.619615080679637, 0.150904783801963, -0.260010480080886, 0.208486803835228,
                    -0.0196800173000042, -0.239450492974754, 0.296090380097894, -0.481614017216307,
                    -0.971026234055273, -0.184495158632302, 0.323977456586704, 0.309476727839069,
                    0.376020268246474, 0.120300007514913, -0.0398688724089595, -0.604606953120197,
                    -0.310697226683914, 0.905995593797587, 0.155847547692019, -0.0308784040572618,
                    -0.262231635184975, 0.367427483852254, 0.582192452543918, -0.0419534375562876,
                    0.625651797741796, -0.148864029448816, 0.522322438865956, 0.127182852428874,
                    0.711980147721761, -0.249584877646328, 0.304494172412319, 1.01513330161015,
                    0.0336508330632883, 0.72771349606063, 0.315027688112281, 0.36420846323233,
                    0.582195385210071, -0.149308226805259, 1.08475062228797, -0.0712645452408935,
                    0.0980013665632669, 0.62599639211146, -0.307192739766245, 1.06616919726002,
                    0.320981472178431, 0.430455695184946, 0.113581175659419, 0.385582347780455,
                    -0.5305153237789, 0.262556900946294, -0.600535120515417, 0.635689775504659,
                    -0.738681845215961, -0.0469111623040703, 0.0403927942321047,
                    -0.645178695492092, 0.289742920101816, -1.09708467444021, -0.187339471360133,
                    0.128506729575264, 0.381207849388403, -0.565380948852408, -0.541745942939341,
                    -0.324112307874032, 0.414252354436115, 0.0760495374605341, -0.554687045262824,
                    -0.646588511386076, 0.799124102094311, -0.0584795500125779, 0.0911981382665203,
                    -0.403471251009762, -0.0806761629594961, -0.291819061323963,
                    0.498336908444486, 0.473962802530991, -0.638504480972443, -0.134775287114085,
                    0.037030793845986, -0.301326073357728, 0.169300546692284, 0.0351141517443523,
                    0.183516722247869, 0.383715970788671, -0.350326352856882, 0.348717276578697,
                    -0.36311966006437, 0.342547168180977, -0.191857760093512, 0.368044526012397,
                    -0.245254556081676, 0.0486276929055123, 0.882891031149379, 0.0957070807524474,
                    0.0450285514989446, -0.216590971103216, 0.0785963874622837, 0.235503169320868,
                    0.226026875500168, -0.522633502372429, 0.609287907078652, 1.21334270065782,
                    1.00922760169993, 0.753562191498929, 0.220981855872765, -0.0016624131126548,
                    0.335931917276923, 0.988449944128451, 1.38104349513725, -0.65635666718336,
                    0.668574214639165, 0.965660697663452, 0.517708304707838, 1.50669610339253,
                    -0.17832054338956, 0.734886602653688, -0.504819283134623, -0.497289250055154,
                    0.104438425548153, 0.189944855170239, -0.598876481879415, 0.786028480369964,
                    -0.978982288922267, 0.0188621449415791, -0.355245329672385, -0.307289318970589,
                    0.449856746317238, -0.816188608437251, 0.031656307220544, 0.45021299734816,
                    -0.623414766970654, 0.0415410087829002, -0.129352348088498, 0.202666573502654,
                    0.395848051661602, -0.229181839591092, 0.227163134662248, -0.07644946632109,
                    0.963158905897103, 0.950335187840864, 0.350722298193091, -0.0207250378260539,
                    0.335415384932882, 0.215958155899557, -0.0365474437209059, 1.22718786151663,
                    0.192597089769805, -0.168798835442105, -0.28256993569579, -0.108257170270832,
                    0.135194449524874, 0.563862096972684, -0.404990747851072, -0.0340094403979379,
                    -0.354194520892833, 0.547340223342133, 0.0085440748380563, -0.148681563628307,
                    -0.157776293699669, 0.175800369056986, 0.400412697717456, 0.403816089685023,
                    -0.104354038620971, 0.51132353963493, -0.0992849631009562, 0.415734092088996,
                    0.753510185751423, 0.840692310380029, 0.651205924261717, 0.256472226285723,
                    0.193577111444199, 0.043259604726245, 1.09202874033127, 0.133058346651467,
                    0.25748221269518, -0.601036185301808, 0.838452775655526, 0.680883938799282,
                    0.772357279953098, -0.261478504754308, 0.30706011627791, 0.169481195505028,
                    0.60565708373413, 0.48510363314773, -0.257685733766237, 0.33718762600688,
                    -0.105682419529614, -0.998281762426987, -0.123324350853238, -1.25922952097031,
                    0.418996848373553, -0.241491859172811, -0.977713918349218, -0.175377808991698,
                    -0.866360478544863, 0.867995085566264, 0.231452273284162, -0.677304307407015,
                    -0.185912467275443, -0.18980659065096, -0.0945376454077689, 0.426979318651007,
                    -0.681992864010152, 0.571384018679055, 0.50391574547001, -0.427255282308664,
                    -0.0654615882244798, 0.434063288648678, 1.1598586156922, 0.0830190295177999,
                    -0.185001197937027, 0.121353735742076, -0.887176056920568, 0.190121750404728,
                    0.201443864203753, 0.549188944046477, 0.204225218240737, 0.242068111901779,
                    0.218447165180919, 0.102050964690875, 0.14524122060439, 0.356942098898173,
                    -0.271650013649074, 0.290739383702316, 0.122423242899401, 0.26496159765415,
                    0.107541328511786, 0.765237958760569, 0.389493179150918, 0.401786179465006,
                    0.230912286949376, 0.944378464943014, 0.763996884243781, -0.423752277364599,
                    -0.265405258101879, -0.290476494035403, 0.430321601244685, 1.24607479527733,
                    0.447565724558046, 0.243989480399193, 0.0170418186606867, 0.46506625880059,
                    1.30243103298073, -0.0823337249984206, 0.126157195649455, 0.39932066471323,
                    0.942655033520765, 1.17589335337471, 0.205013486250607, -0.902065934732462,
                    0.32026915264592, 0.526452902601562, 1.66562152649663, 0.155816768830618,
                    -0.428857742810209, 0.0495986084275728, 0.901600771861325, -0.796146917347812,
                    -0.43374916882982, -1.474886979042, 0.148321626306142, 0.259800744327919,
                    -0.621626994947342, -0.403995164077514, -0.511836099941007, 0.671207050091056,
                    0.357946076305439, -0.478399114410272, 0.185702641684827, 0.154293868362316,
                    0.0441924958332436, 0.623496495439412, -0.831403205765099, 0.202438376363432,
                    0.020467306080775, -0.581245998096193, 0.749281783592238, 0.409123311173726,
                    0.619279835567283, -0.0649607099362868, -0.596589710324622, 1.08714757032419,
                    -0.291992939094171, 0.415424955625383, -0.078559512328269, 0.978006381533248,
                    0.566304601008483, -0.588058622008413, -0.0713323044154208, -0.598497798101903,
                    -0.368393974516408, 0.104802937532526, -0.325607476347221, 0.254018193641627,
                    -0.12569134177977, -0.669082503271028, 0.505845583133834, 1.19215134907941,
                    0.185037910172241, 0.385739472780443, 0.304015803348449, 0.0344533324345093,
                    0.953050518339176, 0.536662844666466, 0.473598728195299, 0.341218072452733,
                    1.32753403644551, 1.39974038332095, 0.167481512483505, 0.673559204375437,
                    -0.338557731341211, 0.86933251257031, 1.7748094765202, 0.139995449875414,
                    0.238734370516454, -0.327565931224228, 0.694311993743582, 0.861526446389338,
                    1.2128806672676, -0.23937482802242, -0.249467626489576, 0.932435087147411,
                    1.22050242491215, 0.823396057425327, 0.557629211228402, -0.518374720382752,
                    0.887120052821984, -1.87438508342905, -0.978939733364996, -1.0230834371892,
                    0.147042316856613, 0.598217735022502, -1.3760030377373, -0.260321913788175,
                    -0.69339254881541, 0.284083391288618, -0.0392237071224031, -0.504165822083762,
                    0.3196250049524, -0.182816888442687, -0.0784571208469629, 0.118860206731483,
                    0.124030592041532, 0.19805728761196, 0.344691249875666, -0.106420982711285,
                    -0.0457347626299622, 0.437212514601943, -0.113570035361363, 0.0275228349568692,
                    -0.526852839135842, 0.562238796796436, 0.107733357730611, -0.344774059676695,
                    -0.397158628003524, -0.0016865466780871, 0.597738130640352, 0.294976426344735,
                    -0.025714367819262, -1.22518541050482, -0.560931615049128, -0.770666426374419,
                    0.051262349944591, 0.480909858155599, 0.0721700715337192, 0.0816072712353844,
                    0.624929978882165, 0.851319115488498, 0.355398706193341, -0.735961798475159,
                    0.261730658184547, 0.531447275557525, 0.926622091438916, -0.318806789329244,
                    -0.184792142751661, -0.579025338249721, 0.389442598746228, 1.38666559555071,
                    0.676002987410186, 0.260626825016971, -1.04537609797794, 0.298705484365769,
                    2.3352816038954, 0.715268171570983, -0.122337583480771, -0.687895029834664,
                    0.850573203487215, 1.69441855318015, 1.1928225112423, -0.298782977088422,
                    -0.298694890745054, 1.22384404188, 0.566449116411383, 1.239797649617,
                    0.0129919526933406, -0.191402611685675, 0.461208357916206, -2.06863314117516,
                    -0.62167584400197, -0.710659372641459, -1.4576982117172, -0.160302416080157,
                    -1.02662898518752, -0.125766642310396, -0.487559699620675, -0.0565049923696677,
                    -0.659043134512501, -0.775224975023655, 0.0269770091786086, -0.148397924716477,
                    -0.346505847365165, -0.0487561032366273, -1.17145515145881, 0.041059079230763,
                    0.596314621354111, 0.179608108313118, 0.119098775322527, -0.129101849037117,
                    0.193164919373601, 0.193460065121485, 0.497152980108297, 0.387862179041122,
                    0.150334392871335, 0.424054566830705, 0.335698642949382, -0.358635006176684,
                    0.069696153189085, -0.847033672422106, -0.57185144137344, -0.388579301620211,
                    -0.443560289157205, -0.773934308772567, -0.0115121537880212,
                    -0.119896430928122, -0.073519607053953, -0.651002363589697, 0.0345510360508285,
                    -0.0968276678977706, 0.20389112270982, -0.688578589364624, -0.911386437268789,
                    0.111412381196942, 0.435513151427802, 0.551778746847788, -0.342770350118627,
                    -0.685471168670148, -0.252195664388012, 0.4672128034809, 0.877346036663425,
                    -0.356102025108896, -0.750808018610257, 0.803224444910114, 1.47783181135571,
                    0.663734921048615, -0.0959995166698207, -0.668733203532781, 0.360805612668315,
                    1.15961978524948, 0.022487242313332, 0.112364187932144, -0.183937183737612,
                    0.141034051766781, 1.34032898304194, 0.0129467170297062, -0.19814035342543,
                    -0.315455495325182, 0.11289652654699, -2.70134969721778, -0.44695385268318,
                    -1.34927392826739, -0.833020237827895, 1.76457636796113, -0.956131990783994,
                    -0.487619704905171, -0.589145853194992, -0.216266029443261, 0.936411587928177,
                    -1.66073816472254, 0.390849286646926, 0.0172244832120612, -0.290857128428343,
                    1.49485048450541, -1.63012039660599, 0.400627695726143, 0.548753433094475,
                    -0.379648103058457, 0.254331051264131, -0.684779528870283, -0.271981715193647,
                    0.502477800247391, 0.151971713169473, -0.132186044287437, -0.242586156564092,
                    -0.588407311461822, 0.0705206626594475, -0.0255913687860243,
                    0.284204345954038, -0.976050823801059, -0.888084658725631, -0.135048778057637,
                    -1.18709968502978, 0.179486791544718, -0.665760355778194, -0.0654203114225993,
                    -0.107187239803835, -1.17664569831351, -0.214959665520723, -0.527171037531839,
                    0.258492570994008, -0.702657297315653, -0.539197463446072, 0.459867936755046,
                    0.346525414569775, -0.384013065835575, -1.02775090833282, 0.21883178281422,
                    0.133769716417419, -0.190429244628989, 0.185936245556377, -0.941570023788557,
                    0.147087947051062, 0.19919223830928, 0.390397875372056, 0.300160760417498,
                    -0.150888044000084, -0.948523502494481, 0.282017406541854, 1.22602983622724,
                    -0.374849712650473, -0.125025891921752, 0.0347640863024594, -0.586381376023951,
                    1.86582406482423, -0.802942714962683, 0.262158884723959, 0.229192287668212,
                    -0.320591904722136, 1.72307680273625, 0.385815696033944, -1.82849900147697,
                    -2.46263910320681, 0.283131454583526, 1.09112980351138, 0.317007017669852,
                    -1.94208075441865, -0.650434565421648, 1.51284318931963, 0.743669863653337,
                    -0.0982253216816957, -0.972276875803161, 0.187624762469483, 1.68524704227019,
                    -0.0326723441919254, 0.633446917113758, -0.450006108032923, 0.140250990807857,
                    -0.430869020814811, -1.15388323886421, 0.232876770249504, -0.668722568329116,
                    -0.246126154205569, 0.285203966524449, -0.281814925900216, 0.549544470857796,
                    -0.532200401146213, 0.308780783291326, 0.357019944431954, -0.571429808275827,
                    0.104061534139305, -1.13972207773719, -0.556047717699471, 1.04284459143105,
                    -0.77087482303356, 0.250979820574173, -0.0134153285397787, -0.187263493774509,
                    -0.391629539229598, -0.579956938409188, 0.603911632621049, -0.547206009859053,
                    -1.07000413316646, -0.763249959269871, -0.759658341190234, 1.3083682283924,
                    0.260382110066615, -0.480550761165375, -0.335331097434789, -1.04292220407834,
                    1.4762213346294, -0.151487151006563, -0.207885615729737, -0.854871724206077,
                    -1.31390491734944, 1.4304661319631, -0.0616286044427596, -0.150663991495321,
                    -0.511571696048662, -1.08892972083375, 1.63910788094984, 0.836514970622744,
                    -0.915212313139932, -0.351461578495787, -0.348104618604, 0.216905666170703,
                    0.84936099802708, -0.174079405070302, -0.478987804575695, -0.798285097795953,
                    -0.218463470935149, 0.620958910941241, 0.667509752252357, -0.44412251046478,
                    -1.62858029026992, -0.0252319877233163, 0.0282954566884541, 0.361721992610995,
                    -1.20018899400091, -0.821446791065826, -0.364940203815593, -0.834974490622628,
                    0.238977715485589, -1.14478667172054, -1.40577336397994, -0.761293914134285,
                    -1.30708857380226, 0.531516544095751, 0.206826755146494, 1.58691147384838,
                    0.000383719381602532, -1.85757505284266, -2.01828649160197, 0.735967260135595,
                    1.68193347142225, 0.0401979638276941, -1.81170445670102, -1.14531376699949,
                    1.36876596778782, 0.92162254142432, 0.422266081043887, -0.982886020094384,
                    -1.44789553092525, 1.35648603314576, 0.277715931911145, 0.0829708502168387,
                    -0.64913670856717, -0.69976469017044, 0.424274617343713, -0.258809245056373,
                    0.612793732231438, -0.638861892901616, -0.783088887836954, 0.456702910591844,
                    -0.166766929263501, 0.893143008801644, -1.16728852692867, -0.604460564256897,
                    0.605891923283963, -0.374503338317033, 0.346398373123066, -1.1280285281053,
                    -1.26385578608216, 0.677420603107947, -0.707984980042515, 0.780280411856523,
                    -1.01673468144425, -0.716986975367519, -0.833916511606503, -1.17039331932678,
                    0.674117081145937, -0.107528952426179, -0.59784800329615, -0.739468094849855,
                    -0.524918758801875, 1.62479372689702, 0.141000346633817, -0.164765124738153,
                    -0.788633701149166, -0.958524954012736, 1.16522182524207, 0.19723561157744,
                    0.274793098630482, -1.20748934251449, -1.75569089858068, 0.734900905317716,
                    -0.221032547584173, -0.319039576101687, -0.0676671794955155,
                    -0.895748525784491, 0.353658777349921, 0.160108040727792, -0.676134842874035,
                    -0.025631896491177, -0.0418570507329104, 0.0327423355940937,
                    0.389060491448305, 0.0880645273982168, -0.512058919954399, -0.5531099021368,
                    -0.0323440371488175, 0.379971461005752, 0.961104055079886, -1.00086675782489,
                    -1.59557975451214, -0.200277716286536, 0.0143694772837648, 0.432395223425374,
                    -1.00836730792102, -1.14882339024112, -0.769456652977169, -0.537087461900271,
                    -0.0859878558796368, -0.998507716441196, -1.01885163502563, -0.670546162119507,
                    -1.59188076222745, 0.207926484768535, 0.0408172111936551, 0.895706342209307,
                    0.245080174490028, -1.62559294451236, -0.952181304031313, 0.199702712752696,
                    1.68764361926049, 1.30705659081244, -1.23688847750938, -1.39163107863658,
                    0.520728742853208, 1.26602700409223, 0.443511236409108, -0.590389994529247,
                    -1.00064895745706, 0.622263034949952, 0.124837727000151, 0.470430097100433,
                    -0.445746509813858, -0.379289479133092, 0.119571930950087, 0.392265954629463,
                    0.536830906412092, 0.0666988929046539, 0.234573174161568, 0.537788604951686,
                    0.359024163754131, 0.750909094535416, -0.673715732350402, -0.743285128343729,
                    -0.0470992736235019, -0.9731298852622, 0.274645285641242, 0.0920990271083326,
                    0.120927915537919, 0.388940729002002, -0.674096779030114, 0.639297597561182,
                    0.846751569381766, -0.947001022347933, -0.665828858226718, -0.810291357819571,
                    0.160056678762857, -0.892431586318409, -0.231058107929925, -0.735757198528944,
                    -0.998934142092175, 0.374159823341159, 0.702737671704512, -0.847322437486156,
                    -0.518993586405003, -0.84599182932769, 0.774841978269365, -1.41138945607877,
                    -0.62095452533643, -0.335651724783133, -1.21085014722976, 1.10767566386517,
                    0.0674405098931317, -0.736101149048892, 0.11466923192765, -1.37642184496627,
                    0.425358291356304, -0.0424871592839509, 0.382162443516421, 0.654282657190265,
                    -0.555484179253704, 0.948663170273193, -1.09012429295205, 0.909247110873991,
                    -0.599277617001547, -0.294192360011866, -0.080304717578457, -0.294787943060072,
                    0.835513534612659, -0.316057093532208, -0.337213969080273, 1.22891345929565,
                    -1.2239831029855, 1.33419085465752, -0.867218942125874, -0.266167508756876,
                    -0.726153943417444, -0.771099028361689, 1.05336498519918, -0.301518799484739,
                    -0.139877530858815, -0.715542459657144, -0.676977646209053, 0.853087159075699,
                    -0.465845615008874, 1.28839536552271, 0.66136966764671, -0.660674568084531,
                    -1.71594774426754, -0.619244791669075, 1.25629062419402, 0.36998961733542,
                    -1.3223491315822, -1.54573570600633, -0.330242813498472, 0.792094997451116,
                    0.251608160534403, -0.60571546159639, 0.140038596632053, -0.29544163700595,
                    -0.485116110714639, 0.168831353121391, 0.222566082145988, -0.670075130386408,
                    -0.0672000801862964, 0.0550491885782423, 0.741195631159286, 0.0414808391533938,
                    -1.04459139689095, 0.809797300714933, -0.426869442681792, 0.841176757100612,
                    -0.887861730040427, -0.23275817811259, 1.12928050574615, 0.0281366616009009,
                    0.455638852833403, -0.112956050654802, -1.01477923101314, 0.640560770927095,
                    -0.426350449555539, 0.371764927595791, 0.862972768374746, -0.860778067870712,
                    -0.448161188677517, -0.496392372360052, 0.770093660557876, 0.366858436398161,
                    -0.492059578991132, -0.724312645489931, -0.946660461029042, 0.0539808290120902,
                    -0.50195392782345, 0.450334484876549, -0.615433111760893, -0.666510460476677,
                    0.690515159615671, -0.501012590396371, 0.193458563007937, -0.949257217383519,
                    -0.699612571970244, 0.402305565499972, -0.979020279464308, -0.778750007021358,
                    -0.803699362185824, -1.26291553899739, -0.184665517003452, -0.219009462678727,
                    -0.087441945629422, -0.329945387691509, -1.58138006650198, 0.298166574574111,
                    -0.174714118543264, 1.15403680897919, -0.978914517741515, -1.13137251317627,
                    0.428554637589315, -0.237954838913865, 0.547753243672162, -0.839013469130492,
                    -1.23170762897502, 0.288662955791807, -0.462339830637745, 0.755322659394096,
                    0.00154908585673752, -0.245019017202029, -0.355801956997936,
                    -0.756804220687613, -0.147331254476455, -0.567186146776265, -0.32797831193025,
                    0.196969959210776, -1.22048984608726, 0.76373209471609, -0.243727562539789,
                    0.50454577784359, 0.0899084024996739, 0.526663399988351, -0.688119277386079,
                    0.29444020476123, 0.707212867674009, 0.141156764669148, 0.0306561937793361,
                    -1.49605803646043, -0.108680055771572, 0.556660286483254, -0.328935602699602,
                    0.11133996079775, 0.892693025117353, -0.0955586573996356, -0.713009700306224,
                    0.387091349005544, -0.351017499405486, -0.13326162665183, 0.206214600310091,
                    -0.314873567555641, 0.651896875695115, 0.262750102575723, -0.0952231759251276,
                    0.168983516882292, 0.430734394388365, 0.261159311598163, -0.324780843551565,
                    0.442356861144276, 0.29213764380637, -0.281155017484801, 0.722910121071948,
                    0.439013025326801, 0.372499715386908, -0.0841641137346301, -0.40979973764454,
                    1.47985248575085, 1.3240614058019, 0.794391443976963, -0.754983641819081,
                    -0.487960053311092, 0.925089793831859, -0.199495744219277, 0.4274376152662,
                    -0.681839420955356, -0.298602411480126, 0.983297661939304, -0.986503488721553,
                    0.62926560201268, -0.0674857532167152, -0.340705859273783, 0.298696463215038,
                    -0.658239149222965, 1.21161449612313, -0.11512579353672, -1.00669396949848,
                    0.312109192559518, 0.255545458785747, 0.219134561875393, 0.285309734235312,
                    -0.847689229026689, -0.0421054301220067, -0.867134791453736,
                    0.272287761919699, 0.390598000166847, -1.38080377699387, 0.0945645110262774,
                    -0.382432909999068, -0.221123744299604, -0.905052610609772, -0.534700544089837,
                    0.686535255060676, -0.946621896313673, 0.559545134622056, -0.643361569161812,
                    -0.132677899235513, 1.39527464209661, -0.506969964175937, 0.452396616658376,
                    -0.178290354713596, 0.0840302907705708, 0.47700061281779, -0.560050802177841,
                    0.028508662689319, 0.0179945233878414, -0.0562169624294682, 0.2014808996596,
                    -2.17966046680502, -0.550318150957196, 0.102366692521059, 0.109201624543974,
                    -0.245743036604466, -1.0063128112243, -0.256828899082109, 0.187284255767139,
                    0.647382499807292, -0.0641562812595282, 0.0931994767543216, -0.0626074781873172,
                    -0.1692042826582, 0.366686781781252, -0.18662501460642, 0.593975031700211,
                    0.545274896561352, 0.125812616907922, 0.241453691820863, -0.0816386711796326,
                    0.872345557045364, 0.0466605093542682, 0.350960630408423, 0.309536158206766,
                    -0.541222865682585, 0.2374244038255, -0.334069241143936, 0.0450375040866788,
                    0.154495798703479, 0.151393036290074, -0.936647346933427, 0.430224930572962,
                    0.00120439774227066, 0.109585620085147, 0.320640436230775, 0.495718110060096,
                    0.358766989732231, -0.144686738654002, -0.1899442274855, 0.598782954919926,
                    -0.169827936729476, 1.31339380848908, -0.651268167247849, -0.956744516160167,
                    0.926712400676223, -0.590229305444524, 0.293919560963388, -0.171020186398152,
                    -0.75136993848348, 0.823025563889701, -0.488225036661476, 0.946768468196211,
                    -0.73336229839947, -0.357126724103045, 0.136246657816883, 0.509183970730234,
                    0.602238236436989, -0.174852151833363, -0.526022836099706, 0.350573016159586,
                    -0.569946302287605, 1.62294806095977, -0.159403724973487, -0.420025549872106,
                    0.345010633797472, -0.789903533289845, 0.421760656905377, -0.10107210951216,
                    -0.823039633032865, -0.282766149774545, -1.09120554407924, 0.582154345027163,
                    0.808555265372632, -0.810453859708594, 0.56730494463983, -0.462453581918965,
                    0.475426032288738, -0.0827693953941662, -0.193234749280989, -0.185851769381202,
                    -0.78212823992078, 0.99619911365744, -0.64822449943432, -0.0559658862081696,
                    -0.396599240422802, -1.13579593750811, 0.724646192571129, -0.0894059370756395,
                    -0.0404048480818293, 0.122051633955522, -2.12654714185047, -0.0711268279671344,
                    0.571224699130966, 0.251940421663926, -0.315206401048869, 0.706924428650254,
                    -1.22593992279232, 0.68747583870649, 0.740137059070765, 0.236227319366009,
                    0.907763719521085, -0.384257443601427, 0.14653810467442, 0.265934136276033,
                    -0.0106523749138603, 0.334679796723222, -0.276826881651544, 1.01850646609859,
                    0.918692843396625, -0.216893210825305, 0.789598178370611, -0.556222278337827,
                    0.847645968890077, 0.745390829431631, 0.0297624741635332, 0.666464824850536,
                    -0.00165604058327062, 0.530014841082192, -0.0309700771159362,
                    -0.143119817559324, 0.309050209064281, -0.483172459409985, 0.26578856402358,
                    -0.222938742077476, 0.041786255362922, 0.165809507700814, 0.0861229755603574,
                    0.2713978597601, 0.146326666363892, 0.271652941799971, -0.237612178695329,
                    1.17139882659981, -0.0144567450864201, -0.228254194445045, 0.378656736984567,
                    -0.548706929062504, 0.696485127530423, -0.253306132607349, -0.195204389731875,
                    0.245510432378467, -1.06952400539912, 0.762041654154103, -0.225029600007957,
                    -0.147955706311864, 0.152555166916727, 0.12053012046465, 0.214434964944132,
                    -0.161779826345366, -0.388833968162145, -0.0748618197098381,
                    -0.180038481429193, 0.494050370630418, -0.886884692414782, -1.11408277614576,
                    -0.413635760614397, 0.228097767346782, 0.581845686470729, -0.578441747174566,
                    -0.887776475639267, 0.59217277754999, -0.225019052627061, 0.308861862361709,
                    0.377036847909636, -0.524015064638955, -0.21138291510897, -0.784379233144747,
                    -0.14760027456532, 0.691544696859893, 0.100484841658552, -0.747765642829307,
                    -0.84085401184808, 0.839492053883163, -0.127329958377913, 0.084939027530907,
                    -0.675616388190861, -1.05091915938075, 0.777682926994601, -0.0678415996853801,
                    -0.283553539377873, -1.09836013025293, -0.590273269713322, -0.0292909848068916,
                    0.677788445981946, 0.66235561585786, 0.395218605301721, -0.436794075619915,
                    0.483510824650409, 1.3482328871267, 0.771458577852441, 0.263196818716978,
                    0.513814832297996, -0.0306233812717664, 0.755906188947932, 0.154784354580863,
                    -0.340990190543494, 0.379664128855751, 0.0999193255988757, 1.31871774072735,
                    0.6381392622686, -0.39165402565277, 0.392383801406766, 0.455286006793415,
                    1.14536636241749, 0.290892498045486, -0.185069789173639, -0.167709382661471,
                    0.554644977415109, 0.36792696503602, -0.0137802404365746, -0.447734037092147,
                    0.125751247907232, 1.48687714892205, -0.549159278311265, 0.637060386419915,
                    -0.0294239449273164, -0.0591758916981462, -0.0949440634886685,
                    -0.0808456340803332, 0.32145471689567, 0.404323146948769, -0.27885051632299,
                    -1.20918057429625, -0.0493572954966055, -0.47617397366721, -0.196516984872517,
                    0.257797919077689, -0.335178418562293, 0.239563303486079, -0.509543517566656,
                    0.420651287464636, -1.12734649723645, 0.661519028870032, 0.218812383831192,
                    0.47538564569071, 0.226377052284363, -0.257012349074174, 0.521966444017551,
                    0.0747806716231025, -0.0417044929199895, 0.254118729281799, 0.665634666394268,
                    0.350636795073143, -0.685005093850123, -0.818884073312464, -0.00407778572576944,
                    0.437891196484841, 0.0329027765284877, -0.0889202097093573, -0.658628330410468,
                    0.60602591191369, -0.834984872016577, 0.454344080252923, 0.0151431402284459,
                    -0.615263000972963, -0.559800509618952, -0.279680985050632, 0.895378156401478,
                    0.25878766979838, -0.210069411456224, -0.294265940706566, -0.432575908468891,
                    0.492450183171371, -0.121365189713458, -0.188258386513733, -0.752843393480532,
                    -0.405664977745388, 0.590625419371577, -0.363758572532475, 0.25137397716255,
                    -0.980202021879222, -0.658214021949021, -0.151326151605417, 0.929644871279805,
                    1.04872710162254, -0.111547898652532, -0.673562473678516, -0.750835330399923,
                    1.05253343152975, 1.25988150935131, 0.554827442079016, 0.406657497989116,
                    0.0633653825964302, 0.709060595415849, 0.64576148801137, -0.00835481456638831,
                    0.310965141166932, -0.0780302199137566, 1.25268841435252, 1.19376791042096,
                    0.158316276744358, -0.722496868346894, -0.435764250169854, 1.04477320479333,
                    0.984427542128286, -0.528745453696501, 0.12662669293691, 0.236080506715382,
                    0.43427067634849, 0.740816023631421, -0.416615949745591, -0.158181641197899,
                    0.499757633286379, 0.262240883981235, 0.390987983970847, -0.292382388001585,
                    -0.648214440798185, -0.965921334492805, 0.44654621599193, -0.329734394388097,
                    0.175298859338881, -0.757544094749036, -0.536738655395908, 0.229261638063377,
                    0.12415651066185, 0.253792115809773, 0.531286487551355, -0.150593253625232,
                    -0.130467982492775, 0.233290631879788, -0.0880121331964006, 0.0802702518355574,
                    -0.0131438276243898, -0.806857201328597, -0.31798254104076, 0.387704243097176,
                    0.411926721862685, 0.839201879945222, 0.216197084478836, 0.275500515709231,
                    0.358896581717971, 0.92385605466778, 0.141735120483872, 0.161741317718793,
                    -0.158942659193455, 0.517457214948536, 0.0591879636453712, -0.0415560079922273,
                    -0.30532088739929, -0.0168741266568559, 0.666940511434007, 0.559189429448679,
                    1.3184372194041, -0.0465498465115672, -0.630649194153364, 0.396962784343629,
                    0.472006054558747, 1.58603979017786, 0.557808272219461, 0.19134350776874,
                    -0.492473067396361, -0.934300831429823, 1.0986276450343, 0.520662704978747,
                    0.0405102010792218, -0.961835396464835, -0.596031128116877, 0.203728864261763,
                    -0.213145647760712, 0.610540532187023, -1.36098797373927, -1.28819002748384,
                    0.847819945052549, 1.23972429385302, 1.49220860637747, -0.32241082911915,
                    -0.383453843529665, -0.0819587615318085, 0.961906607074127, 0.887924346613003,
                    0.408384080238733, -0.275871705254033, -0.203102275282191, 1.13532278808472,
                    1.84255258950989, -0.0947172114621537, -0.864717824980389, -0.79011502830889,
                    1.24328852018527, 1.40163632178926, -0.153763036706493, -0.628222186128411,
                    0.0909144131719302, 0.968197211747992, 1.11484946471079, -0.629795939315841,
                    -0.482838127683252, 0.658038826679189, 0.800037141122516, 0.377166596112525,
                    -0.616068763553727, -1.12917445241278, 1.49959310448276, 0.283976049378737,
                    0.631632590559918, -0.259871047476906, -0.166884333009719, 0.479357357454774,
                    0.220495010477663, 0.709698339450088, 0.216267925549204, 0.464590038099502,
                    0.884183460550958, 0.642945936762459, 0.259268050821657, 0.689573694110976,
                    0.328859512000278, 0.238918473484773, -0.248754309460082, -0.458650298881802,
                    0.362856518060056, -0.232558909801597, -0.329989101046749, 0.0550692490176345,
                    0.0171620383405053, 0.509753795987223, -0.0978986390201584, 0.815370631569351,
                    -0.29006295927627, -0.02547605566809, 0.503394224154647, 0.210052563537033,
                    -0.122714112356334, -0.410737790015446, 0.122721756943855, 0.170473230439541,
                    0.31058580650091, 0.543546798349468, -0.295589833830121, 0.209340936271717,
                    -0.0660987520669778, -0.279508934001356, 1.08932308682157, 0.608427381506544,
                    0.202142906933096, -0.182598335445277, -0.341399047396806, 1.08262818397733,
                    0.388154727214101, -0.223456604863037, -0.905214205164924, -1.13258335329921,
                    0.741796222330038, 0.299054487749203, -0.0775185718926628, -0.726039662367555,
                    -0.834408659638584, 1.39372219808789, 0.535109722549237, 0.856291672932039,
                    -0.678849499399746, -0.352963226508606, 1.49809246607958, 1.01280776354002,
                    0.353337779461268, 0.55265519449344, 0.594367171896394, -0.778087915829104,
                    0.857897195446355, 0.365619652057742, 0.44890430629593, 0.302273512271242,
                    -0.717296820417033, -0.212720679619594, 1.63038891789067, 0.192671494964628,
                    0.183068798180045, 0.0928425984252499, 0.153876397947215, 1.76113774605253,
                    -0.474004372159101, -0.732593292190543, 0.572628505980305, 0.464479686271289,
                    0.852684917382455, 0.00884811351951953, -0.277630218160908, 0.167552165573113,
                    0.375229609263386, 0.728024826784766, -0.610682795600724, -0.74976002695216,
                    0.286273269796681, 0.114342660920784, 0.298025885613267, -0.777476430354767,
                    -0.0496378213474884, -0.218656526935035, -0.713470322558834,
                    0.421321967778455, -0.739482472244504, -0.27448184729262, 0.1145041434132,
                    -0.517428213779453, 0.0575091243133421, -0.131729938109856, -0.213316420504469,
                    0.0938542532421991, -0.764400270050373, -0.310495843131999, -0.123285128526419,
                    -0.0283522143434567, 0.609777699999597, 0.238586791516326, 0.0331739430308846,
                    -0.360837864963315, 0.164745240194776, 0.443124318590606, 0.158822483611516,
                    0.311532580364105, -1.18014523753558, -0.0985910435529718, -0.560235411286556,
                    0.464740714241572, 0.161692545431382, -0.225061708449317, 0.249468262493579,
                    -0.118424182282837, -0.121141018836425, 0.53099470919505, 0.349479741450864,
                    -0.206099486213219, -0.0934951541821931, 0.176005681374505, -0.0135511951543197,
                    -0.150576965298228, -0.663544217600844, -0.0135268292264357,
                    0.980909991279503, -0.635660782095407, -0.435437788675813, -0.49322591291444,
                    1.04820741274816, 1.00354322177165, -0.467462222800795, -0.185700764469914,
                    -0.422817740705227, 0.178738189684537, 0.738536890382022, -0.270420003574737,
                    -0.797962758959954, -1.24901314156597, -0.0131828779317971, 0.94874566062785,
                    0.696969100064875, 0.132907397821356, 0.93259971056429, 0.996603644577988,
                    0.96418428412909, 0.612612493416911, 0.46951415268046, 0.312656480200374,
                    -0.37913440806426, 0.0812740156543725, 0.684347713306423, 0.0639350227312718,
                    -0.152902517955241, 0.346564741980188, 0.0832943267284918, 1.42817452233597,
                    -0.437027907665964, -0.0035338070459261, 0.463866040686014, -0.39887629179708,
                    0.55000644397793, 0.0410712741374084, -0.302374272432147, -1.00559346523585,
                    0.257356942373548, 0.937504422399829, -0.883548573537695, -0.537674036905599,
                    0.318669222924711, -0.537708441409592, -0.00658119175179395,
                    -0.695340860396584, -0.65655062117672, -0.260436667819903, -0.27368612591286,
                    0.0656253089147589, -0.587008850778005, -0.97538904247763, -0.0943361902816946,
                    -0.0703330966743465, -0.229162015468642, -0.262589991475664,
                    -0.40081980934427, 0.10225207291009, -0.346182780827149, -0.661505440973187,
                    -0.201382684672382, -0.137122611208653, 0.75668697511955, 0.303229893175143,
                    -0.360635461658785, -0.307633077547877, 0.117717872438035, -0.427699992536674,
                    0.0994300907971422, 0.0306553206467781, -0.124931218750618, 0.156479426569831,
                    -0.121265875858342, 0.143268421954872, 0.497708246586053, -1.18695654015664,
                    -0.231280150198728, -0.38312272522532, -0.0594979402059721, 0.203344935712142,
                    -0.0256917483971222, -0.00841347561376438, -0.305419752682781,
                    0.192855790505035, -0.645012268469469, -0.29258711776896, 0.450360641202593,
                    0.455618026888601, -0.412034005425539, -0.671249630120423, -0.682853315040376,
                    -0.671594183360838, 0.0494436683972594, 0.303788667348584, -0.11781736711935,
                    -0.53599824578825, -1.02208167388119, -0.647055093120144, 0.672756338951046,
                    0.432024158611592, -0.641954330265946, -0.906978287465719, -0.394398924431925,
                    1.28083647935091, 0.93133517843935, -0.113466409745683, 0.246442636179907,
                    0.573712103075422, 0.328920512124256, 0.120177208730093, 0.660722344062947,
                    -0.130632785192418, 0.235602152620856, 0.360853522526241, 0.447099088485466,
                    0.218374589339959, -0.00497129358996062, 0.690267721951329, -0.182200912603605,
                    0.598821747008485, -0.459863607864776, -1.18449340660061, 0.334177169765269,
                    0.0175984705487053, -0.0177737829966385, 0.305418715778075, -1.11685284196363,
                    0.189365917231214, -0.36004114974719, 0.502597958261244, -0.365719694500115,
                    -0.935849742963612, 0.509626696348647, -0.228079212279624, -0.152945520405907,
                    -0.0709201831796804, 0.368813142994327, 0.996422784298574, -0.46004611227994,
                    -0.226792376225048, -0.217640989925371, -0.30614568956932, -0.732577143996735,
                    -0.363879783489722, -0.120873924204064, 0.270203333701858, 0.0120663802331222,
                    -0.242361850118101, -0.241474512399947, -0.315214051386843, -0.366312635550372,
                    -0.427888892196293, 0.0254853638338943, -0.0373911044989939,
                    -0.85951484665935, -0.176847275660435, -0.0382061230323387, -0.687705652033998,
                    -0.917190059694178, -0.0975409703900592, -0.188674087689066,
                    -0.342750277719986, -0.0728829171344501, -0.143349193101699,
                    0.210924750714606, -0.104333415186943, -0.190171691851404, 0.136066814861834,
                    -0.506308153684358, 0.261936780485096, -0.31497955204614, 0.269440547667884,
                    -0.336984543195825, -0.373433009188773, -0.937277063019102, -1.17164346433719,
                    0.61284380947653, 0.0135405857665489, -0.492317748620451, -1.36125134092924,
                    -0.57074259368261, 0.11729566464787, -1.59396490245172, 0.0110732467933587,
                    -0.510536182387919, -1.10001666966488, -0.982716578890077, -0.18771407181481,
                    0.856546413530367, 0.333716036051945, -1.19359716068593, -1.52303016970274,
                    0.854118510235367, 2.19497512685112, -0.520951658906001, 0.202762427183927,
                    -0.529828236207038, 0.725487837019926, 1.6428966614171, -0.0540291925966201,
                    0.164747598259453, -0.437880519073047, 0.591098970832746, 0.692078530808004,
                    0.239134845145897, -0.0143167004218014, -0.303930536111992, 0.41156832274055,
                    0.585843064397466, 0.571498656569575, 0.0177715227948454, -0.364493221332609,
                    0.298728478781186, 0.590730959742075, 0.982677536494784, 0.0367678293235087,
                    -0.986024829839039, -0.339415736186879, 0.0357538086110765, 1.00511945092666,
                    -0.257770782520453, -0.486263177894614, 0.469604769170104, 0.0412943109502066,
                    0.434764355352211, -0.146893286150345, -0.217562831425582, 0.0124705298027266,
                    -0.251261560462012, 0.180473589410303, 0.543925618684028, -0.814318789824824,
                    -0.879684988480449, -0.576928050386101, 0.125128256231631, -0.298408990763022,
                    -0.75479900687121, -0.62741800118215, -0.565968163140866, -0.203653750986455,
                    0.0865278995766232, -0.359031655849051, -1.14332163912951, -0.316369494478494,
                    -0.49817236986862, 0.229572227663311, 0.763758306601503, -0.693449084226017,
                    -0.777132297849761, -0.269037858959229, 0.0439959222209438, 0.519296136169572,
                    -1.38135268101826, -1.00908230433492, -0.122737175913244, -0.181282992508756,
                    0.453317152719066, -1.54246253455622, -0.852295291542333, -1.31998207543469,
                    -0.17502520858207, 1.15983850284842, -1.05170631242152, -0.315787307806346,
                    -0.992943047531823, 0.302585212788244, 1.25223465275648, -0.73516953762018,
                    0.096632957613167, -1.05055999362859, -0.267661527930642, 0.0529779327606974,
                    -0.768045387016599, 1.11927790949909, -0.470462509984847, -0.526566922532766,
                    0.619735065038569, -0.492685642038367, 0.496432184118589, -0.261689808495577,
                    -1.16213551896578, 0.18220399945595, 0.773042454076744, 1.49158765319083,
                    -0.178818641356081, 0.0399994097710703, -0.24794647015203, 0.830586640633356,
                    1.04613702862122, -0.298640797918645, -0.845157448884947, -0.164859887070151,
                    1.31663981032318, 1.03940350527913, 0.0291965527402906, -0.995639200008583,
                    -1.57673810800181, 1.10821172680545, 0.970532069275492, 0.394060547781673,
                    -0.871902742220072, -1.32237510262679, 0.54668040684615, 0.627231262668573,
                    0.648977009379778, -0.110004784441884, -0.867855774465838, 0.405349304007984,
                    0.366933042904719, 1.01219988425422, -0.433068397251986, -0.615921814779384,
                    0.198594352242509, 0.692579429174714, 0.0894916244331763, -0.410118538866442,
                    -0.820586599927026, 0.387229444446035, -0.323899071586425, 0.325789044391315,
                    -0.158104823415613, -1.11670800482142, -0.211823742951907, 0.14037664762057,
                    0.0408916797567596, -0.187925147470281, -0.391838755519927, -0.46169082237688,
                    -0.274391052828861, -0.207502918133262, -0.846400754719525, -0.0958952772959236,
                    -1.1418468944963, -0.478790544932622, -0.349563316802987, -1.03424508195678,
                    0.896611648924349, -0.455890696830352, -0.534242548247562, -0.879368386007039,
                    -0.887622367766767, 0.843575283258247, -1.16204989306968, 0.151350851113989,
                    -0.523647582273229, -0.567894389034876, 0.754923904135978, -1.16354844181637,
                    -0.689148536299854, -1.03817163349248, -0.469981109763921, 1.25520364955638,
                    -1.34672773833501, -0.67204046247852, -1.10334663674166, -0.00440965450547992,
                    0.308857988817309, -2.25315011862486, 0.227410523616129, -1.17576312978874,
                    -0.174127671525639, -0.944616330070659, -0.96133252535742, 0.786500208607038,
                    -1.34218475180993, -0.239684054305229, 0.194425586043096, -1.24261458011529,
                    0.47399237466206, -0.588065836630756, -0.818239058740687, -0.147665069586154,
                    -0.455977348752867, 1.21629961235001, 0.182379082629598, -0.0824607523425698,
                    -1.06166333853162, 2.61247952727695, 1.07940930954629, 1.02395362213806,
                    -1.42038634245628, -0.808362964066632, 1.70322272019612, 0.576360645470658,
                    0.976483777002288, -0.746072017835615, -0.720161729149609, 1.77835097917077,
                    0.700996392240228, 0.685673647812581, -0.0474293016568765, -0.641578015123415,
                    1.60991326095342, 0.749615122961529, 0.554976059877286, -0.268476849959001,
                    -0.164087105542606, 1.01558407271665, 1.30925620282304, 1.12835243476244,
                    -0.627933991233078, -0.429147257259558, 0.64434279761403, 1.21775123358066,
                    0.732479773387399, 0.234469437939993, 0.74695345644399, 1.31438042470805,
                    -0.0122882479316218, 0.921850366352572, -0.037818754098339, 0.206490336920098,
                    0.998911535149006, -0.329559206656882, -0.216112503151948, -0.602421246283967,
                    0.218063599905165, 0.968981635787401, -0.0623956244022281, 0.00264184384707344,
                    -0.456321831423111, 1.14538917703989, 0.2887265777442, 0.0210675994227927,
                    -0.22307699602135, -0.494716781057832, 0.655517561472536, -0.0961014316038415,
                    0.0630802783044828, -0.352564184213056, -0.860213496806324, 0.243307588809614,
                    -0.201480952364591, 0.349230157187463, -0.601216896763173, -0.90671046618084,
                    -0.0577550120819978, -1.55743355227025, -0.595426581649633, -1.34049361012017,
                    -0.929881580645144, 0.00963113002821103, -0.948524535469289,
                    0.253644007616202, -0.999157315366022, -0.147304223136558, -0.87436338579333,
                    -1.10461634831756, 0.45074969414361, -1.11271111610395, -0.496164903839173,
                    0.253372405811859, -0.836076875787711, 0.603185849301241, -1.62350195624571,
                    -0.401853892893143, -1.0207655805572, 0.0835889723630316, 0.736241481822276,
                    -1.10852468550472, -1.34111234611228, -0.845403183580444, 0.395735177447278,
                    1.16617663172525, 0.436316678678814, -0.0612144318511412, -1.1901569490015,
                    2.99793836781646, 1.72008370262785, 0.316212003831762, -0.262554976659601,
                    -0.248073087888088, 3.18363687255845, 0.491844442611978, 0.449766036811848,
                    -0.125042945900759, -1.30135709045473, 2.62378817593384, 0.516419945546659,
                    1.14450122253073, -0.00454490922760442, -0.569982374495567, 2.49119285302874,
                    0.481666197701571, 0.762550166442488, 0.340609639713932, -0.709206971926548,
                    1.6920948132875, 0.681252652431084, 1.02184647316707, 0.363305028381319,
                    1.14213769139182, 2.25080346615837, 0.50356487774431, 0.897650444587789,
                    0.841266518091118, -0.161004732949976, 1.46533571591084, 0.313636377901474,
                    0.329838435262984, 1.28710447011679, -0.722270846408861, 1.16957601089624,
                    -0.231895436665582, 0.254945815831293, 0.744447440922073, 0.695644813092585,
                    1.21759058405968, -0.288011700681551, -0.411084607481004, 0.548343184133472,
                    1.18098522712354, 1.25381695590704, -0.170284766316144, -0.0717001418724119,
                    0.0127863815490343, 0.654118691419939, -0.373066772129107, -8.81371868501901e-05,
                    -0.586923734259681, -0.269028643803153, 0.134111229778981, -0.571808915740917,
                    -0.0585121070431285, -0.911159651170185, -0.670342145324241,
                    -0.0595093438786449, 0.112093146148425, -0.394627591466332, -0.77811524658234,
                    -0.467781956171933, -0.0814553201188204, -0.489184014116546,
                    -0.00603526907923424, -0.745718459567602, -0.54164712909261,
                    0.501250482070219, -0.673891978971767, -0.158294995441508, -0.630854963248181,
                    -0.94793908921099, 0.162778243028342, -0.382905895732563, 0.0766649910713397,
                    -0.312356986204198, -0.468779954468175, -0.190061269937374, -0.0278234100711359,
                    0.75586894330786, -1.31545236411743, -0.588445726293145, -0.206217665340986,
                    0.439865999874257, 0.749750805358985, 0.432401023458954, -0.204417951723315,
                    -1.60941875044665, 1.99031392001027, 1.26806034603928, -0.151348760853132,
                    -0.540373312504641, -0.824311025496144, 2.54616321135842, 2.0810521372951,
                    -0.400041313898154, -0.669185696184769, -1.19489380831941, 2.99923463744079,
                    1.90780177188975, 0.553558147180901, -0.620630708162167, -0.635559724076131,
                    2.57189738717076, 1.24536707371055, 0.882250667753852, 0.12572029117274,
                    -0.658359734661431, 1.65471135059286, 0.589629864114897, 0.977697732663484,
                    1.46852198116764, 0.503521215329656, 0.976685394907037, 1.13468148783218,
                    0.513193179686864, 1.37102710896616, 0.447983082306436, 1.06413498528607,
                    1.01419270641577, 0.0456363380915259, 0.972596229368804, 0.167993980077654,
                    1.39148980239352, -0.133005931688302, 0.0577767661823718, 0.483580157243097,
                    0.598460793902007, 0.852571115107594, -0.632800157510189, -0.0079091626534642,
                    0.180946417556783, 0.149434981640056, 1.00844377137054, 0.238099275278586,
                    0.129811289031151, 0.552434464411817, 0.15627776045126, -0.364912371769179,
                    -1.03555153116056, -0.501279710511676, 0.343849026824557, -0.0571923589997994,
                    0.161136590533232, -0.118741796830733, -0.441037161631484, 0.0875778383957311,
                    -0.432580607771257, 0.531645461503672, -0.019248304570965, -0.391445887189874,
                    -0.00323421440947413, 0.460119964245592, -1.12790944628748, -0.286299825311374,
                    -0.797865499321479, -0.35382257735148, 0.486545671130088, -0.47824235961744,
                    -0.501616994094633, -0.51427898034429, -1.01958936142191, 0.802274949431578,
                    0.526321393644764, -0.278530648361823, 0.144607785911177, -0.893575023942219,
                    -0.134023046947306, 0.321064290252849, 0.520005022149661, 0.131366015842249,
                    -0.766596466001687, 0.212142249558696, -0.538536880740594), .Dim = c(5L, 520L)))



    expect_error(recover_alpha(out, n_message = -5), "n_message must be a positive integer")
    expect_error(recover_alpha(out, n_message = "aa"), "n_message must be a positive integer")
    expect_error(recover_alpha(out, n_message = TRUE), "n_message must be a positive integer")

    expect_error(recover_alpha(out, constraint = "adfads"), "constraint must be either \"unconstrained\", \"overall\", \"resolution\", or \"predicted\"")
    expect_error(recover_alpha(out, constraint = 1), "constraint must be either \"unconstrained\", \"overall\", \"resolution\", or \"predicted\"")
    expect_error(recover_alpha(out, constraint = "predicted"), constraint = "constraint = \"predicted\" is not currently supported -- developer note: add W_pred to function call to enable this in future results")


    class(out) <- "mcmc_mra"
    expect_error(recover_alpha(out), "out must be of class \"mcmc_mra_integrated\" which is the output of the function mcmc_mra_integrated()")

})


# unscale_estimates ------------------------------------------------------------

test_that("unscale_estimates", {
    set.seed(111)
    y <- rnorm(10)
    X <- matrix(1:20, 10, 2)
    params <- list(
        n_mcmc    = 5,
        n_adapt   = 5,
        n_thin    = 1,
        n_message = 5
    )
    locs <- matrix(1:20, 10, 2)

    priors <- list(
        mu_beta = rep(0, ncol(X)),
        Sigma_beta = diag(ncol(X)),
        alpha_sigma2 = 1,
        beta_sigma2 = 1,
        alpha_tau2 = 1,
        beta_tau2 = 1
    )
    inits <- list(
        beta = rep(0, ncol(X)),
        sigma2 = 1,
        alpha = rep(0, 520), ## use M = 2 and n_coarse_grid = 4
        tau2  = 1:2
    )

    out <- mcmc_mra(y, X, locs, params, priors, inits = inits, M = 2, n_coarse_grid = 4)
    expect_equal(
        unscale_beta(out),
        structure(c(-5.95373792478096, -0.516903430166984, -2.10956374385706,
                    -0.711234895826659, -4.6983160529007, -0.34740352232307, 0.0145838680998062,
                    -0.101566142662089, -0.00604333507590154, -0.261082656840876), .Dim = c(5L, 2L)))

    expect_equal(
        unscale_sigma2(out),
        c(1.953453347397, 1.07518264163489, 1.24130939364163, 1.0427615248997,
          1.35807561190342))

    out <- mcmc_mra_integrated(y, X, locs, params, priors, inits = inits, M = 2, n_coarse_grid = 4)
    expect_equal(
        unscale_beta(out),
        structure(c(-0.669013535873582, -0.669013535873582, -0.669013535873582,
                    -0.669013535873582, -0.669013535873582, 0, 0, 0, 0, 0), .Dim = c(5L, 2L)))
    expect_equal(
        unscale_sigma2(out),
        c(0.894666104064748, 0.893605112882278, 0.869138154172264, 0.717748384620361,
          0.602167601575069))

    class(out) <- "aaa"
    expect_error(unscale_beta(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
    expect_error(unscale_sigma2(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
})


