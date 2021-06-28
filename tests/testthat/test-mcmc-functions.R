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
        alpha_sigma = 1,
        beta_sigma = 1,
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

    priors$alpha_sigma <- rep(1, 2)
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma in priors must be a positive numeric value.")
    priors$alpha_sigma <- -1
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma in priors must be a positive numeric value.")
    priors$alpha_sigma <- "aa"
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma in priors must be a positive numeric value.")
    priors$alpha_sigma <- 0
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma in priors must be a positive numeric value.")
    priors$alpha_sigma <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter alpha_sigma in priors must be a positive numeric value.")
    priors$alpha_sigma <- 1

    priors$beta_sigma <- rep(1, 2)
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma in priors must be a positive numeric value.")
    priors$beta_sigma <- -1
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma in priors must be a positive numeric value.")
    priors$beta_sigma <- "aa"
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma in priors must be a positive numeric value.")
    priors$beta_sigma <- 0
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma in priors must be a positive numeric value.")
    priors$beta_sigma <- NA
    expect_error(mcmc_mra_integrated(y, X, locs, params, priors), "If specified, the parameter beta_sigma in priors must be a positive numeric value.")
    priors$beta_sigma <- 1

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

    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_sigma = "TRUE")), "If specified, sample_sigma must be TRUE or FALSE")
    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_sigma = NA)), "If specified, sample_sigma must be TRUE or FALSE")
    expect_error(mcmc_mra_integrated(y, X, locs, params, config = list(sample_sigma = 3)), "If specified, sample_sigma must be TRUE or FALSE")


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
        structure(c(-92.5571186738049, -151.095926637375, -117.466089660476,
                    -30.2685748577414, -346.718350268487, -92.1224463809536, -152.76099512319,
                    -117.77660483886, -29.0377487712637, -346.722966870522, -91.7439528057472,
                    -152.81279696138, -118.041665041661, -30.3553389491394, -347.39539720861,
                    -92.5040446623465, -152.82284843874, -118.385915121575, -29.9809025218088,
                    -347.107269266711, -92.5996435272952, -152.8562559597, -118.148453814749,
                    -30.1810371741872, -347.8229323295, -92.6687485909039, -152.912155590648,
                    -118.306155617692, -29.8751991814534, -348.168018096208, -93.2364992628706,
                    -153.293617319634, -117.661141089327, -30.6655560120206, -347.878932289819,
                    -94.1369909452893, -152.900819188572, -117.731576949571, -31.0453744914782,
                    -347.99629036776, -94.106442804318, -153.312146404338, -118.373761650414,
                    -30.5622425760047, -348.604888487889, -94.397309854117, -153.869296777732,
                    -119.005976951128, -30.0025323454183, -348.977360657177, -93.5871503458946,
                    -152.435059305804, -117.584817782964, -30.7984353229188, -348.609160497612,
                    -93.5882164227514, -151.482270480518, -116.920116191821, -31.298190714304,
                    -348.113908551794, -94.2077274500073, -152.521887841438, -117.903705866978,
                    -32.0198296885061, -348.537171454087, -94.4579165473642, -151.725155561725,
                    -116.968188718744, -32.1093405252799, -348.607562001371, -93.3306509999713,
                    -152.28098124747, -117.032665385647, -30.1330309917355, -346.400315265161,
                    -92.6332350846113, -151.942583530685, -117.511955921028, -29.4100153179587,
                    -346.922099460581, -93.20887329997, -152.70175016561, -117.15169503438,
                    -30.6143613664282, -347.156841683271, -93.3185142449063, -153.032930139752,
                    -117.98981156366, -30.0801648503645, -346.941146909175, -92.8661748966044,
                    -152.808497294209, -118.395864871188, -29.9442410958716, -347.77255277127,
                    -93.4666050429784, -152.67070613502, -117.915410417357, -30.1915046462778,
                    -347.468527158723, -94.0487473968152, -152.991341750722, -118.211621618979,
                    -30.7410384532203, -347.915780701446, -94.4912649417034, -153.34172766878,
                    -118.656859317722, -31.4221727458085, -348.223262602481, -94.0043158760031,
                    -152.944534970304, -118.227656066715, -30.4064876303905, -347.407895098249,
                    -93.7005232250487, -152.715953974735, -118.364100886819, -30.1305888025497,
                    -348.665830747635, -94.0010536402063, -152.811303281105, -117.706436877065,
                    -31.1061725807496, -348.646337238406, -94.4947409904472, -152.325874759845,
                    -118.05460766521, -31.5409169186747, -347.865360367361, -94.7839250933935,
                    -153.02227612161, -117.412814248983, -31.5130883333914, -348.380336089704,
                    -94.575404239923, -152.25171595773, -117.841528241454, -31.7032978915093,
                    -348.025186413733, -93.5380304940109, -152.205470701388, -118.06902920851,
                    -31.3156195753736, -347.081067364322, -93.0979812367722, -152.556869201623,
                    -117.120699747226, -30.3802828254381, -346.893785178949, -93.7537006483559,
                    -152.270034892744, -117.597095187212, -30.4286060180712, -346.997100897986,
                    -93.6471742181251, -152.389803570047, -117.725843300755, -30.100341838624,
                    -346.933734459268, -93.549553708346, -152.676339960715, -117.934134223347,
                    -30.3632475645679, -347.384520275119, -93.7313488896701, -152.929784316285,
                    -118.009771657053, -30.8614400811037, -347.890550935865, -94.3278197222212,
                    -153.554285115868, -118.066204447853, -30.5965996812432, -347.490972977416,
                    -94.5514941058422, -152.750650047766, -118.663279670845, -30.1784691644656,
                    -347.878988006423, -94.0538783175738, -152.450906237652, -117.952963542407,
                    -30.0542049429378, -347.438396941287, -94.3663232342725, -152.448387766153,
                    -118.191283727594, -30.6323507651803, -348.119228780788, -94.9244435452419,
                    -152.212742969889, -118.46907920831, -31.468971048469, -347.871072102792,
                    -94.4822604082243, -152.456310712097, -118.169963219528, -31.5132035512287,
                    -347.90784041498, -94.8819740796419, -152.403798029277, -117.916210215755,
                    -32.4875891840439, -347.965228818548, -95.4977883388265, -153.2057130609,
                    -118.446869714803, -32.756703270792, -347.521946639348, -93.706149804711,
                    -152.375119674142, -118.429950598002, -29.7827353449148, -346.83020022431,
                    -93.3453545308877, -152.714510438979, -117.896847267533, -30.4585735180351,
                    -346.221454573977, -93.0172411692128, -152.863026752263, -117.462009290125,
                    -30.1018729136923, -347.868634381482, -93.8247642167063, -152.921514974294,
                    -116.57708521341, -30.823638640917, -347.412580787982, -94.2704618657903,
                    -152.900980679266, -117.784675567146, -31.2812876614171, -347.137749177756,
                    -93.1925676065335, -153.188540175394, -118.323784871888, -31.0769224808183,
                    -347.51686996739, -93.6317631085901, -152.630100294921, -118.251833350694,
                    -31.1375732754864, -347.488526826176, -94.1949719958927, -152.797696853932,
                    -118.357159796141, -30.523735419539, -347.37364977134, -93.4953922401089,
                    -151.879362950438, -118.555547668544, -29.9684974768163, -347.771836768,
                    -94.242098434051, -152.382208563925, -119.251895693037, -30.5971012592379,
                    -348.20979781175, -94.9253009950015, -151.799326162275, -118.781381741845,
                    -31.7163392159527, -348.142097771923, -95.116360741491, -152.5431881159,
                    -118.135706484026, -31.94376077956, -348.233268020969, -94.9526221794632,
                    -151.960423718837, -117.95190898911, -31.9262497659884, -348.68818570033,
                    -95.476040551346, -152.21892429772, -117.438762726745, -31.7971159572529,
                    -348.306367763273, -93.6617497358909, -153.168300699281, -117.441919013839,
                    -28.9223555909944, -346.871811655731, -93.5972341388608, -152.678863187022,
                    -117.191828512996, -30.3806103841729, -346.520114241627, -94.1108100111948,
                    -152.566283015785, -117.846318287796, -30.528475394816, -346.943399316046,
                    -93.9501292843511, -153.001533228238, -117.066734286243, -31.061346596814,
                    -347.546993857195, -93.6956174150396, -153.264973851703, -116.846443114676,
                    -31.1771621097408, -347.332698431414, -93.0206556914473, -152.30575792596,
                    -117.700833969793, -31.3509145259761, -347.234465170926, -93.1094957629042,
                    -152.687992662174, -117.973269057793, -31.6199672004365, -347.716991750086,
                    -93.7306390695553, -153.5897315142, -118.438417723859, -30.7904117863963,
                    -347.454762432263, -93.3114855102783, -152.134771579586, -117.556404573299,
                    -30.704191003123, -347.069814157959, -93.7418973386303, -152.027440682272,
                    -117.689253588295, -30.6826041774034, -346.749406743771, -94.3635084477434,
                    -152.86414225757, -118.574058266088, -30.9954000976388, -347.533760598998,
                    -94.2321527745501, -152.950665688575, -118.078346213221, -30.7534467208113,
                    -348.327926746438, -95.203245633457, -152.294466448642, -117.981748521475,
                    -31.2169847627361, -348.180856910005, -94.9487449176093, -153.205610592077,
                    -118.73496520117, -31.3971933497027, -349.169314591394, -94.2317820892293,
                    -153.734920948864, -117.648807969011, -30.1528665111485, -346.704154404137,
                    -94.2556314812002, -152.205627405823, -117.857589122852, -30.3644043963574,
                    -346.772918302071, -93.6294010096917, -152.86183337878, -117.411731372187,
                    -30.7890274843056, -347.406509396416, -93.4254431340381, -152.513240508001,
                    -117.329266017543, -30.97894497928, -347.657507455302, -93.2232543047343,
                    -152.400376295126, -117.016982791255, -30.9148313652103, -346.536864184269,
                    -92.9811690711902, -153.412726321244, -117.170236583055, -30.7352573933799,
                    -346.929557584924, -94.2355538906206, -153.658511757231, -117.663794545098,
                    -30.8222572559081, -348.432612209164, -93.6709478832287, -153.485120550288,
                    -117.941726198621, -30.8014342069455, -347.903321171142, -93.5285767692955,
                    -153.050872960528, -117.501080690567, -30.6814858930364, -346.978089892975,
                    -94.531198061517, -152.641297869172, -117.285328388731, -30.4719268110903,
                    -347.416643682308, -93.8522380496532, -152.608781472061, -117.691629109309,
                    -30.9836199258945, -348.216095007879, -94.250577385715, -153.158267164764,
                    -117.871186539834, -30.8911927684359, -348.662015112577, -94.4239495269001,
                    -153.184465990272, -118.004497654624, -31.1312698815767, -348.392745300751,
                    -94.0071973918583, -153.777534056813, -118.005539885623, -30.2893826758745,
                    -349.204775962202, -92.6538880872137, -152.881092877064, -118.452314103202,
                    -30.8529741678834, -346.708452532466, -93.5005650091935, -152.553775363307,
                    -117.524134828167, -29.9620672625444, -346.652181559687, -94.4688212149932,
                    -152.513851991172, -117.583239521464, -30.7794373265001, -347.646781732896,
                    -93.7738064512408, -153.37704835708, -118.123764872565, -30.4728194641736,
                    -347.387644292099, -93.6323964418061, -153.337895564703, -117.254807123955,
                    -31.2489692357136, -347.362846102902, -94.1105753900257, -153.603306021788,
                    -117.591734885869, -31.2193521243439, -347.194944528891, -94.9899656015638,
                    -153.808627265249, -117.223432672207, -30.7776988445894, -348.102753747747,
                    -94.2712951755912, -154.079205111465, -116.9180419729, -30.9958658428503,
                    -348.584928806245, -94.037921535209, -152.666220237725, -117.515798107271,
                    -30.4899264662505, -347.974163253631, -94.3833689304373, -152.931253530618,
                    -117.934193667011, -30.984912617924, -348.474468593076, -94.1836467578111,
                    -152.350353042162, -118.342232213717, -30.8153258786029, -348.15897304939,
                    -94.0733515648989, -153.939828871895, -117.899161893194, -31.4445363954011,
                    -347.985522008173, -94.5026759270635, -153.648584635623, -118.345622966083,
                    -30.4265764806987, -348.913131394397, -93.9925086978916, -153.895257035927,
                    -117.774738540376, -30.4293897409586, -349.043925558926, -93.4758003097593,
                    -152.558135556554, -118.173065404049, -30.2950341129854, -347.17636599553,
                    -93.4446338765059, -152.912385305698, -117.918897939145, -30.3912319332617,
                    -346.680940484026, -93.3180816659154, -152.155773713896, -117.601722988669,
                    -30.459136754795, -346.971427249181, -93.8335053290777, -152.995178085287,
                    -117.622245318817, -30.917065488966, -347.326076371494, -93.5600591958495,
                    -152.395160054486, -117.100920712709, -31.8557637374231, -347.684522956245,
                    -93.4188859478002, -152.537393366655, -117.387173984317, -31.9116733496511,
                    -348.019747369257, -93.964500200215, -152.891931507505, -117.640596112398,
                    -31.0697442325028, -347.70528943516, -94.6826886168424, -153.296756021442,
                    -117.091218917471, -30.6893813383955, -348.408967624162, -94.3103266031987,
                    -152.929266054496, -117.615280664507, -30.7454550877891, -348.313122042429,
                    -94.1580730120757, -153.22103308333, -117.790141865977, -31.2160604952956,
                    -348.737596216959, -95.0146126450161, -152.771600977197, -117.63599007262,
                    -31.7898952756954, -348.996154452793, -94.2761343933877, -153.264411161614,
                    -118.095885435673, -30.9659693857034, -348.101397993832, -94.3404096868306,
                    -153.94367432907, -118.35311351363, -31.5283312391793, -348.886135054314,
                    -93.9872591126169, -154.205831165764, -117.90639421572, -30.6426952933715,
                    -349.338503268929, -94.1773250839022, -152.277789749538, -117.998000573534,
                    -30.4552163090555, -347.428905546544, -93.5146852838224, -153.341670545812,
                    -117.748506754183, -30.5462849860032, -347.199102628732, -93.6674808906083,
                    -153.011843119563, -117.782527542031, -30.6287014365625, -347.578447362443,
                    -93.9666607093209, -152.785892558701, -118.300470530279, -30.206249366216,
                    -348.032243957723, -93.3019078367543, -152.651573861649, -117.797882980146,
                    -31.4938294712027, -347.838639025834, -93.4903030796757, -152.361225340971,
                    -117.895282027373, -31.5588052223771, -348.158828024029, -93.7602579259653,
                    -152.15507773637, -117.587828873602, -31.0049915165488, -348.007350208717,
                    -93.9630343621607, -152.379944149935, -116.636402204539, -30.5278251694701,
                    -348.185396453306, -93.8308161015924, -152.912646841071, -117.62320122291,
                    -30.7731879959967, -348.551980767615, -93.5763906982318, -153.463658694983,
                    -117.627037984892, -31.1484179858439, -348.637439759955, -93.6658850856336,
                    -153.514700189579, -117.384171340238, -31.4609826273791, -348.577408978771,
                    -94.5848987100896, -153.084011452441, -117.732243195406, -31.4052024856073,
                    -348.249712151724, -93.977606178448, -153.637095291205, -118.081597423911,
                    -31.4384925885304, -348.665621438042, -93.9994350037871, -153.992040591279,
                    -118.631512620535, -30.9089806816399, -348.766615090811, -94.0538892706922,
                    -152.780287560809, -117.738350294311, -29.7949999842308, -348.707085062758,
                    -93.8953143216533, -152.755829681718, -117.820079164996, -30.2454089113154,
                    -347.947452117424, -94.3934864280695, -152.717804264724, -118.008553251829,
                    -30.9938746047888, -348.330493434718, -93.8863806500143, -153.543267794639,
                    -117.969838084388, -30.6198785053434, -347.995011011757, -93.4951669427832,
                    -152.328550133823, -117.768130771193, -30.6937293695009, -348.028330745257,
                    -93.8483515391992, -152.60876383677, -117.861461650821, -30.933251175682,
                    -348.436065749893, -93.7799054002288, -152.316933910318, -118.244266809258,
                    -30.4590789510567, -348.334872967557, -94.395702652009, -153.284382578856,
                    -117.395374775777, -30.4504255011144, -348.506279824445, -94.5833796296614,
                    -153.124445514894, -117.468382647942, -30.7605666827013, -348.56871769708,
                    -94.204999049181, -153.680643601976, -117.125154021013, -30.6638086560082,
                    -348.991567704653, -94.457384333834, -153.513260674335, -117.598993908615,
                    -31.7666305998557, -348.324419769385, -94.2524125128222, -153.132569900889,
                    -118.169025381961, -31.238842295634, -348.410807328711, -94.7743199294981,
                    -153.895301191056, -117.68822536423, -31.2634811154688, -348.892197866426,
                    -94.2525361345336, -153.019981648269, -118.370022447394, -30.5987846362483,
                    -349.417739689694, -93.8397695048726, -152.928554324524, -117.666632237526,
                    -30.0684742581498, -348.013162956096, -94.0171769070164, -152.785315101087,
                    -118.428802230489, -29.987376183263, -348.358450206687, -94.6528448410646,
                    -152.897645264265, -118.568448872977, -30.8722639916502, -347.973022416248,
                    -93.6756295372108, -153.345498012165, -118.140943403349, -31.0831929118049,
                    -348.717376957469, -93.3559906464248, -152.757661955213, -118.775799510458,
                    -30.8202935542833, -348.597611865979, -93.6257426305428, -153.06906945639,
                    -118.597950320764, -31.0565000265246, -348.501408498101, -94.2502810757072,
                    -152.904340293197, -118.604658972285, -31.562237234292, -347.947276716356,
                    -94.1588440383627, -152.896221821701, -117.861913752355, -31.2861671022398,
                    -347.834839556221, -94.2886343062655, -152.869653379222, -117.68860037303,
                    -31.7341072223584, -348.740214575188, -94.1813483237276, -153.279599976043,
                    -117.629338664068, -31.370302658724, -348.972078433978, -94.3474165456486,
                    -153.038052006123, -117.756727596984, -31.2014914496214, -349.099947269549,
                    -94.439279067959, -153.173621019984, -118.336780192335, -31.0923969421756,
                    -349.116533419653, -94.3372881225154, -153.138033133976, -118.321705051207,
                    -31.2407687313845, -349.191740608592, -94.1006898328786, -153.017660855245,
                    -118.331273220996, -30.4878878938976, -349.079681620855, -94.1374235450744,
                    -153.191765078427, -118.51847301835, -31.4628402416365, -348.157138366429,
                    -94.0194690832404, -152.871345119977, -118.718468460771, -31.3788089767291,
                    -347.858000544591, -93.4377786109768, -152.426730169223, -118.06647353877,
                    -30.722008014937, -348.47634832762, -94.4231324859761, -152.882840399046,
                    -118.750984630399, -30.8452784321413, -347.675195470002, -94.1851963997101,
                    -153.007359440192, -118.537082484681, -30.7472025684968, -348.694317697486,
                    -93.6459971169579, -153.289836743536, -118.754784573414, -30.8069313560489,
                    -348.403674133359, -93.8761054898816, -153.392740188889, -118.316677560032,
                    -30.7718235188802, -348.788203385486, -94.061414213219, -152.83162667673,
                    -117.99073928739, -31.8386774853688, -348.622406126916, -93.5523324402712,
                    -153.138979152439, -118.177927970372, -31.3644556580939, -348.835774166822,
                    -93.8993626214312, -153.106459761393, -118.076997240211, -31.1826316509329,
                    -349.234381824689, -94.3959108748889, -152.741954690268, -118.846279847974,
                    -30.6113777386246, -348.083844917953, -93.7042034988561, -152.534580707482,
                    -118.179300680375, -31.3542398610169, -348.813985025677, -94.1187963764902,
                    -153.099037720797, -117.748462706132, -31.0723742451312, -349.290610054954,
                    -93.3215999556696, -152.760660370963, -117.809838491613, -31.0663992766417,
                    -348.409763210746, -93.5819801796913, -152.306294245884, -118.705960842814,
                    -31.5162986023824, -347.345215747273, -93.0159491481985, -152.165893194525,
                    -118.698102461941, -30.5919626571502, -347.970682899004, -93.224083893578,
                    -152.730144969129, -118.124903200172, -30.5785437795759, -348.180054402207,
                    -93.268235644994, -152.718100591955, -117.728652191867, -30.9573311264571,
                    -347.71813776941, -93.7911137894453, -152.800068347506, -118.707857697637,
                    -31.5496301696341, -348.464391562099, -93.6448325886443, -153.341991102074,
                    -118.91239215313, -31.4694608598999, -348.859491522567, -93.6287983177923,
                    -152.662773969172, -118.465580945886, -31.4278514160825, -348.249987184264,
                    -93.9885954678046, -153.009511643227, -118.68982459849, -31.440925122874,
                    -349.372990628541, -94.0201614626025, -153.082496602103, -117.609770237988,
                    -31.7112807639204, -349.445837348851, -94.5149731478094, -153.326913422402,
                    -117.379789348612, -31.3092696236943, -348.670790810806, -94.3885755938555,
                    -152.445377356505, -118.623303875116, -30.069364753849, -348.856361028696,
                    -94.1485954657239, -152.764446897144, -117.959318818528, -30.5311568990971,
                    -348.981760125009, -93.9588676071444, -152.59020821865, -117.703057122074,
                    -30.9942972371238, -348.778066230952, -94.0805858521375, -152.882371350008,
                    -117.723919704237, -30.8896403032028, -348.49352339084, -93.8143946449072,
                    -152.211351271902, -118.490361622272, -31.7537127644641, -347.198662503112,
                    -93.3499849323329, -152.028238390196, -118.445270913274, -31.1506799060713,
                    -347.560262034284, -93.3577018247597, -151.891574460063, -118.472436506048,
                    -30.329017289604, -348.924381511923, -93.0248351048636, -151.928349278571,
                    -118.293833195956, -31.1375550608435, -348.270336816424, -93.4731073303893,
                    -152.67174622379, -118.598217436541, -32.2074320368365, -348.443011994366,
                    -92.4558738891444, -152.144559521086, -119.852962177215, -31.8263250866051,
                    -348.519010416296, -93.4224074976529, -152.44144782634, -118.798020959411,
                    -30.996163373, -348.307298958629, -93.8182524871698, -152.809431106556,
                    -118.608697789746, -31.4906635439172, -348.748551406158, -93.7898714677998,
                    -153.355804311948, -117.60152657117, -32.0726138546277, -349.49550645986,
                    -93.6628338093796, -153.410180753727, -118.031841278258, -31.4730870484153,
                    -349.005574461541, -93.7098810461039, -152.467888169436, -117.685028151683,
                    -30.3475778213604, -350.054804738725, -93.3567633666173, -153.068748608072,
                    -117.559023564358, -30.3389237713901, -348.698833036988, -93.3576230239866,
                    -153.042460862576, -118.135236474203, -30.2134902696012, -347.397635542451,
                    -93.7070419221857, -152.980360116233, -116.81727898202, -31.858352931826,
                    -347.744807550045, 69.7752849175067, 113.232096654962, 87.0023766134501,
                    21.7461154940217, 255.204226552191, 68.8814993473878, 112.633100910214,
                    86.8762152698827, 21.6809359427671, 255.059187493185, 68.6551565066297,
                    113.59015841142, 87.1652804858403, 21.9098624305583, 254.833079573775,
                    68.7295525147179, 112.931199522362, 87.4260574416111, 21.3628140359476,
                    254.789683339626, 69.0264145818632, 112.713785808872, 87.402151429758,
                    22.6848655395471, 254.557095339589, 68.7427642822949, 113.148597392924,
                    87.0454265488931, 23.3356931781075, 255.883696573479, 68.5830286202713,
                    113.187453995248, 86.2327927534229, 22.8645167686344, 255.751457265094,
                    69.0262615417214, 112.779393091349, 86.4417741718197, 23.1032172326954,
                    255.064705067458, 68.6268130327518, 112.801274868921, 86.2053132365006,
                    23.2823398554052, 256.635483139504, 69.154002653756, 112.628101406807,
                    86.2405426921775, 22.6637899817422, 257.058827903609, 69.4504919494537,
                    112.606429740973, 85.7348563115391, 22.4718435982741, 255.318369141821,
                    69.4079045916923, 113.595895651647, 85.4044106810365, 22.1279937817863,
                    254.987711724423, 68.8001024417125, 112.832913994696, 84.8180140352139,
                    21.9539006338474, 255.335815450183, 69.2261837268336, 113.485592616643,
                    85.3047447739661, 21.9020457758448, 254.302214818203, 70.0286153795313,
                    112.741206591171, 85.3176727407407, 21.8176550992384, 253.906521660217,
                    70.3109551603661, 112.952128088596, 84.6559070656605, 20.5318654778329,
                    254.79341396275, 69.5992398701324, 112.679938545988, 85.1798798730133,
                    20.6701587596629, 253.321904084947, 68.6164238615423, 112.483143370145,
                    85.0019106157728, 21.069076373926, 254.362772314663, 68.9212513896476,
                    113.835244053375, 86.52931210759, 21.8082598645266, 254.625308888204,
                    68.3176677851772, 113.118994211597, 87.3009635639692, 21.6203667245636,
                    255.056607242647, 68.2799880089693, 113.362898145408, 87.9299065274187,
                    21.413811775714, 254.745229240814, 68.7791811474869, 112.728178799143,
                    87.1815448787912, 21.9388277126885, 254.486806397076, 69.57421103721,
                    112.521108944668, 86.9227590935414, 22.9015979173312, 254.635877992661,
                    69.1487997532414, 112.825724141587, 86.6789557998876, 23.7342327514955,
                    255.59595752953, 68.6174723753525, 112.717026202146, 85.682551147808,
                    23.1265947547277, 255.489995679338, 68.8973541582386, 113.39559979484,
                    86.6850395116035, 23.8480976544768, 255.61580699997, 68.8873286507514,
                    112.717318349981, 86.4356300973564, 23.1064654245427, 256.429350621719,
                    68.6575770331711, 112.512531811721, 86.0853594566321, 22.5230588412254,
                    256.840670898838, 69.6757787215651, 112.416821957574, 85.5510807782268,
                    22.7891432288888, 255.823756447675, 70.1498308767198, 113.291671236264,
                    85.2767711892153, 22.2548407206922, 255.414516017877, 69.8094342826785,
                    113.438077619518, 85.3990199166369, 22.1013794620135, 254.605800258504,
                    69.7080269435996, 113.532363079838, 86.0311978296435, 22.0302427666523,
                    255.370000266647, 69.6355047993065, 113.591822910688, 85.9695150907869,
                    21.4941460725068, 254.567034474015, 70.4447239548029, 113.348519230918,
                    84.1814187489658, 20.8030252429442, 253.960575236642, 69.9830987039329,
                    112.433526355123, 84.1087012889817, 21.8144667449903, 253.637463205838,
                    69.1677621548147, 111.862880273354, 84.3273664834175, 21.7666695546331,
                    254.333845380765, 69.094997470708, 113.60088490514, 86.6285570886064,
                    22.1279491808675, 254.959381455861, 68.6361297910765, 113.54421184438,
                    87.7720841113183, 21.3699938979497, 255.274528007168, 68.6721718979179,
                    113.261161298472, 87.8702298935243, 21.5712636362748, 254.439006584955,
                    68.6434731239237, 113.748789194688, 87.7100509701604, 21.6521219100518,
                    255.36123014799, 69.7472973548198, 112.508599334961, 87.1531937439353,
                    22.6139629861724, 256.18960715627, 69.0991399744518, 112.152176183499,
                    87.1311528431972, 22.8051412738727, 256.469516189067, 69.1953423623541,
                    112.432104314156, 86.9324482394386, 22.9243204643421, 255.925483439841,
                    68.3700314404269, 112.861430058824, 86.2246535650907, 22.8280301549319,
                    256.358005506271, 68.7651285221974, 112.389785174954, 86.6065833321195,
                    22.5085950681413, 256.480787392973, 68.9251267327782, 112.467993038413,
                    86.1639250259952, 22.4507344395374, 257.213343522557, 68.7484320383314,
                    112.152110929983, 85.7928126702086, 22.1801982130836, 255.758917046538,
                    69.1055797745922, 112.755174193251, 85.8616063180229, 22.6966040793417,
                    255.522679540601, 68.8693470299805, 113.348214936907, 85.8456900295562,
                    22.7485531364933, 255.515976804644, 69.4288012541402, 112.823290792179,
                    85.6391430969344, 23.2531678551829, 255.88575568726, 69.1882874527181,
                    112.766651824928, 86.2433436618263, 21.7842654654226, 254.496990372143,
                    69.1180203590101, 113.208582257596, 85.0063480803255, 21.5658597792921,
                    254.322681128413, 69.4021396553958, 112.440822918618, 84.434115971237,
                    22.5919499559069, 254.087711032119, 68.8695992898458, 112.358085680794,
                    84.5980260065737, 22.9596380781645, 252.878459318185, 68.3921096565028,
                    112.871293679753, 87.3456021948729, 21.481726909911, 254.520420427944,
                    68.4370602124828, 113.356118702907, 87.1007568618209, 21.9162303499353,
                    254.50905443619, 68.5699734530592, 112.769624187779, 87.630574077719,
                    21.5187169291768, 255.167150902488, 68.8758925273343, 112.868510712834,
                    87.3598429448527, 21.5142067603883, 254.752732528548, 69.289034091416,
                    112.07382578137, 87.0378113770609, 21.7995223172893, 254.681368202915,
                    69.1522909849061, 111.8448122619, 87.0518688132963, 22.3231098564185,
                    255.341749027672, 68.6831530488855, 111.51147353313, 87.4313076948014,
                    22.7378093203919, 256.294264737628, 67.9678917648079, 112.44365702135,
                    86.1550463336239, 22.8722070680568, 257.018372795283, 68.7150396837,
                    112.380219912607, 86.2994976769428, 22.5451576728058, 256.514388502968,
                    68.8639181694688, 112.40383699799, 85.9297598243679, 22.1435980736813,
                    256.874398102964, 69.0874909816211, 111.909715063677, 86.2015805696842,
                    22.1840387650339, 255.648843400276, 68.9487264327114, 112.600502269528,
                    85.7087144199847, 22.170925525751, 256.036654723091, 68.5810517016067,
                    112.920013660588, 85.8646653690745, 22.6289984627778, 255.054709533508,
                    69.6197243014517, 112.050743589305, 86.5239755556611, 22.5726581148495,
                    254.923529396081, 69.6707221310724, 112.471513974668, 85.8906878629606,
                    22.091005213287, 254.348491770956, 68.4795650891373, 112.943006336416,
                    85.0066704544444, 22.6240841203868, 254.6747130304, 69.0990157448257,
                    113.074126411902, 83.6824412717374, 22.8816777513884, 254.335252296499,
                    68.7726532958244, 113.477923066003, 84.0026232894484, 22.844500956102,
                    253.674030778314, 68.5608475911948, 112.410468886134, 86.2772991171518,
                    21.4793315760148, 254.93245616732, 68.8796771933927, 112.505101150645,
                    86.156242144814, 22.0967726595471, 254.182418188926, 68.5294329342745,
                    112.898779356706, 86.9773381575852, 21.3801537181084, 255.052013683247,
                    68.9873136965812, 112.809059273518, 86.7878775352841, 21.4264423377225,
                    255.239722028222, 69.2341912304071, 112.597399641614, 86.8084338901211,
                    21.7368567969717, 255.368696160216, 68.9921076026789, 112.673349151685,
                    87.3198942440206, 22.004895501021, 255.360105592857, 68.8775402216298,
                    111.669738325683, 87.0157083715137, 22.8468769901007, 255.773517031234,
                    68.5673635886865, 111.825350027683, 86.3667788356022, 22.8148151235038,
                    255.958699863814, 68.4703869548414, 111.893129667, 86.4378694730489,
                    22.5993586578284, 256.460948984257, 68.3894230576382, 111.621153929131,
                    85.9962827576855, 21.8177284155697, 256.415934237502, 69.0149162759525,
                    112.081203026336, 86.0302088300894, 23.3976289708962, 255.201564308681,
                    68.4166209741774, 111.926478604612, 85.7723635293898, 23.014107237145,
                    255.89919764485, 69.0572172513616, 111.994518062444, 85.54626604361,
                    22.4873605227143, 255.893446219213, 69.8681936118366, 111.872808505286,
                    86.0642214424403, 22.4371898727298, 255.460050875918, 69.5389226434193,
                    111.88813811346, 86.1578789192592, 22.9627729917221, 254.772173784497,
                    68.7741433594816, 111.929309668876, 86.1108322983486, 22.5516317227244,
                    254.568832863631, 69.4785516274709, 112.215374932347, 85.5278184571651,
                    22.8335641226299, 253.961294167418, 69.2818450104104, 112.957388511452,
                    85.3783704040249, 22.8280058602075, 252.779152208426, 68.4425431014122,
                    112.7956118158, 85.9517015453481, 23.1215301894809, 253.958718130244,
                    68.4401711792392, 112.62259082355, 86.4687475915945, 22.481564730501,
                    253.857134717206, 68.3950774477579, 112.614998158964, 86.3227110213842,
                    22.203530174241, 254.817664740634, 68.8604416172187, 112.457922214703,
                    86.3205607241218, 22.2337338327927, 254.643849311044, 69.4057411546591,
                    113.129755079047, 86.354183239565, 22.3205248449929, 255.799386137459,
                    69.4666375072295, 113.354180713261, 86.4932636038604, 23.2525087032204,
                    254.874711803221, 69.0220807234889, 112.650134778995, 86.4231416940719,
                    22.7924388888712, 255.320512141983, 68.7531514534245, 112.980942081193,
                    85.8281763627463, 22.6787849839222, 256.301144425259, 68.8880532120238,
                    112.471006090377, 85.9145766899935, 22.3569488352826, 255.883948275235,
                    69.2397851335435, 111.804640222101, 86.0370792856076, 22.7606441852095,
                    255.223877131975, 69.5384977901607, 112.074677633721, 86.6651106975233,
                    22.3679548747046, 255.719703068553, 69.2842142875463, 111.595795228737,
                    86.7905545433714, 21.905022335464, 255.741214636392, 69.5025874560386,
                    111.917474398123, 85.9979382299189, 21.9234835117948, 255.86512298664,
                    69.79726309228, 112.247088838204, 86.3965059185127, 22.1807891184713,
                    255.344198657079, 69.3295244027147, 112.470698874407, 86.0242923109723,
                    22.690297232237, 255.655667055589, 69.5236790738213, 112.664480705874,
                    85.7216716254899, 21.8921723987322, 254.663783503588, 69.0945191453518,
                    111.433466847823, 85.899499712391, 22.6339152421072, 253.684553106031,
                    69.2472158565146, 111.473112601981, 85.9507591575723, 22.7304107310608,
                    253.697957455335, 68.3129962573257, 112.668399425505, 86.4133664954366,
                    21.6962055524492, 254.088661119662, 68.7376081338805, 112.859728161774,
                    86.335639898298, 22.6417688999322, 254.720678582386, 68.3515330004621,
                    113.023367736043, 86.3467822642107, 22.4472737864617, 254.357882467917,
                    68.5719326827016, 113.27881636698, 86.2501522302173, 22.7703231730613,
                    255.121760180477, 68.6507931930883, 112.79387968162, 86.7964248453814,
                    23.1248135184492, 255.43303363856, 69.540394150888, 112.847586942598,
                    86.5989117339251, 23.7648413715674, 255.087620521771, 69.0626854941113,
                    112.38642148834, 86.301468149336, 23.0560384968422, 255.047527016016,
                    69.1134211240392, 112.143531031398, 86.4528388783756, 22.3468688978627,
                    256.369951670048, 69.0571570794805, 112.960240486643, 86.035348611525,
                    22.8098149488125, 256.289641234734, 69.4073793078101, 112.191066009713,
                    85.5228621502836, 23.6334923360311, 255.928150703919, 68.9141165183808,
                    112.228812735661, 86.442505184006, 22.8065640480406, 255.862109150169,
                    69.5836556066031, 112.306512742107, 86.2923831203769, 23.0019321522742,
                    255.69229408956, 69.4104857298323, 112.357692290743, 86.1036227810096,
                    22.7985654258459, 256.055670616728, 69.5033326001109, 112.327310051848,
                    86.254869247355, 21.9425674691835, 256.161293384279, 69.5783857291147,
                    112.48053022865, 85.9885978251832, 21.8124266488952, 255.327255670863,
                    69.387613763529, 111.961498904457, 85.7064923011851, 22.460387922981,
                    254.801780332979, 68.9625824891754, 112.112579584262, 85.891450028086,
                    20.9925655891523, 255.176931348665, 68.9652577490366, 111.436538924416,
                    85.2877985072582, 22.4370716669423, 255.542153317171, 69.046347402274,
                    112.9215628756, 85.8376834043555, 22.3850780681337, 255.409655287364,
                    69.6245484136731, 112.81202045462, 85.6113648473253, 23.2086469650433,
                    255.042516418903, 69.3020423470433, 112.884857113631, 86.0107021446536,
                    23.5033837307559, 254.894113234028, 69.9934969993222, 113.217233243464,
                    86.0583021591087, 23.0084315567863, 254.894749039818, 68.509613677552,
                    112.723377722318, 86.3652309776267, 23.362925906818, 254.842894491922,
                    69.4606130156307, 113.0507558998, 85.9214504661909, 22.9131955148748,
                    255.34239810336, 69.5682647565678, 112.811351968105, 86.0466290091432,
                    22.8599060788085, 255.704705350249, 68.288307903295, 112.437253935827,
                    85.9798434898511, 23.3728920753185, 254.720908646267, 68.9177756986564,
                    112.630482372806, 86.6115599912503, 23.2184808359496, 255.967528769902,
                    69.6017258038984, 113.011208295849, 86.1262928146752, 21.6566240822656,
                    256.118313454194, 68.9757505572262, 112.953604611319, 86.5643959260688,
                    22.2772748579256, 255.709544806802, 69.6180021336087, 113.083303383056,
                    86.4040435482882, 22.9109434188316, 255.571256356872, 69.4796556397258,
                    112.544086983095, 86.4425757292046, 22.1193999598504, 255.898355515275,
                    70.2264215425161, 112.382816098633, 86.4198168052465, 21.5204327386043,
                    255.575216566251, 69.6718323579454, 112.768860043003, 86.1869620116717,
                    21.5173143978799, 255.233179882631, 69.3115916025643, 112.23784507651,
                    86.1191850260099, 21.5696313883164, 255.558218099834, 68.9061865564449,
                    111.921931960545, 85.5357994131726, 21.7324816948096, 255.035257876625,
                    69.1467416259289, 112.047276200921, 86.0835295575026, 21.8834253508726,
                    254.074276551679, 69.3436919506019, 113.416784293529, 84.9988462594549,
                    22.7715866506946, 255.624596622476, 69.4398767734668, 113.806561337658,
                    86.1758345137295, 23.3993348494889, 255.563899286543, 69.6561121899796,
                    112.351626348494, 86.0207859271023, 23.4814511949458, 254.21512401673,
                    69.2456175284265, 112.770578600525, 85.703020024417, 23.5556238360699,
                    254.635801540296, 69.702510559057, 112.904979981059, 85.8146303419234,
                    23.0092013818925, 254.591549837911, 68.9368852523133, 112.018083720777,
                    85.9325211073732, 22.8094973969082, 255.281458693113, 68.5175727866672,
                    112.64498766487, 85.9494555805134, 23.208347323894, 256.049165700038,
                    68.4736205420158, 112.451227739736, 85.8293561882538, 22.7804155604392,
                    255.79577501685, 68.7623435582394, 113.075027839624, 86.3529701591058,
                    22.9682773690439, 255.879888760733, 69.5568663252979, 113.079582162929,
                    86.4665701486575, 23.4839403606116, 255.836730841841, 68.884944436965,
                    112.335212274944, 86.0867267528107, 22.2830838588065, 255.146750134144,
                    69.0247458457024, 113.386600682665, 85.6670005228798, 22.3146594584466,
                    256.692456827479, 69.528669301116, 112.76851546821, 86.0614603064544,
                    22.4351486847221, 255.817292198225, 69.6863073794195, 112.381092195588,
                    86.4165045729492, 21.0206026648678, 256.189531569263, 68.9905717856538,
                    111.636141435783, 86.5122812583529, 21.1653123915199, 255.775396721227,
                    68.1510488358307, 111.758580226925, 85.7326967093329, 22.4893634612405,
                    255.259377912606, 68.7213411451226, 112.048354447367, 85.6872620937278,
                    21.7507532124972, 254.61229039897, 68.7089908470018, 111.77670805893,
                    85.969000448766, 22.1734132466426, 254.588123977231, 68.503518353963,
                    112.140113194093, 85.2287571672756, 22.9393951658097, 254.927210278014,
                    69.2644323182359, 112.433007780684, 86.0188572042321, 22.6774311159575,
                    254.496719960522, 70.1526519550381, 112.811381724125, 86.27684960378,
                    23.7581612459714, 253.959694181693, 69.0304214325432, 112.635582621747,
                    86.0263419232357, 23.2482643098552, 254.703262346168, 68.5484601831405,
                    112.455727180365, 85.308776458829, 23.5331767367187, 254.414625060992,
                    68.6789186892832, 111.926627531302, 85.6117779736152, 23.1432433965431,
                    254.489306358797, 68.2495553540691, 111.794211608261, 85.6547994900858,
                    22.6635734129712, 256.046662532293, 68.9536614190896, 112.398359872466,
                    86.0830472957659, 21.509056076586, 256.127692082003, 69.6101788782365,
                    112.569327584375, 85.7532407965104, 22.9404123874454, 255.316186649399,
                    69.0131542334947, 112.658615680327, 86.040091596776, 23.4896502314863,
                    256.577340323698, 69.0937046130925, 112.168138164758, 86.5816692504674,
                    23.0490569242928, 255.961770156751, 69.5629035606421, 112.610393705591,
                    85.4826467709995, 22.5601597331626, 255.829293092093, 69.3093357692697,
                    112.439621333141, 86.0250278229107, 22.6053438972026, 255.810267958332,
                    69.4237064279244, 111.44440600077, 85.9063221452904, 22.5206408008442,
                    256.174786396461, 68.60995488337, 111.619661561518, 86.3380104453577,
                    21.7689408027459, 255.715071587849, 67.5875160628947, 111.574828955126,
                    85.9703858642646, 22.1271004282272, 255.224475491354, 68.3842218490166,
                    112.098713340184, 85.9495315948265, 21.7464721515857, 254.753827320058,
                    68.8038960488085, 111.893290741976, 85.6487594643083, 22.7140968284483,
                    255.897542104913, 68.9996384079331, 112.363995602232, 85.668369078303,
                    22.399592763304, 255.078294406082, 69.0140308176888, 112.359539321354,
                    86.4971182042167, 22.8964708935004, 254.577039217744, 68.9500584802133,
                    112.583597966305, 86.3588350348091, 23.2793700618159, 254.437123534264,
                    68.6606257681442, 112.672071078306, 86.6962103990174, 23.0528256242295,
                    254.390834852749, 68.3223153368224, 112.637293696562, 85.3923618099002,
                    22.8698075126039, 254.64666463507, 68.5428338485933, 112.312833936717,
                    86.0512434128181, 22.2778623646704, 254.403521635467, 68.6857614249551,
                    111.429108736937, 85.7547548358359, 22.3668830071096, 256.055591173332,
                    69.6355979365832, 112.428195870614, 85.9821639400317, 21.5202851873674,
                    256.160880519997, 69.4734731192785, 112.431012643306, 85.108162774702,
                    22.8377792975588, 256.26013308978, 69.5761288839622, 113.027009542539,
                    86.133160759112, 23.3116616243672, 256.19022473716, 68.7484972739203,
                    112.323496411347, 85.9373841051459, 23.2839158824381, 255.460166526621,
                    69.2430749306283, 112.049192306687, 85.6048712724244, 22.7296223067825,
                    255.512338452788, 69.2592635051474, 112.811894382511, 85.9297637745315,
                    22.6021455646415, 255.044191447695, 69.1518937899455, 112.584757086884,
                    86.0641584633875, 21.868972448771, 255.856134184098, 69.3819473067615,
                    112.833844091327, 86.1152964959307, 22.0632660098261, 255.800729085454,
                    68.1953589901362, 111.818968403466, 85.9853459488249, 22.0146668438557,
                    255.802184118022, 68.7000446328651, 112.194136084216, 85.0830412838098,
                    22.0839974777955, 255.607899607649, 69.2188847015242, 112.266387668212,
                    86.1324030256822, 23.0829556511873, 256.230724818832, 68.0076526056336,
                    112.114385648672, 86.3592238535067, 23.2533964534454, 255.630573252289,
                    69.271841132465, 112.146731314587, 85.8575672029426, 23.1929854751723,
                    254.603429829262, 68.864897747767, 112.248860768069, 86.176383059118,
                    22.1354791577891, 254.700198986887, 68.9374393943973, 112.374010935999,
                    86.6147806769729, 22.9497181863654, 254.703782786876, 68.7656552100103,
                    111.985828730658, 85.6516415584757, 22.1508362739756, 254.909132081157,
                    68.8580195085147, 111.605853531196, 85.9237537359385, 21.8752972665894,
                    254.062181406841, 69.0615328898936, 111.94702424009, 86.1681697580142,
                    21.825362428786, 254.888592353037, 69.1238839388365, 112.822671571365,
                    86.2073492954062, 22.6340941076135, 254.548972019385, 69.5188295187454,
                    112.822034409185, 85.7956517228593, 22.7831357879237, 255.869256076189,
                    68.912926633539, 113.215514763401, 86.0370281678266, 23.1867788675285,
                    256.682020442339, 68.6793812083464, 112.511531010593, 86.3886758613407,
                    22.8244891828488, 256.258967582514, 68.9389054071822, 112.435170685881,
                    85.9526436009773, 23.2316512562974, 256.014785993528, 68.2623162925038,
                    112.25300765473, 86.368411874931, 22.6995573102051, 255.320366842486,
                    68.8133473728476, 112.483090255605, 86.1857419352841, 22.6316157171689,
                    255.76985373988, 69.8677625489797, 111.80820731926, 85.5466050358447,
                    22.4538474091364, 255.547489915784, 69.8681140633606, 111.907364220436,
                    86.0163258617646, 22.6454495482973, 255.385769753722, 69.2114764689363,
                    112.459254450969, 85.8384398332596, 23.2593915858919, 256.852972821615,
                    69.7740545469435, 111.520999373332, 86.1471014027612, 22.9596346219383,
                    256.228823734372, 68.4189647760392, 112.552708399283, 85.785520490855,
                    23.2740590079665, 254.491651544995, 69.5166248655645, 112.525425929926,
                    85.921628090959, 22.9092155134179, 254.961769056387, 68.6159827204641,
                    112.204273394593, 86.5560255467116, 22.1516080689664, 254.574123584135,
                    68.6158353280309, 111.826795203834, 86.5508734129887, 21.9251089489078,
                    253.959980989414, 68.7161498656797, 111.767290385973, 85.6853087457008,
                    21.9088810573356, 254.686364073104, 69.0336659010086, 112.223883475812,
                    86.166058146642, 22.3376360211993, 255.074852694987, 68.9352362851429,
                    112.293774430853, 86.2094096435269, 22.1891732568399, 254.61635174626,
                    69.9574308318977, 113.070155024852, 85.9905975585748, 22.3249523494252,
                    254.727598086342, 69.2245305559686, 113.042026620357, 85.6333733644101,
                    23.1247910579871, 255.515117438067, 69.4701668898927, 112.677096605677,
                    86.2051957410164, 23.5472407088568, 255.405684749901, 68.9031616996068,
                    112.82122745895, 86.5116156081121, 22.6849192652921, 255.172081659088,
                    69.3352759531882, 112.517167136578, 86.2554168412775, 22.0684510964382,
                    256.222017428818, 69.4755536553628, 112.203426646094, 86.7980039335422,
                    22.230053786989, 256.152776787903, 68.9263721543416, 112.261037288966,
                    86.1567351620592, 22.7713006275546, 256.423530384056, 69.3887786614085,
                    111.604529099644, 86.4989695871186, 23.6749972294141, 255.9510981101,
                    69.2361547019296, 111.982570605334, 86.9719535770407, 23.786178802349,
                    255.966056578704, 69.6899479130497, 112.025314436198, 86.3651655331973,
                    23.4354974707676, 256.005432071184, 69.402554802601, 111.650503330426,
                    86.0404092601297, 23.6098294226628, 255.18154897801, 69.4802821957693,
                    111.899307854336, 86.1438975635464, 21.6003188874058, 254.474657804872,
                    69.4045060902536, 112.147938669013, 85.0533877010683, 21.82650189066,
                    254.509897430777, 69.1350522693383, 111.917533187982, 85.4810095901558,
                    21.2671630420204, 254.379066591573, 68.7729631458017, 111.638154651249,
                    86.2279224641636, 22.3110707313432, 253.920688671393, 69.1667039355754,
                    111.947635579794, 86.0955762412547, 22.7959276568157, 254.629720162124,
                    68.5791806782263, 112.414524248262, 85.7924582344399, 22.1375957551481,
                    255.18073127113, 68.4351591076318, 112.248813482996, 85.5204856191904,
                    22.1600705492287, 254.697860802964, 69.5979454179554, 111.489426949676,
                    86.1591369101107, 21.9861786547594, 255.187198704417, 70.3561366430911,
                    112.231678320942, 86.219680590295, 22.9491873529751, 255.898106164205,
                    69.2876567596279, 111.939176240102, 86.2268500434944, 23.0218487415962,
                    255.9928656528, 69.1893938897145, 111.522046930253, 86.7911278510848,
                    22.7891659283258, 255.065246380615, 69.6501298642732, 112.149470983562,
                    86.0062283794377, 22.6506096700881, 256.56219323969, 69.1500217972711,
                    111.847196199412, 86.4034701646163, 22.3817662108205, 255.740874133218,
                    69.1612600567963, 112.422434936757, 86.6197988954224, 23.3641235808054,
                    255.709701454914, 69.7122387056248, 112.582587412143, 87.2981367820714,
                    23.6264342601829, 256.578258492688, 69.2345037137483, 111.657344817524,
                    87.0736202366634, 23.779492880046, 255.571059195479, 69.0963548384299,
                    111.725627268724, 86.6273146701206, 22.9660278509947, 256.437079554301,
                    68.4628769268051, 112.512149504278, 86.1084717978213, 23.6431033508564,
                    255.952541012137, 69.8725271670123, 111.311602002957, 85.9181293953202,
                    21.9691128972138, 253.718726427282, 69.7369285878823, 111.122132577911,
                    85.8482554498, 22.4320357345948, 254.080621929879, 69.4584248818234,
                    111.582673264222, 85.269293767113, 21.9171874252397, 253.650717114466,
                    69.2361124704404, 111.080189129978, 85.8170395007446, 22.1900950592271,
                    253.286321782, 69.0126669457685, 111.443108217318, 85.3140757762565,
                    21.9748454972217, 254.432001845705, 69.0100802242735, 112.359504909664,
                    85.6155965522727, 21.8291742115101, 254.409202330288, 69.5318550594242,
                    112.38590959125, 85.6972852560671, 22.8828814993839, 254.310423978667,
                    69.3847089466893, 112.313389713873, 86.2091150648594, 21.8114931019795,
                    254.828269595248, 70.2584835294883, 112.354706708959, 86.1909275072646,
                    22.7044747324433, 254.980221032243, 70.0709516197363, 112.070729513581,
                    86.3103709241125, 22.8593665790836, 254.618155103714, 69.6596043520292,
                    112.387552221656, 86.0385122881725, 22.5380931747164, 255.760118373098,
                    69.8999967621992, 113.001785019059, 85.836794136659, 22.7045882726235,
                    256.738060138002, 69.7901540192721, 111.779452390926, 86.7296856799922,
                    23.3572529015767, 255.741655465749, 69.2032580724752, 112.525785740699,
                    86.0292311903605, 22.7263851121385, 257.446820501116, 69.6697908686485,
                    113.474799688129, 86.6260783525167, 21.6645947374701, 256.160680811467,
                    69.0774125705067, 111.499805964845, 86.5800214863683, 23.6044238878003,
                    255.800798840218, 69.1012337749259, 111.381206188257, 85.9009649733081,
                    23.5909069423667, 256.259072393753, 68.3666684073255, 111.821722599163,
                    85.2113055009459, 24.0718322902707, 256.043414760465, 69.3159746645286,
                    111.207008239025, 85.897082580506, 22.804577751264, 253.880290812104,
                    70.0299400427514, 111.508544871734, 86.0604295113071, 22.6351748524692,
                    254.259915735213, 69.0588178869372, 112.171386314909, 85.8115045157824,
                    21.894432135997, 253.486128657999, 69.3850843774723, 112.001283386592,
                    85.705204072151, 21.2566853204157, 253.187288281103, 69.2453736849122,
                    112.288924058285, 85.7845933961743, 22.0947735803465, 254.100727072435,
                    69.2514606893192, 112.551455199357, 85.7866379911488, 21.9225588701838,
                    255.148434683939, 69.3218570656357, 112.456359110495, 86.1422385891333,
                    21.9294589360553, 254.202238998525, 69.7282816503952, 112.206781414225,
                    86.6927294545882, 22.2933038545969, 255.057918669916, 69.7147177864229,
                    113.051731064063, 86.549505566054, 22.035822867777, 254.851960114693,
                    69.7511998741909, 112.022709927017, 87.2412607866066, 22.7285352456525,
                    254.972260813613, 69.3897330549255, 112.661199527769, 86.2058386405926,
                    23.0368492672996, 255.51860547242, 69.8636584794898, 112.278890648487,
                    86.4512811906018, 22.4699703093489, 255.121146179971, 69.3322536215855,
                    111.877703965385, 86.6445644609418, 23.3925345431662, 256.006914962764,
                    69.9238355030256, 111.672001236713, 86.4499502632777, 23.1824143274105,
                    256.874162804362, 69.4633554029261, 112.357685558014, 86.9025342403765,
                    22.7239420174088, 256.738302065105, 68.6282937319809, 112.083464559939,
                    86.2484615075327, 23.4211830729909, 256.364583988412, 68.0157103600209,
                    111.855987925418, 86.1563180375251, 23.985836942165, 256.736324487453,
                    68.9678742695656, 112.262019385731, 85.7935891303893, 24.2506268689351,
                    257.238494015061, 69.2506463794678, 111.420916576801, 85.4928588033814,
                    21.9686349699773, 254.425693980397, 69.5963747187638, 110.989293925885,
                    85.9166872940112, 22.6159378434212, 254.766627959268, 69.2305117405026,
                    112.233491284936, 85.4413826758384, 22.0386110221454, 254.812956596339,
                    69.1288493921258, 112.433606090407, 86.0111985723453, 21.3756396874233,
                    254.12768989311, 69.3104440245934, 112.999981703478, 85.816161972383,
                    23.140866630054, 253.283861901693, 69.8509553417326, 112.936200660335,
                    86.4268875646408, 22.7080136411116, 253.70354543491, 69.3411569658801,
                    113.307456286767, 86.2144500479505, 22.7175526821507, 254.251031977644,
                    69.3488119288936, 113.058463598587, 86.863864064445, 22.8439235769351,
                    255.082571927371, 69.2501336934764, 112.816726926459, 86.6805831134926,
                    22.934007949715, 255.57726592615, 69.0719405821241, 112.840437573647,
                    86.4218083822771, 23.4293704499118, 256.070286137022, 69.9380380817073,
                    112.098635987538, 85.9978986819706, 24.0718391103225, 256.382510016436,
                    69.7729719167585, 112.434201962202, 86.4981371146527, 24.1431982661055,
                    255.639226639859, 69.8281978668615, 112.266721241828, 86.5562693965774,
                    24.3369995139037, 256.197335943542, 69.3238308840605, 112.574159878663,
                    86.2011349696137, 24.4424809952485, 257.182584813202, 68.8518059934481,
                    113.618955363404, 86.9261229093828, 24.7512099308755, 257.554797020047,
                    68.5765398257987, 112.596250922988, 86.5837938233244, 25.1705271921982,
                    258.113001927396, 67.9887857549017, 111.838025848488, 86.615090128769,
                    24.3657304498025, 257.431622137361, 68.7268849984558, 112.537101221595,
                    86.6813306642779, 24.1321648989061, 258.102967545067, 68.8390273890043,
                    110.549501469084, 84.4687882561983, 22.8425861214141, 254.508128818085,
                    69.0316163837444, 111.491384717355, 85.2023740711363, 22.4956568128519,
                    255.750049466748, 68.5841010463138, 112.468013004427, 84.7046656481844,
                    22.1011824802222, 254.177178884669, 68.6231334362957, 112.44639659131,
                    84.8249747255269, 22.1124102011529, 255.309131121222, 68.6708685468095,
                    112.699332843949, 86.2221219047427, 22.8187468223587, 254.514498214761,
                    69.6651073465831, 113.468308687664, 86.6336193662203, 22.4447149885268,
                    254.001139999206, 69.2939748697433, 113.779468229476, 86.139195100773,
                    22.4080451721837, 254.54758863051, 69.5882520013619, 112.98444731973,
                    86.9206005858477, 23.1650200855593, 253.491809905631, 69.194695198047,
                    112.862948520523, 87.4423487171709, 23.4437514712923, 255.496222016644,
                    69.6187699773486, 113.375580124678, 85.9504189536386, 24.851749351824,
                    255.90963553672, 70.1726724148418, 113.790648670002, 86.8280311559439,
                    24.1273559972389, 256.844272393982, 70.3254413261604, 113.662348336229,
                    86.2136665894022, 24.4248461228033, 256.613100070943, 69.9360695641835,
                    113.611569777036, 86.5259462304858, 23.9300088063979, 256.612853780014,
                    69.6526918556188, 112.894678672333, 86.5986311214627, 24.1082985403636,
                    256.533904630626, 69.1369974488895, 112.155625718158, 86.8945250471484,
                    24.724512944867, 256.075385530915, 69.4571113367866, 111.934294952288,
                    86.7200737331805, 25.5724626178311, 257.489757552112, 68.6346665822081,
                    112.483140314474, 86.3095322796336, 26.0953002438016, 256.581352890986,
                    68.4611811260839, 112.422874043969, 87.3050412145553, 23.9844964559488,
                    256.826272178091), .Dim = c(5L, 520L)))

    expect_equal(
        recover_alpha(out, n_message = 1, constraint = "resolution"),
        structure(c(-0.0614319107745018, 0.956210025182244, 0.202364767222321,
                    -0.0575119519976397, -0.0834761834007622, 0.869432113981247,
                    0.629523292690088, -0.899800838297153, -0.90326416680287, -0.120747341207533,
                    0.644605955162547, -0.0174618044175077, 0.130816579484645, -1.03667920236796,
                    -0.414987914462699, 0.513392547020629, 0.39519102538641, 0.0267497526988905,
                    -0.749129126881797, -1.06350350183627, 0.22260579705835, 0.55636084761359,
                    -0.00388584278618964, -1.06743423005675, -0.925404475678278,
                    -0.230209496203571, 1.02804664815753, 0.447297123630676, -1.54003078232213,
                    -0.691218161913582, -0.082437064696677, -0.0422712312287459,
                    -0.476705937784857, -1.17046851952018, -0.612067314025239, -0.728625798046721,
                    0.386354199060278, -0.466392458033397, -0.814945561752083, -0.467737129956902,
                    0.223826685152858, 0.272034168236061, 0.168238016609166, -1.10555195453874,
                    0.101861100337004, -0.748470114285936, 0.139753084565974, 0.585583538301364,
                    -1.56773150212447, 1.71268327031214, 0.21236743315302, -0.433240831228296,
                    0.747040275650875, -1.30671163965893, 2.37384813162774, 1.34067058340366,
                    -0.543091203075264, 0.354943296476222, -0.874790505578531, 2.5932021070864,
                    0.802187889628186, -0.302530444678958, 0.715073471871904, -0.313911933007134,
                    2.20548472812848, 1.00534058955773, -0.4344481212662, 0.40820567009024,
                    0.469558455920321, 2.07445283761029, -0.0207923196719833, 0.968483139129745,
                    0.263251054956584, -0.466191130238926, -0.212513857189236, 0.696473336087053,
                    1.07558322998058, -0.416795649348558, -0.732276899022125, 0.00623861775582668,
                    0.443369297493433, 0.16852575704398, 0.282346507856971, -0.99506052632978,
                    -0.345977178010457, 0.137507651013465, 0.168177295907874, -0.278988336088048,
                    -0.68914093146752, -0.911136054510763, -0.179752677704101, 0.492887887394076,
                    -0.73458213698575, -0.784874406392703, -0.788834574951181, 0.0778044129520481,
                    0.841388070093153, 0.35314863346548, -1.18303916523526, -0.461069825678635,
                    0.376716437394009, 0.363385230627742, -0.0349381378642875, -1.06181276153242,
                    -0.138308035278783, -0.294467732415455, 0.47114783783546, 0.0976001548297063,
                    -1.01925342717885, -0.731639778724968, -0.194745180289772, 0.97241209505647,
                    -0.343288505783107, -0.793885224980755, 0.010883850469952, -1.0108718147419,
                    -0.0721542864509033, 0.352887841598374, -1.16816481637019, 1.46565536999285,
                    0.555800851094006, -0.392363256454978, 0.593903126524367, -1.38312295698204,
                    1.91709230503889, 1.29585703256737, -0.551361488851278, 0.0492275854081186,
                    -0.820021566904302, 2.11634454527285, 0.842537768188095, -0.433760002759414,
                    0.740862462623397, -0.1633681062915, 2.15625898018615, 1.13358305313201,
                    -0.990323181117681, 0.589412344636742, -0.025841709836854, 2.28489534483381,
                    0.837703844403762, 0.138664359679908, -0.143947373933656, -0.948050518945777,
                    0.919005765976124, 0.4379242923915, 0.410581657937605, -1.13194538522708,
                    -0.81749040136701, 0.628824621074614, -0.342360864173003, -0.701733029288114,
                    -1.08700861717784, -0.785484348285493, -0.604557424389924, 0.0235350370715537,
                    0.0958072977507314, -0.136679216431503, -1.22349063273015, -0.936477077682156,
                    0.120314712826968, 0.181818616845291, 0.103060343544655, -0.111404282650454,
                    -0.676274033374852, 0.60820632786951, -0.331870356109974, -0.0620787540937613,
                    -0.94439724586428, 0.112621371315555, 1.5320550464055, -0.298953903956118,
                    -0.14819827731975, -0.512112197498439, 0.464399495671103, 1.28891581346281,
                    0.254766531065115, -0.74245080559615, -1.18887770285856, 0.360704259160627,
                    0.760731586832208, 1.08282853855142, -0.204559576553891, -1.20020318055528,
                    0.985858804342058, -0.499501548929885, 1.18391649552124, 0.0990640452109233,
                    -1.90194184324449, 1.06976403658931, 0.354031948402536, -0.101956093068139,
                    0.401956447613827, -1.82417053923119, 1.52189502959442, 0.122908926736358,
                    -0.571917960974829, 0.128615176824781, -1.23516624931864, 1.88786882220148,
                    -0.0546577328475735, -0.948074096124657, 0.997687297470009, -0.501361822609965,
                    2.17322032660971, 1.20964204292719, -1.05353172056927, 1.03036000648655,
                    -0.924900508558821, 2.20566139938936, 0.189871784456045, -0.0578435709665541,
                    0.169842907644238, -0.110037783151469, -0.0790531809732897, 0.279341790779512,
                    -0.113773764777306, -0.605767602699302, -0.577165151645417, 0.560425049172323,
                    0.740605941611193, -0.792647156997745, -0.222584739524507, -0.309453140803001,
                    -0.185880846724302, 0.124974513008993, -0.701321266742504, 0.0559019159304626,
                    0.0940185923841881, 0.108069412102566, -0.0301951173535429, 0.122543451680791,
                    -0.536722055719167, 0.587358497388067, -0.260848259559793, 0.194315527891405,
                    -0.241064413668084, -0.0759356123805333, 0.137907312255336, -0.112701613465447,
                    1.32276797541078, 0.114105633886602, 0.023698440390433, -0.190933395063993,
                    -0.209304606079627, 0.899829456324341, 0.265562301224776, -0.211881860196627,
                    -0.0975321124577846, 0.177689006932667, 1.08915370211341, 0.558965611431056,
                    0.0551403980152827, 0.174989735141253, 0.332337133176509, 0.710948182315946,
                    0.0610667626531551, 0.340399930432682, -0.750557106863027, 0.827711179617268,
                    0.394376877181202, 0.359724337243833, 0.664459751043069, -1.04212495819012,
                    0.41125784477498, 0.307385567153659, -0.230344737679609, 0.342611537774076,
                    -0.352524005943337, 1.40077340548589, -0.114641422152204, -0.41753415308316,
                    0.417295230650154, -0.786043445702418, 1.61191991497358, 0.197815882937164,
                    0.158909836075949, 0.312307166351218, -0.41813901526659, 1.82373790062192,
                    0.161012938104761, 0.647389302340486, -0.261375059078773, -0.0597874473125586,
                    -0.347829718836124, -0.0777511508695738, -0.775515791481105,
                    -0.0385737593660451, 0.218425815474419, 0.881896443874183, 0.810078830328536,
                    -0.474017148855374, -0.416127849130248, 0.0141126448941549, -0.433289005743433,
                    0.859389100628135, 0.171074504100808, 0.510712133485299, 0.898991361165031,
                    -0.0944979682019493, 0.292784404369542, -0.0114452590939962,
                    0.124640451398022, 0.845612204774966, -0.172681187166035, -0.843153730695633,
                    0.615792547330486, 0.132566679073136, -0.195148618121053, 0.165848356213388,
                    -0.392934204669871, 0.04807351190275, -0.127446911122263, 0.185936625814577,
                    -0.230442314552334, -0.0432310689928386, 0.490299792400087, -0.186149978200547,
                    -0.347510693239826, 0.0899900691827895, 0.345298370604041, 0.303260391005409,
                    0.248581074918491, -0.141613822227654, -1.06452498108473, 0.405080454571021,
                    -0.173358421164622, -0.245370860010837, -0.725293131642331, 0.201943711683839,
                    0.477908757123998, -0.201314455253133, 0.358530537342517, 0.0901481085295899,
                    0.164959761299031, 0.372584000837946, 0.492858399085399, 0.553134961209025,
                    -0.116804527981628, 1.86914589529724, 0.512549878349489, 0.236437454566044,
                    0.00201463278847314, 0.167536361569091, 1.2489152604665, -0.00751352264859406,
                    0.804492828049526, 0.251044564387712, -0.287955290129844, 1.41619018266763,
                    0.553180723916938, 0.322602252403442, 0.268206129358703, -0.563522280700525,
                    0.531722697051805, 0.651556323744842, -0.183861847000927, 0.302117247490401,
                    0.198334707223054, 0.769678788227566, 1.01400140934848, 0.423790713457691,
                    -0.323090037416932, 0.210665895935108, 0.733743749630577, 1.00662090308242,
                    -0.0965023970899281, 0.316800278929747, 0.189481185574991, 0.307811863277351,
                    0.5596201661686, -0.825770509281853, 0.272177621056983, -0.379210920991881,
                    -0.46000284904315, -0.630487593433401, -0.486625303474369, -0.000448543212542063,
                    0.748421423644658, -0.949017504297473, -0.114299527720647, 0.335564446942051,
                    0.0734217682673943, 0.204148557715101, -0.328273826285141, 0.414343386884127,
                    0.682449924561453, 0.672634439480476, -0.165738565990608, -0.343127228443109,
                    -0.122102129349116, 0.241499880421472, -0.0285631668460269, 0.175620344668516,
                    -0.0194982515168078, 0.275733632356008, 0.457203123698235, 0.765345534273877,
                    0.463178408865559, 0.07070164464713, 0.416439624965193, 0.473247301938699,
                    0.164063862541525, 0.451724257012998, 0.980774035096715, 0.695333997648476,
                    0.0164895389733601, 1.04007979638669, 0.455423403621097, 1.32540175137865,
                    -0.172034585383244, 0.675008707642121, 0.771101923145181, -0.26082678929788,
                    1.33889987573448, 0.269830073545478, 0.605220959720022, 0.146202017519045,
                    0.371479556163038, 1.38014391881757, 0.712105821098902, 0.176258341822219,
                    -0.51540700381662, 0.489323589412706, -0.265639770074358, 0.763404779019936,
                    -0.891641465557854, -0.279703161659029, 0.473975296593906, -0.102856709098319,
                    0.258725991541922, -0.418934481279244, 0.0344701981480284, 0.0413120359177555,
                    0.744645831204224, -0.500815075699023, 0.550842440235883, 0.241699435026021,
                    0.0528727123889325, 0.726067676634187, -0.701494539000095, -0.4041989705836,
                    0.521687412822871, -0.159157505653113, 0.644196932324633, 0.11247226769845,
                    -0.234349642766993, 0.187036254915768, 0.00594408873470798, -0.248289096219423,
                    0.270580756243618, -0.508336009130204, -1.0419240894775, -0.137664711887851,
                    0.282073147860075, 0.347790743912981, 0.338296895928849, 0.053935696072017,
                    -0.0104362384590502, -0.667971930273723, -0.274534373024522,
                    0.897929033260695, 0.183801799148625, -0.0392505413440034, -0.301284122096774,
                    0.393467222050333, 0.5811617609086, -0.00408200150985749, 0.611824173046671,
                    -0.129005424889044, 0.548724395617569, 0.122367221300514, 0.738218142906476,
                    -0.260571526858143, 0.355896967563552, 1.03900185763527, 0.02687450183479,
                    0.748865895647498, 0.306430596947848, 0.428099634947046, 0.605863354296133,
                    -0.156406815590358, 1.10359583998229, -0.0781786818858023, 0.1662416945166,
                    0.649707818635392, -0.31417642545513, 1.08410701343928, 0.314855426883994,
                    0.500111233335417, 0.123811105205277, 0.405956487698276, -0.561715637723353,
                    0.263475817956191, -0.551474894079092, 0.646266096941304, -0.71901955385389,
                    -0.0797078302730938, 0.0416729377558767, -0.596948095669632,
                    0.301478243285828, -1.07951130835504, -0.224325796564642, 0.130562587685716,
                    0.427611710794935, -0.549626882999235, -0.529958580758787, -0.371483419317713,
                    0.41710896250963, 0.119946596314293, -0.541956085410021, -0.640369687956337,
                    0.740305356807312, -0.0523967272625754, 0.126284695004273, -0.391557096013258,
                    -0.0836504338935171, -0.374215583672935, 0.516566597578972, 0.487378839350072,
                    -0.619000985572541, -0.13405332611547, -0.0431306953933586, -0.272111941057844,
                    0.128375624987228, 0.0529841560549471, 0.198488508769373, 0.441321418933171,
                    -0.37644629109792, 0.280223888353873, -0.351043142022689, 0.349944863760527,
                    -0.0528083313069487, 0.316960153654051, -0.239747094863162, 0.065522623647098,
                    0.881510979589002, 0.174077461638149, 0.0202515250769437, -0.157741574379685,
                    0.102332780529593, 0.225129752096016, 0.263522804476281, -0.538145374631995,
                    0.689855484504619, 1.24012905805409, 0.995832638668055, 0.775705084971662,
                    0.210348492912431, 0.0851538390045334, 0.363767520541643, 0.976427615669806,
                    1.40089452888924, -0.663587076289588, 0.752281611752181, 0.993542679067843,
                    0.506780660933799, 1.52567294145348, -0.184285817154532, 0.817228154880979,
                    -0.492157790484555, -0.48174313013466, 0.074284809745734, 0.191991498623764,
                    -0.547207591878873, 0.799068980902717, -0.964445901648361, -0.0123992523691641,
                    -0.352887383426719, -0.256051539978586, 0.463709178309642, -0.80398208484808,
                    -0.0020655163483525, 0.453220544429541, -0.573102234384521, 0.0564271201702695,
                    -0.121265844720043, 0.165012295686381, 0.399914183658694, -0.180431538790231,
                    0.240551330030456, -0.0734263617226958, 0.922892141784246, 0.95657400519687,
                    0.396268822028674, -0.00931951193786063, 0.341190127898528, 0.179650672851096,
                    -0.028976620969928, 1.26743816350441, 0.201025620151512, -0.150417227772863,
                    -0.291403559809339, -0.109548971638887, 0.158512387129036, 0.573447477358258,
                    -0.388116064726105, 0.0342736854697137, -0.390218423952106, 0.569737571832206,
                    0.0394559598875617, -0.173456283404143, -0.07293728930415, 0.123495273078902,
                    0.498411002783087, 0.440291775777041, -0.154304190149233, 0.556724991237616,
                    -0.128150276665309, 0.547624297850462, 0.79082812511183, 0.802557384227647,
                    0.681455978646625, 0.241786091183727, 0.312775511096078, 0.0787397880173444,
                    1.06653461537774, 0.152943026870219, 0.249926021152191, -0.493336035265202,
                    0.8714588808202, 0.664215820600504, 0.791796465733057, -0.266887368034034,
                    0.404491572308075, 0.201580109892348, 0.591880166556251, 0.504245297389161,
                    -0.262225100735378, 0.430849347442631, -0.0909683718133181, -0.986553929903693,
                    -0.151323486905426, -1.25636645379035, 0.473705512964415, -0.226420102936117,
                    -0.966983304076237, -0.203751291113178, -0.863263028145411, 0.922734179053748,
                    0.247200065920453, -0.668674469786609, -0.21489776098008, -0.186256337680078,
                    -0.0396794561334843, 0.443529071201112, -0.676663840887386, 0.542126606101704,
                    0.508077300498357, -0.372010194310548, -0.0484609419320634, 0.435076548212606,
                    1.13241198369883, 0.0877435587339477, -0.12859253807818, 0.139474109167537,
                    -0.891563654706346, 0.169142647729473, 0.205428009408806, 0.609782281125575,
                    0.224013547510168, 0.226757725017485, 0.214195426817582, 0.100184544507343,
                    0.214038775130888, 0.381646425296225, -0.311910217633582, 0.32176183458165,
                    0.110720863478093, 0.363102707246526, 0.145143354061048, 0.708409357931593,
                    0.449718826854053, 0.390572386218693, 0.37423331823868, 0.990890206737959,
                    0.723470559573364, -0.381461464727295, -0.27231105108892, -0.149806063589807,
                    0.474854162615514, 1.22031451466444, 0.47138366925742, 0.239494860522072,
                    0.13810287413267, 0.504604411598649, 1.28357776774486, -0.0643890802548412,
                    0.123882742936416, 0.511707679826486, 0.979264753258576, 1.16051425114784,
                    0.223892851447033, -0.9043754204328, 0.424925708221394, 0.561861559209433,
                    1.65188653691837, 0.174825737584854, -0.431101706997765, 0.150810764470862,
                    0.918009666204711, -0.787240153926675, -0.459219478963064, -1.47144187155496,
                    0.206039636631004, 0.27658542828587, -0.613598595311863, -0.429243265515993,
                    -0.508217564215215, 0.729358791097837, 0.375463306079126, -0.472145924972047,
                    0.161114185312897, 0.158228327512205, 0.103328538893606, 0.642060955404958,
                    -0.827816713769352, 0.17949492975788, 0.0247726120097838, -0.520282796074014,
                    0.769226047432909, 0.409211819443037, 0.599996587048508, -0.060447110461169,
                    -0.532340019754912, 1.10929514091626, -0.296340634968637, 0.403108290879288,
                    -0.0745249817658191, 1.04823164158836, 0.592126224942916, -0.598486528918109,
                    -0.0723214650689101, -0.596302812369998, -0.288066170294314,
                    0.137494630262267, -0.342757377550335, 0.265692747764263, -0.125618082389678,
                    -0.571595313282614, 0.547999349335925, 1.17228114745498, 0.2067626455447,
                    0.384835385739166, 0.419553601773913, 0.0794851061185113, 0.940459061988051,
                    0.560237127444793, 0.476669608045228, 0.462859284959961, 1.37115453223197,
                    1.39099012706202, 0.189794212928348, 0.68056523269378, -0.220601500465449,
                    0.909865334046572, 1.76283201554963, 0.159645354276222, 0.241110237478622,
                    -0.215496167415324, 0.732797958117061, 0.849266410526617, 1.23200532878559,
                    -0.238685490892564, -0.141872031414437, 0.969952422956652, 1.20845347556931,
                    0.84240193414098, 0.557746171659943, -0.413056529222473, 0.904848004133441,
                    -1.86742102532449, -1.00210342649285, -1.01922971773352, 0.207335942298187,
                    0.616358590670814, -1.36978000627593, -0.282882070805499, -0.689395423723312,
                    0.345097208382185, -0.0202517246617049, -0.499397793617852, 0.298448021093293,
                    -0.178553146143557, -0.015886080952896, 0.139106800229683, 0.126705839147753,
                    0.179412618567994, 0.349302860112914, -0.0411989962613006, -0.0236703840953112,
                    0.437314492089982, -0.127996286012092, 0.032512867024991, -0.457451198700454,
                    0.586942817905083, 0.105069570598147, -0.352789206713908, -0.391713236315951,
                    0.0740436616872842, 0.626308633437532, 0.290121163608823, -0.0250329208843141,
                    -1.2185303882278, -0.476134200238363, -0.736389364370666, 0.0477301468611486,
                    0.491656850821528, 0.0838869092353036, 0.178519225686813, 0.663813458597389,
                    0.848896727649105, 0.374798779344246, -0.725090824100718, 0.36995566109124,
                    0.575251051944917, 0.931138167335071, -0.299259499051288, -0.175953851672617,
                    -0.467733512651193, 0.431531366395774, 1.38261250153818, 0.695784902131209,
                    0.267211782620151, -0.933285295584199, 0.339192157379557, 2.32723531738242,
                    0.734485782796639, -0.118255027617707, -0.577554815784339, 0.889857183817242,
                    1.68478392227588, 1.21178601117083, -0.296208970811861, -0.190357021102528,
                    1.26250142817861, 0.556297293786052, 1.25868164931171, 0.0148974615200927,
                    -0.0842557893514595, 0.479842461991922, -2.06287076174387, -0.643136456175057,
                    -0.706540446653491, -1.39554916286855, -0.141223611342753, -1.02149734557296,
                    -0.146418491745749, -0.483307196608393, 0.00653386969008807,
                    -0.639059883552932, -0.771304329154468, 0.00806235596684246,
                    -0.143886149174975, -0.281593534373741, -0.0273705501752772,
                    -1.16921066098107, 0.0250270838815823, 0.601201981732132, 0.247560171714269,
                    0.142461411094956, -0.128793907337425, 0.181402981291342, 0.198849591294476,
                    0.569557656456226, 0.413895811633211, 0.148780225475981, 0.418055446576091,
                    0.341800625391727, -0.280138487960897, 0.0991754572493733, -0.849830826346356,
                    -0.570868339365177, -0.381316426858859, -0.357340596946017, -0.740386324205332,
                    -0.0142150297735952, -0.111628361534827, -0.0644748090274163,
                    -0.556047013747957, 0.071759192599643, -0.0989417103676047, 0.218067872935706,
                    -0.679555187536394, -0.808733954511183, 0.151444377141786, 0.434507373318517,
                    0.568571388170596, -0.334826253875391, -0.578568564458067, -0.211751539036442,
                    0.463280894351755, 0.895396093488898, -0.349689069806658, -0.642033079453483,
                    0.843265567501646, 1.47131185110159, 0.682210047439014, -0.0912041245341726,
                    -0.559870783219992, 0.400311510420892, 1.15153940664671, 0.0411149695191853,
                    0.115982811193362, -0.0756683365800086, 0.180204894471473, 1.33155709541828,
                    0.031629339689232, -0.195114793727996, -0.207671089020465, 0.131992082703078,
                    -2.69615825664895, -0.467520146569427, -1.34502337277176, -0.769905578765787,
                    1.78413337617542, -0.951511489070633, -0.507291680504352, -0.58476366811685,
                    -0.152185760185063, 0.956908251684865, -1.65719973977711, 0.373051503220779,
                    0.0218679796639663, -0.224769841752817, 1.51680021642423, -1.62804626939175,
                    0.385820974288833, 0.55378996228005, -0.310378825262376, 0.278298030269127,
                    -0.684340062638839, -0.282572100959612, 0.508056530824561, 0.225740196498435,
                    -0.105597474971233, -0.243649826792243, -0.59360980941932, 0.0768307993573103,
                    0.0540401270985171, 0.313969442313038, -0.978127133353695, -0.887102646517235,
                    -0.127799081707877, -1.10047019900449, 0.212714207111986, -0.668128460262864,
                    -0.0582548785795893, -0.099011162013781, -1.08260842835614, -0.178590499597306,
                    -0.529496165437514, 0.270738787891275, -0.694423558471911, -0.438670489423487,
                    0.498539862453725, 0.344032177843019, -0.368616598172821, -1.02024917144479,
                    0.323722952144266, 0.173384331563796, -0.194578048812566, 0.203086790553527,
                    -0.935242648336413, 0.254331876846834, 0.238920032683325, 0.38437660867887,
                    0.318165870682208, -0.145820609509286, -0.840457821587194, 0.32154505191869,
                    1.21863480058381, -0.356460053104627, -0.12094635696134, 0.142854778802302,
                    -0.547032131946253, 1.85774060319153, -0.784406573854938, 0.265711432077865,
                    0.337129837085342, -0.310398359943605, 1.74467176325231, 0.437930515015779,
                    -1.83723254043903, -2.44272158784318, 0.293376621955304, 1.11285749455499,
                    0.36902958796091, -1.95103867828853, -0.630318014059924, 1.52313352787147,
                    0.765584506142638, -0.0465638327719375, -0.98161211297969, 0.208021483941366,
                    1.69545772554409, -0.010718457382076, 0.684137830816667, -0.45971081642832,
                    0.160722512085471, -0.421000003523574, -1.13246980634399, 0.281553037943894,
                    -0.678526592752235, -0.226309572410742, 0.294313578149456, -0.261965049427502,
                    0.594713380279245, -0.541565997960696, 0.326790315260411, 0.364985231404063,
                    -0.554200583919076, 0.144161873386253, -1.14805250043362, -0.540953007455983,
                    1.04935674202301, -0.756827440961139, 0.284787213438619, -0.020257158652548,
                    -0.175547855859179, -0.386619575462873, -0.568828247423262, 0.631030958727195,
                    -0.552538946057815, -1.06140986678449, -0.759458513391365, -0.750780200101239,
                    1.32921367811045, 0.256344643234947, -0.47450866553524, -0.332306899977283,
                    -1.03551733886614, 1.49176751484467, -0.154519902862894, -0.203850092359005,
                    -0.85217451474162, -1.30730247671318, 1.44188885652559, -0.0639450248820665,
                    -0.14815008416457, -0.508863451147448, -1.08264943659253, 1.64752608494007,
                    0.834666796811838, -0.913776530327027, -0.34852810066036, -0.341851143604089,
                    0.223253780719261, 0.847789411993951, -0.173329040555842, -0.475729313466843,
                    -0.791908498628082, -0.213468995208757, 0.619530560047863, 0.667888979532137,
                    -0.440534437687603, -1.62203402167984, -0.0210729247422989, 0.0269282106143649,
                    0.361945732341582, -1.19633985143102, -0.81475306882799, -0.361252632446082,
                    -0.836323099011373, 0.23916336767067, -1.14079449167011, -1.39899609041878,
                    -0.757816676421328, -1.30843446872865, 0.531704199785622, 0.216968677302219,
                    1.60837370384918, 0.0525907871439415, -1.86608420682195, -1.9985680122364,
                    0.746218878987293, 1.70360694155355, 0.0924893667591675, -1.82050945224695,
                    -1.12527834978465, 1.37918113276066, 0.943684891034543, 0.474537063560291,
                    -0.992229099059202, -1.42729343909647, 1.36695872705792, 0.300249517346089,
                    0.134705835069923, -0.659111572253451, -0.678563429679713, 0.434561374429308,
                    -0.236372710786185, 0.662962712030577, -0.649203661023598, -0.762120195670235,
                    0.466197441188967, -0.145859956713551, 0.939873130343386, -1.17725087004023,
                    -0.585343260521796, 0.614166022027661, -0.356712923798952, 0.387723088409103,
                    -1.13681236940492, -1.24829682518998, 0.683981804034062, -0.694200749173774,
                    0.814482925046804, -1.02359681291215, -0.705529038252138, -0.829190216762242,
                    -1.15993276952611, 0.700822216842255, -0.112648464090448, -0.589822937672011,
                    -0.736127918371778, -0.516817891713281, 1.64466456998784, 0.137253634006591,
                    -0.159268627601364, -0.786049764141836, -0.951790940073678, 1.17959219139809,
                    0.194491243141961, 0.278343665810866, -1.20513015647973, -1.7495687260997,
                    0.74520469464062, -0.223100883035045, -0.316969160280934, -0.0651731321730722,
                    -0.889763588727931, 0.361142550040398, 0.158451525954916, -0.675091766314637,
                    -0.0227981989286832, -0.0357535091215837, 0.0383739995822054,
                    0.387622258078075, 0.0885006108367179, -0.508804997210802, -0.546779848221512,
                    -0.027867787489015, 0.37862524054519, 0.961267632745177, -0.997210173246003,
                    -1.58901127014813, -0.196482574330275, 0.0130446986057606, 0.43250156320542,
                    -1.0044001329577, -1.14206576567562, -0.766030239493006, -0.538420145940975,
                    -0.0858422948376472, -0.994372498939754, -1.0119908101566, -0.667279258052986,
                    -1.59322394374153, 0.208116143968425, 0.0507978140306307, 0.916824601536661,
                    0.297295155903555, -1.63358187179242, -0.932978798565699, 0.20990693361432,
                    1.70908522913802, 1.35972158188881, -1.24529830321625, -1.37192653224653,
                    0.531374751458713, 1.28815470447332, 0.49690728991979, -0.599647213629957,
                    -0.979873989305588, 0.633241205436349, 0.148519298072678, 0.5242391601858,
                    -0.456256409101996, -0.356526742467793, 0.130882717307202, 0.417158101170038,
                    0.590365451439965, 0.055073051487831, 0.258312796912747, 0.548908032148603,
                    0.385851380795579, 0.80380886364577, -0.68336361821919, -0.7119224363272,
                    -0.0380238970553819, -0.953888654946589, 0.318911172676877, 0.082118559123515,
                    0.137493808094405, 0.395672988489736, -0.661258201970419, 0.674270406455577,
                    0.840048226411916, -0.936468938429925, -0.661835019939076, -0.801462947537154,
                    0.185684539330197, -0.896967852594443, -0.224506546073854, -0.733498170383768,
                    -0.992603378589443, 0.3917222441045, 0.699652168401727, -0.842954177448803,
                    -0.517381398323721, -0.84068367837984, 0.786602630384351, -1.41351913003593,
                    -0.618354692943342, -0.333990174394302, -1.20568284898509, 1.11561395777727,
                    0.0658844724112058, -0.734927036871596, 0.116744292822489, -1.37098809502202,
                    0.430939725066111, -0.0437584763550731, 0.382392467687879, 0.656935999516975,
                    -0.549638478775115, 0.952881689187876, -1.09130290527884, 0.909034425853037,
                    -0.596010699166499, -0.287920769528597, -0.0768210005301455,
                    -0.295981461873367, 0.835246194786748, -0.312239925634117, -0.330573978676334,
                    1.23203230088649, -1.22523606752182, 1.33408333529815, -0.86299118705594,
                    -0.259260041984341, -0.72319790685367, -0.772413196418263, 1.05346557826536,
                    -0.297072501890739, -0.132829954347493, -0.71264539862446, -0.678328611759056,
                    0.853322919932879, -0.456249949456151, 1.30884630365139, 0.713142553088375,
                    -0.667722370301043, -1.69776325357995, -0.60930613847222, 1.27713763385071,
                    0.422747143705578, -1.32993729228311, -1.5269304113358, -0.319256335441736,
                    0.813420268443792, 0.306447337889551, -0.614481533964195, 0.160069094336137,
                    -0.283217901825424, -0.454305917277928, 0.229572780395301, 0.214208832122054,
                    -0.628704545658053, -0.0533358363890102, 0.0945403308301422,
                    0.808178633387598, 0.0336309175807585, -0.983002963331188, 0.824448206677971,
                    -0.385014966255817, 0.908058762956387, -0.89515939341689, -0.168043980390365,
                    1.14250834519109, 0.0597442444081011, 0.506652232118796, -0.116596757423515,
                    -0.978841083153668, 0.64882952433075, -0.409295282931218, 0.408410439719461,
                    0.86145253107955, -0.844798856884864, -0.445903417996817, -0.490708621583941,
                    0.793364737901157, 0.363621730167559, -0.488778741038999, -0.724222735837543,
                    -0.94357483535515, 0.0669711561628219, -0.503883289160456, 0.453159633539542,
                    -0.61548887506595, -0.663509932324409, 0.697686687333686, -0.50214537651761,
                    0.194764953221039, -0.94865745075748, -0.695807452334854, 0.40641286611023,
                    -0.97977510312505, -0.778953830633242, -0.802208058666764, -1.25817847497834,
                    -0.181980367741062, -0.21970356592, -0.0885263527005691, -0.327507694861907,
                    -1.57580614666855, 0.300343835636085, -0.175525498545682, 1.15278730112351,
                    -0.97557127936247, -1.12510189601866, 0.430675895487212, -0.238951116980708,
                    0.546840511024868, -0.834896054962883, -1.22489520895986, 0.290903426527507,
                    -0.463519223197459, 0.754952988972036, 0.00622946568859106, -0.237834341603985,
                    -0.353420126850438, -0.758124278312749, -0.14720268474619, -0.562210226654372,
                    -0.320603874064545, 0.199438201568739, -1.22188539329161, 0.764149124991889,
                    -0.234859822041045, 0.523933323200943, 0.140254551011367, 0.521097081569522,
                    -0.671573605417535, 0.303408453099507, 0.727383087356081, 0.192909816774488,
                    0.0245272510170196, -1.47875639256696, -0.0989216305413549, 0.578123800768054,
                    -0.273915377245174, 0.104293262560539, 0.911347586502444, -0.0823700105403589,
                    -0.67519101527175, 0.454010405692117, -0.357941398079703, -0.0750976855032377,
                    0.22460544041067, -0.261494454298841, 0.735476503897502, 0.253824346247498,
                    9.07692475209387e-05, 0.191592033784332, 0.490133526304795, 0.350350426862519,
                    -0.336796928096511, 0.537788361169277, 0.313011361568243, -0.239036256005221,
                    0.774105154604342, 0.444989164881292, 0.401600108552969, -0.0736169705727434,
                    -0.39295931833766, 1.51654007434593, 1.33235747303463, 0.806995161004693,
                    -0.755562029075747, -0.487533263872876, 0.941592538286017, -0.196129511324848,
                    0.433288360495624, -0.685940817881288, -0.301274951213088, 0.987253944668055,
                    -0.986765938102934, 0.63161070847346, -0.0700106707669477, -0.34090264294575,
                    0.298524294439829, -0.657956435189696, 1.21161889942033, -0.115823817716972,
                    -1.00437838140382, 0.310743423789432, 0.255909091436607, 0.216923172127395,
                    0.286162427887177, -0.843553762411574, -0.0432308284350995, -0.867073683984017,
                    0.26917344096637, 0.392860887334564, -1.37536147925253, 0.0942486294877085,
                    -0.382809436548527, -0.223911950989901, -0.901501682044284, -0.528276005826541,
                    0.687118837517744, -0.947422717368234, 0.557780722106088, -0.638732698622306,
                    -0.125523502314863, 1.39661459554253, -0.50811823371199, 0.451809617288518,
                    -0.172889924782737, 0.0916746684980794, 0.478863183541307, -0.561441925369472,
                    0.0288749886317419, 0.0237956063350282, -0.0483259009054393,
                    0.203606735473869, -2.18117608503063, -0.549431390676205, 0.110406000229858,
                    0.126743102822338, -0.198230528845329, -1.00983502136452, -0.242678017768839,
                    0.19459243017954, 0.666365309204906, -0.0152679733537013, 0.0888848825875357,
                    -0.0474064306447133, -0.160352732486871, 0.390628755329004, -0.13231250229083,
                    0.589926737071202, 0.570404627000414, 0.139587545353635, 0.280549165587068,
                    -0.0167290784709451, 0.864486433540179, 0.103900644508087, 0.3705435194247,
                    0.367292644856178, -0.465602420498016, 0.223644717922134, -0.245153379927245,
                    0.075187209730899, 0.210559917035269, 0.230485811196814, -0.960546007402392,
                    0.49809248645883, 0.0179546129550658, 0.139812692083304, 0.389138857921694,
                    0.498077904440322, 0.339469761950653, -0.141109149365548, -0.191634628170299,
                    0.629919383918008, -0.153582004586355, 1.30568036877179, -0.663082868062048,
                    -0.976219544448576, 0.924633644546105, -0.586123872320066, 0.296758600437357,
                    -0.18182059950982, -0.764441652967811, 0.811686623064389, -0.48761865467673,
                    0.949296718026062, -0.738606784056373, -0.360557435285145, 0.125795940999581,
                    0.511346429652271, 0.600815742711262, -0.176571791341104, -0.524504286302808,
                    0.342300207368169, -0.568080769476261, 1.6174162433355, -0.159049116805932,
                    -0.415978633268786, 0.339505541293107, -0.788952105534705, 0.415387376642656,
                    -0.0988618759622284, -0.817404366800723, -0.285665121533214,
                    -1.09116055678949, 0.577129759569743, 0.812523983740391, -0.803623018458687,
                    0.566493945180895, -0.463135792384499, 0.472656320923619, -0.0773226859888609,
                    -0.185498497623044, -0.18517857911624, -0.783349981420869, 0.995618873179779,
                    -0.641733112933451, -0.0476185093601202, -0.394996577195087,
                    -1.13737648515044, 0.725683165646217, -0.0823790382572724, -0.0317504791152032,
                    0.124098328406291, -2.12830732593161, -0.0692499033476963, 0.579166707279711,
                    0.266194501737239, -0.27190333410914, 0.706238710702763, -1.21523399839936,
                    0.694486353777279, 0.754782709380464, 0.280761275415927, 0.907927909131317,
                    -0.374440565994519, 0.151188499532367, 0.281804347155578, 0.0375706418838604,
                    0.337459742391694, -0.268590660363344, 1.03069580327883, 0.944466305913551,
                    -0.167198435907693, 0.783287907473941, -0.5261608690804, 0.866394632822846,
                    0.781093157169522, 0.0894842634609461, 0.634749122975343, 0.0439659026867218,
                    0.543083780673612, -0.0103821617987876, -0.067024800421251, 0.321855603838003,
                    -0.503612769352145, 0.282842363171198, -0.238138644881666, 0.0960041295470262,
                    0.201833117245975, 0.0536222489483862, 0.25639088772067, 0.0848009274731538,
                    0.262516906752438, -0.220667177532597, 1.1718673471814, -0.0515266264224294,
                    -0.287922223713274, 0.331921987668096, -0.561177099197366, 0.709434891096237,
                    -0.272337932660207, -0.219278563872486, 0.20832137646002, -1.07296719245632,
                    0.769249115255548, -0.229128478742524, -0.150104937041505, 0.129110848957964,
                    0.125092051573326, 0.210704521388578, -0.163070482699396, -0.385691562488923,
                    -0.0906314761664646, -0.176053870136478, 0.481930263758017, -0.886809547482926,
                    -1.10918439239528, -0.42335895182876, 0.229931850110084, 0.570023289443526,
                    -0.576187026685204, -0.881555466287296, 0.587208864127007, -0.224731793409056,
                    0.30069471880995, 0.381703849808247, -0.516487755827882, -0.212984713798562,
                    -0.785130499479784, -0.151309881521343, 0.698242559065761, 0.109097233283954,
                    -0.7472044990445, -0.84232995101646, 0.839490829589749, -0.119238450432675,
                    0.0942935365796913, -0.673788190599282, -1.05286830122685, 0.780167809182416,
                    -0.0590533727079503, -0.273828870825113, -1.0959485462113, -0.592457655676562,
                    -0.0255839442687886, 0.686564647685628, 0.672930727508884, 0.433081342381485,
                    -0.435493208595197, 0.491660838962389, 1.35637437016823, 0.780934078781499,
                    0.300918251542925, 0.51669195809658, -0.0254990640375183, 0.762673239438755,
                    0.160742420215207, -0.30429306753382, 0.387006847328891, 0.0952961270787114,
                    1.32991285238066, 0.64872324092422, -0.353231689358267, 0.392707987666967,
                    0.469628788700277, 1.16046513278712, 0.302375608202766, -0.142938325576381,
                    -0.174626010177732, 0.57905855450403, 0.380609468941344, -0.0328542840570947,
                    -0.411824365799788, 0.144092464434159, 1.49838059553397, -0.565163934966534,
                    0.581643307337337, -0.0301536340533914, -0.0531692791160268,
                    -0.0497787589837912, -0.162243041594763, 0.23599903684871, 0.339761827772918,
                    -0.329413900930547, -1.19397117327804, -0.118848177173447, -0.534077444192611,
                    -0.276330700217983, 0.219149780471007, -0.34493480290908, 0.232694590464973,
                    -0.510297145075675, 0.372454682548209, -1.11478377744702, 0.671390358252552,
                    0.228208019517751, 0.492619659108733, 0.196783073372103, -0.24165991465631,
                    0.514625304022019, 0.0759428141978447, -0.0313902797872601, 0.232944814827619,
                    0.671592829681828, 0.327659856403528, -0.686023186344414, -0.812700869955513,
                    -0.0167318883783594, 0.440004229147256, 0.012273719059209, -0.0868537081913985,
                    -0.651805251733386, 0.600394219894611, -0.834963638550299, 0.442232095587144,
                    0.020889846811599, -0.606818007961131, -0.560993935061106, -0.280815160362131,
                    0.89147780775707, 0.267373899498324, -0.200237914428747, -0.292920955072503,
                    -0.434557514986686, 0.494250251373131, -0.110976633683492, -0.177524787405211,
                    -0.750105993999568, -0.408220672468595, 0.595822158819885, -0.352512298326644,
                    0.262539104822622, -0.97684216168696, -0.661057853811577, -0.144566836880443,
                    0.939889985156071, 1.05672285558705, -0.0789841874279489, -0.671851280542818,
                    -0.742215529015056, 1.06254559648949, 1.26660468544205, 0.586619357063938,
                    0.409358225438268, 0.0705189580447723, 0.718810849422255, 0.649393409372536,
                    0.0214123399490802, 0.315603030205324, -0.0731929349244638, 1.26299015029659,
                    1.1916601462454, 0.183834373455596, -0.71621854952781, -0.428407196245729,
                    1.05364240150396, 0.98338783079943, -0.50387879801093, 0.127117003376689,
                    0.247907875195551, 0.431287472393535, 0.72605135845366, -0.403231389808894,
                    -0.154954189143538, 0.504687086233872, 0.228203211389825, 0.344958842833329,
                    -0.315009400235397, -0.672581420230635, -0.977404070106275, 0.393083691422248,
                    -0.406706531189279, 0.11003065032321, -0.826708124645556, -0.562741545018014,
                    0.239619719356568, 0.107309749561743, 0.211816210324116, 0.522714331024368,
                    -0.125377222892212, -0.0655388076489487, 0.298330547569577, -0.11064279508247,
                    0.146805450548314, 0.0159803286882436, -0.77682950368569, -0.271723602453534,
                    0.353758890869585, 0.443027570684599, 0.805500445274902, 0.212174023381124,
                    0.288042793689783, 0.332740754704929, 0.925368191182777, 0.0894554811977173,
                    0.155692590751215, -0.154025881309707, 0.505045841347368, 0.0576386355712302,
                    -0.0763854031420124, -0.304038215839057, -0.0104310176978828,
                    0.663225185024999, 0.558008246295294, 1.30268583037153, -0.0388827533471954,
                    -0.621051106621316, 0.398077587490548, 0.47018099266235, 1.58445991898608,
                    0.569320066111459, 0.202878512132791, -0.489198242614322, -0.937061447761039,
                    1.10453275136112, 0.534292917336387, 0.0530934637897715, -0.957418842128448,
                    -0.599479326820784, 0.213471556904778, -0.198583608172726, 0.623577647519824,
                    -1.35605737667495, -1.29198144255977, 0.859194109270504, 1.25167126879629,
                    1.4988975805292, -0.294374347699545, -0.382321858945964, -0.0714029472023583,
                    0.973818416018673, 0.893713874716994, 0.435499441702149, -0.274295002892906,
                    -0.193069377084186, 1.1472428527083, 1.84650679742906, -0.0696557282555545,
                    -0.862488033630719, -0.780653319245289, 1.25532704326426, 1.40309659183483,
                    -0.13217425755596, -0.62588237442661, 0.10080841119796, 0.979743925112928,
                    1.11473110929938, -0.612962387090107, -0.482155180520735, 0.66708082772152,
                    0.808210512595622, 0.370639448185919, -0.607941124127763, -1.13428279690547,
                    1.50605072076209, 0.286173017179664, 0.614763989524924, -0.265815745969604,
                    -0.181522121611771, 0.480969415969042, 0.237133032150922, 0.71230549613302,
                    0.20809161127022, 0.472152102801054, 0.907038835829184, 0.716126732846845,
                    0.33512765956899, 0.707280423957798, 0.415107616234025, 0.302861773664461,
                    -0.173786183770318, -0.373734078707137, 0.368510239095826, -0.151597766333445,
                    -0.346806437396765, 0.0490916090507199, 0.0312140361111588, 0.478197744592421,
                    -0.110425624571569, 0.699150667377012, -0.328847330181418, -0.0378884676511504,
                    0.480543508965042, 0.182816724403743, -0.216218744427437, -0.431303539494792,
                    0.117992231282784, 0.163659254457428, 0.302996045505658, 0.500952736125612,
                    -0.294814571455163, 0.215218210997044, -0.0644303269204443, -0.281444265486932,
                    1.08054005037161, 0.620554582128847, 0.214112150765978, -0.176505195788309,
                    -0.34362332006166, 1.09005533062509, 0.404318367633891, -0.209329434695071,
                    -0.898991249064522, -1.13637095185548, 0.755453757937062, 0.3171129476994,
                    -0.0624912398818367, -0.719316266638288, -0.839093728848155,
                    1.4102169586653, 0.553919354624071, 0.871654628521931, -0.671834122645805,
                    -0.358045441082687, 1.51571295135957, 1.02649176548991, 0.35961941981018,
                    0.577085566014375, 0.59447523045035, -0.765073172388639, 0.871665226771739,
                    0.371411406306898, 0.472475872568282, 0.302517818320709, -0.704336326425249,
                    -0.198741007155348, 1.63532403004379, 0.214446133060306, 0.183433560265016,
                    0.105925253565118, 0.168261976271509, 1.76525073791918, -0.455062387959288,
                    -0.732425102253004, 0.586343734030805, 0.479567258468137, 0.85691191674232,
                    0.0238708627980486, -0.278178831921664, 0.182299267642435, 0.392287809312876,
                    0.734416806678794, -0.601051827469433, -0.75085099706272, 0.303845239659864,
                    0.134857713864989, 0.307374198863251, -0.76814498238727, -0.0481999371049966,
                    -0.199239391332156, -0.680086669332383, 0.445879542514668, -0.719957783802045,
                    -0.256937626745319, 0.142775051283053, -0.476214135388972, 0.10102237706036,
                    -0.100491546467453, -0.182601992304356, 0.11222465675462, -0.777618906123536,
                    -0.297339300261569, -0.116702244952052, -0.0424982965093221,
                    0.551345512997045, 0.151999992384447, -0.0049897028158199, -0.383083809941326,
                    0.112791045863105, 0.340114156371044, 0.0849386800455534, 0.278961401190813,
                    -1.20718155316419, -0.119009601257446, -0.626730362315755, 0.442736883719732,
                    0.154561728948096, -0.23941348350391, 0.252621339682051, -0.138372238698736,
                    -0.111388892265495, 0.541133778204582, 0.351118884613101, -0.204761186349138,
                    -0.0814720650011509, 0.196520271557603, 0.00407480292264495,
                    -0.143035717708983, -0.666709934753413, 0.0106966697999553, 1.00386709835385,
                    -0.617683681672361, -0.426637324343062, -0.49870634863391, 1.07301054166089,
                    1.02717357683835, -0.449426283197511, -0.176462068624204, -0.429240005796672,
                    0.203696518470451, 0.762345287068285, -0.252395584107127, -0.788570621332184,
                    -1.2557833011024, 0.0118096531772096, 0.964082660238404, 0.703333292726171,
                    0.154590464602677, 0.931547595280968, 1.01213156656675, 0.979680925938908,
                    0.618773229787621, 0.490480046987869, 0.311584181585616, -0.363422728962263,
                    0.0971190311127543, 0.690229207873969, 0.0834585414749824, -0.15408575698963,
                    0.362757931380656, 0.0997308722854768, 1.43400410826469, -0.419646137301925,
                    -0.00501700753366663, 0.481003197683378, -0.38151649477058, 0.556527825055412,
                    0.0557557666546984, -0.30432889432732, -0.986934255960904, 0.276127536380699,
                    0.946110331953967, -0.871665156490536, -0.539903640069561, 0.339736837928854,
                    -0.517505406152168, 0.00634912308109392, -0.684590648409767,
                    -0.658167906575365, -0.236480823853498, -0.252155892586757, 0.0828373630787951,
                    -0.575083214927616, -0.975037174121482, -0.0681590508762611,
                    -0.0564331897828936, -0.212299709714948, -0.256593013254019,
                    -0.406162821435963, 0.119312710729218, -0.364754851994348, -0.666149675121858,
                    -0.22093193778295, -0.164008920963155, 0.742871163861821, 0.25789514338615,
                    -0.390484053432779, -0.342640175888022, 0.0977079375398873, -0.463119165228647,
                    0.0807992881846076, 0.00756880797538884, -0.134523479947234,
                    0.177863971392775, -0.138505829689166, 0.163974387232344, 0.503884635972263,
                    -1.18209441066199, -0.213073536151157, -0.363729808767914, -0.024533307583738,
                    0.22520511618292, -0.0164179372149533, -0.0109034279424662, -0.267145966809323,
                    0.228468673805565, -0.621112861724122, -0.280813567491407, 0.44285760188292,
                    0.497776436695034, -0.380519934569037, -0.649131907445195, -0.670660915136864,
                    -0.680134230710991, 0.0858117110656167, 0.333486123889912, -0.0967024606118869,
                    -0.523959459916, -1.0308350702059, -0.613512198570589, 0.701741541868302,
                    0.452698521821155, -0.630031990012469, -0.915784286084644, -0.362000145163847,
                    1.29766683446729, 0.937985379757968, -0.0938134749114283, 0.24425053015316,
                    0.591569446451473, 0.345957032991322, 0.126782712707261, 0.679807769256627,
                    -0.132930931258059, 0.253767263679038, 0.378320724624672, 0.453699632296455,
                    0.236346361368618, -0.00751351256037935, 0.709108988294702, -0.16404512115534,
                    0.605624223233086, -0.443486521918999, -1.18745653788424, 0.354158170952019,
                    0.0367429468271609, -0.0103507535434062, 0.319868748948437, -1.12040991134545,
                    0.211050880365974, -0.339579806037136, 0.511177920665901, -0.353251698787716,
                    -0.940105278726968, 0.533710133167006, -0.206060147480883, -0.142970160086534,
                    -0.0601166634018853, 0.363718899154058, 1.02365219824679, -0.436090622564819,
                    -0.216895215911428, -0.208807504772864, -0.312912606597436, -0.700987876647247,
                    -0.342851361573423, -0.115082989673994, 0.275605109829172, 0.002511332284584,
                    -0.214095227960172, -0.226563501972862, -0.316148066983203, -0.369959145047261,
                    -0.439875798314965, 0.0457626161993403, -0.020300966062706, -0.865573293898194,
                    -0.186393733787668, -0.0401816237273778, -0.67158681490875, -0.887609053416696,
                    -0.0943908861927127, -0.18976698243037, -0.33939468061584, -0.043330825982224,
                    -0.1003967504804, 0.237185186146107, -0.0900626742932218, -0.199187593768066,
                    0.188882001248288, -0.459792205444785, 0.296219435260724, -0.295062680878516,
                    0.255123910597142, -0.280068963137438, -0.331017673328432, -0.908523340392335,
                    -1.15462440123483, 0.59931591026541, 0.066441615501958, -0.454785219621385,
                    -1.33584474773882, -0.554898504808506, 0.105109946902367, -1.54921658450326,
                    0.0457334440739032, -0.486904581871869, -1.08521496224813, -0.99396185317854,
                    -0.14726764429858, 0.889996169054569, 0.356599799698657, -1.17926106346533,
                    -1.53392460942358, 0.892779422358217, 2.21309267170048, -0.513970751630012,
                    0.220952739442595, -0.533054292513867, 0.7453668341482, 1.6612485458385,
                    -0.0470186581734637, 0.182498698103132, -0.441266479673551, 0.61134912630888,
                    0.710910011413404, 0.246247545644067, 0.00258435765027798, -0.307654895560297,
                    0.432594086409352, 0.60541800613106, 0.578855402223397, 0.033476291343419,
                    -0.368763257777069, 0.320989097122549, 0.611331932612046, 0.990465834578067,
                    0.0510383877596894, -0.991079818723577, -0.315399530875425, 0.0576650484925949,
                    1.01343500126922, -0.245035769685842, -0.492404404362055, 0.495956524298663,
                    0.0647507015438293, 0.443258359063122, -0.135730899741816, -0.225300068683964,
                    0.041759637449843, -0.226147496889098, 0.187620599607854, 0.553149697236648,
                    -0.824741824129497, -0.84688305633756, -0.550683481409152, 0.127122546035189,
                    -0.292289389981079, -0.770135937916677, -0.591098608741049, -0.539844111445319,
                    -0.201169703733676, 0.0916061274652691, -0.378856900668339, -1.10669209549062,
                    -0.287058717650297, -0.491933417747212, 0.237749126901747, 0.745784312608379,
                    -0.652446108432002, -0.738779585023536, -0.250024987663316, 0.0585933956253371,
                    0.506200344170608, -1.33156587609616, -0.962803604838456, -0.0890884170415163,
                    -0.163796800008669, 0.441053397967565, -1.48654665201693, -0.80656092911201,
                    -1.28435271237538, -0.156321233800384, 1.1472411314269, -0.995607279603092,
                    -0.269558855902517, -0.96235808323852, 0.324594967757662, 1.23433903637621,
                    -0.679750551880517, 0.137497684122792, -1.02360716571832, -0.248982869422861,
                    0.0384039501975053, -0.719396077347284, 1.15724526474307, -0.445343076692531,
                    -0.50959266624929, 0.606547294306807, -0.447839155820645, 0.533136050537763,
                    -0.237344481267201, -1.14585127504409, 0.169571953250802, 0.816179983673152,
                    1.5107580481494, -0.171536655258748, 0.0571663119985431, -0.252046572495196,
                    0.852116133188872, 1.06555902002216, -0.291297771979941, -0.8283298450541,
                    -0.169155168025881, 1.33857056030237, 1.05933539943246, 0.0366795308673389,
                    -0.979462608221439, -1.58143732941764, 1.13096274134332, 0.991243591315765,
                    0.401784055613255, -0.85663237062866, -1.3277127690003, 0.570699910073756,
                    0.649004496181021, 0.657034876239841, -0.0958123653511223, -0.874107398035221,
                    0.431116788537054, 0.390059295249756, 1.02059982141703, -0.420029286739265,
                    -0.623438958730461, 0.226612623043621, 0.717360623288698, 0.0980297184690428,
                    -0.398231604357107, -0.829877044173116, 0.418002773909421, -0.297099266440142,
                    0.333991631342215, -0.147323981468816, -1.12855905648352, -0.177813781887693,
                    0.169867106858334, 0.0485581276107467, -0.17790439598312, -0.407373166899582,
                    -0.423990509228105, -0.239835557065874, -0.197295848262179, -0.835093427898244,
                    -0.116256399164293, -1.09946343877832, -0.442817586792245, -0.33500466085696,
                    -1.02197982474823, 0.87759105370997, -0.410373425641296, -0.495705095708786,
                    -0.859263030356438, -0.872286723005033, 0.82871154414633, -1.11200773274004,
                    0.192986884068091, -0.499275300656215, -0.548477768011111, 0.742333381242751,
                    -1.11100023331488, -0.64266431592975, -1.00905620672707, -0.446893086432567,
                    1.23738152366356, -1.29052612543003, -0.627510499810654, -1.07516675709466,
                    0.0171862943857093, 0.29180006004561, -2.19975580988569, 0.269002928653435,
                    -1.14915200253675, -0.154694531182912, -0.959971049796188, -0.911673217320697,
                    0.826140839292322, -1.31663677430553, -0.221551639047789, 0.180125805785686,
                    -1.19546190209022, 0.512686863313917, -0.563033051898756, -0.800696680818533,
                    -0.161478996925005, -0.410072158320133, 1.23627126098276, 0.189901107681898,
                    -0.0659779617294589, -1.06644230809871, 2.63525825781343, 1.09964310161661,
                    1.03155022726709, -1.40417052112079, -0.813358803580066, 1.72641505749601,
                    0.59712322806692, 0.984236455229897, -0.730364684161572, -0.725601308022732,
                    1.8023790204771, 0.722562411068452, 0.693670088499069, -0.0324215947672712,
                    -0.647707798925659, 1.63521215647197, 0.772269309682116, 0.563295784259118,
                    -0.25428721436495, -0.171183800310956, 1.04260003152947, 1.33329554465834,
                    1.13704067215005, -0.614591915547663, -0.437532538693986, 0.673523311603731,
                    1.24349356197379, 0.741535621725966, 0.247034837149243, 0.736897112218983,
                    1.34615640310091, 0.0155252558691359, 0.93130916209001, -0.0258271510697909,
                    0.194334019902598, 1.0336758046856, -0.299197239532987, -0.20585065835138,
                    -0.590546009914561, 0.203475059299372, 1.00707007910754, -0.0289207450863671,
                    0.0150308830855437, -0.4436802738783, 1.128621656783, 0.330403162515353,
                    0.0568668247268533, -0.207476809129105, -0.480825101494759, 0.638693700823353,
                    -0.051197335791926, 0.101329098744571, -0.333195209040063, -0.84415913728796,
                    0.22759064737496, -0.153328786114287, 0.38994807455029, -0.578485088964655,
                    -0.888174982609939, -0.0730232967309377, -1.50688912829452, -0.552727595402018,
                    -1.31511263427068, -0.909643410047181, -0.00669237338742334,
                    -0.896581831396787, 0.296108756204291, -0.973172695487833, -0.127224242269932,
                    -0.890625987200053, -1.05343510258291, 0.492086814427125, -1.0869260716696,
                    -0.476890999185372, 0.237895979537683, -0.786477930452492, 0.64349412313868,
                    -1.59807339151561, -0.383274006481116, -1.0356082832164, 0.131788699840882,
                    0.77598045065119, -1.08331963605532, -1.32290187144537, -0.859913138984695,
                    0.44316054105397, 1.18668739057233, 0.444004162609282, -0.0451487835553053,
                    -1.19539791563399, 3.02155272957654, 1.74086264822009, 0.323980695098349,
                    -0.246709419186743, -0.253542616508405, 3.2076686999894, 0.513163068063648,
                    0.457700725834613, -0.10961373104729, -1.30729056136575, 2.64865809373126,
                    0.53855572936186, 1.15269107481435, 0.0103185775509758, -0.576627569577013,
                    2.5173249318097, 0.504904350489312, 0.771086519025886, 0.354825980871681,
                    -0.716827062406132, 1.71991175451274, 0.705887252380109, 1.03082391284737,
                    0.376879185538428, 1.13326674902294, 2.28071531410527, 0.529900151648178,
                    0.907188710785448, 0.854307501644968, -0.171398066989269, 1.49772151655796,
                    0.341986292256593, 0.340153338094488, 1.299849404805, -0.734400177448713,
                    1.20475870633078, -0.201226410546809, 0.266478912198714, 0.757294474797163,
                    0.681748899061091, 1.25580318993556, -0.254828871586824, -0.397597552085152,
                    0.561835171329875, 1.16568866961883, 1.29514730000705, -0.134790656953982,
                    -0.0556172536546029, 0.0273849607399654, 0.638335904275458, -0.328800696543311,
                    0.0375905638397569, -0.567953379286458, -0.252892365065065, 0.118454739246559,
                    -0.524898099664739, -0.0188998678876544, -0.889525875383491,
                    -0.652593724718287, -0.0750956527861177, 0.16099321568268, -0.353653839534303,
                    -0.754426919738137, -0.448894702123198, -0.0973081180126272,
                    -0.439169475803936, 0.0352668204598103, -0.720979442229023, -0.522454379263735,
                    0.485449999559549, -0.623886566361904, -0.117311941703747, -0.605739097310888,
                    -0.928936478505733, 0.147332561685602, -0.333550395640117, 0.117181366652858,
                    -0.28718079867722, -0.450077203438099, -0.205145918600728, 0.020798512926433,
                    0.796083087216743, -1.29029856521649, -0.569936566547028, -0.221090901905264,
                    0.488037172789006, 0.770532487800892, 0.440172758949515, -0.188549354911743,
                    -1.61489315230659, 2.01434644728803, 1.2891129520095, -0.143492773842965,
                    -0.524701767237246, -0.830018862658456, 2.57061390417135, 2.10264932717172,
                    -0.392013779591267, -0.653885214786243, -1.20107338950063, 3.02452236118319,
                    1.93022211003131, 0.561850074041971, -0.605830023978314, -0.642457159521669,
                    2.5984399477579, 1.26889511441371, 0.890909061731087, 0.139958375677054,
                    -0.666227264597929, 1.68291922985429, 0.61455549533224, 0.986844635194537,
                    1.48221920956792, 0.494436151595053, 1.00694953089646, 1.16129574076527,
                    0.522998053596794, 1.38430655226659, 0.437466363408014, 1.09680766602861,
                    1.04277456013821, 0.0563657911390507, 0.985696347544859, 0.155922221156246,
                    1.42685790798464, -0.102224537831205, 0.0698453485904622, 0.496856134001717,
                    0.584891567001208, 0.890820055721321, -0.599706855839401, 0.00603403514570289,
                    0.19482719595058, 0.134694973416131, 1.04960988172044, 0.27341495736124,
                    0.146085244603931, 0.56730883541934, 0.140923520397443, -0.32099332546224,
                    -0.998191896012031, -0.482483929782525, 0.359992782481214, -0.0727322837054771,
                    0.207461543496464, -0.0796572744944797, -0.419894129995946, 0.104998455441091,
                    -0.448189651289404, 0.579790458039213, 0.021033388141376, -0.368446348790506,
                    0.0151354606229575, 0.444419067429109, -1.07869947921861, -0.245513021376155,
                    -0.773698242616931, -0.335021423258297, 0.470904821955109, -0.428771993282652,
                    -0.46084036457728, -0.4895157659451, -1.00074832414325, 0.786853781134219,
                    0.575517112930214, -0.237970617564201, 0.169614306464688, -0.874855676526181,
                    -0.149200021034453, 0.369825581740628, 0.560392109502914, 0.156446175559132,
                    -0.747982212420567, 0.197117144232983, -0.490070648530633),
                  .Dim = c(5L, 520L)))



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
        structure(c(-6.11420734746121, -6.90052639635516, -7.07895510496375,
                    -4.46692474389574, -6.13571374777994, -0.354422064854689, -0.403503410506908,
                    -0.424436653906349, -0.258765608581638, -0.365787493819805), .Dim = c(5L, 2L)))

    expect_equal(
        unscale_sigma2(out),
        c(2.02176705952465, 0.798706141860752, 1.01269594272481, 0.535766835394645,
          0.889899599677774))

    out <- mcmc_mra_integrated(y, X, locs, params, priors, inits = inits, M = 2, n_coarse_grid = 4)
    expect_equal(
        unscale_beta(out),
        structure(c(-0.669013535873582, -0.669013535873582, -0.669013535873582,
                    -0.669013535873582, -0.669013535873582, 0, 0, 0, 0, 0), .Dim = c(5L, 2L)))
    expect_equal(
        unscale_sigma2(out),
        c(0.615786367308992, 0.615786367308992, 0.615786367308992, 0.615786367308992,
          0.615786367308992))

    class(out) <- "aaa"
    expect_error(unscale_beta(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
    expect_error(unscale_sigma2(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
})


