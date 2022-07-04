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

    params$n_mcmc <- 10
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

    params$n_adapt <- 10
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
                                         0, 0), dim = c(50L, 3L)),
             Sigma_tune = structure(c(2.97876608285749,
                                      0.412648437574117, 0.432599330050228, 0.412648437574117, 3.87407264820285,
                                      0.587518320732799, 0.432599330050228, 0.587518320732799, 4.70641828672913
             ), dim = c(3L, 3L)),
             Sigma_tune_chol = structure(c(1.72591021865492,
                                           0, 0, 0.239090326434078, 1.95369098477945, 0, 0.250649961611196,
                                           0.270047998225292, 2.13791182281402), dim = c(3L, 3L))))
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
                                         0, 0), dim = c(50L, 3L)),
             Sigma_tune = structure(c(2.97876608285749,
                                      0.412648437574117, 0.432599330050228, 0.412648437574117, 3.87407264820285,
                                      0.587518320732799, 0.432599330050228, 0.587518320732799, 4.70641828672913
             ), dim = c(3L, 3L)),
             Sigma_tune_chol = structure(c(1.72591021865492,
                                           0, 0, 0.239090326434078, 1.95369098477945, 0, 0.250649961611196,
                                           0.270047998225292, 2.13791182281402), dim = c(3L, 3L))))

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
        structure(c(-55.7413786355926, -120.397997808002, -60.5073216554941,
                    40.2287465535353, -391.52340363644, -55.3071576106442, -122.06358370527,
                    -60.8181089915441, 41.4601513129057, -391.526798263586, -54.9296863741474,
                    -122.116270095828, -61.0835728999231, 40.1437395471594, -392.196851077704,
                    -55.6916168505625, -122.127248091819, -61.4280435091201, 40.5199697179713,
                    -391.905384235902, -55.7901911047978, -122.161131564182, -61.190095935052,
                    40.3221659849157, -392.617192170157, -55.8635555549914, -122.216548960483,
                    -61.345927054694, 40.6304866315159, -392.958756995231, -56.4363063530628,
                    -122.596423051482, -60.6971236398901, 39.8418154467378, -392.6677611016,
                    -57.3412518194609, -122.201538968865, -60.7619747774417, 39.4618417222475,
                    -392.786002151808, -57.3136943073188, -122.611427917883, -61.397726312415,
                    39.9426183334513, -393.398335778229, -57.6060035567961, -123.168455878075,
                    -62.0237382460533, 40.4983262403619, -393.776345949422, -56.7961736828601,
                    -121.735264447158, -60.5973404032608, 39.6977984710337, -393.414095390959,
                    -56.7970038511444, -120.784074139994, -59.9286986812442, 39.1937857359592,
                    -392.924066148986, -57.4161435785065, -121.825171683805, -60.909711245295,
                    38.4689700923321, -393.351106673855, -57.6660881012975, -121.029305009045,
                    -59.9729311242236, 38.3777796358401, -393.42345922259, -56.5144596938562,
                    -121.582535006642, -60.0736252229986, 40.3637117466486, -391.206590608003,
                    -55.8173752434951, -121.24480497185, -60.5533285258018, 41.0872850269738,
                    -391.727086402727, -56.3937905872642, -122.005181333795, -60.1937860689367,
                    39.8841017985187, -391.959256928297, -56.5049496636514, -122.337780360058,
                    -61.0326467477304, 40.4201702135741, -391.739778106998, -56.0554384768086,
                    -122.114331084741, -61.4388914075504, 40.5588103286877, -392.566478425026,
                    -56.6606712675465, -121.976204908192, -60.957100057498, 40.3149781747688,
                    -392.257655203355, -57.2491008291326, -122.29464579614, -61.2494000057945,
                    39.76817389636, -392.701816453718, -57.6969889710255, -122.641799764267,
                    -61.6881055887692, 39.0872435271544, -393.010121852504, -57.2131158081545,
                    -122.242500630303, -61.251390929922, 40.1000202984715, -392.199539892977,
                    -56.9103294931197, -122.013944620863, -61.3808974892973, 40.3708922513825,
                    -393.464404440682, -57.2106425200307, -122.110955945347, -60.7176609536074,
                    39.3896937650642, -393.451999029011, -57.7036638101613, -121.627796754553,
                    -61.0618271348651, 38.9499788573248, -392.676963045824, -57.992214496565,
                    -122.326174541915, -60.4175054890317, 38.9742143977396, -393.196086930595,
                    -57.7833312192903, -121.556731010002, -60.845007674096, 38.7821426498925,
                    -392.843045636403, -56.7210564549941, -121.505839367115, -61.1093041207737,
                    39.1799868835695, -391.88985308275, -56.2811046614376, -121.858181456847,
                    -60.1616351417276, 40.1158125667916, -391.7012109666, -56.9370606515037,
                    -121.573214769566, -60.6393436783314, 40.0685327365464, -391.801689869026,
                    -56.831162837391, -121.695539839476, -60.7699380785708, 40.3986115878442,
                    -391.733803382586, -56.7355588083469, -121.984450524882, -60.980073424283,
                    40.1389370895812, -392.178203807032, -56.9231898113335, -122.238248203191,
                    -61.0559538574352, 39.6465416159737, -392.676362809342, -57.5296361074082,
                    -122.859124697033, -61.1083995333637, 39.9169791643099, -392.270299880068,
                    -57.7609759844901, -122.048948944002, -61.6963536771427, 40.3362594352182,
                    -392.659032675467, -57.2662962457858, -121.745424622488, -60.9758954597101,
                    40.4560685361592, -392.22609588679, -57.5784657919016, -121.743572926083,
                    -61.2057499028299, 39.8703405900925, -392.917231675397, -58.1350466397889,
                    -121.511235844667, -61.4774368220503, 39.025883710902, -392.678608000355,
                    -57.6912859756921, -121.758643946906, -61.1743432921329, 38.975234025561,
                    -392.722680900143, -58.0898657647473, -121.70917082115, -60.9182708953376,
                    37.9964856422404, -392.784835061214, -58.7051083198235, -122.512710350092,
                    -61.4478743671992, 37.7251954611625, -392.343876245246, -56.8882955691105,
                    -121.673359656923, -61.4688801029411, 40.7112459013008, -391.642856249102,
                    -56.5273221065341, -122.013989476164, -60.9366928373805, 40.0357624923904,
                    -391.032796092964, -56.1986518505785, -122.165266566712, -60.5038820477298,
                    40.3931701211333, -392.67708030965, -57.0049181718763, -122.228319390715,
                    -59.6224414485488, 39.6725165922504, -392.215952780478, -57.4493919433347,
                    -122.213389113138, -60.8351303726086, 39.2174702635679, -391.932674157906,
                    -56.3780726054856, -122.503948132302, -61.3787054232504, 39.4328186016168,
                    -392.297950993304, -56.8373552197226, -121.93939191923, -61.30358061344,
                    39.3848253158972, -392.253939438399, -57.4129411677633, -122.09146237739,
                    -61.3930825474689, 40.0019103815829, -392.139943349336, -56.7146395845286,
                    -121.165986773616, -61.5760002930151, 40.548668555753, -392.552940191632,
                    -57.4578413737373, -121.671947297901, -62.2616742604603, 39.9070548694529,
                    -393.008000854072, -58.1369563682681, -121.096126602408, -61.7847418724234,
                    38.7760845436551, -392.953736191798, -58.3250454119958, -121.846733393779,
                    -61.1355890114669, 38.5401641659439, -393.053885142929, -58.1595539682482,
                    -121.268685942023, -60.9500757423912, 38.5524104992738, -393.514089739042,
                    -58.6821818298643, -121.529528321149, -60.4362365717097, 38.6790648709246,
                    -393.134691115631, -56.8431934924598, -122.463172944714, -60.4785870465444,
                    41.5696456782836, -391.689652492694, -56.778424571016, -121.974987283305,
                    -60.2294731353409, 40.1116139842586, -391.336952576963, -57.290971202693,
                    -121.865778943152, -60.8863590221112, 39.9640667462759, -391.757880546192,
                    -57.1265769768168, -122.308373756515, -60.1118239266407, 39.4308567417701,
                    -392.356886520366, -56.8656872445314, -122.584841224696, -59.9051795062624,
                    39.3206400641795, -392.129640970767, -56.1966620001608, -121.640964764727,
                    -60.7785015487909, 39.1610686946324, -392.004755905792, -56.3281717621967,
                    -122.011662699231, -61.0443147685162, 38.9217456540037, -392.453841110885,
                    -56.9679233656307, -122.870545397558, -61.4763702872744, 39.7575025643557,
                    -392.193816929559, -56.5423448481881, -121.403244230028, -60.5695818351292,
                    39.8243977182443, -391.842032286666, -56.9618238220601, -121.307786194936,
                    -60.6900889277035, 39.8217391899102, -391.551474174942, -57.5750992360719,
                    -122.159567131468, -61.5694800797976, 39.4913591072263, -392.353958216651,
                    -57.4392787270497, -122.257450901962, -61.0719693665193, 39.7227314595445,
                    -393.158012290528, -58.4082551544691, -121.607993248204, -60.9749984826105,
                    39.2534756983471, -393.01592634787, -58.152916983382, -122.522163149989,
                    -61.7282153303861, 39.0707901376129, -394.006451003858, -57.4127771623812,
                    -123.025173606098, -60.6822379398436, 40.3369316820385, -391.528182554896,
                    -57.4368206861074, -121.496276580972, -60.8916132753648, 40.1256136565728,
                    -391.596608057322, -56.8109997683702, -122.1546127733, -60.4472888053975,
                    39.7016761947304, -392.229254693365, -56.6069343061311, -121.814981229489,
                    -60.3686154032532, 39.5142035857403, -392.477644560827, -56.40690088613,
                    -121.740237188517, -60.0740165355918, 39.5904966942067, -391.341359307486,
                    -56.2151682295834, -122.799595845356, -60.274733755626, 39.8024930761115,
                    -391.701849441469, -57.5376731955874, -122.975726815376, -60.7641524019836,
                    39.756820246743, -393.134582801652, -56.9492708596448, -122.706017923021,
                    -60.9758639801964, 39.7828700574124, -392.594884587282, -56.7732595104739,
                    -122.275531162885, -60.49893002737, 39.8613333677432, -391.744715947439,
                    -57.7527109293117, -121.90904366151, -60.2702921034023, 40.0259426032409,
                    -392.234294618311, -57.0598933937711, -121.906549801472, -60.6754977408629,
                    39.4904715862888, -393.053094083358, -57.453796215868, -122.472911066973,
                    -60.8571248507405, 39.5723753415354, -393.506473111461, -57.625757803891,
                    -122.506973194034, -60.9924574633719, 39.3277500298224, -393.239897056946,
                    -57.2085627896463, -123.103060907557, -60.9945664671716, 40.1679264968118,
                    -394.052792409369, -55.8342403451934, -122.166329428183, -61.4819118304808,
                    39.6344010894691, -391.539006392527, -56.6815375648364, -121.837215915867,
                    -60.5535994118298, 40.5252787085221, -391.483328797646, -57.6524036979367,
                    -121.793085075127, -60.612070627521, 39.7077292541117, -392.479426842096,
                    -56.9682489814125, -122.648313278795, -61.1500142701836, 40.0136668260889,
                    -392.22338837038, -56.8408148336886, -122.606104605566, -60.2839517830592,
                    39.2480247625376, -392.18196111966, -57.3643209886408, -122.842869938855,
                    -60.6144281969126, 39.318896556351, -391.96412326862, -58.2522589647493,
                    -122.937081174793, -60.2349380259446, 39.8190141476524, -392.804713313324,
                    -57.4908996116367, -123.119864582451, -59.9457632987439, 39.6263672306684,
                    -393.281928924691, -57.2491227354529, -121.817165126546, -60.5036494028121,
                    40.0627357991879, -392.763620510632, -57.5839193077747, -122.203321986351,
                    -60.8981042810281, 39.4895627486883, -393.324808486894, -57.3779456480065,
                    -121.663611791363, -61.3113829006442, 39.6328434409623, -393.024662798763,
                    -57.2696373119026, -123.271343734168, -60.8774635432193, 39.0004464405458,
                    -392.849117627945, -57.7003152860741, -122.984915898157, -61.3292473474853,
                    40.0181654211569, -393.774194535259, -57.190624548492, -123.232777826485,
                    -60.7606083470565, 40.0155343783855, -393.903686731769, -56.6548894549038,
                    -121.840152000119, -61.19896403139, 40.1899474948187, -392.012852186997,
                    -56.6234001532473, -122.191028705712, -60.9443032879474, 40.09359206011,
                    -391.518640711352, -56.4954005554619, -121.425221309155, -60.6257775082322,
                    40.0249994841223, -391.81236836863, -57.0039761294853, -122.241007759687,
                    -60.6430066671999, 39.5641785215922, -392.175208329291, -56.7290358758157,
                    -121.624222825937, -60.129249363516, 38.6383427486633, -392.520102685267,
                    -56.5936616078684, -121.721498740753, -60.4308218400916, 38.6149297524634,
                    -392.815226321366, -57.1687659113176, -122.001005697577, -60.6851104202997,
                    39.5219416675945, -392.454420608275, -57.881447554133, -122.409328914603,
                    -60.0655062468979, 39.9173367741305, -393.202001082612, -57.4859641848568,
                    -122.207312876731, -60.5215138567257, 39.7343404410406, -393.194270148619,
                    -57.3495388706531, -122.598018999049, -60.7075827705102, 39.1956813600119,
                    -393.643208764506, -58.2071800421475, -122.144980387017, -60.5919103764103,
                    38.6268229837141, -393.886697270697, -57.4661203020432, -122.626236696765,
                    -61.0703786565365, 39.4574826266922, -392.984569583799, -57.5347572482781,
                    -123.297456522613, -61.3354797737867, 38.9017095015526, -393.759879068954,
                    -57.1826019076194, -123.55653542416, -60.8913526728113, 39.7900259901191,
                    -394.208707199429, -57.3554739846146, -121.559959129097, -61.0206933794756,
                    40.0275292636325, -392.270110033558, -56.692369800454, -122.621339554488,
                    -60.7707104188249, 39.9365471696956, -392.041328989158, -56.843936888702,
                    -122.28592734223, -60.8037478170418, 39.854608934565, -392.422734545917,
                    -57.1407298693731, -122.051139759028, -61.3205160930725, 40.2793080544455,
                    -392.879012891288, -56.4752071306908, -121.911970316099, -60.8187034966285,
                    39.0020735676659, -392.683353427027, -56.6601622294836, -121.612048163518,
                    -60.9143614437433, 38.9477005318394, -392.99217148055, -56.9203812035697,
                    -121.407455764501, -60.5916792133593, 39.5097530573285, -392.838456917664,
                    -57.1237256702909, -121.73241833012, -59.5961307950848, 39.9456937112319,
                    -393.066330109017, -57.0312432740228, -122.435126005992, -60.5585684955294,
                    39.5897477100166, -393.494753401824, -56.7974751424876, -122.960105788321,
                    -60.5660700871251, 39.2060127220491, -393.56532531253, -56.8758904248612,
                    -122.9317324795, -60.3556795480481, 38.9304630368855, -393.477787864323,
                    -57.7840733863102, -122.466182796342, -60.7189590956987, 39.0053326450422,
                    -393.138234622652, -57.172028361569, -123.003363009292, -61.0679864051791,
                    38.9807551766558, -393.546158835272, -57.1929999767468, -123.352850382368,
                    -61.6181519249724, 39.5142183087557, -393.643721694827, -57.2315623110534,
                    -122.065110247245, -60.7583264201527, 40.6854229703303, -393.551985971907,
                    -57.0726812126323, -122.039608714183, -60.8395753938078, 40.2350397742928,
                    -392.793165661315, -57.5702378542138, -121.999777351128, -61.0271306048749,
                    39.4868410644438, -393.177705755174, -57.0624311977852, -122.823946244428,
                    -60.9872181956851, 39.8618937567479, -392.843953203553, -56.6715184449334,
                    -121.612309398754, -60.7843309107288, 39.7905832436454, -392.878286131171,
                    -57.0269279884332, -121.904882231476, -60.8754987734228, 39.55165924184,
                    -393.287924808833, -56.9636876696692, -121.640597720208, -61.2486857984522,
                    40.0113928583189, -393.193587200574, -57.5972232290888, -122.672040802689,
                    -60.370568053928, 39.9735313867751, -393.390521885768, -57.7991153839552,
                    -122.580517742677, -60.4218598078979, 39.595221095069, -393.482130230155,
                    -57.4099141294825, -123.117008103142, -60.1000305769294, 39.6840354214756,
                    -393.898502704665, -57.6510960418328, -122.905166436477, -60.5914170778744,
                    38.6143246354957, -393.217126272206, -57.4415028209437, -122.505207818508,
                    -61.1634611599771, 39.1628515494147, -393.297674325297, -57.9649214513549,
                    -123.263608734873, -60.6780598244143, 39.1497350834986, -393.774974365764,
                    -57.443466075289, -122.385439045051, -61.3585929223476, 39.8188432862622,
                    -394.298318174013, -57.0172728342646, -122.217073971811, -60.6843715802961,
                    40.4096003472372, -392.860947652639, -57.1945354877953, -122.07396644943,
                    -61.4460300022078, 40.4904477796773, -393.206964792221, -57.8299771087983,
                    -122.186995904957, -61.5846616700414, 39.6050673664532, -392.822928778908,
                    -56.8526596399475, -122.637232259663, -61.155640793162, 39.3933104335198,
                    -393.569209084721, -56.533470364084, -122.055505715989, -61.7883623182187,
                    39.6543711798199, -393.451917757705, -56.8045136956906, -122.378944660579,
                    -61.6070955474577, 39.4126284251564, -393.359748023195, -57.4310526996724,
                    -122.233247063753, -61.6077903715584, 38.8921692388225, -392.813223995383,
                    -57.3408593396966, -122.247515696846, -60.8585477655906, 39.1391097138827,
                    -392.714476150565, -57.4752675133853, -122.251250202908, -60.676753552138,
                    38.6658714814425, -393.634491059954, -57.3672510838691, -122.66604880261,
                    -60.6169057639135, 39.0303433996731, -393.866273056034, -57.5304143718395,
                    -122.415065906216, -60.7488137450954, 39.2067999223947, -393.989855234977,
                    -57.6255559648784, -122.545480067365, -61.3281715109165, 39.3179895528001,
                    -394.002157986931, -57.5252517779446, -122.506900276757, -61.3116476577074,
                    39.1735265317395, -394.074863727809, -57.2893131603284, -122.384915710685,
                    -61.3205108812342, 39.9285806842707, -393.961533970457, -57.3149019121101,
                    -122.483849985506, -61.5344871491014, 39.013136657023, -393.007076961373,
                    -57.1969009182517, -122.164301192907, -61.7339311789991, 39.0966322259338,
                    -392.708654283632, -56.6151675722517, -121.721774049746, -61.080822212452,
                    39.7522744399789, -393.328414745018, -57.6005903637596, -122.181904537781,
                    -61.7636184735282, 39.6269665948594, -392.529369532828, -57.3629626003119,
                    -122.313365766679, -61.5472909596813, 39.7215119577641, -393.551414226128,
                    -56.8242535866909, -122.606468634255, -61.7616341505536, 39.6556008259355,
                    -393.264920004044, -57.0544460546855, -122.723394175067, -61.3191778673301,
                    39.6809440279807, -393.655259683753, -57.2374575087442, -122.177248359083,
                    -60.9891666722403, 38.6036241220346, -393.496337622943, -56.7313031903822,
                    -122.500440675854, -61.167389938668, 39.0684368122197, -393.714205788517,
                    -57.0812045913093, -122.472622026384, -61.0561016525799, 39.2492639502708,
                    -394.109059041026, -57.5791281561082, -122.111208511343, -61.8333137976925,
                    39.8091735368757, -392.9654948238, -56.8892592967923, -121.903497936512,
                    -61.1684014220712, 39.0630256386966, -393.696585213559, -57.305149251981,
                    -122.467084845282, -60.7377696931311, 39.3447355352536, -394.1728491156,
                    -56.5085763418343, -122.12810039772, -60.7990383908738, 39.3510832722287,
                    -393.291618656019, -56.7594801163953, -121.601073246905, -61.7208011740702,
                    38.9581162855464, -392.196593096523, -56.1934505791543, -121.461937350297,
                    -61.7123627587017, 39.8817188369853, -392.822778257291, -56.4016177581491,
                    -122.028949638867, -61.1379885364757, 39.8935684523209, -393.033585907272,
                    -56.4458818915141, -122.021572692385, -60.7399330258904, 39.5121486550444,
                    -392.573838947424, -56.9689845266772, -122.110553863224, -61.7166453696072,
                    38.9157859923211, -393.323051851429, -56.8229806370227, -122.661983148095,
                    -61.9179364526954, 38.9900573432354, -393.721982653306, -56.8070891877848,
                    -121.994229570258, -61.4671738138149, 39.0239785088586, -393.117087731594,
                    -57.1668184867288, -122.352966121652, -61.686971598534, 39.0026955359712,
                    -394.244815392686, -57.1998288157404, -122.436922328845, -60.6013383401105,
                    38.7249549349476, -394.321072736508, -57.696451986989, -122.687931989123,
                    -60.3667807177889, 39.1218477798834, -393.547235859219, -57.5715491247276,
                    -121.81029924669, -61.6111483718121, 40.3553878752187, -393.735775282438,
                    -57.3329716038394, -122.130955820326, -60.9479895300128, 39.8898575488961,
                    -393.862647342767, -57.1442832695774, -121.957172318022, -60.6920418226131,
                    39.4250985727072, -393.659443721165, -57.266538807691, -122.249389450354,
                    -60.7129747547825, 39.5292289848555, -393.374998316409, -56.9919146570275,
                    -121.507559212113, -61.5046081885277, 38.7198735065272, -392.050760597633,
                    -56.5275250198695, -121.325875269598, -61.4589237145293, 39.3220777479825,
                    -392.413080874074, -56.5353006442936, -121.192233002286, -61.4848880467973,
                    40.1419879074305, -393.778644578486, -56.2025576113574, -121.2338833561,
                    -61.3044496806452, 39.3305906443094, -393.126775671258, -56.6510297838166,
                    -121.984217813725, -61.6063345158352, 38.2565201000864, -393.302364314291,
                    -55.6340480057008, -121.465954697645, -62.8579092584376, 38.6319694570553,
                    -393.381968231907, -56.6008593455163, -121.773169720062, -61.7992008242204,
                    39.4552499179232, -393.174328954797, -56.9971430442111, -122.151746010074,
                    -61.6056974350214, 38.9534518602508, -393.619583031722, -56.9698682721362,
                    -122.707572650356, -60.5941986421416, 38.3646984081079, -394.369746576239,
                    -56.8442663122098, -122.76874513869, -61.0212897437762, 38.9584986367558,
                    -393.882027797454, -56.8927029110782, -121.830794419199, -60.6737101080896,
                    40.0787495678393, -394.933479581675, -56.540822927838, -122.433981082211,
                    -60.5477764713646, 40.0837200819703, -393.578989976183, -56.5426037045587,
                    -122.40874311205, -61.1241425273301, 40.2070994532871, -392.278544299328,
                    -56.8925087402485, -122.347010291143, -59.806259533856, 38.5613765736472,
                    -392.625999391268, 42.6869110681662, 90.58316291876, 45.0674050273316,
                    -30.1350993272486, 288.028949744069, 41.7928656945848, 89.9845522754147,
                    44.9420490820312, -30.2008069565364, 287.883095845726, 41.5661872913255,
                    90.9425623850892, 45.2328804302578, -29.9725526641008, 287.655925976849,
                    41.6405220693725, 90.2854835972499, 45.496522464351, -30.5195993918803,
                    287.612346996875, 41.937975907048, 90.0713471439413, 45.4763535698987,
                    -29.19577078859, 287.381944782851, 41.6559371669728, 90.5111603363973,
                    45.1239834237718, -28.5404597737009, 288.714505429144, 41.4987037402463,
                    90.5563423065157, 44.315911724064, -29.0045007460017, 288.592111449798,
                    41.9449987365894, 90.1550304040984, 44.5294802083905, -28.7572001010647,
                    287.917681820937, 41.5484631157573, 90.1831874448529, 44.2971633095141,
                    -28.5695461161564, 289.50127504039, 42.0779787043469, 90.0152205318336,
                    44.3357933910217, -29.1807981161338, 289.936525691644, 42.3761247809406,
                    89.9974153206779, 43.8327295770738, -29.367105357191, 288.206326699921,
                    42.334613479959, 90.9894227885759, 43.5042901399137, -29.7068968382915,
                    287.88403302752, 41.7274429671185, 90.2278563858655, 42.9194814126206,
                    -29.8782428924949, 288.238658006219, 42.1538441488178, 90.8811074405852,
                    43.4075152592343, -29.9283592458683, 287.209943775801, 42.9563970347626,
                    90.1367508090117, 43.4215141631173, -30.0117378681988, 286.81776205363,
                    43.2387495524162, 90.3474320219313, 42.7605856347322, -31.2969982772634,
                    287.707041484055, 42.5270056713769, 90.0749536044773, 43.2851377584547,
                    -31.1584648411706, 286.237001451083, 41.5441642750934, 89.8779747306645,
                    43.1074653241399, -30.7594615710952, 287.278567415082, 41.8331373437695,
                    91.1859252157712, 44.5935351232045, -30.0724268787103, 287.45084691942,
                    41.2291098914129, 90.4698780697313, 45.3658366421158, -30.2612320574175,
                    287.880762705316, 41.190744461205, 90.7143746263269, 45.9964075257827,
                    -30.4692771811319, 287.567196440248, 41.6894977015703, 90.0810657150687,
                    45.2511378624526, -29.945361147999, 287.307101522487, 42.4847525734318,
                    89.8769459315135, 44.996343616345, -28.9817446874438, 287.456953809532,
                    42.0610819621151, 90.1869633256704, 44.757305313742, -28.1445722611565,
                    288.422880468558, 41.5325876557315, 90.0854908797256, 43.7656451482817,
                    -28.7438875036118, 288.328172624308, 41.8162405398621, 90.7717108458899,
                    44.7731885776612, -28.012250860606, 288.468291174814, 41.8095656543089,
                    90.1002996399998, 44.5282235809819, -28.7441870585059, 289.296051782862,
                    41.5822222704515, 89.9009910311636, 44.1813882146165, -29.3198705253758,
                    289.720014803959, 42.6019922771682, 89.8091324347331, 43.6495704171098,
                    -29.0482249195528, 288.713610030843, 43.076984184627, 90.686324676176,
                    43.3770689229053, -29.5787386577339, 288.312679813031, 42.7370865486458,
                    90.8338628321555, 43.5007721047117, -29.7297554752226, 287.510277665915,
                    42.6358860289146, 90.9284209426551, 44.1342004856647, -29.7994358047071,
                    288.279104189993, 42.5633949509662, 90.9876368069334, 44.0735903035767,
                    -30.3347640529955, 287.479399175351, 43.3725596744773, 90.743871754592,
                    42.2863551483631, -31.0255494540739, 286.875120041006, 42.9108613020652,
                    89.8284362367372, 42.2142416678669, -30.0140023573926, 286.553332682524,
                    42.0954771806726, 89.257527935904, 42.4332180147103, -30.0617827345758,
                    287.250338215467, 42.0075870761784, 90.950978270396, 44.6913246434219,
                    -29.751297445592, 287.787116894964, 41.5478894624649, 90.8941382684017,
                    45.8350185977778, -30.5108826305886, 288.099736274653, 41.5824082618459,
                    90.6108924744173, 45.9341678807548, -30.3128894537809, 287.259578391536,
                    41.5524717971451, 91.0987393157391, 45.7778723608985, -30.2357885432495,
                    288.176703976079, 42.6553034115648, 89.8602598472634, 45.2253871974933,
                    -29.2757772314508, 289.001842916055, 42.0316064684256, 89.5209270203004,
                    45.2156957371548, -29.1032221054239, 289.282391762321, 42.108959773398,
                    89.7997244780725, 45.0152066907802, -28.9494375266453, 288.758129921145,
                    41.2894874881165, 90.2387909826766, 44.3138563086517, -29.0318118302513,
                    289.210087346056, 41.688894834518, 89.775227488604, 44.700929369566,
                    -29.3388510108074, 289.349972056695, 41.8512810721481, 89.8594190239791,
                    44.2617062365797, -29.3883742975012, 290.096466610694, 41.6758548840692,
                    89.5472545995609, 43.8925714685317, -29.6538166801905, 288.652866120835,
                    42.0336182963838, 90.1521651261725, 43.9626759669478, -29.1343197286067,
                    288.424682427542, 41.7976031752217, 90.7457020935597, 43.9479310147375,
                    -29.0805851363237, 288.423776159753, 42.3570347557891, 90.2204083106204,
                    43.7425558348864, -28.575096345803, 288.7976003143, 42.116371599136,
                    90.1629528055779, 44.3478640351903, -30.0437288007722, 287.411560056365,
                    42.0459160958658, 90.6039646410986, 43.1117995171679, -30.2621960257986,
                    287.238998646542, 42.3298735011484, 89.8354472003598, 42.5402316706467,
                    -29.2362758665644, 287.005053933658, 41.7972413557358, 89.7522874263099,
                    42.7044855131996, -28.8687117420583, 285.796273245188, 41.3062328474039,
                    90.2209661891034, 45.4067473572453, -30.3944496976921, 287.35288044655,
                    41.349830822402, 90.7051159996597, 45.1609161874663, -29.9622772657216,
                    287.337424797776, 41.479404111926, 90.1171958655671, 45.6890691510796,
                    -30.3660194923611, 287.986788833804, 41.8189195250607, 90.231678456792,
                    45.4363229262482, -30.4173342502602, 287.553732343295, 42.2650375390544,
                    89.4535127439829, 45.1317679373937, -30.1817371213236, 287.461587419906,
                    42.1472890267024, 89.2392524934397, 45.1562236531164, -29.6646119546914,
                    288.123027240326, 41.6542504497371, 88.9067371994828, 45.5333281842765,
                    -29.1828129690619, 289.107250885929, 40.9154112958168, 89.8381196298561,
                    44.2546467879439, -28.9975542594326, 289.868311881103, 41.6422575423654,
                    89.7702209677708, 44.3973043889579, -29.2950234377187, 289.389221068793,
                    41.7927011316098, 89.8004957375852, 44.0307510729522, -29.6884465356489,
                    289.764412812371, 42.0167759478075, 89.3095623463416, 44.303275264584,
                    -29.6446201146341, 288.549660139216, 41.878048220992, 90.0011817350727,
                    43.8106954984898, -29.6560371615255, 288.944754046299, 41.5101519931634,
                    90.3200686231971, 43.9673947348481, -29.1972316823201, 287.967582032548,
                    42.5484589310731, 89.4493676629512, 44.6278624929027, -29.2535875787227,
                    287.839434941564, 42.5990393265795, 89.8684040990788, 43.9958299688323,
                    -29.7357417199133, 287.26624366134, 41.4074943247152, 90.338258084359,
                    43.1129037296431, -29.2033443142918, 287.593525712787, 42.026650726958,
                    90.4681295282396, 41.7894527470227, -28.9463508108002, 287.254622980842,
                    41.7001306218841, 90.871256358258, 42.1100345786235, -28.983870392886,
                    286.593636823079, 41.4778569485087, 89.760395752175, 44.3378077238078,
                    -30.3914440046783, 287.773730422449, 41.795446111913, 89.8539197264955,
                    44.2143047252033, -29.7754682454013, 287.018593417864, 41.4411629528398,
                    90.2442713795112, 45.0305702191451, -30.4984372017034, 287.875304945465,
                    41.9628325664009, 90.1844235480214, 44.8723734856062, -30.5381347199204,
                    288.031932203405, 42.2757847033834, 90.005797393262, 44.9225092946613,
                    -30.3218711109338, 288.120607786558, 42.0610078854613, 90.1032942402821,
                    45.4399264547417, -30.0367541122235, 288.12080993499, 41.8858912900366,
                    89.0890709599504, 45.1278790592284, -29.049667706796, 288.582404851272,
                    41.536190215258, 89.2393885357034, 44.4750757495797, -29.0362975825341,
                    288.821298782704, 41.4126806142653, 89.2999661894731, 44.5484004695759,
                    -29.229694631936, 289.348395944103, 41.3218977419807, 89.0265145641166,
                    44.1029651345232, -30.0025012945203, 289.31722068223, 41.9465285445387,
                    89.4881102822751, 44.1342564842763, -28.4239843578949, 288.11276815753,
                    41.3474843477648, 89.3323032885972, 43.874794133812, -28.8079306792465,
                    288.815902812759, 41.9873058544453, 89.3979472170362, 43.6490745057765,
                    -29.3352134040117, 288.813131991575, 42.7974811416097, 89.2732421944096,
                    44.1684849818091, -29.3865516222751, 288.381204041387, 42.4674434138308,
                    89.2855218088646, 44.2638067569419, -28.8625463466333, 287.69389343394,
                    41.702012387298, 89.3240410363988, 44.2181603811436, -29.2752507155106,
                    287.490643500323, 42.40595114861, 89.6081780757559, 43.6361040945229,
                    -28.9945245199867, 286.883015276099, 42.208999940567, 90.3491835263634,
                    43.4871332896392, -29.0007298043452, 285.700785731375, 41.3636790639325,
                    90.1469013287232, 44.0140196225546, -28.742379040047, 286.813905648088,
                    41.3617368675218, 89.9734489413528, 44.5270979192902, -29.3795247731303,
                    286.708899754026, 41.3331714807551, 89.9729307603178, 44.3872864919823,
                    -29.6787506916149, 287.66036499363, 41.8498886409379, 89.8410070995768,
                    44.4059474033775, -29.699204852913, 287.461274785352, 42.456967003498,
                    90.5467012707487, 44.4650558160034, -29.6569784110996, 288.589951724108,
                    42.5136745302964, 90.7920288778029, 44.6024033542257, -28.6815025726092,
                    287.70415741999, 41.9722313915003, 90.0284846598047, 44.5103613410617,
                    -29.0404621319228, 288.192072517082, 41.698488019892, 90.3796806800066,
                    43.9506798176796, -29.1310792010627, 289.202438216892, 41.8258285852817,
                    89.8862422315439, 44.0459247091988, -29.4465516521194, 288.792540384409,
                    42.1748508254122, 89.2242245533488, 44.155505369531, -29.051040224186,
                    288.143406236645, 42.4723238403893, 89.4912740558415, 44.7704936381112,
                    -29.4475719343107, 288.645710112372, 42.2166451219455, 89.0080780886007,
                    44.8914397662014, -29.9119793873085, 288.669046363607, 42.4336906734561,
                    89.3248776807745, 44.0997485690207, -29.8948026386145, 288.793135243669,
                    42.727069208256, 89.6493866710537, 44.5009368388971, -29.6400379028949,
                    288.271500351667, 42.2581337312639, 89.8681784749675, 44.1312699336676,
                    -29.1336092547446, 288.581829960817, 42.4513056579585, 90.057967589207,
                    43.8305372064306, -29.9345209385568, 287.588772609847, 42.0214582898028,
                    88.8241499226873, 44.0095398181278, -29.194792663383, 286.608622696772,
                    42.1738037997857, 88.8623612110515, 44.0613508878821, -29.0993362709926,
                    286.621530870735, 41.23782909929, 90.0214829804305, 44.4814617927042,
                    -30.16365705181, 286.961184381246, 41.663463106413, 90.2125868401882,
                    44.4011737670137, -29.2161805435712, 287.594878299538, 41.2802156562023,
                    90.3755939920411, 44.4075918468955, -29.405889113729, 287.237698749501,
                    41.5357705422402, 90.6489199378306, 44.3370659625426, -29.114017981909,
                    287.983905343621, 41.6363651851587, 90.1719381740508, 44.8933467192312,
                    -28.7671529211596, 288.313650539289, 42.4976390431145, 90.2102316226593,
                    44.6784898351164, -28.0685009375703, 287.99719281356, 42.0177519022652,
                    89.7735048550502, 44.4307380273211, -28.769046351318, 287.993997165673,
                    42.0504003560178, 89.5782153300704, 44.6246582452711, -29.431851826486,
                    289.302875398399, 41.9908984849009, 90.4156454348781, 44.2033690383655,
                    -28.9678170970848, 289.214651719633, 42.3415728200823, 89.6315594014803,
                    43.650308411005, -28.1796472744914, 288.864941400699, 41.8482020232544,
                    89.6555086741728, 44.5417821127796, -29.010995059611, 288.802201315901,
                    42.5175863029663, 89.7258197928059, 44.3863001275848, -28.8102238629999,
                    288.632896529666, 42.3425730460434, 89.769295574044, 44.2027395320356,
                    -29.0141765049268, 288.992900450679, 42.4335569882751, 89.7311848109127,
                    44.3595414277268, -29.8748064838857, 289.095171834903, 42.5068695727869,
                    89.8772802193026, 44.0972839769567, -30.0103596020206, 288.257896841483,
                    42.3147025992611, 89.3525523948434, 43.8176088137207, -29.3668885954683,
                    287.729693625354, 41.8887250184322, 89.4997732479086, 44.0038675024248,
                    -30.8377370508497, 288.102926081289, 41.890927504242, 88.8218066618508,
                    43.4007707068364, -29.3947257691747, 288.467170449681, 41.9738549931141,
                    90.2766653490388, 43.9141173502359, -29.4726510716125, 288.297837837122,
                    42.5528871178161, 90.1672872395483, 43.6862451145403, -28.6490358008894,
                    287.935210707525, 42.2322944620449, 90.2405959201595, 44.0815159719863,
                    -28.3564303953604, 287.79976710129, 42.9400989346373, 90.5824927853456,
                    44.1470017835686, -28.8550567422037, 287.793115189436, 41.470480783936,
                    90.099511557619, 44.4644654140323, -28.4956122913193, 287.755352306977,
                    42.4003667154222, 90.4364167357801, 44.0339761985003, -28.911813790255,
                    288.276663590445, 42.4892880979688, 90.2230043411131, 44.2190912597581,
                    -28.9145877420526, 288.619840696421, 41.2098665630043, 89.9010850015038,
                    44.2170949315741, -28.3828888447944, 287.654079324579, 41.8627053152585,
                    90.1155381738845, 44.8077642395258, -28.6167816500214, 288.941538951914,
                    42.5386645972813, 90.4577290433109, 44.2301076709238, -30.2052284519669,
                    289.077980889631, 41.9100425105204, 90.3804438550332, 44.6368949210317,
                    -29.5460549939527, 288.659322027873, 42.5532217951464, 90.5073113927473,
                    44.4858789828079, -28.8898568704504, 288.519282549076, 42.4127466026254,
                    89.959915023884, 44.5386432066272, -29.6837516444878, 288.844781703321,
                    43.1569418193137, 89.7876640280894, 44.5262717035499, -30.2927078706177,
                    288.515557669589, 42.5998451801886, 90.1634521556037, 44.298940303163,
                    -30.3052736472832, 288.167789916887, 42.2376939884267, 89.6246285004219,
                    44.2336818690993, -30.2596924572743, 288.48824625983, 41.8310589380833,
                    89.3037023124449, 43.6513404926772, -30.1009470051415, 287.962236791495,
                    42.0710204190176, 89.4266371405084, 44.1994369515888, -29.9519173154232,
                    286.999759801185, 42.2730431770147, 90.7737413739902, 43.0851725326135,
                    -29.0840553983523, 288.523926721542, 42.3696169116649, 91.1639285189929,
                    44.2625739406606, -28.4559035048756, 288.466640306857, 42.5862728481769,
                    89.7097887862927, 44.1093880678382, -28.3727860370684, 287.124599595688,
                    42.1746514243107, 90.1294520020353, 43.7993056271081, -28.2942130967227,
                    287.551352973477, 42.6429016006929, 90.2774840027551, 43.9250934988794,
                    -28.8387926837745, 287.516104791112, 41.870998817816, 89.402448582405,
                    44.0621708023672, -29.0265957464881, 288.22122005332, 41.4571753825643,
                    90.0597376879263, 44.1135919052601, -28.616054799633, 288.994635308103,
                    41.4264182742101, 89.9019770964702, 44.0153677827236, -29.0596291315573,
                    288.749916529905, 41.7179817918876, 90.5069881923595, 44.4639461082752,
                    -28.9296354609552, 288.840674186639, 42.4950431492657, 90.4699780939071,
                    44.4788485021524, -28.4044554123126, 288.784649645165, 41.8223984002184,
                    89.7529691114326, 44.1053994247224, -29.5158665079971, 288.093697492713,
                    41.9657867410278, 90.8254292443138, 43.7354957372542, -29.4486113582419,
                    289.644558825195, 42.4647113419673, 90.1945389888277, 44.1619801162959,
                    -29.3499518952852, 288.769611028291, 42.6170603133031, 89.7861890009582,
                    44.529606216877, -30.7888461694498, 289.13598130969, 41.9175165682941,
                    89.0261285421649, 44.6305565311519, -30.6597890834041, 288.712826424295,
                    41.0754843396667, 89.1382979676859, 43.852048198159, -29.3446386575586,
                    288.189988311088, 41.6443069634683, 89.422097827821, 43.8066047002317,
                    -30.0879924354191, 287.538699916327, 41.631285089424, 89.1477527883573,
                    44.0882093519417, -29.6673886161764, 287.512577678369, 41.4343243041633,
                    89.4985147807036, 43.3245626139398, -28.9145634872203, 287.834277005141,
                    42.195542282275, 89.792090202834, 44.1160062308471, -29.175960254507,
                    287.406184075133, 43.0842684387435, 90.1720520859443, 44.3774193097685,
                    -28.0938982688121, 286.873650176766, 41.9626676640921, 89.9994894375022,
                    44.1339982089774, -28.6011812528486, 287.622882791503, 41.4826396468584,
                    89.8272858623403, 43.428017336553, -28.3133923586588, 287.342038375981,
                    41.619559984195, 89.3064572438511, 43.7432473062159, -28.702387041116,
                    287.424129599457, 41.1988099135257, 89.1824344259859, 43.7931294626158,
                    -29.1855615220059, 288.985492049856, 41.9052595522963, 89.7810919811429,
                    44.1921189780393, -30.3591794384914, 289.07145924299, 42.5410282653175,
                    89.9295766253688, 43.7800669657625, -28.9386761156273, 288.27523324242,
                    41.9393261420028, 90.0409916744137, 44.0279451911345, -28.3165200297803,
                    289.522108189505, 42.0433264044161, 89.6301439095047, 44.6027065833988,
                    -28.6545746639653, 288.885190827143, 42.5112730634106, 90.1002438531206,
                    43.5690451785343, -29.1438024537613, 288.785732672704, 42.2464774181498,
                    89.8895279648115, 44.1506394784614, -29.1614965368625, 288.776762024257,
                    42.3529503940812, 88.8506797977995, 44.0371523517983, -29.2922062589984,
                    289.126642785251, 41.5345327539685, 89.0002033283145, 44.4666801123085,
                    -30.0654260986705, 288.653120226957, 40.5092450913905, 88.9431860074199,
                    44.0956770971767, -29.7157370785381, 288.153049704858, 41.3044829945027,
                    89.4611540399138, 44.0728005493361, -30.1002777585103, 287.677422374996,
                    41.7235341646407, 89.2533206102162, 43.7711361732452, -29.1342204457967,
                    288.819010441464, 41.9315950680819, 89.7231608594529, 43.7723101185221,
                    -29.4532497762745, 287.990700373785, 41.9463081017415, 89.7194298801879,
                    44.6025997310928, -28.9558380661798, 287.491131933673, 41.8830075612497,
                    89.9451277378808, 44.4677064056706, -28.5717938321355, 287.354387375654,
                    41.5947808510321, 90.0365958960223, 44.8107393555807, -28.7964910836324,
                    287.312395890847, 41.258879429349, 90.0068567440673, 43.5145641860133,
                    -28.9772114506619, 287.573619255375, 41.4841357538527, 89.6889214266825,
                    44.1812032998666, -29.566459820038, 287.33664174261, 41.6304697984286,
                    88.8095729574323, 43.8843106787422, -29.478026157329, 288.990309679508,
                    42.5775991722472, 89.8089329656244, 44.0921185516459, -30.3320483313646,
                    289.092341989429, 42.4070649911624, 89.8194037789846, 43.1835290510567,
                    -29.0010598392375, 289.175188155073, 42.5088840151651, 90.4635806225234,
                    44.1997620286924, -28.4384004773916, 289.089174476977, 41.6954165340047,
                    89.8352460781956, 44.0413043754371, -28.3932648376928, 288.411709556663,
                    42.1951997954109, 89.585987073797, 43.7691181033189, -28.9879167189257,
                    288.530006149873, 42.1954338754268, 90.2736969598696, 44.1036346939845,
                    -29.19148358033, 288.04573057046, 42.0766117822858, 89.9737356576555,
                    44.2165140257445, -29.9696052458805, 288.8128114519, 42.3014889803894,
                    90.1932400232121, 44.25465700673, -29.7906023616179, 288.73267971298,
                    41.1130005841986, 89.1696964053541, 44.1152207701869, -29.8415642528944,
                    288.724226879512, 41.6167562521359, 89.5417685821456, 43.2091065534952,
                    -29.7731017336618, 288.525827759988, 42.1352436611123, 89.6129918437765,
                    44.2570552948075, -28.7744443986328, 289.147081072408, 40.9404393518271,
                    89.4735892748889, 44.4697600006235, -28.5988635525799, 288.546631712342,
                    42.2049345634514, 89.5065160989765, 43.9695318727516, -28.6588525595538,
                    287.520666769271, 41.7986452212667, 89.6098941160021, 44.291288353136,
                    -29.7154912356909, 287.61964460335, 41.8723003208371, 89.7371105720271,
                    44.7341664705098, -28.899920225161, 287.626188032235, 41.7022751282514,
                    89.35190993102, 43.7767213415923, -29.6970315911398, 287.834856102149,
                    41.7969825131045, 88.9755347391008, 44.0548680346983, -29.9701996161593,
                    286.99065898773, 42.0019621502001, 89.3206401693718, 44.3047381743666,
                    -30.0164235870575, 287.816802869026, 42.0641721331988, 90.1982336992665,
                    44.3412804159053, -29.2070961074218, 287.477437538075, 42.4583550496007,
                    90.2087251716765, 43.9287210532202, -29.04654678816, 288.799622499371,
                    41.8588246664774, 90.6424490099416, 44.1905653171006, -28.5988624191791,
                    289.627938041152, 41.6388754872571, 89.9950175089138, 44.5787842968474,
                    -28.919878747747, 289.24838829474, 41.8957435708406, 89.9142803743476,
                    44.1627426200205, -28.5586805573901, 289.025725155731, 41.1990598693954,
                    89.646698731141, 44.5501288610717, -29.1607723949443, 288.295100281249,
                    41.7320621314485, 89.8161416610927, 44.3299677885218, -29.2464960539357,
                    288.697581337632, 42.782119032359, 89.1302161138746, 43.6800963651036,
                    -29.4237084425851, 288.457231133631, 42.7807022800308, 89.2349060344121,
                    44.1451035825546, -29.2256561305626, 288.295501108842, 42.1240612468832,
                    89.7900115655051, 43.9649048670045, -28.608624203271, 289.762691360967,
                    42.6867816906048, 88.8531495338568, 44.2726162315148, -28.9071490418197,
                    289.138496006209, 41.3522749234845, 89.9113692363241, 43.9012232221772,
                    -28.578040435805, 287.410184017444, 42.4501870857642, 89.8844375345012,
                    44.0385638021843, -28.942597266301, 287.881120023979, 41.5500691760005,
                    89.5639925943324, 44.6754248885761, -29.6996331647498, 287.495000023229,
                    41.5507565591614, 89.1875543819334, 44.6739025534385, -29.9252897307405,
                    286.882871295297, 41.6522415150878, 89.1292712959841, 43.8129254097574,
                    -29.940436145585, 287.611422710518, 41.9711668355604, 89.5868236878235,
                    44.2989072551639, -29.5103754441287, 288.001708374418, 41.8741113500521,
                    89.6565380351943, 44.3480897092094, -29.657370047487, 287.543981343895,
                    42.8981972652339, 90.4296493916444, 44.1366813722705, -29.5206725285887,
                    287.65530869052, 42.1685760829162, 90.4017321453686, 43.7884624955155,
                    -28.716982880197, 288.447888924206, 42.4230955437384, 90.0516648827872,
                    44.3783675241302, -28.2888078768129, 288.352567194785, 41.8662892023064,
                    90.210214585859, 44.7038451069404, -29.1560828464024, 288.155639866779,
                    42.2859666660542, 89.8793448629048, 44.4385726212903, -29.7903173944016,
                    289.211350070425, 42.3846361426724, 89.5119927514928, 44.9360115231215,
                    -29.6550018778657, 289.079244512581, 41.819338779316, 89.54619931947,
                    44.272570650131, -29.1250520807352, 289.312149278208, 42.2897295148458,
                    88.8958516984725, 44.6201984861269, -28.2143134750339, 288.843334229069,
                    42.1421159261791, 89.2899675019163, 45.0970253185542, -28.0958417001364,
                    288.864137460264, 42.5982600452352, 89.3410184219107, 44.4906678489478,
                    -28.4415771318092, 288.906974449295, 42.3117924960497, 88.9695925218993,
                    44.1658364435202, -28.2651057296281, 288.084491000595, 42.4138636717119,
                    89.2570751346672, 44.2635338991709, -30.2519066581641, 287.394845794572,
                    42.3382649370844, 89.5058202661424, 43.1740638029739, -30.0255705160021,
                    287.430655449593, 42.069167166654, 89.2756058563344, 43.603736811921,
                    -30.5846200399972, 287.300885458009, 41.7076090389392, 88.9963916178657,
                    44.3536372265047, -29.5403271390869, 286.84390957306, 42.1020284492846,
                    89.3057786293676, 44.2250848653929, -29.0550631447749, 287.554484705475,
                    41.5152546975265, 89.7718593740479, 43.9264436400286, -29.7130927161676,
                    288.106988172915, 41.3719627390738, 89.603817392275, 43.6597045433504,
                    -29.6906803088934, 287.625602393908, 42.5355784529237, 88.8391215880457,
                    44.3049100000036, -29.8655910012976, 288.117387971053, 43.2949087974008,
                    89.5704569609676, 44.3746663791863, -28.9055341026452, 288.834905920022,
                    42.2268269750588, 89.2811651562903, 44.3865186205971, -28.8327404166371,
                    288.933170090272, 42.1275863842134, 88.8640888247894, 44.949604212256,
                    -29.0692183262006, 288.0114997243, 42.5708895017337, 89.4700605950484,
                    44.148174468347, -29.2232326570498, 289.491232005268, 42.044625478467,
                    89.1314616737022, 44.5212326369007, -29.5119228727134, 288.629379872493,
                    42.0508144521447, 89.6967683535964, 44.7297650347231, -28.534079448821,
                    288.584841847216, 42.6083372895126, 89.8627063516328, 45.420094344528,
                    -28.2736791461965, 289.457650616007, 42.1369616110512, 88.9539778372292,
                    45.197884837038, -28.1094806020747, 288.461161751103, 42.0018214750882,
                    89.0312282460964, 44.7523596180981, -28.9173190452482, 289.332560612036,
                    41.3695507313046, 89.8215629324597, 44.2337362034887, -28.2378438396815,
                    288.850152429945, 42.8062026005639, 88.6683614097796, 44.0406595689658,
                    -29.8833918890623, 286.639999905118, 42.6707053817477, 88.4788074131715,
                    43.9716605888075, -29.4204324847475, 287.002296181414, 42.3923932755814,
                    88.9391261741453, 43.3944124480623, -29.9352333920221, 286.573137220636,
                    42.170335400835, 88.4361621004192, 43.9446335637557, -29.6623238592365,
                    286.209731672997, 41.9471534387596, 88.7981274131989, 43.444793440823,
                    -29.8777141644747, 287.356523578243, 41.9447472217713, 89.7127582419452,
                    43.7499615178073, -30.0238265487119, 287.334868123704, 42.4664874660145,
                    89.73613185987, 43.835722391543, -28.971120501405, 287.237314576386,
                    42.3188988674745, 89.6588913512545, 44.3519188976841, -30.0444883304927,
                    287.75677471085, 43.1912897673557, 89.694140719575, 44.3377476320943,
                    -29.1550318634849, 287.911029221114, 42.9996519787319, 89.4069594171454,
                    44.4577360347425, -29.0058504711018, 287.549984855087, 42.5839217820953,
                    89.7168123127182, 44.1828062867809, -29.3306708781334, 288.689362847607,
                    42.8147542205934, 90.3118521998969, 43.9719046139618, -29.1728723970194,
                    289.648954574665, 42.6939842557843, 89.0698226777325, 44.851045844413,
                    -28.5316023880644, 288.628570344571, 42.1003997329611, 89.8043812258105,
                    44.1516568788661, -29.1742758034187, 290.324032221383, 42.5692269115593,
                    90.7601650843855, 44.7500742151473, -30.2317718698497, 289.042964794505,
                    41.9799038976559, 88.7943947613333, 44.7047285466702, -28.2858372353112,
                    288.689615281874, 42.0056564875712, 88.6818596641482, 44.0261134434258,
                    -28.2954853673695, 289.151740273151, 41.271985490717, 89.1252727148614,
                    43.3366265865803, -27.8127272328093, 288.937825933568, 42.2496426953755,
                    88.5628443438997, 44.0216316268108, -29.0482428426523, 286.802249004378,
                    42.9636445440724, 88.8641502989007, 44.1857251108364, -29.2177000146996,
                    287.182161138531, 41.9925768403935, 89.526477480548, 43.9382528157959,
                    -29.9585809132269, 286.408906073709, 42.3188753191642, 89.3554661959351,
                    43.8340292163382, -30.5966120040507, 286.110765104719, 42.1791052152751,
                    89.6416314895136, 43.9159964017566, -29.7590545861628, 287.02497377482,
                    42.1849357604288, 89.9018919381921, 43.9209576476558, -29.9321940369032,
                    288.073428625559, 42.2547261422719, 89.8034893060336, 44.2795994085175,
                    -29.926816016284, 287.127888889439, 42.6599696541103, 89.5494370461459,
                    44.8329144356885, -29.5653436225176, 287.984051079098, 42.6442803038019,
                    90.3889569255547, 44.6916313330266, -29.8262835776579, 287.778058247382,
                    42.6771936293793, 89.3541505445801, 45.3832546577624, -29.1381374873573,
                    287.896628591824, 42.3110416777691, 89.9844389225496, 44.3452091272422,
                    -28.8347421012004, 288.437839957462, 42.7784293305806, 89.5898129001156,
                    44.5853821492529, -29.4084055321829, 288.029269420557, 42.2402945763274,
                    89.1765243600442, 44.7725865232393, -28.4933292820007, 288.902343876505,
                    42.8275393096528, 88.9637913201072, 44.5759755060267, -28.708920557485,
                    289.762637567775, 42.3663380592327, 89.6500810658328, 45.0280666327485,
                    -29.1674262331136, 289.626404376376, 41.5319447502563, 89.379239441623,
                    44.3739405433209, -28.468139330185, 289.254711086597, 40.9201261634033,
                    89.1548619394247, 44.2818388240819, -27.9015647543453, 289.628287332552,
                    41.8727290025952, 89.5626028287539, 43.9191395115073, -27.6357122000307,
                    290.131448237165, 42.1842705371361, 88.7760600574368, 43.6186801691212,
                    -29.8844471583265, 287.348049676064, 42.529988945879, 88.3441135287779,
                    44.0431672068023, -29.2372597627718, 287.68919971302, 42.1640837175571,
                    89.5876148010596, 43.5691364512267, -29.8148481658529, 287.735923926077,
                    42.0622997446792, 89.7865657604702, 44.1407537804874, -30.4782894762455,
                    287.051163178478, 42.2436176402524, 90.3511726843322, 43.9479115294517,
                    -28.7138361427333, 286.207856213463, 42.783588021674, 90.2848566566297,
                    44.5610374001677, -29.1478941081661, 286.627958814678, 42.27283779101,
                    90.6526624295439, 44.3509515524041, -29.140144742216, 287.175614592777,
                    42.2789424289535, 90.3992284318595, 45.0023135696646, -29.0163035012765,
                    288.006848424896, 42.1778957662214, 90.1521001123269, 44.8200872042973,
                    -28.9295910259724, 288.500350080644, 41.9963441036594, 90.1695048840615,
                    44.560922502648, -28.4384056180071, 288.990594880762, 42.8581942552914,
                    89.4199847692768, 44.1349878248502, -27.800697930469, 289.297724292411,
                    42.6881334863859, 89.7464724632007, 44.631842040248, -27.7347582731838,
                    288.546796203899, 42.7385224379141, 89.5705338884016, 44.6863564002427,
                    -27.5462737615351, 289.096826792161, 42.2306562368247, 89.8724412442723,
                    44.329000383186, -27.44467225096, 290.076567459037, 41.7569837854926,
                    90.9155694152958, 45.0528344369342, -27.1372903541449, 290.446595337488,
                    41.4812541120621, 89.8934921314107, 44.7099497272462, -26.7177313503607,
                    291.004628722112, 40.8935205044811, 89.1365101039165, 44.7409953879722,
                    -27.5218225543609, 290.323724318071, 41.6317163107707, 89.8364274209589,
                    44.8071399354407, -27.754891088401, 290.995456193176, 41.772617604047,
                    87.904276203224, 42.5952233943216, -29.0106420633883, 287.43066595906,
                    41.9651726561615, 88.8457907049987, 43.3294229816432, -29.357717428449,
                    288.672768053031, 41.5175654214072, 89.8216366303271, 42.8328973287908,
                    -29.7525144526853, 287.100225748793, 41.5563983121041, 89.7987379652421,
                    42.954867081451, -29.7418471682704, 288.232585797571, 41.6037484465877,
                    90.0497736697814, 44.3540120837664, -29.0363991893347, 287.438342094317,
                    42.5973085544504, 90.8160988103736, 44.7676479902984, -29.4117629043429,
                    286.925222649755, 42.2250659136255, 91.1237417754796, 44.2752409584568,
                    -29.4503347454545, 287.471579322834, 42.5176431000113, 90.3242716922536,
                    45.0582080303674, -28.6959443501187, 286.415116716921, 42.1216469500531,
                    90.1973832588169, 45.5806756878268, -28.4205348398918, 288.417875260666,
                    42.5424720619725, 90.7036867811677, 44.0882883302818, -27.0165461703629,
                    288.828204303, 43.0924233951718, 91.1114665907629, 44.9641981948466,
                    -27.7454681902188, 289.758016705061, 43.2408360089426, 90.9753466602826,
                    44.3472091865875, -27.4527938765916, 289.520549961924, 42.8473410625725,
                    90.9174763274371, 44.6567020575438, -27.9521304650399, 289.513889577836,
                    42.5607870969512, 90.1955873529093, 44.727269001786, -27.7772059991158,
                    289.430071284493, 42.0431863217333, 89.4541638438763, 45.0218174474881,
                    -27.1627081559249, 288.96886704886, 42.3624049219412, 89.2323865969743,
                    44.8466015263249, -26.315195860045, 290.382237134637, 41.5396441785647,
                    89.7816359143923, 44.4356673548052, -25.7921955464756, 289.473683325902,
                    41.3660855804197, 89.7217849436109, 45.4310133877225, -27.9027794328434,
                    289.718681719603), .Dim = c(5L, 520L)))

    expect_equal(
        recover_alpha(out, n_message = 1, constraint = "resolution"),
        structure(c(-0.0590065558990034, 0.898836915156284, 0.232153019528937,
                    -0.0506496962755136, -0.140352558040234, 0.872336720886778, 0.571905351822149,
                    -0.870851914923833, -0.896256906575246, -0.177197895653137, 0.648420266185582,
                    -0.0752972407563952, 0.158086187864285, -1.02947505657346, -0.470486760381959,
                    0.518426227876489, 0.337722879065705, 0.0515255177905942, -0.741869959549348,
                    -1.11732159198215, 0.228925782439831, 0.500631482578026, 0.0176661186615661,
                    -1.06053002132884, -0.976564655830771, -0.222974057093154, 0.976146671578622,
                    0.465025948909556, -1.53408755558652, -0.738707123658486, -0.0753116579732875,
                    -0.0882287473753074, -0.463321579674982, -1.16591931477333, -0.655532238257734,
                    -0.72314503084209, 0.347489758159433, -0.457814101352199, -0.811685873250866,
                    -0.508098057187595, 0.226356511767477, 0.240220936544418, 0.171801490378243,
                    -1.10300541680832, 0.0626565155728258, -0.749468927637281, 0.114169022560805,
                    0.584418388811628, -1.56528956440448, 1.67257696893085, 0.207973045120241,
                    -0.453776302776266, 0.741883668399907, -1.30399430328609, 2.3314492083085,
                    1.33348550507512, -0.559865179087254, 0.346770719847655, -0.871675480701583,
                    2.54807729889237, 0.793062052074163, -0.316826974611811, 0.704918235580863,
                    -0.310456581736105, 2.15805749484218, 0.995227485188167, -0.447518171917565,
                    0.397075132573207, 0.473202291316312, 2.02573162615005, -0.0188462168413821,
                    0.911354859941543, 0.293878636236968, -0.459473879272196, -0.269816052083314,
                    0.698947490881153, 1.01793795371567, -0.387006740034224, -0.725321519732887,
                    -0.050737825034048, 0.446873942815898, 0.110105535207023, 0.310430644770122,
                    -0.987714516698503, -0.402205071091913, 0.142474397100528, 0.109337658250762,
                    -0.253482610839967, -0.681471783990759, -0.965931298961664, -0.173061841459003,
                    0.43506791520873, -0.712430843133063, -0.777364174085051, -0.841008063657171,
                    0.086065338063662, 0.78737502150414, 0.371398789583225, -1.17666289906987,
                    -0.50891160653461, 0.385376451373105, 0.316277099684733, -0.021092244974966,
                    -1.05736806234233, -0.180852918867345, -0.28768066425647, 0.432325262901621,
                    0.106387392795398, -1.01657010446385, -0.770053051456841, -0.191637654216066,
                    0.941420902943587, -0.340011291427118, -0.791947237636009, -0.0262626752104609,
                    -1.01200369351315, -0.0965577692352042, 0.350985526530678, -1.1661028775614,
                    1.42693997391405, 0.550801578714157, -0.411611633077783, 0.58777103108838,
                    -1.38052791002065, 1.87512664464455, 1.28782202313548, -0.566851415453819,
                    0.0400216991466777, -0.816849180099723, 2.07079627737471, 0.832458438118465,
                    -0.446805565920396, 0.729699868255352, -0.159760913318536, 2.10782329992071,
                    1.12248268197297, -1.00216675247878, 0.577306505813397, -0.0220093902523502,
                    2.23488015523597, 0.838642643085137, 0.0822979083890587, -0.111641791834245,
                    -0.941716401181964, 0.860951998693025, 0.439465556667244, 0.353166994801526,
                    -1.10045038965666, -0.810739405350574, 0.570899492360184, -0.339597494719214,
                    -0.761093566336328, -1.05723631254807, -0.777928981817496, -0.662198463549402,
                    0.0281728590075829, 0.0341570875480386, -0.109667511295697, -1.21492945253064,
                    -0.993438583455855, 0.127530401245423, 0.119120779363435, 0.126357676463272,
                    -0.102312975726591, -0.731170782757602, 0.618663738905241, -0.391094470736846,
                    -0.0428041414429003, -0.936790338999572, 0.0634615818709392,
                    1.54452170226318, -0.348593288071534, -0.133236458075089, -0.507942193867734,
                    0.423939939288914, 1.2988157791477, 0.216439995299698, -0.733003317900454,
                    -1.18778678690569, 0.327103504593623, 0.764976675054356, 1.05390305948784,
                    -0.201899115596973, -1.1997430309994, 0.953605955273673, -0.501138504355197,
                    1.16212619508189, 0.0954748159518886, -1.90066906009457, 1.03412093968771,
                    0.347596135015607, -0.118520718602198, 0.393692873650366, -1.82174201344313,
                    1.48069497532583, 0.113032569819779, -0.584809751554914, 0.117258898247599,
                    -1.23179396647515, 1.84120189949323, -0.0667138349516279, -0.958626320837503,
                    0.984503881295666, -0.497393108374581, 2.12246829628768, 1.19653336393719,
                    -1.06294682080147, 1.01633562190182, -0.920654577862905, 2.15277272232154,
                    0.189200813484632, -0.112399982532281, 0.204637077213135, -0.104503676091696,
                    -0.13798715913498, 0.279330524789742, -0.170060152639234, -0.571654416822412,
                    -0.571006031669555, 0.501395783607649, 0.741975687716263, -0.852604210075356,
                    -0.190086358384121, -0.301889860338292, -0.245330475786574, 0.128579996972121,
                    -0.767024095440128, 0.0853733736363438, 0.10394749225398, 0.0475564218675117,
                    -0.0231186130399976, 0.0504481041811111, -0.511969464323158,
                    0.600045954497915, -0.322137642634658, 0.208200069582318, -0.31159284106166,
                    -0.0553374311193693, 0.148703021130103, -0.166113461538828, 1.34361703138458,
                    0.0602085940782686, 0.040978596710886, -0.187395357980847, -0.245834581365273,
                    0.915930506672737, 0.229643596293471, -0.200501427783138, -0.100481924109715,
                    0.154411666974511, 1.09476351863358, 0.534371723450967, 0.0566467686110741,
                    0.172528647422411, 0.309716114082335, 0.707722964498515, 0.0437991483253271,
                    0.333548441510402, -0.750416588537234, 0.797307091317521, 0.385146208485111,
                    0.347396302801883, 0.652483058781002, -1.03965096804623, 0.370733307680212,
                    0.294407064486393, -0.239305123172418, 0.327839299784728, -0.348604501340049,
                    1.35160606726893, -0.129801464606913, -0.424390598004599, 0.401104822922264,
                    -0.781393995148392, 1.55690307363034, 0.18164631913254, 0.15306033141901,
                    0.295523267835051, -0.41320225781152, 1.76583908935623, 0.158072492560194,
                    0.596372906769897, -0.223411318269768, -0.0556783638078002, -0.407548620467452,
                    -0.080036254128629, -0.828733215149896, -0.000908562537460966,
                    0.223213911904651, 0.822088117492996, 0.809200227373793, -0.532495607618724,
                    -0.379491272461962, 0.0207223801281486, -0.493904227105958, 0.860716758427714,
                    0.102064402995829, 0.544385509756353, 0.909926446029772, -0.158687025326088,
                    0.295751629211878, -0.097798934243972, 0.15218623402933, 0.867175285279188,
                    -0.244181377809838, -0.824734839403618, 0.524734258836396, 0.152328323718621,
                    -0.172878029602998, 0.101752455650541, -0.349263398091239, -0.00845257842018299,
                    -0.103367703957588, 0.18344013168948, -0.255660528853127, -0.0148237502844779,
                    0.463497479004914, -0.168764572593318, -0.36197982493951, 0.0897665791316911,
                    0.350616715985069, 0.286996637429066, 0.24741715277446, -0.149109028957497,
                    -1.06907477997169, 0.397437391095423, -0.183716655657989, -0.258717264977008,
                    -0.726016744298221, 0.179116011563536, 0.463625616173942, -0.207833967614409,
                    0.340511069168201, 0.0935555209546237, 0.123633093729467, 0.354937058258855,
                    0.489093127040611, 0.533569387940304, -0.111622233175069, 1.81468484364223,
                    0.493113877162273, 0.234373789759939, -0.0180074450818779, 0.173309187931665,
                    1.18666607477718, -0.0277534926665055, 0.803215859267041, 0.230907661111473,
                    -0.282040398764764, 1.35039926730599, 0.547315461580581, 0.277326900948061,
                    0.309637985589234, -0.561517234306592, 0.471308296639876, 0.646246225078585,
                    -0.230950299480277, 0.344064531770869, 0.200609154480929, 0.709808871382279,
                    1.01007451779526, 0.372063162330381, -0.280379814118191, 0.213818923812426,
                    0.67473270691363, 1.00560072114169, -0.158873139974062, 0.35971851915491,
                    0.195794864302371, 0.24849623225461, 0.565473350833827, -0.908840355485104,
                    0.312314211192813, -0.356110772747911, -0.526786036948863, -0.587844085586255,
                    -0.532274557284055, 0.0468519177923099, 0.760453990846486, -0.976686048245909,
                    -0.0207223234515084, 0.354820452972238, 0.145710115673921, 0.146899840321083,
                    -0.276080537772174, 0.450797919367858, 0.688574602337763, 0.705743353415187,
                    -0.21504417348612, -0.281990180734397, -0.1291199613336, 0.239065685683158,
                    -0.0426037830362702, 0.162258786513945, 0.00897839968126846,
                    0.257351391775188, 0.455821065938114, 0.737994793490486, 0.464231233903831,
                    0.0556713989038542, 0.393827735973503, 0.473620793474009, 0.136874660181206,
                    0.458421234612189, 0.933280653645994, 0.671443871880314, 0.018972013384456,
                    1.01463128710492, 0.463052838759012, 1.26030073619819, -0.196731635080681,
                    0.678652734127095, 0.746906495535768, -0.253482120632867, 1.26517194134917,
                    0.244715728508226, 0.609303222814443, 0.122597283956793, 0.378514647003698,
                    1.30291916979084, 0.702760578330214, 0.138537135488356, -0.471022460214215,
                    0.488955197895251, -0.327294152797208, 0.754241641424073, -0.929774949227831,
                    -0.233721301175592, 0.473126914043547, -0.163102606963463, 0.250225470334527,
                    -0.457889295128993, 0.0838180266886752, 0.0387319312913235, 0.688431658423298,
                    -0.506606004345628, 0.510737111381641, 0.295505945352829, 0.0439302277976878,
                    0.681402199881973, -0.69096517754204, -0.453603824473639, 0.570574763539526,
                    -0.172975179940096, 0.612736985745215, 0.16118379790646, -0.270327631771263,
                    0.220327447880493, -0.0291526742868768, -0.235583121356569, 0.301798033810599,
                    -0.471239401571626, -0.936897156317002, -0.199299619277269, 0.338799167851533,
                    0.29899358652176, 0.388845227749769, 0.147252574834056, -0.0477028686647998,
                    -0.584633731889681, -0.320442237766684, 0.908073283989395, 0.142685265838111,
                    -0.0270987561393383, -0.251827578049898, 0.360121264431967, 0.582610582610222,
                    -0.0587101051229411, 0.630613395582458, -0.156479364857148, 0.514832346903432,
                    0.129280283130555, 0.700280051582865, -0.245873289122414, 0.287381370229724,
                    1.00839723592739, 0.0365521535862676, 0.718022061879424, 0.317724396738271,
                    0.343377941819568, 0.575515627537868, -0.146331782346056, 1.07588945045876,
                    -0.0692373594362437, 0.0759049069018829, 0.619301803272037, -0.304296693839866,
                    1.05762514407894, 0.322701139249403, 0.407955835929755, 0.11080377679616,
                    0.376201703791729, -0.515975723469239, 0.261213980063587, -0.615777744073313,
                    0.632769409302227, -0.747789015857947, -0.0314600447217117, 0.0389534557471336,
                    -0.660193214569631, 0.286356939290101, -1.10536248676261, -0.169432205613873,
                    0.126880554977788, 0.366680262113334, -0.570321428178147, -0.547617782142133,
                    -0.29986214202097, 0.412568609103857, 0.0621014534339395, -0.558575447709671,
                    -0.64976371830117, 0.828797723338084, -0.0610957304884323, 0.0801627879427542,
                    -0.406876807760199, -0.0796226732225414, -0.252143303682772,
                    0.492153883999947, 0.47025510670619, -0.643196516560806, -0.134874675251602,
                    0.0726408564303824, -0.309817830647376, 0.182139527685905, 0.0311702617597405,
                    0.177875617759963, 0.358160272139287, -0.340798522536545, 0.368226638395242,
                    -0.366247880446991, 0.340740018595696, -0.250950770248892, 0.385117740390712,
                    -0.248742640803528, 0.0438568727432482, 0.885039534495061, 0.0613077374054001,
                    0.0535053210313947, -0.236442087138983, 0.071949198906168, 0.240605578044363,
                    0.209038075001558, -0.517575283781156, 0.583016497220058, 1.20584047247814,
                    1.01507267616671, 0.743422740307039, 0.224254696812267, -0.0297795333277691,
                    0.328084300150948, 0.993526338743383, 1.37159010113282, -0.654305967939024,
                    0.641539486072105, 0.957786704855067, 0.522262559440886, 1.49753845629829,
                    -0.176725093854088, 0.708323499347614, -0.508337845333699, -0.504516813511898,
                    0.118872223353719, 0.188293858408883, -0.615216640238657, 0.782374000417491,
                    -0.985784303658924, 0.03397652438656, -0.356973058431947, -0.32355181912871,
                    0.445915716249857, -0.822003321241091, 0.0482919559487414, 0.448332353829045,
                    -0.63952359453134, 0.0372643263180805, -0.133413688120019, 0.22173598506987,
                    0.393727740684994, -0.245068093270902, 0.223607150612594, -0.0784096515949919,
                    0.983664846834415, 0.9476082821663, 0.335346458004807, -0.0235916529553037,
                    0.332255082859234, 0.234114213648127, -0.0396278705527493, 1.21344856058582,
                    0.190632987170176, -0.176590065855834, -0.277500426782353, -0.107918200140119,
                    0.126488539275016, 0.560938961715607, -0.410189761900114, -0.0636462911345461,
                    -0.341854958978416, 0.537936440248274, -0.000492135159134932,
                    -0.136240987997198, -0.197680459836874, 0.19340879550461, 0.368333580515838,
                    0.393435600468138, -0.0830473911065255, 0.487488140756682, -0.0897573881029814,
                    0.373968373103651, 0.742970190975171, 0.856668472917278, 0.635698763863729,
                    0.260977059039931, 0.155580277297162, 0.0332248557091361, 1.10268309239312,
                    0.12299372361133, 0.259498173111638, -0.635656407382498, 0.829118285704794,
                    0.687813743705924, 0.762701974211211, -0.260171183716153, 0.275697885213702,
                    0.160401423128519, 0.611347407366637, 0.475628139160619, -0.25666975142326,
                    0.307032903910056, -0.109806227418346, -1.00378179406644, -0.109676937331926,
                    -1.26111186141065, 0.401481496992943, -0.245729820300937, -0.982772530378799,
                    -0.161440855347223, -0.868300415220915, 0.850408589950689, 0.227005295901364,
                    -0.681421991544749, -0.171460929044088, -0.191854949581426, -0.112296614450713,
                    0.422310082471341, -0.684591486855723, 0.58626990958993, 0.501725795386562,
                    -0.445367545091386, -0.0701819487823627, 0.433593610959356, 1.1741613695674,
                    0.0807044900912786, -0.203834595720451, 0.116305741056181, -0.884453554854982,
                    0.201834463943271, 0.199423003667249, 0.528561140325934, 0.198617831115342,
                    0.250079600566778, 0.22178527392731, 0.102096687292544, 0.121981819596328,
                    0.349655668656482, -0.253623456176779, 0.276182019227292, 0.125973440841165,
                    0.232865784719706, 0.0966615141834382, 0.789255195265486, 0.359520693288118,
                    0.405267540728715, 0.185625354588979, 0.931232402439605, 0.780971538316415,
                    -0.445974438837766, -0.263572645822052, -0.334912570562011, 0.417726412232327,
                    1.25683667680098, 0.434729989754885, 0.244833181774283, -0.0216218536347412,
                    0.453843922950398, 1.31025648725398, -0.0917161944253735, 0.126233573454442,
                    0.363201070781315, 0.932279212337818, 1.18219150265748, 0.195385774403348,
                    -0.90191929291214, 0.286629901149894, 0.516422068084751, 1.67120843773298,
                    0.146203239816316, -0.428712566302863, 0.0170597765070681, 0.896985871544416,
                    -0.800360836778342, -0.421177679696541, -1.47694306725415, 0.129702226246224,
                    0.255074164930164, -0.625441712599475, -0.391460680636499, -0.513937418238214,
                    0.65239788792465, 0.353006394175168, -0.481397903362961, 0.198050300906459,
                    0.15211096263107, 0.0249642057460733, 0.618263571226578, -0.833148995195796,
                    0.214158236562753, 0.0181907157086698, -0.601216428920566, 0.74367355621456,
                    0.409080906401915, 0.629386305291121, -0.0672811514000102, -0.617807394677413,
                    1.08089915362561, -0.289810129508805, 0.422278620752422, -0.0807077602692061,
                    0.954706276869985, 0.558961370443768, -0.583145618936776, -0.0701591616894461,
                    -0.600024451875242, -0.394889105872835, 0.0954283097985638, -0.31816641321382,
                    0.248597274877113, -0.126505390953326, -0.700768429700048, 0.493865013521344,
                    1.20046910725638, 0.173520960047313, 0.385290878352492, 0.267081516883991,
                    0.0217049005108834, 0.957972135721906, 0.523957590711802, 0.47178782089091,
                    0.302491083763641, 1.31513936903517, 1.40299117217485, 0.155701836346054,
                    0.67036521635103, -0.376211295817541, 0.857827010456845, 1.77951073910374,
                    0.129737843936368, 0.237174997052904, -0.363487418593053, 0.683396369266546,
                    0.866376873131816, 1.20302112427146, -0.240317136719426, -0.284003975645874,
                    0.921798177250395, 1.22527468571881, 0.813658676507259, 0.55690211653507,
                    -0.552197242515945, 0.882125739154219, -1.87771209244204, -0.967407163348014,
                    -1.02526804318606, 0.127508630240186, 0.593103961199461, -1.37899058785909,
                    -0.249040081152145, -0.69561889185195, 0.264280928403906, -0.0445759546704778,
                    -0.506482786714162, 0.330309759928525, -0.185122243306637, -0.0988317192489774,
                    0.113145655334137, 0.122687250828807, 0.20759670890385, 0.342278185590601,
                    -0.127744468784726, -0.0519659706932316, 0.437075553106563, -0.106020435587141,
                    0.0249804470320214, -0.549619643261508, 0.5552445827459, 0.108871495624268,
                    -0.340351724976756, -0.399883664022752, -0.0265463457761257,
                    0.589617164306304, 0.296979190631916, -0.0256935353596646, -1.22840976396219,
                    -0.588661882890477, -0.780444249276343, 0.0522746463700656, 0.475692050515761,
                    0.0669627319292374, 0.0501218080098624, 0.613932644535737, 0.851713218736627,
                    0.34525977569723, -0.740705122020103, 0.226889870602207, 0.519020774063137,
                    0.923504493795608, -0.329048650185769, -0.18868867482206, -0.614668863251723,
                    0.377495997378333, 1.38766388102386, 0.665644896113065, 0.257550085965704,
                    -1.08127040426174, 0.287216103852273, 2.33815998437228, 0.705259436544509,
                    -0.124515157917074, -0.723271472102397, 0.839428938603731, 1.69804858749176,
                    1.18300703809263, -0.300412385211217, -0.333457026370269, 1.21287977111453,
                    0.570328560831722, 1.23005857817293, 0.0116078011544403, -0.225794997042669,
                    0.455954091186101, -2.07141269866175, -0.61093145570095, -0.712930759244642,
                    -1.47787740867366, -0.165684371012759, -1.02912049438675, -0.115391120631728,
                    -0.489873792424618, -0.0769973966531268, -0.664684117250488,
                    -0.777163153198693, 0.0365471158839341, -0.150797030923599, -0.367650001909197,
                    -0.0547979289572709, -1.17262880072038, 0.0492625490392129, 0.593786697224857,
                    0.157425997427168, 0.112490935818244, -0.129402086967723, 0.199295092353793,
                    0.190749055051811, 0.473486733200474, 0.380485914025542, 0.150838332076368,
                    0.427319807660638, 0.332713487792148, -0.384277025956152, 0.0613275554124471,
                    -0.84608605263573, -0.572145781639108, -0.39201768489383, -0.471640966896331,
                    -0.783466566720392, -0.0108094303702906, -0.124000629051523,
                    -0.0776387862607066, -0.68178479330146, 0.0240001241009793, -0.096541990143546,
                    0.196572641831779, -0.692639620229983, -0.944488067248614, 0.100056263291265,
                    0.435146912313531, 0.543043007283018, -0.346386106224742, -0.719819051371019,
                    -0.263671520339216, 0.468194373709043, 0.86794394306628, -0.359140887324656,
                    -0.785711711929451, 0.791863291249655, 1.48001575083962, 0.654131151826732,
                    -0.0984442929935199, -0.703661043272206, 0.349597828972264, 1.16253167086142,
                    0.012832698871307, 0.110350589392425, -0.218680549386136, 0.129922414249023,
                    1.34356502118546, 0.00328235674311372, -0.199936305173992, -0.350047993496617,
                    0.107509994999617, -2.70386985148446, -0.436628779246121, -1.35158938933986,
                    -0.853531737761813, 1.75905757155391, -0.958392741810723, -0.477713946171264,
                    -0.591505388725238, -0.237109832353923, 0.930623685167859, -1.66250875435682,
                    0.399865967702567, 0.0147754304749981, -0.312384632926012, 1.48864655530144,
                    -1.63123323630786, 0.408201873522955, 0.546164916632051, -0.402242659096117,
                    0.247548992151877, -0.685173809067706, -0.266479333096242, 0.49968922764026,
                    0.127897660421354, -0.139720452955203, -0.242355919484766, -0.585604515965827,
                    0.0674544715085972, -0.0515527241218763, 0.275759444100856, -0.975469771611479,
                    -0.888443894705688, -0.138473623843623, -1.21526767865353, 0.170055092337378,
                    -0.66519505605433, -0.0690064748834516, -0.110957202852688, -1.20710764610857,
                    -0.22527760333486, -0.526758914000627, 0.252197515245939, -0.70642316168275,
                    -0.571632883496392, 0.448896734620007, 0.346910807733806, -0.3919935889013,
                    -1.03121750766974, 0.185089100385142, 0.122530165604275, -0.189318949606943,
                    0.17702547174099, -0.944588201135815, 0.112643202567511, 0.187920643870996,
                    0.39236179705231, 0.290811055595668, -0.153437114414373, -0.983211360196236,
                    0.270803327763417, 1.228627366782, -0.384384284139514, -0.127210149983256,
                    7.30973495421949e-05, -0.59754423413969, 1.86874084916337, -0.812542180745481,
                    0.260168779844832, 0.194550544105709, -0.324143211169229, 1.71302036153677,
                    0.354945056102139, -1.82686180460232, -2.46760527400639, 0.279572322519556,
                    1.08100491524661, 0.286192645523386, -1.94035763553853, -0.655444338969943,
                    1.50928713890579, 0.733455075268509, -0.128821465884812, -0.970405126583785,
                    0.182554399412311, 1.68174312849209, -0.042887175807266, 0.60343691015953,
                    -0.447975966066309, 0.135170892050098, -0.434226327242953, -1.16379859698032,
                    0.204089002307683, -0.666611393675339, -0.251033401394778, 0.282134107410883,
                    -0.290927776249113, 0.522873614654969, -0.530172645444395, 0.304336388629991,
                    0.35436159883605, -0.579237535393197, 0.0804220822790853, -1.13795776210765,
                    -0.559744996153427, 1.04067590140824, -0.777135570916695, 0.231019978425849,
                    -0.0120476541888195, -0.190087702950507, -0.393322974287685,
                    -0.584815242178234, 0.587775575076961, -0.546239923698323, -1.07201771844223,
                    -0.764578159891158, -0.76344800174234, 1.29575303598989, 0.261007174884298,
                    -0.481903551098782, -0.336443889815058, -1.04601914048292, 1.46654646763653,
                    -0.151120730851602, -0.208727288251453, -0.855907295060859, -1.31662950205729,
                    1.42306829692177, -0.0614386449729523, -0.151133894630846, -0.512632448544423,
                    -1.09151066860792, 1.63337121468075, 0.836599095847703, -0.915440611615594,
                    -0.352609836771222, -0.350682354276518, 0.212321727370835, 0.849392853748469,
                    -0.174179585314704, -0.480249999448091, -0.800930150754514, -0.222283830715128,
                    0.620973873680555, 0.667451737982759, -0.445495852361717, -1.63131317497911,
                    -0.0285727127136113, 0.0283126272249206, 0.3616542928699, -1.20164868085706,
                    -0.824254745144248, -0.368004960153939, -0.834948733244602, 0.238881743796242,
                    -1.14629320595299, -1.40862355372744, -0.76423362383278, -1.30705639519462,
                    0.531399484636239, 0.203283274295757, 1.5769234797285, -0.0305431882455593,
                    -1.85602377808883, -2.0232090596709, 0.732397220826869, 1.67183003620529,
                    0.00922163147316724, -1.8100440461026, -1.15030655382265, 1.36516086234582,
                    0.911317896130473, 0.391302027475632, -0.981024032793414, -1.45301674771997,
                    1.35288764867207, 0.267201583626438, 0.0523247417288815, -0.647029206574757,
                    -0.705027376285017, 0.420776471037868, -0.269227637493831, 0.583111291403725,
                    -0.636586266862281, -0.788286135895959, 0.453508985312737, -0.176382395076445,
                    0.865557659974485, -1.16508074963963, -0.609189222701303, 0.603155435697403,
                    -0.382552921444123, 0.322110715763586, -1.12613101147116, -1.2676790175363,
                    0.675266313626906, -0.714101192805078, 0.760176394766205, -1.01536206026807,
                    -0.719748739128676, -0.835499926157098, -1.17491782219866, 0.658283943345225,
                    -0.10662343325086, -0.5997117600146, -0.740646469212322, -0.528332500336649,
                    1.6127590740868, 0.141543034707723, -0.165968236596342, -0.789608306856564,
                    -0.961301517925612, 1.15621025153845, 0.197519847489474, 0.274090774078616,
                    -1.20842251033233, -1.75818676849278, 0.728118933415601, -0.220913214416555,
                    -0.319379314456512, -0.0686656078950989, -0.898189048766618,
                    0.348430552293962, 0.160138601084867, -0.67624965497356, -0.0267537238830187,
                    -0.0443642569906366, 0.0285475449043702, 0.389056970924742, 0.0880502994647259,
                    -0.513323904419678, -0.555734440611609, -0.0358804526889287,
                    0.379967322573947, 0.961097892309311, -1.0022649017464, -1.59832540160151,
                    -0.203414775228936, 0.0143802689773338, 0.432346110105158, -1.00986649237876,
                    -1.15166417805314, -0.772370487232848, -0.537059538719973, -0.0860910117795299,
                    -1.00006109802948, -1.02174406044254, -0.67336082512827, -1.59184216225194,
                    0.207788337556337, 0.0373081144381899, 0.885902236251432, 0.21414642389513,
                    -1.62423672725691, -0.956990050631759, 0.196130273768631, 1.67764740609122,
                    1.27585659461568, -1.2353832162641, -1.39654866750271, 0.517032795325974,
                    1.25564099477702, 0.411873607161027, -0.588581707475441, -1.0058079886216,
                    0.618476662606895, 0.113718203148796, 0.438502163988261, -0.44348425705337,
                    -0.384941660052789, 0.115728961776256, 0.380637557125375, 0.50512036839865,
                    0.0693749429767365, 0.228682773742605, 0.534273454538322, 0.346760012230966,
                    0.72027630433422, -0.671378207346947, -0.751585568126529, -0.0500386635727921,
                    -0.981788812091004, 0.248823473810546, 0.0943443795501082, 0.116822690468581,
                    0.386812163198755, -0.679726796192138, 0.618962166609322, 0.848071344254549,
                    -0.949536880167756, -0.667136417447551, -0.81400111123472, 0.144998855293721,
                    -0.891690904751357, -0.232534673944969, -0.736584475135004, -1.00149838081569,
                    0.363481116177326, 0.703093604181674, -0.848216013820974, -0.519667674505286,
                    -0.848091537167679, 0.767287175319296, -1.41128095372298, -0.62137930086368,
                    -0.336375790963302, -1.21289195532296, 1.10218541008771, 0.0674130873718468,
                    -0.736173062698157, 0.113791266079119, -1.3785999129358, 0.421158820042024,
                    -0.0425648558280898, 0.38228545980968, 0.653207018693564, -0.557870207162669,
                    0.945232587250501, -1.09019665299846, 0.90941135430964, -0.600555388627413,
                    -0.296792607633517, -0.0832981700471294, -0.294826730546959,
                    0.835610238980308, -0.317512158542769, -0.339998346481373, 1.22615619844223,
                    -1.22398089164577, 1.33417141977085, -0.868804467573909, -0.269084633387166,
                    -0.728792802170318, -0.771062485102789, 1.05323559396339, -0.303173225535602,
                    -0.142863829536214, -0.718132904700241, -0.676921948294961, 0.852892932727201,
                    -0.469256985334908, 1.27896725496543, 0.630695319605621, -0.659662452570814,
                    -1.72053382706309, -0.622759463969345, 1.24659932199965, 0.338737344866303,
                    -1.32115300137056, -1.55044549695279, -0.334062686879236, 0.781971342539745,
                    0.219149626929834, -0.604111814485918, 0.135093458536858, -0.299388505011962,
                    -0.499428547786515, 0.133652111245169, 0.224577381324906, -0.68123218826608,
                    -0.0713666751921664, 0.0368482183355248, 0.703268154823235, 0.0438425344865152,
                    -1.06160529320104, 0.805705341596962, -0.44584988798043, 0.803872646595437,
                    -0.885594199871008, -0.250895062358097, 1.12578432653157, 0.0143872571001609,
                    0.426312002880394, -0.112159287645056, -1.02499216371896, 0.638460381259392,
                    -0.433343128535093, 0.350964702413791, 0.862945360432491, -0.865222323559351,
                    -0.448852168619624, -0.498512627255842, 0.756709642787399, 0.367239936452108,
                    -0.492672652143824, -0.724461730256337, -0.947694213192932, 0.0459132790950889,
                    -0.501922070332029, 0.449864632885976, -0.615603515903757, -0.667526681046152,
                    0.685476482269507, -0.501191327238814, 0.193427275561277, -0.949668260298154,
                    -0.701006158595277, 0.398880796666049, -0.979280108107673, -0.778396164013998,
                    -0.804413092460948, -1.26475945186314, -0.187314340223054, -0.219251026564478,
                    -0.0869273981964227, -0.33097037686187, -1.58363865631625, 0.295831957145452,
                    -0.174883554058965, 1.15448829004879, -0.980229916294689, -1.13397855983818,
                    0.426305086721214, -0.238035701353766, 0.548001415420771, -0.840572288068358,
                    -1.23458211922333, 0.286403282561594, -0.462339532767345, 0.755326720089471,
                    -0.000184340434913111, -0.2480760518397, -0.358095851977566,
                    -0.756743879360073, -0.147532002321299, -0.569010517875682, -0.331127657948638,
                    0.194652145774029, -1.22039789545286, 0.763416954101018, -0.246937904548986,
                    0.495756854270375, 0.0600713814503138, 0.52714739937926, -0.692358988208241,
                    0.291185198228646, 0.69799563759539, 0.110480553054231, 0.031319690375966,
                    -1.50044839047209, -0.112194792806974, 0.546613008073603, -0.361542326337286,
                    0.112309211351402, 0.888078071361747, -0.0995907899765598, -0.730427783930047,
                    0.348932654092692, -0.349094286490526, -0.149516992399526, 0.201256768412122,
                    -0.339820063529249, 0.606241746190591, 0.2660710179467, -0.122556023605682,
                    0.163431590721473, 0.403664113939495, 0.213683754656671, -0.320768108549515,
                    0.414507088611288, 0.287336911204704, -0.298757995797587, 0.691994832439292,
                    0.43707435944728, 0.363046390644314, -0.086183102099568, -0.415781019045546,
                    1.46003453230793, 1.3212459220793, 0.790440313522822, -0.754496493892873,
                    -0.487422204777239, 0.915899718900647, -0.200894471064174, 0.425866506675362,
                    -0.680747099431812, -0.297036705933341, 0.980128864781221, -0.98693475419671,
                    0.628924131074768, -0.0669331537753237, -0.340243694942888, 0.29758887564202,
                    -0.658834627980752, 1.21203013070556, -0.115161764741856, -1.0073663744841,
                    0.311587871006651, 0.254953867473148, 0.220138587513873, 0.284768811030915,
                    -0.849234636101368, -0.0427418652424763, -0.867594086339352,
                    0.27341761140724, 0.389602810871821, -1.38300214875306, 0.0935549986323565,
                    -0.382715865678023, -0.220244782711632, -0.906452625152866, -0.537391403106227,
                    0.68512479474623, -0.946737420918225, 0.559985575502054, -0.645092955083669,
                    -0.135728401514786, 1.39353665601652, -0.50695046263391, 0.452384870385828,
                    -0.18025534441685, 0.0807431131163732, 0.475041377988418, -0.55993822965138,
                    0.0281461424854683, 0.015909262389556, -0.0596216671417409, 0.199411799492424,
                    -2.17950065405697, -0.550868598592274, 0.099402043335175, 0.101480194492858,
                    -0.273903540165747, -1.00653642515056, -0.260571594686354, 0.184503980981702,
                    0.639041083662011, -0.093165110662909, 0.0932040830138305, -0.0666044387647275,
                    -0.172244013149054, 0.356050426608768, -0.218459101281582, 0.594108059979504,
                    0.538306730490604, 0.12187385428345, 0.223892913054911, -0.118908863718076,
                    0.874493815968236, 0.0300366207640934, 0.346040620355012, 0.283371834138297,
                    -0.585083007783595, 0.242169909050773, -0.360573392036827, 0.0383052477850754,
                    0.129865886998314, 0.10616660118859, -0.929512730225355, 0.409184118302591,
                    -0.00242568964151246, 0.10011092239839, 0.286843549056385, 0.493781578620599,
                    0.363655477265553, -0.1443014830827, -0.186252798946157, 0.585151317892581,
                    -0.174698081931922, 1.31516972712311, -0.647578232352146, -0.947713331644714,
                    0.927645240435112, -0.591652218770577, 0.293073786052304, -0.168082034810396,
                    -0.745487386035638, 0.827846676778677, -0.488903223167483, 0.946462635009937,
                    -0.732037846823033, -0.355155146744423, 0.140545103542024, 0.508003649265675,
                    0.603269507566552, -0.174596670037602, -0.526235626653101, 0.35365652152683,
                    -0.570998064945826, 1.62506483650571, -0.159822527037846, -0.421492488574849,
                    0.34664455002158, -0.790624601598012, 0.42388252019478, -0.102086939692754,
                    -0.825338264164206, -0.282422686726363, -1.09159311257898, 0.583648419853688,
                    0.806997181038298, -0.813362375087475, 0.566660152761926, -0.462571363535289,
                    0.476072408742425, -0.0847711158563129, -0.196584231433114, -0.187174345247342,
                    -0.782047579642295, 0.996070147340902, -0.650534385050832, -0.0596023552428484,
                    -0.398335198653143, -1.13558530284141, 0.72395905321909, -0.0918723592271533,
                    -0.0441824385538467, 0.120121381930346, -2.12627222656539, -0.0721005099983927,
                    0.568321368362774, 0.24590647120251, -0.340842061080593, 0.705764981169295,
                    -1.22893133822288, 0.684835786242559, 0.734027764071556, 0.209900399970735,
                    0.906256894145145, -0.386991841811863, 0.144598762125753, 0.259496415478399,
                    -0.0389875540051889, 0.332130077863127, -0.27915050304594, 1.01509756215773,
                    0.907652792995304, -0.246755652571068, 0.790793044868707, -0.565482614343836,
                    0.843256694608158, 0.730528721441644, -0.0046275933096922, 0.6756356965181,
                    -0.0157855282755008, 0.527452949418574, -0.0360016640395955,
                    -0.180636207952659, 0.303415080792519, -0.477767441293054, 0.263917197055093,
                    -0.211047472736769, 0.0203523569849722, 0.155544558029476, 0.0946861929454315,
                    0.277410984770484, 0.173846584871889, 0.278899017587094, -0.241415125608029,
                    1.17073680073295, -0.0045594777896838, -0.204603134905778, 0.400381719026967,
                    -0.545379956868032, 0.692630323276973, -0.248488532621622, -0.185252164609551,
                    0.263178647527752, -1.06918314682255, 0.760851171280819, -0.224010315781868,
                    -0.146229057721413, 0.163739228235556, 0.118636053678955, 0.216710812717963,
                    -0.161627577593606, -0.389517364020332, -0.0679388382981756,
                    -0.181752551605086, 0.498360312596674, -0.88725962901816, -1.11589370220335,
                    -0.409890629193001, 0.227112119535747, 0.585592439678322, -0.579528992245798,
                    -0.890397174349367, 0.593567018140646, -0.225447520285002, 0.311190960411267,
                    0.37522107588245, -0.527310153842368, -0.21157250917841, -0.784427926836599,
                    -0.146820318014704, 0.689137170763445, 0.0966823997113373, -0.748937210645323,
                    -0.840643725637563, 0.839028697679282, -0.130136368358421, 0.0808074016769567,
                    -0.677348157991503, -1.05054476897689, 0.776399538085485, -0.0708457196307393,
                    -0.287845137006638, -1.10034582804025, -0.589818971174708, -0.0309744439142037,
                    0.674683155434394, 0.658084489558235, 0.372799048072977, -0.438541978484125,
                    0.481013672190329, 1.34529562540416, 0.767834485287707, 0.240868809809456,
                    0.511492091538798, -0.0322489769291394, 0.753300869315865, 0.153188593382154,
                    -0.362709595175431, 0.375734308869255, 0.101131791890964, 1.31547845973188,
                    0.634402130050276, -0.414015766476751, 0.391112355411451, 0.450581611500496,
                    1.14177838054198, 0.287637791037607, -0.207690544876755, -0.166543762001453,
                    0.546916421779031, 0.365669583209637, -0.00223658940185345, -0.463495703009471,
                    0.119929555025919, 1.48274723694502, -0.543090114262867, 0.663014462204625,
                    -0.0253147755138627, -0.0596883349106747, -0.108581030967088,
                    -0.0596551225362134, 0.354841888490636, 0.435208765374028, -0.264031910761162,
                    -1.215305465775, -0.0325023835240472, -0.455020801978034, -0.158592498317034,
                    0.267704278651848, -0.333053706581353, 0.241198665008397, -0.507992107118714,
                    0.446398319988447, -1.13164348145685, 0.661043750389751, 0.216616752848324,
                    0.470250270310569, 0.241753891017453, -0.261935118600405, 0.526414769983717,
                    0.0743746540031793, -0.0452155211548764, 0.263610979799154, 0.663329762350145,
                    0.358358128796965, -0.685151041874747, -0.821356744318695, 0.000951601673733649,
                    0.436812211255699, 0.0391288860365648, -0.0900636513161146, -0.661706478816853,
                    0.607703873605089, -0.83527683277228, 0.457639687991161, 0.0129329076031297,
                    -0.619111701751123, -0.560136766613283, -0.279539796638915, 0.895985864268141,
                    0.255781467945354, -0.214502981976029, -0.295708272254075, -0.432141120519816,
                    0.491229157099724, -0.124869298016563, -0.193054381353612, -0.754877248759698,
                    -0.405042635407369, 0.588325818542359, -0.367498100156325, 0.246408400582766,
                    -0.982497093386286, -0.65750043186091, -0.15411945805684, 0.926169592083568,
                    1.04557176578324, -0.130842901608588, -0.675323994307604, -0.753709776604268,
                    1.04913504728648, 1.25736132155535, 0.535981287685445, 0.404551083103826,
                    0.0608820843466162, 0.705778497275645, 0.644726028601308, -0.0260276331629825,
                    0.308194823381058, -0.0799386189866595, 1.24943188830426, 1.19520828086743,
                    0.143096030550282, -0.725779669642378, -0.43869668545284, 1.04225345518253,
                    0.985662773768937, -0.542500119369834, 0.125621505634797, 0.232011837753049,
                    0.435056687050206, 0.748000864245419, -0.423207278623011, -0.159367723919416,
                    0.497847011279532, 0.270831853713673, 0.409206073493731, -0.281383949684937,
                    -0.640930094784295, -0.963596299586584, 0.459008787174326, -0.302305883635,
                    0.206377769727226, -0.73786886015381, -0.530154340726796, 0.226001983900005,
                    0.131082433643126, 0.280581160785403, 0.53375012498455, -0.15412726593263,
                    -0.146156236580168, 0.212381612828914, -0.0699119868407934, 0.0613158207576703,
                    -0.0150186332149929, -0.813852238306993, -0.333858906781618,
                    0.40392416563725, 0.402308231678916, 0.852555011730374, 0.217110354109593,
                    0.270330334450112, 0.369100176484466, 0.922641693423117, 0.157902495153081,
                    0.162717745954438, -0.161671904126578, 0.521750196577216, 0.0592943224139333,
                    -0.0314497827680782, -0.306451228103413, -0.0202446498099071,
                    0.667564987313085, 0.559387850451543, 1.32245673489777, -0.0494253613164233,
                    -0.635237189082218, 0.39557172025178, 0.472476674226215, 1.58561608333861,
                    0.553905331641403, 0.186056363202439, -0.4947007134665, -0.933535496496006,
                    1.09589878993216, 0.516198411247728, 0.0348569946536159, -0.964501644975371,
                    -0.595064527055996, 0.199828182268732, -0.217856002446865, 0.604731394943386,
                    -1.36385363518406, -1.28712589833233, 0.843423085670338, 1.23580213105345,
                    1.48953391296817, -0.3390301265147, -0.384884087955328, -0.0856016495781944,
                    0.958007709224944, 0.885658483223764, 0.392295293054133, -0.277442785396175,
                    -0.206627027364419, 1.13145462536065, 1.84108633039229, -0.109622680778358,
                    -0.86648006014704, -0.793545357358568, 1.23944574547798, 1.40114799413593,
                    -0.16665870913215, -0.629947922857667, 0.0872602700715959, 0.964495201616245,
                    1.11475743006879, -0.639874094836898, -0.483863645048416, 0.654677219243808,
                    0.797215726288186, 0.379499735897532, -0.621067071953263, -1.12804480712214,
                    1.49706686081956, 0.282610097276034, 0.637475225613812, -0.257106928406074,
                    -0.162952913329178, 0.478361506933624, 0.215614087773929, 0.708413131652353,
                    0.223324307850618, 0.461673858388252, 0.878735095377607, 0.624831929330782,
                    0.233201129979541, 0.689258785816833, 0.303256299709972, 0.225696006131017,
                    -0.266623523805549, -0.48841426163051, 0.362928934680554, -0.257465856633161,
                    -0.320124982215816, 0.0568872874773376, 0.0101563628661765, 0.518733362736725,
                    -0.0966573785391063, 0.849854445902423, -0.280729089414336, -0.0232281267143328,
                    0.50963745514963, 0.216891412313032, -0.0960198040330624, -0.406388990608122,
                    0.122738947220999, 0.171813673202621, 0.312867973620484, 0.555139501654622,
                    -0.296853553462284, 0.205665781836899, -0.0679119847165737, -0.278732040832807,
                    1.09068260287145, 0.604168836024428, 0.196297295742145, -0.186223164747503,
                    -0.340621513011882, 1.07903498836004, 0.382888975420229, -0.229930411287476,
                    -0.908625145165644, -1.13139402240276, 0.736426216826828, 0.293314716412794,
                    -0.0842391209668136, -0.72957749398509, -0.832994062348519, 1.38754478546417,
                    0.529182479758617, 0.849483044179046, -0.682485163626609, -0.351451029989846,
                    1.49159587639568, 1.00841545218682, 0.350734898441118, 0.538181092357718,
                    0.593409039817118, -0.782617381592132, 0.853490313922009, 0.363217338783777,
                    0.434920078853708, 0.301288085851525, -0.721839313401404, -0.217169559734515,
                    1.62831353186619, 0.179706895863859, 0.182086992337958, 0.0882085766052114,
                    0.149331997977953, 1.75930235879919, -0.485383190769312, -0.73342568456593,
                    0.567736305293181, 0.459726845484141, 0.850583448927129, -0.000423356049168433,
                    -0.278083391444341, 0.162287549244809, 0.369853735064254, 0.724681529089697,
                    -0.617428531707503, -0.749818973229623, 0.280018053164184, 0.108135050798126,
                    0.293519822068816, -0.783376186130937, -0.0503725114234328, -0.225123887279779,
                    -0.722528489815033, 0.411676402678513, -0.748203535317259, -0.279894740966881,
                    0.106572370257936, -0.527596557101191, 0.0410964595171777, -0.14483696833048,
                    -0.222920499965056, 0.0900088336080387, -0.760470638008826, -0.317870321818528,
                    -0.129709047782569, -0.0264804005870474, 0.626852035068424, 0.260470677425673,
                    0.043127127741279, -0.35628882979131, 0.17780669930923, 0.471951534087651,
                    0.176883136883248, 0.320843541287701, -1.17034668636251, -0.0920153585396122,
                    -0.541488716398277, 0.469246157908543, 0.162504894517149, -0.219399419681054,
                    0.249974073426358, -0.113885219726825, -0.124915357300523, 0.525546010420868,
                    0.348458288680177, -0.205688124193244, -0.0987803087155328, 0.169424422125701,
                    -0.0217417270840983, -0.154283047520266, -0.662239409721053,
                    -0.0224890572426943, 0.973748241458082, -0.643702702968675, -0.439691241589287,
                    -0.491426055821847, 1.03922685409663, 0.996241424806328, -0.475408777533094,
                    -0.190139238281773, -0.420827479853784, 0.169795816333675, 0.731205288012376,
                    -0.278316203454796, -0.802466258178583, -1.24695527873951, -0.0220983749977393,
                    0.943897770871729, 0.694237463688157, 0.120088616322363, 0.932140984797428,
                    0.991200628320257, 0.959296847301289, 0.609947370760896, 0.457104731313194,
                    0.312225793092352, -0.384616140459769, 0.0762979392441423, 0.681750128860475,
                    0.0523451414561862, -0.153249687529183, 0.34089367735001, 0.0781612225626702,
                    1.4254981555156, -0.447411441002174, -0.00370266134376607, 0.457850019222576,
                    -0.404265370664348, 0.546871289753774, 0.0321881060204134, -0.302270109599732,
                    -1.01214290595547, 0.251581963777852, 0.933252723497723, -0.890860760806873,
                    -0.537325063647842, 0.311289928955887, -0.543788079093076, -0.0129781668684892,
                    -0.70157779895564, -0.656234839205069, -0.268792389224416, -0.279917010956495,
                    0.0576751894617757, -0.593057861634662, -0.975539738976551, -0.103163563261916,
                    -0.0744649014627043, -0.23740866436323, -0.266332358159616, -0.399503725697087,
                    0.0961328530049741, -0.342094255479822, -0.662495117561335, -0.195941825015524,
                    -0.129719069805219, 0.759200067530742, 0.314001234611879, -0.351638939222919,
                    -0.294793242972133, 0.124673809721997, -0.418548730468046, 0.103219585531832,
                    0.0380860116067332, -0.119331108019765, 0.153272651734994, -0.117292721350395,
                    0.136685772506723, 0.494411293373393, -1.18768535875418, -0.234275706216494,
                    -0.390700853991973, -0.069975394044036, 0.193954638551464, -0.0294611188928968,
                    -0.00674732398474021, -0.318928947352987, 0.182205875045639,
                    -0.655277061978822, -0.297711092233072, 0.4529889185645, 0.441044015433931,
                    -0.421530815540308, -0.68080597256613, -0.688308914022656, -0.668880809834008,
                    0.0367976675129853, 0.294814603014885, -0.126944916214427, -0.541457356825447,
                    -1.01939294779382, -0.658751118128606, 0.663990571638891, 0.423090742463209,
                    -0.647390690192907, -0.904307156154715, -0.405706452558093, 1.27557255818184,
                    0.928408273033316, -0.125039230941411, 0.246455278268229, 0.56751425249783,
                    0.323601612505541, 0.117248252162046, 0.649477549050857, -0.130564212492146,
                    0.229291796948047, 0.355418638199637, 0.444125626207278, 0.207772617944215,
                    -0.00477862373270455, 0.683715239236221, -0.187823774125633,
                    0.595684405617856, -0.469545873138088, -1.18409342462306, 0.327225789823217,
                    0.0117030785007444, -0.0212848657425262, 0.296853233354142, -1.11616313539454,
                    0.181828086218928, -0.36629647331182, 0.498466589709231, -0.373102600914137,
                    -0.934814849099013, 0.501269899196927, -0.234749202197918, -0.157775015621478,
                    -0.0772066012658286, 0.370249260912971, 0.986991383042152, -0.467273637254067,
                    -0.2315310136139, -0.222785716190486, -0.304031851842836, -0.743532584646942,
                    -0.370389580267954, -0.124120925856239, 0.266369669218278, 0.0149364753144141,
                    -0.252300484491798, -0.246448472797795, -0.316128533946827, -0.366365283623423,
                    -0.423922106813222, 0.0179015148958399, -0.0428307885661638,
                    -0.858035745135169, -0.173201588661044, -0.0358971005259434,
                    -0.693847930420162, -0.925964273293658, -0.0994796644645781,
                    -0.187349282514589, -0.341543676903171, -0.0832106998030895,
                    -0.155999104483953, 0.199860161027924, -0.109453265656413, -0.186569580197812,
                    0.118408089884184, -0.520018820690723, 0.247604893212188, -0.323010109180103,
                    0.273888702581729, -0.356143556789959, -0.385993139948283, -0.949357495361731,
                    -1.17881602319764, 0.616930559851056, -0.00427422408364464, -0.503457728509943,
                    -1.37202170878347, -0.577600892524543, 0.120966317883699, -1.60922709416198,
                    0.000741363813517637, -0.52061006546144, -1.10652268136374, -0.979336439901459,
                    -0.201602269398734, 0.846554778276726, 0.323939536542309, -1.19994363013092,
                    -1.51976336466063, 0.840807447950908, 2.18935015275767, -0.524071782259909,
                    0.192107540078382, -0.529400156512423, 0.718607657143593, 1.63720730534925,
                    -0.0571795284748973, 0.15435263248969, -0.43738085305759, 0.584089613846118,
                    0.686256830950612, 0.235904878472084, -0.0242076442036989, -0.30328124142406,
                    0.404291192090056, 0.579814998803641, 0.56811020280383, 0.00859344906783122,
                    -0.363606816207891, 0.291029294146711, 0.584416655507596, 0.979037069685049,
                    0.0284542394431355, -0.984805042663339, -0.347709442333752, 0.0290728750797768,
                    1.0011862535768, -0.265138320578473, -0.484598400043211, 0.460526106639833,
                    0.0341768376712253, 0.43071335553396, -0.153274387105938, -0.215282873527144,
                    0.00241288396588857, -0.258869550797556, 0.177016162275084, 0.538671053839181,
                    -0.811080885441442, -0.890892950808308, -0.584915276251422, 0.123911320898781,
                    -0.302060165929731, -0.749936961853933, -0.639675940215341, -0.57396423394195,
                    -0.205309667175104, 0.0836188065133001, -0.352349767117772, -1.15565771424892,
                    -0.325214545841533, -0.501306372965846, 0.226050111379976, 0.770161651833661,
                    -0.707011396350879, -0.788490055155222, -0.277207269615797, 0.0385331318053659,
                    0.523870430885808, -1.39783629209302, -1.02272138300407, -0.136810504576133,
                    -0.188217876616164, 0.456930788782017, -1.56128759040399, -0.865808698754677,
                    -1.33485493675594, -0.182915404194489, 1.16335549882447, -1.07067846399082,
                    -0.329229431365576, -1.00565024843392, 0.293257431829971, 1.25734860675354,
                    -0.753815285433944, 0.0846227388065017, -1.06188202351251, -0.275668911814364,
                    0.0572553903388666, -0.784550285687004, 1.10805475800635, -0.4810831817496,
                    -0.533923403636322, 0.623638687117051, -0.507971938558569, 0.485554928775684,
                    -0.272012007795311, -1.16923255594463, 0.185953144334434, 0.758304992878891,
                    1.48566600804006, -0.182101770139635, 0.0300025354698903, -0.247174539033722,
                    0.823153308712563, 1.04014517792123, -0.30196309481024, -0.854946685924574,
                    -0.164007170592271, 1.30907004850616, 1.03326901177303, 0.0257889379678993,
                    -1.00502796460361, -1.57571967026519, 1.10036422855717, 0.964178672493574,
                    0.390514507593675, -0.880728238222815, -1.32109854611895, 0.538405885119147,
                    0.620578436774252, 0.645247876007037, -0.118148050093765, -0.866217515446806,
                    0.396490157618842, 0.359896409920538, 1.00828993007394, -0.440460952183031,
                    -0.613797342174024, 0.188987851326772, 0.685068449655034, 0.0855077448629338,
                    -0.41673442162076, -0.81780556862509, 0.37671688739232, -0.331998808655158,
                    0.321965908363815, -0.163946080327605, -1.11301222814947, -0.223384566885585,
                    0.131484161199609, 0.0373428658519437, -0.193101174672194, -0.386849310475696,
                    -0.474418794608454, -0.284751012326836, -0.212056859310636, -0.851652038844435,
                    -0.0891301555192285, -1.15607535654718, -0.489507376188669, -0.355860680640262,
                    -1.03957295605258, 0.902943064101976, -0.471091050564127, -0.545676200088252,
                    -0.887877386425252, -0.894054846258172, 0.848299184780672, -1.17878522817267,
                    0.139011779585445, -0.533956362799213, -0.576078410501282, 0.758561302751588,
                    -1.18127255970339, -0.702658634833483, -1.05027076497971, -0.479744877584864,
                    1.2602667242783, -1.36561255150102, -0.685018931872548, -1.11509344933721,
                    -0.0135890286835831, 0.313777297840431, -2.27115180929206, 0.215212105439235,
                    -1.18692801740119, -0.182448992210652, -0.940130633926856, -0.978150404204428,
                    0.774826959661624, -1.35294932657882, -0.247499546548255, 0.198633332525731,
                    -1.25862920890484, 0.462575395342299, -0.598635263251634, -0.825827219628252,
                    -0.143588061114713, -0.47159237353992, 1.21015150173378, 0.178972116523234,
                    -0.0920072508409646, -1.06062834133746, 2.60462947315494, 1.07318740157663,
                    1.02050551404063, -1.42976268627939, -0.807242133204227, 1.6952338601032,
                    0.569989618724179, 0.972951621625725, -0.75512139933079, -0.718866545922566,
                    1.7700824002813, 0.694398190197433, 0.682014688998429, -0.0560211817035281,
                    -0.640014891093557, 1.60122100337028, 0.742708153401736, 0.551155987585528,
                    -0.276518271983264, -0.162154885576541, 1.00632221591351, 1.30195440980947,
                    1.12435882833288, -0.635377524382534, -0.42673343536805, 0.634367159764025,
                    1.20996115879123, 0.728328345132812, 0.227620819560968, 0.749977374741292,
                    1.30355516715872, -0.0206757398025559, 0.917547942853009, -0.0441373071143687,
                    0.210265062959081, 0.98711673104404, -0.338682229990354, -0.220713746266544,
                    -0.608381638799301, 0.222698436575826, 0.956116971899888, -0.0724079182652417,
                    -0.00279966914293439, -0.462278248523177, 1.150800694059, 0.274703089385156,
                    0.0104135580047569, -0.229800582996635, -0.500937395887121, 0.660921500392277,
                    -0.111173874910293, 0.0517512546748833, -0.360796032446586, -0.867175327197174,
                    0.248208002917238, -0.217650857489559, 0.337223855506807, -0.610812536394768,
                    -0.914656633443741, -0.0531929552185488, -1.57445216227393, -0.607936325547513,
                    -1.35113545326963, -0.938524959141859, 0.0144087690620438, -0.966035418694901,
                    0.241207431808704, -1.01005327048578, -0.155890591369626, -0.869618051964892,
                    -1.12190110972753, 0.438617570588832, -1.12354273800673, -0.504432995851914,
                    0.257907950696989, -0.852860756641604, 0.591331402634893, -1.63420526933004,
                    -0.40984989921958, -1.01640092068078, 0.0672496573932619, 0.724541048149433,
                    -1.11914619132559, -1.34896429931234, -0.841129049504246, 0.379642193003718,
                    1.15987585301559, 0.432827017253683, -0.0704807094820126, -1.18894471893296,
                    2.98981039740966, 1.71370705887315, 0.312680989789737, -0.271675234949385,
                    -0.246772661741034, 3.17536982692309, 0.485314939042723, 0.446152097042827,
                    -0.133883483207569, -1.2998787496972, 2.61524247642902, 0.509658530376555,
                    1.14076365514336, -0.0129961296911745, -0.568233837961657, 2.48222878006004,
                    0.474591140479532, 0.758651575904299, 0.332622630675104, -0.707093296962455,
                    1.68257442807794, 0.673779157793632, 1.0177535021429, 0.35581349064779,
                    1.14471236819821, 2.24059453000891, 0.4956048430777, 0.89332464108179,
                    0.834250013395177, -0.157878639078973, 1.45431768477559, 0.305099245457683,
                    0.325204548641523, 1.28048052599897, -0.718526473731544, 1.15764753989993,
                    -0.241095257737356, 0.249833593814166, 0.738056868437525, 0.700008472041901,
                    1.20467819336602, -0.297923852008466, -0.416971888206902, 0.54194980677255,
                    1.18582739816711, 1.23988857152059, -0.180841328291564, -0.0786234163794113,
                    0.00615287143983778, 0.659089104598991, -0.387961348423914, -0.011247536108101,
                    -0.595003382872594, -0.276147616131738, 0.138952913820617, -0.587593575433146,
                    -0.0701986385299165, -0.920310765498634, -0.678045987127206,
                    -0.0548015032268836, 0.0956197842231177, -0.406672858409138,
                    -0.788089162136171, -0.475918302946411, -0.0767468440641892,
                    -0.506052412951988, -0.018162929946584, -0.75612011089396, -0.549895624752054,
                    0.505908830030108, -0.690780760615922, -0.170334048991208, -0.641417294988699,
                    -0.956107762107536, 0.16732473262249, -0.399599463902433, 0.0647530105175917,
                    -0.322952536101919, -0.476828442398357, -0.185620055898454, -0.0442891759009569,
                    0.744039068251737, -1.3260441418014, -0.596417418679351, -0.201836931621585,
                    0.423541386390497, 0.743373223670915, 0.428870019332724, -0.21355002769775,
                    -1.60811748366871, 1.98204710853031, 1.26160596137106, -0.154921107695071,
                    -0.549371186922343, -0.82292072214716, 2.53775755878647, 2.07444320862774,
                    -0.403696336248061, -0.677926985173571, -1.19332459149689, 2.99055153693121,
                    1.90095887408268, 0.549779366780569, -0.629016163534061, -0.633720717565382,
                    2.56279943766489, 1.23820872416761, 0.878306916243979, 0.117756434999279,
                    -0.656160468537422, 1.64506467555279, 0.582072770448747, 0.973543849086298,
                    1.46100287712721, 0.506166332122746, 0.966363704592283, 1.12664205103624,
                    0.508768251489258, 1.36392519055252, 0.451144489882154, 1.05302552542315,
                    1.00559152431634, 0.0408412406092609, 0.965826083134061, 0.171706992013924,
                    1.37950114427301, -0.142232908809717, 0.0524502885615981, 0.476995581089966,
                    0.602694049319439, 0.839643071513734, -0.642680085744274, -0.0139812760225979,
                    0.174353407866448, 0.154058076782036, 0.994560689457614, 0.227598619026338,
                    0.122808707275777, 0.545633388715004, 0.161071619312452, -0.379705189617226,
                    -1.04661700973352, -0.509292067860887, 0.336672319673454, -0.0524042907718183,
                    0.145535795144468, -0.130278177392967, -0.449992003323537, 0.0799677648682575,
                    -0.427851945724939, 0.515420552930721, -0.0311054361585832, -0.401146940562274,
                    -0.0111838860377986, 0.464810040997406, -1.14451001481737, -0.298289572584395,
                    -0.808039902200122, -0.361925171725403, 0.491178763273155, -0.494950757798136,
                    -0.513601443358397, -0.524699484240635, -1.02769897715515, 0.806825800924372,
                    0.509685549242391, -0.290455195308226, 0.134083008869055, -0.901632603732992,
                    -0.129550073981562, 0.304558723466982, 0.508127811148427, 0.120807738475452,
                    -0.774611102070139, 0.216569102964712, -0.554951970866682), .Dim = c(5L,
                                                                                         520L)))



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
        c(1.04217141145773, 1.04037731093668, 0.99985058907344, 0.781924338781402,
          0.781924338781402))

    class(out) <- "aaa"
    expect_error(unscale_beta(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
    expect_error(unscale_sigma2(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
})


