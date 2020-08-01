context("mcmc functions")

test_that("mcmc_mra", {
    y <- rnorm(10)
    X <- matrix(1:20, 10, 2)
    params <- list(
        n_mcmc    = 500,
        n_adapt   = 500,
        n_thin    = 1,
        n_message = 50
    )
    locs <- matrix(1:20, 10, 2)

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
    expect_error(mcmc_mra(y, X, locs, params), "X must be a numeric matrix with N rows.")
    X <- matrix(0, 5, 3)
    expect_error(mcmc_mra(y, X, locs, params), "X must have the same number of rows as the length of y.")

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

    # mcmc_mra(y, X, locs, params, priors,  M = 2, n_coarse_grid = 4)

    # expect_error(mcmc_mra(y, X, params), 'argument "priors" is missing, with no default')
    # priors <- default_priors_pg_lm(y, X)
    # out <- pg_lm(y, X, params, priors)
    # expect_true(class(out) == "pg_lm")

})

