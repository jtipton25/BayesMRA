context("testing plotting functions")

# plot_fitted_alphas -----------------------------------------------------------

test_that("plot_fitted_alphas", {
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
    p <- plot_fitted_alphas(out)
    expect_true(is.ggplot(p))
    expect_error(plot_fitted_alphas(out, width = -1), "width must be a positive number")
    expect_error(plot_fitted_alphas(out, width = NA), "width must be a positive number")
    expect_error(plot_fitted_alphas(out, width = "aaa"), "width must be a positive number")
    expect_error(plot_fitted_alphas(out, height = -1), "height must be a positive number")
    expect_error(plot_fitted_alphas(out, height = NA), "height must be a positive number")
    expect_error(plot_fitted_alphas(out, height = "aaa"), "height must be a positive number")
    expect_error(plot_fitted_alphas(out, base_size = -1), "base_size must be a positive number")
    expect_error(plot_fitted_alphas(out, base_size = NA), "base_size must be a positive number")
    expect_error(plot_fitted_alphas(out, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_fitted_alphas(out, file = NA), "file must be a character string")
    expect_error(plot_fitted_alphas(out, file = NA), "file must be a character string")
    # check writing to tempfile
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    # suppress the geom_raster warning
    suppressWarnings(plot_fitted_alphas(out, file = file))
    expect_true(file.exists(file))
    file.remove(file)

    class(out) <- NULL
    expect_error(plot_fitted_alphas(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")

})

# plot_fitted_MRA --------------------------------------------------------------

test_that("plot_fitted_MRA", {
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
    preds <- predict_mra(out, new_data)
    p <- plot_fitted_MRA(out, preds)
    expect_true(is.ggplot(p))
    expect_error(plot_fitted_MRA(out, preds, width = -1), "width must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, width = NA), "width must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, width = "aaa"), "width must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, height = -1), "height must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, height = NA), "height must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, height = "aaa"), "height must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, base_size = -1), "base_size must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, base_size = NA), "base_size must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_fitted_MRA(out, preds, file = NA), "file must be a character string")
    expect_error(plot_fitted_MRA(out, preds, file = NA), "file must be a character string")
    # check writing to tempfile
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_fitted_MRA(out, preds, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(out) <- NULL
    expect_error(plot_fitted_MRA(out, preds), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
    class(out) <- "mcmc_mra"
    class(preds) <- NULL
    expect_error(plot_fitted_MRA(out, preds), "preds must be of class \"mcmc_mra_pred\" which is the output of the function predict_mra()")

})

# plot_fitted_process ----------------------------------------------------------

test_that("plot_fitted_process", {
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
    data <- data.frame(Lat = locs[, 1], Lon = locs[, 2], Observed = y)
    p <- plot_fitted_process(out, data)
    expect_true(is.ggplot(p))

    data <- list(Lat = locs[, 1], Lon = locs[, 2], Observed = y)
    expect_error(plot_fitted_process(out, data), "data must be a data.frame with variables Observed, Lat, and Lon")
    data <- data.frame(Lat = locs[, 1], Lon = locs[, 2], Observed = NA)
    expect_error(plot_fitted_process(out, data), "The data.frame data must contain a numeric vector named Observed of the observed values")
    data <- data.frame(Lat = locs[, 1], Lon = NA, Observed = y)
    expect_error(plot_fitted_process(out, data), "The data.frame data must contain a numeric vector named Lon of the longitude locations")
    data <- data.frame(Lat = NA, Lon = locs[, 2], Observed = y)
    expect_error(plot_fitted_process(out, data), "The data.frame data must contain a numeric vector named Lat of the latitude locations")

    expect_error(plot_fitted_process(out, data, width = -1), "width must be a positive number")
    expect_error(plot_fitted_process(out, data, width = NA), "width must be a positive number")
    expect_error(plot_fitted_process(out, data, width = "aaa"), "width must be a positive number")
    expect_error(plot_fitted_process(out, data, height = -1), "height must be a positive number")
    expect_error(plot_fitted_process(out, data, height = NA), "height must be a positive number")
    expect_error(plot_fitted_process(out, data, height = "aaa"), "height must be a positive number")
    expect_error(plot_fitted_process(out, data, base_size = -1), "base_size must be a positive number")
    expect_error(plot_fitted_process(out, data, base_size = NA), "base_size must be a positive number")
    expect_error(plot_fitted_process(out, data, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_fitted_process(out, data, file = NA), "file must be a character string")
    expect_error(plot_fitted_process(out, data, file = NA), "file must be a character string")
    # check writing to tempfile
    data <- data.frame(Lat = locs[, 1], Lon = locs[, 2], Observed = y)
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_fitted_process(out, data, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(out) <- NULL
    expect_error(plot_fitted_process(out, data), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
    class(out) <- "mcmc_mra"

})

# plot_MRA_grid ----------------------------------------------------------------

test_that("plot_MRA_grid", {
    locs <- matrix(rnorm(20), 10, 2)
    MRA <- mra_wendland_2d(locs, M = 2, n_coarse_grid = 4)
    p <- plot_MRA_grid(MRA)
    expect_true(is.ggplot(p))
    expect_error(plot_MRA_grid(MRA, width = -1), "width must be a positive number")
    expect_error(plot_MRA_grid(MRA, width = NA), "width must be a positive number")
    expect_error(plot_MRA_grid(MRA, width = "aaa"), "width must be a positive number")
    expect_error(plot_MRA_grid(MRA, height = -1), "height must be a positive number")
    expect_error(plot_MRA_grid(MRA, height = NA), "height must be a positive number")
    expect_error(plot_MRA_grid(MRA, height = "aaa"), "height must be a positive number")
    expect_error(plot_MRA_grid(MRA, base_size = -1), "base_size must be a positive number")
    expect_error(plot_MRA_grid(MRA, base_size = NA), "base_size must be a positive number")
    expect_error(plot_MRA_grid(MRA, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_MRA_grid(MRA, file = NA), "file must be a character string")
    expect_error(plot_MRA_grid(MRA, file = NA), "file must be a character string")
    # check writing to tempfile
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_MRA_grid(MRA, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(MRA) <- NULL
    expect_error(plot_MRA_grid(MRA), 'MRA must be of class "mra_wendland_2d"')

    })


# plot_posterior_params --------------------------------------------------------

test_that("plot_posterior_params", {
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
    p <- plot_posterior_params(out)
    expect_true(is.ggplot(p))
    p <- plot_posterior_params(out, alpha = 0.5)
    expect_true(is.ggplot(p))

    expect_error(plot_posterior_params(out, alpha = -1), "alpha must be a number between 0 and 1")
    expect_error(plot_posterior_params(out, alpha = 2.5), "alpha must be a number between 0 and 1")
    expect_error(plot_posterior_params(out, alpha = NA), "alpha must be a number between 0 and 1")
    expect_error(plot_posterior_params(out, alpha = "aaa"), "alpha must be a number between 0 and 1")

    expect_error(plot_posterior_params(out, width = -1), "width must be a positive number")
    expect_error(plot_posterior_params(out, width = NA), "width must be a positive number")
    expect_error(plot_posterior_params(out, width = "aaa"), "width must be a positive number")
    expect_error(plot_posterior_params(out, height = -1), "height must be a positive number")
    expect_error(plot_posterior_params(out, height = NA), "height must be a positive number")
    expect_error(plot_posterior_params(out, height = "aaa"), "height must be a positive number")
    expect_error(plot_posterior_params(out, base_size = -1), "base_size must be a positive number")
    expect_error(plot_posterior_params(out, base_size = NA), "base_size must be a positive number")
    expect_error(plot_posterior_params(out, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_posterior_params(out, file = NA), "file must be a character string")
    expect_error(plot_posterior_params(out, file = NA), "file must be a character string")
    # check writing to tempfile
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_posterior_params(out, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(out) <- NULL
    expect_error(plot_posterior_params(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")

})

# plot_predicted_process -------------------------------------------------------

test_that("plot_predicted_process", {
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
    data <- data.frame(Lat = locs[, 1], Lon = locs[, 2], Observed = y)
    new_data <- list(locs_pred = locs, X_pred = X)
    preds <- predict_mra(out, new_data)
    p <- plot_predicted_process(out, data, preds)
    expect_true(is.ggplot(p))
    expect_error(plot_predicted_process(out, data, preds, width = -1), "width must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, width = NA), "width must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, width = "aaa"), "width must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, height = -1), "height must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, height = NA), "height must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, height = "aaa"), "height must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, base_size = -1), "base_size must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, base_size = NA), "base_size must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_predicted_process(out, data, preds, file = NA), "file must be a character string")
    expect_error(plot_predicted_process(out, data, preds, file = NA), "file must be a character string")
    # check writing to tempfile
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_predicted_process(out, data, preds, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(out) <- NULL
    expect_error(plot_predicted_process(out, data, preds), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
    class(out) <- "mcmc_mra"
    class(preds) <- NULL
    expect_error(plot_predicted_process(out, data, preds), "preds must be of class \"mcmc_mra_pred\" which is the output of the function predict_mra()")

})


# plot_sim_vs_fitted -----------------------------------------------------------

test_that("plot_sim_vs_fitted", {
    y <- rnorm(10)
    X <- matrix(1:20, 10, 2)
    beta <- rnorm(ncol(X))
    locs <- matrix(rnorm(20), 10, 2)
    MRA <- mra_wendland_2d(locs, M = 2, n_coarse_grid = 4)
    W <- MRA$W
    alpha <- rnorm(ncol(MRA$W))
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
    p <- plot_sim_vs_fitted(out, X, beta, W, alpha)
    expect_true(is.ggplot(p))

    # check data
    expect_error(plot_sim_vs_fitted(out, X, cbind(beta, beta), W, alpha), "beta must be a numeric vector")
    expect_error(plot_sim_vs_fitted(out, X, rep(NA, length(beta)), W, alpha), "beta must be a numeric vector")
    expect_error(plot_sim_vs_fitted(out, X, rep("aaa", length(beta)), W, alpha), "beta must be a numeric vector")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, cbind(alpha, alpha)), "alpha must be a numeric vector")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, rep(NA, length(alpha))), "alpha must be a numeric vector")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, rep("aaa", length(alpha))), "alpha must be a numeric vector")

    expect_error(plot_sim_vs_fitted(out, cbind(X, 1), beta, W, alpha), "X must have the same number of columns as the length of beta")
    expect_error(plot_sim_vs_fitted(out, X, c(beta, 1), W, alpha), "X must have the same number of columns as the length of beta")

    expect_error(plot_sim_vs_fitted(out, rbind(X, 1), beta, W, alpha), "X must have the same number of rows as the basis function matrix W")
    expect_error(plot_sim_vs_fitted(out, X, beta, cbind(W, 1), alpha), "W must have the same number of columns as the length of alpha")
    expect_error(plot_sim_vs_fitted(out, rbind(X, 1), beta, cbind(W, 1), alpha), "W must have the same number of columns as the length of alpha")
    expect_error(plot_sim_vs_fitted(out, X, beta, rbind(W, rep(1, ncol(W))), alpha), "The basis function matrix W must have the same number of rows as the fitted basis function matrix  W in the MRA object contained within the out object")


    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, width = -1), "width must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, width = NA), "width must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, width = "aaa"), "width must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, height = -1), "height must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, height = NA), "height must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, height = "aaa"), "height must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, base_size = -1), "base_size must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, base_size = NA), "base_size must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, file = NA), "file must be a character string")
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha, file = NA), "file must be a character string")
    # check writing to tempfile
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_sim_vs_fitted(out, X, beta, W, alpha, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(out) <- NULL
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")
    class(out) <- "mcmc_mra"
    class(W) <- "aaa"
    expect_error(plot_sim_vs_fitted(out, X, beta, W, alpha), 'W must be of class "spam" and is typically in an object created by "mra_wendland_2d\\(\\)"')

})

# plot_sparse_matrix

test_that("plot_sparse_matrix", {
    locs <- matrix(rnorm(20), 10, 2)
    MRA <- mra_wendland_2d(locs, M = 2, n_coarse_grid = 4)
    p <- plot_sparse_matrix(MRA$W, base_size = 9, title = "A title")
    expect_true(is.ggplot(p))

    expect_error(plot_sparse_matrix(MRA$W, tile_size = -1), "tile_size must be a positive number")
    expect_error(plot_sparse_matrix(MRA$W, tile_size = NA), "tile_size must be a positive number")
    expect_error(plot_sparse_matrix(MRA$W, tile_size = "aaa"), "tile_size must be a positive number")
    expect_error(plot_sparse_matrix(MRA$W, tile_size = TRUE), "tile_size must be a positive number")

    expect_error(plot_sparse_matrix(MRA$W, base_size = -1), "base_size must be a positive number")
    expect_error(plot_sparse_matrix(MRA$W, base_size = NA), "base_size must be a positive number")
    expect_error(plot_sparse_matrix(MRA$W, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_sparse_matrix(MRA$W, base_size = TRUE), "base_size must be a positive number")

    expect_error(plot_sparse_matrix(MRA$W, title = -1), "title must be a character string")
    expect_error(plot_sparse_matrix(MRA$W, title = NA), "title must be a character string")
    expect_error(plot_sparse_matrix(MRA$W, title = TRUE), "title must be a character string")
})


# plot_test_data ---------------------------------------------------------------

test_that("plot_test_data", {
    data("code_test")
    code_test <- code_test %>%
        mutate(Observed = MaskTemp, Truth = TrueTemp)
    p <- plot_test_data(code_test)
    expect_true(is.ggplot(p))

    expect_error(plot_test_data(select(code_test, - Observed)), "The data.frame data must contain a numeric vector named Observed of the observed values")
    expect_error(plot_test_data(select(code_test, - Truth)), "The data.frame data must contain a numeric vector named Truth of the full values")
    expect_error(plot_test_data(select(code_test, - Lat)), "The data.frame data must contain a numeric vector named Lat of the latitude locations")
    expect_error(plot_test_data(select(code_test, - Lon)), "The data.frame data must contain a numeric vector named Lon of the longitude locations")

    expect_error(plot_test_data(code_test, width = -1), "width must be a positive number")
    expect_error(plot_test_data(code_test, width = NA), "width must be a positive number")
    expect_error(plot_test_data(code_test, width = "aaa"), "width must be a positive number")
    expect_error(plot_test_data(code_test, height = -1), "height must be a positive number")
    expect_error(plot_test_data(code_test, height = NA), "height must be a positive number")
    expect_error(plot_test_data(code_test, height = "aaa"), "height must be a positive number")
    expect_error(plot_test_data(code_test, base_size = -1), "base_size must be a positive number")
    expect_error(plot_test_data(code_test, base_size = NA), "base_size must be a positive number")
    expect_error(plot_test_data(code_test, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_test_data(code_test, file = NA), "file must be a character string")
    expect_error(plot_test_data(code_test, file = NA), "file must be a character string")
    # check writing to tempfile
        file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_test_data(code_test, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(code_test) <- NULL
    expect_error(plot_test_data(code_test), "data must be a data.frame with variables Observed, Lat, and Lon")
})

# plot_trace -------------------------------------------------------------------

test_that("plot_trace", {
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
    p <- plot_trace(out)
    expect_true(is.ggplot(p))

    expect_error(plot_trace(out, width = -1), "width must be a positive number")
    expect_error(plot_trace(out, width = NA), "width must be a positive number")
    expect_error(plot_trace(out, width = "aaa"), "width must be a positive number")
    expect_error(plot_trace(out, height = -1), "height must be a positive number")
    expect_error(plot_trace(out, height = NA), "height must be a positive number")
    expect_error(plot_trace(out, height = "aaa"), "height must be a positive number")
    expect_error(plot_trace(out, base_size = -1), "base_size must be a positive number")
    expect_error(plot_trace(out, base_size = NA), "base_size must be a positive number")
    expect_error(plot_trace(out, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_trace(out, file = NA), "file must be a character string")
    expect_error(plot_trace(out, file = NA), "file must be a character string")
    # check writing to tempfile
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_trace(out, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(out) <- NULL
    expect_error(plot_trace(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")

})

# plot_trace_alpha -------------------------------------------------------------

test_that("plot_trace_alpha", {
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
    p <- plot_trace_alpha(out)
    expect_true(is.ggplot(p))
    p <- plot_trace_alpha(out, alpha = 0.5)
    expect_true(is.ggplot(p))

    expect_error(plot_trace_alpha(out, alpha = -1), "alpha must be a number between 0 and 1")
    expect_error(plot_trace_alpha(out, alpha = 2.5), "alpha must be a number between 0 and 1")
    expect_error(plot_trace_alpha(out, alpha = NA), "alpha must be a number between 0 and 1")
    expect_error(plot_trace_alpha(out, alpha = "aaa"), "alpha must be a number between 0 and 1")

    expect_error(plot_trace_alpha(out, width = -1), "width must be a positive number")
    expect_error(plot_trace_alpha(out, width = NA), "width must be a positive number")
    expect_error(plot_trace_alpha(out, width = "aaa"), "width must be a positive number")
    expect_error(plot_trace_alpha(out, height = -1), "height must be a positive number")
    expect_error(plot_trace_alpha(out, height = NA), "height must be a positive number")
    expect_error(plot_trace_alpha(out, height = "aaa"), "height must be a positive number")
    expect_error(plot_trace_alpha(out, base_size = -1), "base_size must be a positive number")
    expect_error(plot_trace_alpha(out, base_size = NA), "base_size must be a positive number")
    expect_error(plot_trace_alpha(out, base_size = "aaa"), "base_size must be a positive number")
    expect_error(plot_trace_alpha(out, file = NA), "file must be a character string")
    expect_error(plot_trace_alpha(out, file = NA), "file must be a character string")
    # check writing to tempfile
    file <- paste(tempfile(), "png", sep = ".")
    expect_false(file.exists(file))
    plot_trace_alpha(out, file = file)
    expect_true(file.exists(file))
    file.remove(file)

    class(out) <- NULL
    expect_error(plot_trace_alpha(out), "out must be of class \"mcmc_mra\" or \"mcmc_mra_integrated\"")

})
