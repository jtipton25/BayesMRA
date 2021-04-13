library(spam)
library(tidyverse)
library(patchwork)
library(BayesMRA)

# setup image directory
if (!dir.exists(here::here("images"))) {
    dir.create(here::here("images"))
}
if (!dir.exists(here::here("images", "test-data"))) {
    dir.create(here::here("images", "test-data"))
}

# load the data
data("code_test", package = "BayesMRA")

locs <- cbind(code_test$Lat, code_test$Lon)

## create a mask in the upper corner
code_test$MaskTemp[code_test$Lat > 0.60 & code_test$Lon > 0.60] <- NA

# Plot the data ----------------------------------------------------------------

dat_plot <- data.frame(
    Lat = code_test$Lat,
    Lon = code_test$Lon,
    Observed = code_test$MaskTemp,
    Truth = code_test$TrueTemp)

plot_test_data(dat_plot)

# Setup the model -------------------------------------------------------------
dat_fit <- code_test %>%
    filter(!is.na(MaskTemp))

y <- dat_fit$MaskTemp
X <- model.matrix(~ Lon + Lat, data = dat_fit)
locs <- cbind(dat_fit$Lat, dat_fit$Lon)

params <- list(
    n_mcmc    = 1000,
    n_adapt   = 1000,
    n_thin    = 1,
    n_message = 50
)

priors <- list(
    alpha_tau2   = 1,
    beta_tau2    = 1,
    alpha_sigma2 = 1,
    beta_sigma2  = 1,
    mu_beta      = rep(0, ncol(X)),
    Sigma_beta   = 5 * diag(ncol(X)))

# Fit unconstrained model -----------------------------------------------------

if (file.exists(here::here("results", "test-fit-resolution-unconstrained.RData"))) {
    load(here::here("results", "test-fit-resolution-unconstrained.RData"))
} else {
    start   <- Sys.time()

    # to do: center and scale X
    out <- mcmc_mra(
        y = y,
        X = X,
        locs = locs,
        params = params,
        priors = priors,
        M = 3,
        n_coarse_grid = 20,
        joint = FALSE,
        constraint = "unconstrained")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "test-fit-resolution-unconstrained.RData"))
    # cell phone notification when done
    pushoverr::pushover(message = "Finished fitting unconstrained resolution first order model")
}

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 224
sum(out$MRA$W %*% out$alpha[iter_idx, ])
sum(out$MRA$W[, out$MRA$dims_idx == 1] %*% out$alpha[iter_idx, out$MRA$dims_idx == 1])
sum(out$MRA$W[, out$MRA$dims_idx == 2] %*% out$alpha[iter_idx, out$MRA$dims_idx == 2])
sum(out$MRA$W[, out$MRA$dims_idx == 3] %*% out$alpha[iter_idx, out$MRA$dims_idx == 3])

# Trace plots of parameters ---------------------------------------------------

plot_trace(out)

# Trace plots for spatial random effects  -------------------------------------
plot_trace_alpha(out)

# Plot Posterior distributions-------------------------------------------------
plot_posterior_params(out)

# Plot fitted process ----------------------------------------------------------

plot_fitted_process(out, dat_plot)

# Predict spatial process ------------------------------------------------------

new_data <- list(
    locs_pred = cbind(code_test$Lat, code_test$Lon),
    X_pred    = model.matrix(~ Lon + Lat, data = code_test)
)

preds <- predict_mra(out, new_data)

plot_predicted_process(out, dat_plot, preds)

# Plot fitted spatial random effects alpha -------------------------------------

plot_fitted_alphas(out)

# Plot fitted MRA resolutions --------------------------------------------------

plot_fitted_MRA(out, preds)












# Fit constrained model --------------------------------------------------------

if (file.exists(here::here("results", "test-fit-resolution-constrained.RData"))) {
    load(here::here("results", "test-fit-resolution-constrained.RData"))
} else {
    start   <- Sys.time()
    out <- mcmc_mra(
        y = y,
        X = X,
        locs = locs,
        params = params,
        priors = priors,
        M = 3,
        n_coarse_grid = 20,
        joint = FALSE,
        constraint = "resolution")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "test-fit-resolution-constrained.RData"))
    # cell phone notification when done
    pushoverr::pushover(message = "Finished fitting constrained resolution first order model")
}

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 211
sum(out$MRA$W %*% out$alpha[iter_idx, ])
sum(out$MRA$W[, out$MRA$dims_idx == 1] %*% out$alpha[iter_idx, out$MRA$dims_idx == 1])
sum(out$MRA$W[, out$MRA$dims_idx == 2] %*% out$alpha[iter_idx, out$MRA$dims_idx == 2])
sum(out$MRA$W[, out$MRA$dims_idx == 3] %*% out$alpha[iter_idx, out$MRA$dims_idx == 3])

# Trace plots of parameters ----------------------------------------------------

plot_trace(out)

# Trace plots for spatial random effects  --------------------------------------

plot_trace_alpha(out)

# Plot Posterior distributions--------------------------------------------------

plot_posterior_params(out)

# Plot fitted process ----------------------------------------------------------

plot_fitted_process(out, dat_plot)

# Predict spatial process ------------------------------------------------------

new_data <- list(
    locs_pred = cbind(code_test$Lat, code_test$Lon),
    X_pred    = model.matrix(~ Lon + Lat, data = code_test)
)

preds <- predict_mra(out, new_data)

plot_predicted_process(out, dat_plot, preds)

# Plot fitted spatial random effects alpha -------------------------------------

plot_fitted_alphas(out)

# Plot fitted MRA resolutions --------------------------------------------------

plot_fitted_MRA(out, preds)


