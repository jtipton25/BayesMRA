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

# clean up this data object
# load(here::here("data", "SmallTestData.RData"))
data("SmallTestData", package = "BayesMRA")

code.test$MaskTemp = code.test$MaskedData
code.test$TrueTemp = code.test$FullData

locs <- cbind(code.test$Lat, code.test$Lon)

## create a mask in the upper corner
code.test$MaskTemp[code.test$Lat > 0.60 & code.test$Lon > 0.60] <- NA

# Plot the data ----------------------------------------------------------------

dat_plot <- data.frame(
    Lat = code.test$Lat,
    Lon = code.test$Lon,
    Observed = code.test$MaskTemp,
    Truth = code.test$TrueTemp)

plot_test_data(dat_plot, file = here::here("images", "test-data", "test-data.png"))

# Setup the model --------------------------------------------------------------
dat_fit <- code.test %>%
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

# Fit unconstrained model ------------------------------------------------------

if (file.exists(here::here("results", "test-fit-unconstrained.RData"))) {
    load(here::here("results", "test-fit-unconstrained.RData"))
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
        joint = TRUE,
        constraint = "unconstrained")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "test-fit-unconstrained.RData"))
    # cell phone notification when done
    pushoverr::pushover(message = "Finished fitting unconstrained first order model")
}

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 224
sum(out$MRA$W %*% out$alpha[iter_idx, ])
sum(out$MRA$W[, out$MRA$dims_idx == 1] %*% out$alpha[iter_idx, out$MRA$dims_idx == 1])
sum(out$MRA$W[, out$MRA$dims_idx == 2] %*% out$alpha[iter_idx, out$MRA$dims_idx == 2])
sum(out$MRA$W[, out$MRA$dims_idx == 3] %*% out$alpha[iter_idx, out$MRA$dims_idx == 3])

# Trace plots of parameters ----------------------------------------------------

plot_trace(out, file = here::here("images", "test-data", "test-data-trace-plots-no-constraint.png"))

# Trace plots for spatial random effects  --------------------------------------

plot_trace_alpha(out, file = here::here("images", "test-data", "test-data-trace-alpha-no-constraint.png"))

# Plot Posterior distributions--------------------------------------------------

plot_posterior_params(out, file = here::here("images", "test-data", "test-data-posterior-no-constraint.png"))

# Plot fitted process ----------------------------------------------------------

plot_fitted_process(out, dat_plot, file = here::here("images", "test-data", "test-data-fitted-process-no-constraint.png"))

# Predict spatial process ------------------------------------------------------

new_data <- list(
    locs_pred = cbind(code.test$Lat, code.test$Lon),
    X_pred    = model.matrix(~ Lon + Lat, data = code.test)
)

preds <- predict_mra(out, new_data)

plot_predicted_process(out, dat_plot, preds, file = here::here("images", "test-data", "test-data-predicted-process-no-constraint.png"))

# Plot fitted spatial random effects alpha -------------------------------------

plot_fitted_alphas(out, file = here::here("images", "test-data", "test-data-fitted-alpha-no-constraint.png"))

# Plot fitted MRA resolutions --------------------------------------------------

plot_fitted_MRA(out, preds, file = here::here("images", "test-data", "test-data-fitted-MRA-no-constraint.png"))












# Fit constrained model --------------------------------------------------------

if (file.exists(here::here("results", "test-fit-constrained.RData"))) {
    load(here::here("results", "test-fit-constrained.RData"))
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
        joint = TRUE,
        constraint = "resolution")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "test-fit-constrained.RData"))
    # cell phone notification when done
    pushoverr::pushover(message = "Finished fitting constrained first order model")
}

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 211
sum(out$MRA$W %*% out$alpha[iter_idx, ])
sum(out$MRA$W[, out$MRA$dims_idx == 1] %*% out$alpha[iter_idx, out$MRA$dims_idx == 1])
sum(out$MRA$W[, out$MRA$dims_idx == 2] %*% out$alpha[iter_idx, out$MRA$dims_idx == 2])
sum(out$MRA$W[, out$MRA$dims_idx == 3] %*% out$alpha[iter_idx, out$MRA$dims_idx == 3])

# Trace plots of parameters ----------------------------------------------------

plot_trace(out, file = here::here("images", "test-data", "test-data-trace-plots-constrained.png"))

# Trace plots for spatial random effects  --------------------------------------

plot_trace_alpha(out, file = here::here("images", "test-data", "test-data-trace-alpha-constrained.png"))

# Plot Posterior distributions--------------------------------------------------

plot_posterior_params(out, file = here::here("images", "test-data", "test-data-posterior-constrained.png"))

# Plot fitted process ----------------------------------------------------------

plot_fitted_process(out, dat_plot, file = here::here("images", "test-data", "test-data-fitted-process-constrained.png"))

# Predict spatial process ------------------------------------------------------

new_data <- list(
    locs_pred = cbind(code.test$Lat, code.test$Lon),
    X_pred    = model.matrix(~ Lon + Lat, data = code.test)
)

preds <- predict_mra(out, new_data)

plot_predicted_process(out, dat_plot, preds, file = here::here("images", "test-data", "test-data-predicted-process-constrained.png"))

# Plot fitted spatial random effects alpha -------------------------------------

plot_fitted_alphas(out, file = here::here("images", "test-data", "test-data-fitted-alpha-constrained.png"))

# Plot fitted MRA resolutions --------------------------------------------------

plot_fitted_MRA(out, preds, file = here::here("images", "test-data", "test-data-fitted-MRA-constrained.png"))


