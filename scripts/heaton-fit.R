library(spam)
library(tidyverse)
library(patchwork)
library(BayesMRA)
library(scoringRules)

load(here::here("data", "AllSatelliteTemps.RData"))

locs <- cbind(all.sat.temps$Lat, all.sat.temps$Lon)

# Plot the data ----------------------------------------------------------------

dat_plot <- data.frame(
    Lat = all.sat.temps$Lat,
    Lon = all.sat.temps$Lon,
    Observed = all.sat.temps$MaskTemp,
    Truth = all.sat.temps$TrueTemp)

plot_test_data(dat_plot)

# Setup the model --------------------------------------------------------------
dat_fit <- all.sat.temps %>%
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


# Fit model -------------------------------------------------------------------

if (file.exists(here::here("results", "heaton-fit.RData"))) {
    load(here::here("results", "heaton-fit.RData"))
} else {
    start   <- Sys.time()
    # to do: center and scale X
    out <- mcmc_mra(
        y             = y,
        X             = X,
        locs          = locs,
        params        = params,
        priors        = priors,
        M             = 3,
        n_coarse_grid = 40,
        joint         = FALSE,
        constraint    = "resolution")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "heaton-fit.RData"))
    # cell phone notification when done
    pushoverr::pushover(message = "Finished fitting Heaton model")
}

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 224
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
    locs_pred = cbind(all.sat.temps$Lat, all.sat.temps$Lon),
    X_pred    = model.matrix(~ Lon + Lat, data = all.sat.temps)
)

preds <- predict_mra(out, new_data)

plot_predicted_process(out, dat_plot, preds)

# Plot fitted spatial random effects alpha -------------------------------------

plot_fitted_alphas(out)

# Plot fitted MRA resolutions --------------------------------------------------

plot_fitted_MRA(out, preds)

## calculate RMSE
oos <- (!is.na(all.sat.temps$TrueTemp) & is.na(all.sat.temps$MaskTemp))
RMSE <- sqrt(mean((preds$y_pred_mean[oos] - all.sat.temps$TrueTemp[oos])^2))








# Fit finer grid model ---------------------------------------------------------

if (file.exists(here::here("results", "heaton-fit-fine-grid.RData"))) {
    load(here::here("results", "heaton-fit-fine-grid.RData"))
} else {
    start   <- Sys.time()
    # to do: center and scale X
    out <- mcmc_mra(
        y             = y,
        X             = X,
        locs          = locs,
        params        = params,
        priors        = priors,
        M             = 4,
        n_coarse_grid = 40,
        joint         = FALSE,
        constraint    = "resolution")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "heaton-fit-fine-grid.RData"))
    # cell phone notification when done
    pushoverr::pushover(message = "Finished fitting Heaton fine-grid model")
}

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 224
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
    locs_pred = cbind(all.sat.temps$Lat, all.sat.temps$Lon),
    X_pred    = model.matrix(~ Lon + Lat, data = all.sat.temps)
)

preds <- predict_mra(out, new_data)

plot_predicted_process(out, dat_plot, preds)

# Plot fitted spatial random effects alpha -------------------------------------

plot_fitted_alphas(out)

# Plot fitted MRA resolutions --------------------------------------------------

plot_fitted_MRA(out, preds)

## runtime
as.numeric(runtime, units = "mins")

## calculate model metrics
oos <- (!is.na(all.sat.temps$TrueTemp) & is.na(all.sat.temps$MaskTemp))
y_oos <- all.sat.temps$TrueTemp[oos]
# RMSE
RMSE <- sqrt(mean((preds$y_pred_mean[oos] - y_oos)^2))
#MAE
MAE  <- mean(abs(preds$y_pred_mean[oos] - y_oos))
# CRPS
CRPS <- BayesComposition::makeCRPS(preds$y_pred[, oos], y_oos, nrow(preds$y_pred[, oos]))
# Cov
CI <- apply(preds$y_pred[, oos], 2, quantile, prob = c(0.025, 0.975))
in_out <- rep(FALSE, sum(oos))
for (i in 1:sum(oos)) {
    in_out[i] <- (y_oos[i] > CI[1, i]) & (y_oos[i] < CI[2, i])
}
COV <- mean(in_out)
# Int

RMSE
MAE
COV
mean(CRPS)
# Int



## calculate predictive score by distance to observation -- seems strange
locs <- cbind(all.sat.temps$Lat, all.sat.temps$Lon)
in_samp <- !is.na(all.sat.temps$MaskTemp)
locs_oos <- locs[oos, ]
locs_is <- locs[in_samp, ]
D_min <- rep(0, sum(oos))
for (i in 1:sum(oos)) {
    if (i %% 500 == 0) {
        message("On iteration ", i, " out of ", sum(oos))
    }
    D_min[i] <- min(fields::rdist(locs_oos[i, , drop = FALSE], locs_is))
}
nbreaks <- 5
dat_scores <- data.frame(
    RMSE = (preds$y_pred_mean[oos] - y_oos)^2,
    MAE  = abs(preds$y_pred_mean[oos] - y_oos),
    D_min = cut(D_min, breaks = nbreaks)
)

dat_scores %>%
    group_by(D_min) %>%
    summarize(RMSE = sqrt(mean(RMSE)), MAE = mean(MAE)) %>%
    ggplot(aes(x = D_min, y = RMSE, group = 1)) +
    geom_line()
dat_scores %>%
    group_by(D_min) %>%
    summarize(RMSE = sqrt(mean(RMSE)), MAE = mean(MAE)) %>%
    ggplot(aes(x = D_min, y = RMSE, group = 1)) +
    geom_line() +

##  try to plot MSE/MAE on the map...
    geom_point(data = dat_scores, alpha = 0.01)
