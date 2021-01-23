library(spam)
library(tidyverse)
library(patchwork)
library(BayesMRA)

load(here::here("data", "SmallTestData.RData"))

str(code.test)
code.test$MaskTemp = code.test$MaskedData
code.test$FullTemp = code.test$FullData
code.test$TrueTemp = code.test$FullData

locs <- cbind(code.test$Lat, code.test$Lon)

## create a mask in the upper corner
code.test$MaskTemp[code.test$Lat > 0.60 & code.test$Lon > 0.60] <- NA


# Plot the data ---------------------------------------------------------------

dat_plot <- data.frame(
    Lat = code.test$Lat,
    Lon = code.test$Lon,
    Observed = code.test$MaskTemp,
    Truth = code.test$FullTemp)

plot_test_data(dat_plot)


# Setup the model -------------------------------------------------------------
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

M <- 3
n_coarse_grid <- 10

sd_y <- sd(y)
mu_y <- mean(y)

MRA <- mra_wendland_2d(locs, M = 3, n_coarse_grid = 20)
MRA$n_dims
sqrt(MRA$n_dims)
sum(MRA$n_dims)
nrow(code.test)

# Fit unconstrained model -----------------------------------------------------

if (file.exists(here::here("results", "test-fit-unconstrained.RData"))) {
    load(here::here("results", "test-fit-unconstrained.RData"))
} else {
    start   <- Sys.time()

    # to do: center and scale X
    out <- mcmc_mra(
        y = (y - mu_y) / sd_y,
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
}

# pushoverr::pushover(message = "Finished fitting unconstrained first order model")

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

# Plot fitted process

plot_fitted_process(out, dat_plot)

# Predict spatial process -----------------------------------------------------

new_data <- list(
    locs_pred = cbind(code.test$Lat, code.test$Lon),
    X_pred    = model.matrix(~ Lon + Lat, data = code.test)
)

preds <- predict_mra(out, new_data)


# MRA_pred <- mra_wendland_2d_pred(locs_pred, out$MRA)
#
# # y_pred <- matrix(NA, nrow = nrow(code.test), ncol = nrow(out$alpha))
# X_pred <- model.matrix(~ Lon + Lat, data = code.test)
# W_pred <- MRA_pred$W_pred
#
# Xbeta_pred  <- t(X_pred %*% t(out$beta)) * sd_y + mu_y
# Walpha_pred <- t(W_pred %*% t(out$alpha)) * sd_y
# y_pred <- Xbeta_pred + Walpha_pred
#
#
# na_idx <- is.na(code.test$MaskTemp)
#
# Wpred_mean <- apply(Walpha_pred, 2, mean)
# Wpost_mean <- apply(Walpha_post, 2, mean)
# plot(Wpred_mean[!na_idx], Wpost_mean)
# range(Walpha_post)
# range(Walpha_pred)
#
# y_pred_mean <- apply(y_pred, 2, mean)
# y_pred_sd   <- apply(y_pred, 2, sd)

dat_pred <- data.frame(
    y = y_pred_mean,
    sd = y_pred_sd,
    Lat = code.test$Lat,
    Lon = code.test$Lon)

zlims = range(
    range(y_pred_mean),
    range(code.test$TrueTemp, na.rm = TRUE))

p_obs <- code.test %>%
    ggplot(aes(x = Lat, y = Lon, fill = TrueTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)
p_fit <- dat_pred %>%
    ggplot(aes(x = Lat, y = Lon, fill = y)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)
p_obs + p_fit


p_sd <- dat_pred %>%
    ggplot(aes(x = Lat, y = Lon, fill = sd)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")
p_sd


# Plot fitted spatial random effects alpha ------------------------------------

plot_fitted_alphas(out)


# Plot fitted MRA resolutions -------------------------------------------------

plot_fitted_MRA(out, preds)










# Fit constrained model -------------------------------------------------------


if (file.exists(here::here("results", "test-fit-constrained.RData"))) {
    load(here::here("results", "test-fit-constrained.RData"))
} else {
    start   <- Sys.time()

    # to do: center and scale X
    out <- mcmc_mra(
        y = (y - mu_y) / sd_y,
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
}

# pushoverr::pushover(message = "Finished fitting constrained first order model")

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 211
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

# Plot fitted process

plot_fitted_process(out, dat_plot)

# Predict spatial process -----------------------------------------------------

new_data <- list(
    locs_pred = cbind(code.test$Lat, code.test$Lon),
    X_pred    = model.matrix(~ Lon + Lat, data = code.test)
)

preds <- predict_mra(out, new_data)

# locs_pred <- cbind(code.test$Lat, code.test$Lon)
# MRA_pred <- mra_wendland_2d_pred(locs_pred, out$MRA)
#
# # y_pred <- matrix(NA, nrow = nrow(code.test), ncol = nrow(out$alpha))
# X_pred <- model.matrix(~ Lon + Lat, data = code.test)
# W_pred <- MRA_pred$W_pred
#
# Xbeta_pred  <- t(X_pred %*% t(out$beta)) * sd_y + mu_y
# Walpha_pred <- t(W_pred %*% t(out$alpha)) * sd_y
# y_pred <- Xbeta_pred + Walpha_pred
#
#
# na_idx <- is.na(code.test$MaskTemp)
#
# Wpred_mean <- apply(Walpha_pred, 2, mean)
# Wpost_mean <- apply(Walpha_post, 2, mean)
# plot(Wpred_mean[!na_idx], Wpost_mean)
# range(Walpha_post)
# range(Walpha_pred)


y_pred_mean <- apply(y_pred, 2, mean)
y_pred_sd   <- apply(y_pred, 2, sd)

dat_pred <- data.frame(
    y = y_pred_mean,
    sd = y_pred_sd,
    Lat = code.test$Lat,
    Lon = code.test$Lon)


zlims = range(
    range(y_pred_mean),
    range(code.test$TrueTemp, na.rm = TRUE))

p_obs <- code.test %>%
    ggplot(aes(x = Lat, y = Lon, fill = TrueTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)
p_fit <- dat_pred %>%
    ggplot(aes(x = Lat, y = Lon, fill = y)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)
p_obs + p_fit


p_sd <- dat_pred %>%
    ggplot(aes(x = Lat, y = Lon, fill = sd)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")
p_sd


# Plot fitted spatial random effects alpha ------------------------------------

plot_fitted_alphas(out)

# Plot fitted MRA resolutions -------------------------------------------------

plot_fitted_MRA(out, preds)


