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



theme_custom <- function(scale_text = 1) {
    theme(
        plot.title = element_text(
            size  = 24 * scale_text,
            face  = "bold",
            hjust = 0.5
        ),
        axis.text.x  = element_text(size = 14 * scale_text),
        axis.text.y  = element_text(size = 14 * scale_text),
        strip.text.x = element_text(size = 10 * scale_text),
        axis.title.x = element_text(
            size   = 16 * scale_text,
            margin = margin(t = 10, r = 0, b = 0, l = 0)
        ),
        axis.title.y = element_text(
            size   = 16 * scale_text,
            margin = margin(t = 0, r = 10, b = 0, l = 0)
        ),
        legend.text  = element_text(size = 10 * scale_text),
        legend.title = element_text(size = 14 * scale_text)
    )
}


plot_mask <- code.test %>%
    ggplot(aes(x = Lat, y = Lon, fill = MaskTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")
plot_full <- code.test %>%
    ggplot(aes(x = Lat, y = Lon, fill = TrueTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")
plot_mask + plot_full

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

if (file.exists(here::here("results", "test-fit-resolution-unconstrained.RData"))) {
    load(here::here("results", "test-fit-resolution-unconstrained.RData"))
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
        joint = FALSE,
        constraint = "unconstrained")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "test-fit-resolution-unconstrained.RData"))
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

# Posterior distributions ----
## predicted vs. estimated beta
dat_plot <- data.frame(
    beta = c(out$beta),
    parameter = factor(rep(1:ncol(out$beta), each = nrow(out$beta)))
)

p_beta <- dat_plot %>%
    ggplot(aes(x = parameter, y = beta, color = parameter)) +
    geom_boxplot() +
    geom_point(position = "jitter") +
    geom_hline(yintercept = 0, col = "red") +
    ggtitle("Posterior Distribution for beta") +
    theme(legend.position = "none")

## predicted vs. estimated sigma2
p_sigma2 <- data.frame(sigma2 = out$sigma2) %>%
    ggplot(aes(y = sigma2)) +
    geom_boxplot() +
    ylab("sigma2") +
    ggtitle("Posterior Distribution for sigma2")

## predicted vs. estimated tau2
dat_plot <- data.frame(
    tau2 = c(out$tau2),
    resolution = factor(rep(1:M, each = nrow(out$tau2)))
)

p_tau2 <- dat_plot %>%
    ggplot(aes(x = resolution, y = tau2)) +
    geom_boxplot() +
    geom_point(position = "jitter") +
    geom_hline(yintercept = 0, col = "red") +
    ggtitle("Posterior Distribution for tau2")



(p_beta / p_sigma2 / p_tau2)

# Plot fitted process

Xbeta_post <- t(X %*% t(out$beta)) * sd_y + mu_y
Walpha_post <- t(out$MRA$W %*% t(out$alpha)) * sd_y
mu_post <- Xbeta_post + Walpha_post

dat_plot <- data.frame(
    Lat          = locs[, 1],
    Lon          = locs[, 2],
    mean_Xbeta   = apply(Xbeta_post, 2, mean),
    mean_Walpha  = apply(Walpha_post, 2, mean),
    mean_mu      = apply(mu_post, 2, mean))
zlims = range(
    range(dat_plot$mean_mu),
    range(code.test$TrueTemp, na.rm = TRUE))

plot_Xbeta <- dat_plot %>%
    ggplot(aes(x = Lat, y = Lon, fill = mean_Xbeta)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")

plot_Walpha <- dat_plot %>%
    ggplot(aes(x = Lat, y = Lon, fill = mean_Walpha)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")

plot_mu <- dat_plot %>%
    ggplot(aes(x = Lat, y = Lon, fill = mean_mu)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)

plot_mask <- code.test %>%
    ggplot(aes(x = Lat, y = Lon, fill = MaskTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)

(plot_Xbeta + plot_Walpha) / (plot_mu + plot_mask)

# Predict spatial process -----------------------------------------------------
locs_pred <- cbind(code.test$Lat, code.test$Lon)
MRA_pred <- mra_wendland_2d_pred(locs_pred, out$MRA)

# y_pred <- matrix(NA, nrow = nrow(code.test), ncol = nrow(out$alpha))
X_pred <- model.matrix(~ Lon + Lat, data = code.test)
W_pred <- MRA_pred$W_pred

Xbeta_pred  <- t(X_pred %*% t(out$beta)) * sd_y + mu_y
Walpha_pred <- t(W_pred %*% t(out$alpha)) * sd_y
y_pred <- Xbeta_pred + Walpha_pred


na_idx <- is.na(code.test$MaskTemp)

Wpred_mean <- apply(Walpha_pred, 2, mean)
Wpost_mean <- apply(Walpha_post, 2, mean)
plot(Wpred_mean[!na_idx], Wpost_mean)
range(Walpha_post)
range(Walpha_pred)

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

fitted_alphas <- data.frame(
    x = unlist(sapply(1:M, function(i) out$MRA$locs_grid[[i]][, 1])),
    y = unlist(sapply(1:M, function(i) out$MRA$locs_grid[[i]][, 2])),
    res = factor(unlist(sapply(1:M, function(i) rep(i, each = nrow(out$MRA$locs_grid[[i]]))))),
    alpha = unlist(sapply(1:M, function(i) apply(out$alpha, 2, mean)[out$MRA$dims_idx == i]))
) %>%
    ggplot(aes(x = x, y = y, fill = alpha)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Posterior mean spatial random effects") +
    facet_wrap( ~ res, ncol = 2) +
    geom_rect(
        data = NULL,
        aes(xmin = min(locs[, 1]),
            xmax = max(locs[, 1]),
            ymin = min(locs[, 2]),
            ymax = max(locs[, 2])
        ),
        fill  = NA,
        color = "black",
        alpha = 0.5
    ) +
    theme_bw() +
    theme_custom(0.5)
fitted_alphas

# Plot fitted MRA resolutions -------------------------------------------------

W_alpha_res = unlist(sapply(1:M, function(m) W_pred[, out$MRA$dims_idx == m] %*% apply(out$alpha, 2, mean)[out$MRA$dims_idx == m], simplify = "matrix"))
dimnames(W_alpha_res) <- list(
    site = 1:nrow(W_alpha_res),
    res  = paste("resolution", 1:M)
)
dat_locs <- data.frame(
    site = 1:nrow(W_alpha_res),
    Lat  = code.test$Lat,
    Lon  = code.test$Lon
)

estimated_MRA <- merge(
    dat_locs,
    as.data.frame.table(W_alpha_res, responseName = "W_alpha")
) %>%
    ggplot(aes(x = Lon, y = Lat, fill = W_alpha)) +
    geom_raster() +
    facet_wrap( ~ res, ncol = 2) +
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle("Estimated Multiresolution process") +
    scale_fill_viridis_c(option = "B") +
    theme_bw() +
    theme_custom(0.5)

estimated_MRA















# Fit constrained model -------------------------------------------------------


if (file.exists(here::here("results", "test-fit-resolution-constrained.RData"))) {
    load(here::here("results", "test-fit-resolution-constrained.RData"))
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
        joint = FALSE,
        constraint = "resolution")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "test-fit-resolution-constrained.RData"))
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

# Posterior distributions ----
## predicted vs. estimated beta


p_beta <- data.frame(
    beta = c(out$beta),
    parameter = factor(rep(1:ncol(out$beta), each = nrow(out$beta)))
) %>%
    ggplot(aes(x = parameter, y = beta)) +
    geom_boxplot() +
    geom_point(position = "jitter", alpha = 0.125) +
    geom_hline(yintercept = 0, col = "red") +
    ggtitle("Posterior Distribution for beta") +
    theme(legend.position = "none")


## predicted vs. estimated sigma2
p_sigma2 <- data.frame(sigma2 = out$sigma2) %>%
    ggplot(aes(y = sigma2)) +
    geom_boxplot() +
    ylab("sigma2") +
    ggtitle("Posterior Distribution for sigma2")

## predicted vs. estimated tau2

p_tau2 <- data.frame(
    tau2 = c(out$tau2),
    resolution = factor(rep(1:M, each = nrow(out$tau2)))
) %>%
    ggplot(aes(x = resolution, y = tau2)) +
    geom_boxplot() +
    geom_point(position = "jitter", alpha = 0.125) +
    geom_hline(yintercept = 0, col = "red") +
    ggtitle("Posterior Distribution for tau2")

(p_beta + p_sigma2 + p_tau2)

# Plot fitted process

Xbeta_post <- t(X %*% t(out$beta)) * sd_y + mu_y
Walpha_post <- t(out$MRA$W %*% t(out$alpha)) * sd_y
mu_post <- Xbeta_post + Walpha_post

dat_plot <- data.frame(
    Lat          = locs[, 1],
    Lon          = locs[, 2],
    mean_Xbeta   = apply(Xbeta_post, 2, mean),
    mean_Walpha  = apply(Walpha_post, 2, mean),
    mean_mu      = apply(mu_post, 2, mean))
zlims = range(
    range(dat_plot$mean_mu),
    range(code.test$TrueTemp, na.rm = TRUE))

plot_Xbeta <- dat_plot %>%
    ggplot(aes(x = Lat, y = Lon, fill = mean_Xbeta)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")

plot_Walpha <- dat_plot %>%
    ggplot(aes(x = Lat, y = Lon, fill = mean_Walpha)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")

plot_mu <- dat_plot %>%
    ggplot(aes(x = Lat, y = Lon, fill = mean_mu)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)

plot_mask <- code.test %>%
    ggplot(aes(x = Lat, y = Lon, fill = MaskTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)

(plot_Xbeta + plot_Walpha) / (plot_mu + plot_mask)

# Predict spatial process -----------------------------------------------------
locs_pred <- cbind(code.test$Lat, code.test$Lon)
MRA_pred <- mra_wendland_2d_pred(locs_pred, out$MRA)

# y_pred <- matrix(NA, nrow = nrow(code.test), ncol = nrow(out$alpha))
X_pred <- model.matrix(~ Lon + Lat, data = code.test)
W_pred <- MRA_pred$W_pred

Xbeta_pred  <- t(X_pred %*% t(out$beta)) * sd_y + mu_y
Walpha_pred <- t(W_pred %*% t(out$alpha)) * sd_y
y_pred <- Xbeta_pred + Walpha_pred


na_idx <- is.na(code.test$MaskTemp)

Wpred_mean <- apply(Walpha_pred, 2, mean)
Wpost_mean <- apply(Walpha_post, 2, mean)
plot(Wpred_mean[!na_idx], Wpost_mean)
range(Walpha_post)
range(Walpha_pred)


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

fitted_alphas <- data.frame(
    x = unlist(sapply(1:M, function(i) out$MRA$locs_grid[[i]][, 1])),
    y = unlist(sapply(1:M, function(i) out$MRA$locs_grid[[i]][, 2])),
    res = factor(unlist(sapply(1:M, function(i) rep(i, each = nrow(out$MRA$locs_grid[[i]]))))),
    alpha = unlist(sapply(1:M, function(i) apply(out$alpha, 2, mean)[out$MRA$dims_idx == i]))
) %>%
    ggplot(aes(x = x, y = y, fill = alpha)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Posterior mean spatial random effects") +
    facet_wrap( ~ res, ncol = 2) +
    geom_rect(
        data = NULL,
        aes(xmin = min(locs[, 1]),
            xmax = max(locs[, 1]),
            ymin = min(locs[, 2]),
            ymax = max(locs[, 2])
        ),
        fill  = NA,
        color = "black",
        alpha = 0.5
    ) +
    theme_bw() +
    theme_custom(0.5)
fitted_alphas

# Plot fitted MRA resolutions -------------------------------------------------

W_alpha_res = unlist(sapply(1:M, function(m) W_pred[, out$MRA$dims_idx == m] %*% apply(out$alpha, 2, mean)[out$MRA$dims_idx == m], simplify = "matrix"))
dimnames(W_alpha_res) <- list(
    site = 1:nrow(W_alpha_res),
    res  = paste("resolution", 1:M)
)
dat_locs <- data.frame(
    site = 1:nrow(W_alpha_res),
    Lat  = code.test$Lat,
    Lon  = code.test$Lon
)

estimated_MRA <- merge(
    dat_locs,
    as.data.frame.table(W_alpha_res, responseName = "W_alpha")
) %>%
    ggplot(aes(x = Lon, y = Lat, fill = W_alpha)) +
    geom_raster() +
    facet_wrap( ~ res, ncol = 2) +
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle("Estimated Multiresolution process") +
    scale_fill_viridis_c(option = "B") +
    theme_bw() +
    theme_custom(0.5)

estimated_MRA
