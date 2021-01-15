load("~/Downloads/heatoncomparison-master/Data/AllSatelliteTemps.RData")

str(all.sat.temps)
locs <- cbind(all.sat.temps$Lat, all.sat.temps$Lon)


library(tidyverse)
library(patchwork)

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


plot_mask <- all.sat.temps %>%
    ggplot(aes(x = Lat, y = Lon, fill = MaskTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")
plot_full <- all.sat.temps %>%
    ggplot(aes(x = Lat, y = Lon, fill = TrueTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")
plot_mask + plot_full

dat_fit <- all.sat.temps %>%
    filter(!is.na(MaskTemp))
library(BayesMRA)
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


if (file.exists(here::here("results", "heaton-fit.RData"))) {
    load(here::here("results", "heaton-fit.RData"))
} else {
    start   <- Sys.time()
    out <- mcmc_mra(
    y = (y - mu_y) / sd_y,
    X = X,
    locs = locs,
    params = params,
    priors = priors,
    M = 3,
    n_coarse_grid = 10)
    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "heaton-fit.RData"))
}

# Trace plots of parameters ---------------------------------------------------

dat_tau2 <- data.frame(
    tau2      = c(out$tau2),
    iteration = rep(1:nrow(out$tau2), times = M),
    resolution = factor(rep(1:M, each = nrow(out$tau2)))
)

p_tau2 <- ggplot(dat_tau2, aes(x = iteration, y = tau2, group = resolution, color = resolution)) +
    geom_line() +
    scale_color_viridis_d(end = 0.8) +
    ggtitle("Trace plots for tau2") +
    theme_bw() +
    theme_custom(0.7)

dat_lambda <- data.frame(
    lambda      = c(out$lambda),
    iteration = rep(1:nrow(out$lambda), times = M),
    resolution = factor(rep(1:M, each = nrow(out$lambda)))
)

p_lambda <- ggplot(dat_lambda, aes(x = iteration, y = lambda, group = resolution, color = resolution)) +
    geom_line() +
    scale_color_viridis_d(end = 0.8) +
    ggtitle("Trace plots for lambda") +
    theme_bw() +
    theme_custom(0.7)


p_sigma2 <- data.frame(sigma2 = c(out$sigma2), iteration = 1:length(out$sigma2)) %>%
    ggplot(aes(x = iteration, y = sigma2)) +
    geom_line() +
    scale_color_viridis_d(end = 0.8) +
    ggtitle("Trace plots for sigma2") +
    theme_bw() +
    theme_custom(0.7)

dat_beta <- data.frame(
    beta      = c(out$beta),
    iteration = rep(1:nrow(out$beta), times = ncol(X)),
    covariate = factor(rep(1:ncol(X), each = nrow(out$beta)))
)
p_beta  <- ggplot(dat_beta, aes(x = iteration, y = beta, group = covariate, color = covariate)) +
    geom_line() +
    scale_color_viridis_d(end = 0.8) +
    ggtitle("Trace plots for beta") +
    theme_bw() +
    theme_custom(0.7)

(p_tau2 + p_lambda) /  (p_sigma2 + p_beta)

# Trace plots for spatial random effects  -------------------------------------
dat_alpha <- data.frame(
    alpha = c(out$alpha),
    iteration = rep(1:nrow(out$alpha), each = ncol(out$alpha)),
    knot      = factor(unlist(sapply(1:M, function(m) 1:out$MRA$n_dims[m]))),
    resolution = factor(paste("resolution", out$MRA$dims_idx))
)

plot_idx <- sapply(1:M, function(m) sample(out$MRA$n_dims[m], 10))

dat_alpha %>%
    subset(
        (resolution == "resolution 1" & knot %in% plot_idx[, 1]) |
            (resolution == "resolution 2" & knot %in% plot_idx[, 2]) |
            (resolution == "resolution 3" & knot %in% plot_idx[, 3])
    ) %>%
    ggplot(aes(x = iteration, y = alpha, group = knot, color = knot)) +
    geom_line(alpha = 0.5) +
    facet_wrap( ~ resolution, ncol = 1) +
    scale_color_viridis_d(end = 0.8) +
    ggtitle("Trace plots for alpha") +
    theme_bw() +
    theme_custom(0.7) +
    theme(legend.position = "none")


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
    tau2 = c(out$beta),
    resolution = factor(rep(1:M, each = nrow(out$tau2)))
)

p_tau2 <- dat_plot %>%
    ggplot(aes(x = resolution, y = tau2)) +
    geom_boxplot() +
    geom_point(position = "jitter") +
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
    range(all.sat.temps$TrueTemp, na.rm = TRUE))

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

plot_mask <- all.sat.temps %>%
    ggplot(aes(x = Lat, y = Lon, fill = MaskTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B", limits = zlims)

(plot_Xbeta + plot_Walpha) / (plot_mu + plot_mask)

# Predict spatial process -----------------------------------------------------
locs_pred <- cbind(all.sat.temps$Lat, all.sat.temps$Lon)
MRA_pred <- mra_wendland_2d_pred(locs_pred, out$MRA)

# y_pred <- matrix(NA, nrow = nrow(all.sat.temps), ncol = nrow(out$alpha))
X_pred <- model.matrix(~ Lon + Lat, data = all.sat.temps)
W_pred <- do.call(cbind, MRA_pred$W_pred)
# for (k in 1:nrow(out$alpha)) {
#     if (k %% 10 == 0) {
#         message("Generating prediction for iteration ", k, " out of ", nrow(out$alpha), " iterations")
#     }
#     Xbeta_pred <- X_pred %*% out$beta[k, ] * sd_y + mu_y
#     Walpha_pred <- W_pred %*% out$alpha[k, ] * sd_y
#     y_pred[, k] <- Xbeta_pred + Walpha_pred
# }

Xbeta_pred  <- t(X_pred %*% t(out$beta)) * sd_y + mu_y
Walpha_pred <- t(W_pred %*% t(out$alpha)) * sd_y
y_pred <- Xbeta_pred + Walpha_pred

Walpha_post <- t(out$MRA$W %*% t(out$alpha)) * sd_y

y_pred_mean <- apply(y_pred, 2, mean)
y_pred_sd   <- apply(y_pred, 2, sd)

dat_pred <- data.frame(
    y = y_pred_mean,
    sd = y_pred_sd,
    Lat = locs_pred[, 1],
    Lon = locs_pred[, 2])

p_obs <- all.sat.temps %>%
    ggplot(aes(x = Lat, y = Lon, fill = TrueTemp)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")
p_fit <- dat_pred %>%
    ggplot(aes(x = Lat, y = Lon, fill = y)) +
    geom_raster() +
    scale_fill_viridis_c(option = "B")
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

# Plot fitted MRA resoultions -------------------------------------------------

W_alpha_res = unlist(sapply(1:M, function(m) W_pred[, out$MRA$dims_idx == m] %*% apply(out$alpha, 2, mean)[out$MRA$dims_idx == m], simplify = "matrix"))
dimnames(W_alpha_res) <- list(
    site = 1:nrow(W_alpha_res),
    res  = paste("resolution", 1:M)
)
dat_locs <- data.frame(
    site = 1:nrow(W_alpha_res),
    Lat  = all.sat.temps$Lat,
    Lon  = all.sat.temps$Lon
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
