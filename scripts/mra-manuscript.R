library(tidyverse)
library(BayesMRA)
library(mvnfast)
library(tidyverse)
library(patchwork)
library(spam64)
library(fields)
library(spBayes)
library(coda)
library(here)

res <- 100

# setup the project folder directories
if (!dir.exists(here::here("results"))) {
    dir.create(here::here("results"))
}
if (!dir.exists(here::here("images"))) {
    dir.create(here::here("images"))
}

# Simulate data ---------------------------------------------------------------

set.seed(11)

N <- 200^2
## setup the spatial process
locs <- as.matrix(
    expand.grid(
        seq(0, 1, length.out = sqrt(N)),
        seq(0, 1, length.out = sqrt(N))
    )
)
# D <- fields::rdist(locs)

## fixed effects include intercept, elevation, and latitude
# X <- cbind(1, as.vector(mvnfast::rmvn(1, rep(0, N), 3 * exp(-D / 20))), locs[, 2])
X <- cbind(1, rnorm(N))
p <- ncol(X)

beta <- c(2, -1.5)

## MRA spatio-temporal random effect
M <- 4
n_coarse <- 20

MRA    <- mra_wendland_2d(locs, M = M, n_coarse = n_coarse, use_spam = TRUE)

# MRA    <- mra_wendland_2d(locs, M = M, n_max_fine_grid = 2^8, use_spam = TRUE)
W        <- MRA$W
n_dims   <- MRA$n_dims
dims_idx <- MRA$dims_idx

Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(0.999, length(n_dims)), prec_model = "CAR")

tau2         <- 10 * 2^(2 * (1:M - 1))
Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2)

## initialize the random effect
## constraint on each resolution of Walpha
constraints <- make_constraint(MRA, constraint = "resolution", joint = TRUE)
A_constraint <- constraints$A_constraint
a_constraint <- constraints$a_constraint
alpha   <- as.vector(rmvnorm.prec.const(n = 1, mu = rep(0, sum(n_dims)), Q = Q_alpha_tau2, A = A_constraint, a = a_constraint))


## verify that the parameters sum to 0
sum(W %*% alpha)
sum(W[, dims_idx == 1] %*% alpha[dims_idx == 1])
sum(W[, dims_idx == 2] %*% alpha[dims_idx == 2])
sum(W[, dims_idx == 3] %*% alpha[dims_idx == 3])

sigma2 <- runif(1, 0.25, 0.5)

y <- X %*% beta + W %*% alpha + rnorm(N, 0, sqrt(sigma2))


# Plot the simulated data -----------------------------------------------------


#
# if (!file.exists(here::here("images", "MRA-grid-neighbors.png"))) {
#     png(file = here::here("images", "MRA-grid-neighbors.png"), width = 16, height = 9, units = "in", res = res)
#
#     nnn <- sapply(1:M, function(i) {
#         D <- rdist(MRA$locs_grid[[i]])
#         sapply(1:nrow(MRA$locs_grid[[i]]), function(j) {
#             sum((D[j, ]) < MRA$radius[i])}
#         )
#     })
#
#     p <- data.frame(
#         x = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 1])),
#         y = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 2])),
#         resolution = factor(unlist(sapply(1:M, function(i) rep(paste("resolution", i), each = nrow(MRA$locs_grid[[i]]))))),
#         neighbors <- unlist(nnn)
#     )  %>%
#         ggplot(aes(x = x, y = y, fill = neighbors)) +
#         geom_raster() +
#         facet_wrap(~ resolution) +
#         ggtitle("number of neighbors") +
#         scale_fill_viridis_c() +
#         geom_rect(
#             data = NULL,
#             aes(xmin = min(locs[, 1]),
#                 xmax = max(locs[, 1]),
#                 ymin = min(locs[, 2]),
#                 ymax = max(locs[, 2])
#             ),
#             fill  = NA,
#             color = "black",
#             alpha = 0.5
#         ) +
#         theme_bw() +
#         theme_custom()
#
#     print(p)
#     dev.off()
# }
#
#
# ## plot the random effects
# if (!file.exists(here::here("images", "MRA-random-effects.png"))) {
#     png(file = here::here("images", "MRA-random-effects.png"), width = 16, height = 9, units = "in", res = res)
#
    sim_alphas <- data.frame(
        x = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 1])),
        y = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 2])),
        resolution = factor(unlist(sapply(1:M, function(i) rep(paste("resolution", i), each = nrow(MRA$locs_grid[[i]]))))),
        alpha = unlist(sapply(1:M, function(i) alpha[dims_idx == i]))
    ) %>%
        ggplot(aes(x = x, y = y, fill = alpha)) +
        geom_raster() +
        ggtitle("MRA spatial random effects") +
        scale_fill_viridis_c() +
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
        facet_wrap( ~ resolution, ncol = 2) +
        theme_bw()

    plot(sim_alphas)
#     dev.off()
# }
#
#
# ## plot the simulated data
# if (!file.exists(here::here("images", "MRA-simulated-data.png"))) {
#     png(file = here::here("images", "MRA-simulated-data.png"), width = 16, height = 9, units = "in", res = res)
#
#     zlims <- c(
#         min(c(X %*% beta), c(y), c(W %*% alpha)),
#         max(c(X %*% beta), c(y), c(W %*% alpha))
#     )
#
#
#     g_fixed <- data.frame(
#         temp = c(X %*% beta),
#         x    = locs[, 1],
#         y    = locs[, 2]
#     ) %>%
#         ggplot(aes(x = x, y = y, fill = temp)) +
#         geom_raster() +
#         xlab("x") +
#         ylab("y") +
#         ggtitle("Fixed effects") +
#         colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = zlims) +
#         theme_bw() +
#         theme_custom(0.70)
#
#     dat_map <- data.frame(
#         temp = c(y),
#         re   = c(W %*% alpha),
#         x    = locs[, 1],
#         y    = locs[, 2]
#     )
#
#     g_full <- dat_map %>%
#         ggplot(aes(x = x, y = y, fill = temp)) +
#         geom_raster() +
#         xlab("x") +
#         ylab("y") +
#         ggtitle("Observed surface") +
#         colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = zlims) +
#         theme_bw() +
#         theme_custom(0.70)
#
#     g_re <- dat_map %>%
#         ggplot(aes(x = x, y = y, fill = re)) +
#         geom_raster() +
#         xlab("x") +
#         ylab("y") +
#         ggtitle("Spatial process") +
#         colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = zlims) +
#         theme_bw() +
#         theme(legend.position = "none") +
#         theme_custom(0.70)
#
#     print(g_fixed + g_re + g_full + plot_layout(guides = "collect"))
#     dev.off()
# }
#
#
#
# # plot the MRA spatial process
# if (!file.exists(here::here("images", "MRA-resolutions.png"))) {
#     png(file = here::here("images", "MRA-resolutions.png"), width = 16, height = 9, units = "in", res = res)
#
#     W_alpha_res = unlist(sapply(1:M, function(m) W[, dims_idx == m] %*% alpha[dims_idx == m], simplify = "matrix"))
#     dimnames(W_alpha_res) <- list(
#         site = 1:N,
#         res  = paste("resolution", 1:M)
#     )
#     dat_locs <- data.frame(
#         site = 1:N,
#         x    = locs[, 1],
#         y    = locs[, 2]
#     )
#
#     sim_MRA <- rbind(
#         merge(dat_locs, as.data.frame.table(W_alpha_res, responseName = "W_alpha")),
#         data.frame(
#             site    = 1:N,
#             x       = locs[, 1],
#             y       = locs[, 2],
#             W_alpha = W %*% alpha,
#             res     = "full process"
#         )
#     ) %>%
#         ggplot(aes(x = x, y = y, fill = W_alpha)) +
#         geom_raster() +
#         facet_wrap( ~ res, ncol = 2) +
#         xlab("x") +
#         ylab("y") +
#         ggtitle("Multi-resolution process") +
#         colorspace::scale_fill_continuous_diverging("Blue-Red 3") +
#         theme_bw() +
#         theme_custom()
#     print(sim_MRA)
#     dev.off()
# }


# Sample the simulated data --------------------------------------------------

set.seed(111)
n_sites   <- 5000
s         <- sample(N, n_sites)
y_s       <- y[s, ]
y_oos     <- y[ - s, ]
X_s       <- X[s, ]
X_oos     <- X[ - s, ]
locs_s    <- locs[s, ]
locs_oos  <- locs[ - s, ]

y_obs <- y
y_obs[-s] <- NA
dat_plot <- data.frame(
    Lat = locs[, 1],
    Lon = locs[, 2],
    Observed = y_obs,
    Truth = y)

plot_test_data(dat_plot)


if (!file.exists(here::here("images", "sample-data.png"))) {
    png(file = here::here("images", "sample-data.png"), width = 16, height = 9, units = "in", res = res)
    zlims = range(y)
    p <- data.frame(
        x     = locs_s[, 1],
        y     = locs_s[, 2],
        site  = factor(1:n_sites),
        y_s     = y_s
    ) %>%
        ggplot(aes(x = x, y = y, fill = y_s)) +
        geom_raster() +
        colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = zlims) +
        ggtitle("Sampled y") +
        theme_bw() +
        theme(legend.position = "none")
    print(p)
    dev.off()
}



# Fitting the model -----------------------------------------------------------

params <- list(
    n_mcmc    = 500,
    n_adapt   = 500,
    n_thin    = 1,
    n_message = 50
)

priors <- list(
    alpha_tau2   = 1,
    beta_tau2    = 1,
    alpha_sigma2 = 1,
    beta_sigma2  = 1,
    mu_beta      = rep(0, ncol(X_s)),
    Sigma_beta   = 100 * diag(ncol(X_s)))

set.seed(999)

if (file.exists(here::here("results", "fit-mra-sim.RData"))) {
    load(here::here("results", "fit-mra-sim.RData"))
} else {
    start   <- Sys.time()
    out     <- mcmc_mra(
        y             = y_s,
        X             = X_s,
        locs          = locs_s,
        params        = params,
        priors        = priors,
        M             = M,
        n_coarse_grid = n_coarse,
        joint         = FALSE,
        constraint    = "resolution")

    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "fit-mra-sim.RData"))
    pushoverr::pushover(message = "Finished fitting Simulation example model")
}

# The fitting of `mcmc_mra()` took `r format(runtime, nsmall = 2, digits = 2)`.

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 224
sum(out$MRA$W %*% out$alpha[iter_idx, ])
sum(out$MRA$W[, out$MRA$dims_idx == 1] %*% out$alpha[iter_idx, out$MRA$dims_idx == 1])
sum(out$MRA$W[, out$MRA$dims_idx == 2] %*% out$alpha[iter_idx, out$MRA$dims_idx == 2])
sum(out$MRA$W[, out$MRA$dims_idx == 3] %*% out$alpha[iter_idx, out$MRA$dims_idx == 3])


# Examining the MRA MCMC output -----------------------------------------------

plot_trace(out, file = here::here("images", "trace-plots.png"))


plot_trace_alpha(out, file = here::here("images", "trace-plots-alpha.png"))

plot_posterior_params(out, file = here::here("images", "posterior-parameters.png"))


png(file = here::here("images", "estimated-vs-predicted.png"), width = 16, height = 9, units = "in", res = res)
## predicted vs. estimated beta

# undo the center and scaling transformation
beta_post <- unscale_beta(out)
sigma2_post <- unscale_sigma2(out)

p_beta <- data.frame(
    truth = beta,
    mean  = apply(beta_post, 2, mean),
    lower = apply(beta_post, 2, quantile, prob = 0.025),
    upper = apply(beta_post, 2, quantile, prob = 0.975),
    parameter = factor(1:length(beta))
) %>%
    ggplot(aes(x = truth, y = mean, color = parameter)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    ylab("estimate") +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. predicted beta") +
    theme(legend.position = "none")

## predicted vs. estimated sigma2
p_sigma2 <- data.frame(sigma2 = sigma2_post) %>%
    ggplot(aes(y = sigma2)) +
    geom_boxplot() +
    geom_point(data = data.frame(sigma2 = sigma2),
               aes(x = 0, y = sigma2), color = "red", size = 2) +
    ylab("estimate") +
    xlab("sigma2") +
    ggtitle("Estimated vs. predicted sigma2") +
    theme(legend.position = "none",
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())

## predicted vs. estimated tau2

p_tau2 <- data.frame(
    truth = tau2,
    mean  = apply(out$tau2, 2, mean),
    lower = apply(out$tau2, 2, quantile, prob = 0.025),
    upper = apply(out$tau2, 2, quantile, prob = 0.975),
    parameter = factor(1:length(tau2))) %>%
    ggplot(aes(x = truth, y = mean, color = parameter)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    ylab("estimate") +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. predicted tau2") +
    theme(legend.position = "none")


p_alpha <- data.frame(
    truth = alpha,
    mean = apply(out$alpha * sd_y, 2, mean),
    lower = apply(out$alpha * sd_y, 2, quantile, prob = 0.025),
    upper = apply(out$alpha * sd_y, 2, quantile, prob = 0.975),
    parameter = factor(1:length(alpha)),
    res       = factor(dims_idx)) %>%
    ggplot(aes(x = truth, y = mean, color = res)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. predicted alpha \nNon-identifiable but this is expected") +
    ylab("estimate")

print((p_beta + p_sigma2) / (p_tau2 + p_alpha))
dev.off()

W_alpha_res <- list()
W_alpha_post_res <- list()
for (m in 1:M) {
    W_alpha_res[[m]] <- W[s, dims_idx == m] %*% alpha[dims_idx == m]
    W_alpha_post_res[[m]] <- out$MRA$W[, out$MRA$dims_idx == m] %*% t(out$alpha[, out$MRA$dims_idx == m])
}

p_Walpha <- data.frame(
    truth = c(unlist(W_alpha_res),
              Reduce('+', W_alpha_res)),
    mean  = c(unlist(sapply(1:out$MRA$M, function(m) apply(W_alpha_post_res[[m]] * sd_y, 2, mean))),
              apply(Reduce('+', W_alpha_post_res) * sd_y, 2, mean)),
    lower = c(unlist(sapply(1:out$MRA$M, function(m) apply(W_alpha_post_res[[m]] * sd_y, 2, quantile, prob = 0.025))),
              apply(Reduce('+', W_alpha_post_res) * sd_y, 2, quantile, prob = 0.025)),
    upper = c(unlist(sapply(1:out$MRA$M, function(m) apply(W_alpha_post_res[[m]] * sd_y, 2, quantile, prob = 0.975))),
              apply(Reduce('+', W_alpha_post_res) * sd_y, 2, quantile, prob = 0.9755)),
    parameter = factor(1:n_sites),
    res       = factor(c(rep(1:M, each = n_sites), rep("full process", n_sites)))) %>%
    ggplot(aes(x = truth, y = mean, color = res)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.1) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. predicted Walpha \nNon-identifiable but this is expected") +
    ylab("estimate")

p_Walpha


plot_fitted_alphas(out, file = here::here("images", "fitted-spatial-random-effects.png"))

new_data <- list(
    locs_pred = locs,
    X_pred    = X
)

preds <- predict_mra(out, new_data)

plot_fitted_MRA(out, preds, file = here::here("images", "fitted-process.png"))


plot_predicted_process(out, new_data, preds)

plot_predicted_process(out, dat_plot, preds)


plot_sim_vs_fitted(out, X_s, beta, W[s, ], alpha)
plot_sim_vs_fitted(out, X_s, beta, W[s, ], base_size = 22, alpha, file = here::here("images", "sim-vs-fitted-process.png"))

# spBayes fit to subsample data -------------------------

set.seed(999)

knots.pp <- as.matrix(
        expand.grid(
            seq(0.01, 0.99, length.out = 8),
            seq(0.01, 0.99, length.out = 8)
        )
    )

if (file.exists(here::here("results", "fit-pp.RData"))) {
    load(here::here("results", "fit-pp.RData"))
} else {
    start   <- Sys.time()
    fit     <- spLM(y_s ~ X_s - 1,
                    coords = locs_s,
                    knots = knots.pp,
                    n.samples = 1000,
                    cov.model = "matern",
                    priors = list(
                        "beta.Norm" = list(rep(0, ncol(X_s)), diag(1000, ncol(X_s))),
                        "phi.Unif" = c(1/100, 100),
                        "nu.Unif" = c(1/100, 100),
                        "sigma.sq.IG" = c(0.1, 0.1),
                        "tau.sq.IG" = c(0.1, 0.1)
                    ),
                    starting = list(
                        # "beta" = rnorm(ncol(X_s)),
                        # "phi" = runif(1, 0.1, 10),
                        # "nu"  = runif(1, 0.5, 5),
                        # "sigma.sq" = runif(1, 1, 10),
                        # "tau.sq" = runif(1, 1, 10)

                        "beta" = rep(0, ncol(X_s)),
                        "phi" = 3,
                        "nu"  = 1,
                        "sigma.sq" = 10,
                        "tau.sq" = 1
                    ),
                    amcmc = list("n.batch" = 100,
                                 "batch.length" = 50,
                                 "accept.rate" = 0.44),
                    tuning = list(
                        "phi" = 1,
                        "nu"  = 0.25,
                        "sigma.sq" = 0.7,
                        "tau.sq" = 0.5
                    ),
                    n.report = 10
    )
    end     <- Sys.time()
    runtime <- end - start
    save(fit, runtime, file = here::here("results", "fit-pp.RData"))
    pushoverr::pushover(message = "Finished fitting PP example model")
}

if (file.exists(here::here("results", "fit-pp-recover.RData"))) {
    load(here::here("results", "fit-pp-recover.RData"))
} else {
    start   <- Sys.time()
    fit_out <- spRecover(fit, start = 2500, thin = 2, verbose = FALSE)
    end     <- Sys.time()
    runtime_recover <- end - start
    save(fit_out, runtime_recover, file = here::here("results", "fit-pp-recover.RData"))
    pushoverr::pushover(message = "Finished recovering PP example model parameters")
}

if (file.exists(here::here("results", "predict-pp-recover.RData"))) {
    load(here::here("results", "predict-pp-recover.RData"))
} else {
    ## predict using the tuned Matern fit
    start   <- Sys.time()
    preds_pp <- spPredict(
        fit_out,
        start = 2500,
        thin = 10,
        pred.coords = locs,
        pred.covars = X
    )
    end     <- Sys.time()
    runtime_pred <- end - start
    save(preds_pp, runtime_pred, file = here::here("results", "predict-pp-recover.RData"))
    pushoverr::pushover(message = "Finished recovering PP example model predictions")
}


theta.samps <- mcmc.list(fit_out$p.theta.samples, fit_out$p.theta.samples)
plot(theta.samps, density = FALSE)
beta.samps <- mcmc.list(fit_out$p.beta.recover.samples, fit_out$p.beta.recover.samples)
plot(beta.samps, density = FALSE)
round(summary(theta.samps)$quantiles, 3)
round(summary(beta.samps)$quantiles, 3)



# Examining the MRA MCMC output -----------------------------------------------

png(file = here::here("images", "estimated-vs-predicted-pp.png"), width = 16, height = 9, units = "in", res = res)
## predicted vs. estimated beta
dat_plot <- data.frame(
    truth = beta,
    mean  = apply(beta.samps[[1]], 2, mean),
    lower = apply(beta.samps[[1]], 2, quantile, prob = 0.025),
    upper = apply(beta.samps[[1]], 2, quantile, prob = 0.975),
    parameter = factor(1:length(beta))
)

p_beta <- dat_plot %>%
    ggplot(aes(x = truth, y = mean, color = parameter)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    ylab("estimate") +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. predicted beta") +
    theme(legend.position = "none")

## predicted vs. estimated sigma2
p_sigma2 <- data.frame(sigma2 = c(theta.samps[[1]][, "sigma.sq"])) %>%
    ggplot(aes(y = sigma2)) +
    geom_boxplot() +
    geom_point(data = data.frame(sigma2 = sigma2),
               aes(x = 0, y = sigma2), color = "red", size = 2) +
    ylab("estimate") +
    xlab("sigma2") +
    ggtitle("Estimated vs. predicted sigma2") +
    theme(legend.position = "none",
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
print((p_beta + p_sigma2))
dev.off()



png(file = here::here("images", "fitted-process-pp.png"), width = 16, height = 9, units = "in", res = res)
## Estimated predictive process

print(
    data.frame(
        W    = apply(preds_pp$p.y.predictive.samples, 1, mean),
        site = 1:N,
        lat  = locs[, 2],
        lon  = locs[, 1]
    ) %>%
        ggplot(aes(x = lat, y = lon, fill = W)) +
        geom_raster() +
        xlab("Longitude") +
        ylab("Latitude") +
        ggtitle("Estimated Predictive process") +
        colorspace::scale_fill_continuous_diverging("Blue-Red 3") +
        theme_bw()
)
dev.off()



png(file = here::here("images", "sim-vs-fitted-pp-process.png"), width = 16, height = 9, units = "in", res = res)

## plot estimated mean response
Xbeta_post <- t(X_s %*% t(beta.samps[[1]]))
W_post <- t(fit_out$p.w.recover.samples)
mu_post <- Xbeta_post + W_post

dat_plot <- data.frame(
    mean_Xbeta   = apply(Xbeta_post, 2, mean),
    lower_Xbeta  = apply(Xbeta_post, 2, quantile, prob = 0.025),
    upper_Xbeta  = apply(Xbeta_post, 2, quantile, prob = 0.975),
    truth_Xbeta  = X_s %*% beta,

    mean_W  = apply(W_post, 2, mean),
    lower_W = apply(W_post, 2, quantile, prob = 0.025),
    upper_W = apply(W_post, 2, quantile, prob = 0.975),
    truth_W = (W %*% alpha)[s],

    mean_mu      = apply(mu_post, 2, mean),
    lower_mu     = apply(mu_post, 2, quantile, prob = 0.025),
    upper_mu     = apply(mu_post, 2, quantile, prob = 0.975),
    truth_mu     = X_s %*% beta + (W %*% alpha)[s]
)

plot_Xbeta <- dat_plot %>%
    ggplot(aes(x = truth_Xbeta, y = mean_Xbeta)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_errorbar(aes(ymin = lower_Xbeta, ymax = upper_Xbeta)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated fixed effects")

plot_W <- dat_plot %>%
    ggplot(aes(x = truth_W, y = mean_W)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_errorbar(aes(ymin = lower_W, ymax = upper_W)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated spatial process")

plot_mu <- dat_plot %>%
    ggplot(aes(x = truth_mu, y = mean_mu)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_errorbar(aes(ymin = lower_mu, ymax = upper_mu)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated mean response")

print(plot_Xbeta + plot_W + plot_mu)
dev.off()
