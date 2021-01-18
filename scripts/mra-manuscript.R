library(tidyverse)
library(BayesMRA)
library(mvnfast)
# library(MCMCpack)
library(tidyverse)
library(patchwork)
library(spam64)
library(fields)
library(spBayes)
library(coda)
library(here)

# setup the project folder directories
if (!dir.exists(here::here("results"))) {
    dir.create(here::here("results"))
}
if (!dir.exists(here::here("images"))) {
    dir.create(here::here("images"))
}

res <- 100 # can increase before submission to have better resolution figures

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
# knitr::opts_chunk$set(fig.width = 8, fig.height = 4.5, dpi = 200, dev.args = list(type = "cairo"))


# Simulate data ---------------------------------------------------------------

set.seed(11)

# N <- 40^2
N <- 1e4
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

beta <- rnorm(ncol(X))

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
## set up a linear constraint so that each resolution sums to one
# A_constraint <- sapply(1:M, function(i){
#     tmp = rep(0, sum(n_dims))
#     tmp[dims_idx == i] <- 1
#     return(tmp)
# })
# a_constraint <- rep(0, M)

# constraint on full effect Walpha
# A_constraint <- drop(rep(1, N) %*% W)
# a_constraint <- 0

## constraint on each resolution of Walpha
A_constraint_tmp <- sapply(1:M, function(i){
    tmp <- rep(1, N) %*% W[, dims_idx == i]
    return(tmp)
})
A_constraint <- matrix(0, M, sum(n_dims))
for (i in 1:M) {
    A_constraint[i, dims_idx == i] <- A_constraint_tmp[[i]]
}

a_constraint <- rep(0, M)

alpha   <- as.vector(rmvnorm.prec.const(n = 1, mu = rep(0, sum(n_dims)), Q = Q_alpha_tau2, A = A_constraint, a = a_constraint))


## verify that the parameters sum to 0
sum(alpha)
sum(alpha[dims_idx == 1])
sum(alpha[dims_idx == 2])
sum(alpha[dims_idx == 3])

sum(W %*% alpha)
sum(W[, dims_idx == 1] %*% alpha[dims_idx == 1])
sum(W[, dims_idx == 2] %*% alpha[dims_idx == 2])
sum(W[, dims_idx == 3] %*% alpha[dims_idx == 3])

sigma2 <- runif(1, 0.25, 0.5)

y <- X %*% beta + W %*% alpha + rnorm(N, 0, sqrt(sigma2))


# Plot the simulated data -----------------------------------------------------

# plot precision matrices
## plot Q_alpha_tau2 structure
## note: the display function gives a warning about the default cex setting
png(file = here::here("images", "precision-matrix-sparsity.png"), width = 16, height = 9, units = "in", res = res)
layout(matrix(1:4, 2, 2, byrow = TRUE))
for (m in 1:M) {
    display(Q_alpha_tau2[which(dims_idx == m), which(dims_idx == m)])
    title(main = paste("Nonzero elements of \n precision matrix; resolution", m))
}
dev.off()


## Plot the MRA grid
if (!file.exists(here::here("images", "MRA-grid.png"))) {
    png(file = here::here("images", "MRA-grid.png"), width = 16, height = 9, units = "in", res = res)
    p <- data.frame(
        x = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 1])),
        y = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 2])),
        resolution = factor(unlist(sapply(1:M, function(i) rep(i, each = nrow(MRA$locs_grid[[i]])))))
    ) %>%
        ggplot(aes(x = x, y = y, color = resolution)) +
        geom_point(alpha = 0.75, size = 0.75) +
        ggtitle("MRA grid points") +
        theme_bw() +
        scale_colour_viridis_d(end = 0.8) +
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
        theme_custom()
    print(p)
    dev.off()
}



if (!file.exists(here::here("images", "MRA-grid-neighbors.png"))) {
    png(file = here::here("images", "MRA-grid-neighbors.png"), width = 16, height = 9, units = "in", res = res)

    nnn <- sapply(1:M, function(i) {
        D <- rdist(MRA$locs_grid[[i]])
        sapply(1:nrow(MRA$locs_grid[[i]]), function(j) {
            sum((D[j, ]) < MRA$radius[i])}
        )
    })

    p <- data.frame(
        x = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 1])),
        y = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 2])),
        resolution = factor(unlist(sapply(1:M, function(i) rep(paste("resolution", i), each = nrow(MRA$locs_grid[[i]]))))),
        neighbors <- unlist(nnn)
    )  %>%
        ggplot(aes(x = x, y = y, fill = neighbors)) +
        geom_raster() +
        facet_wrap(~ resolution) +
        ggtitle("number of neighbors") +
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
        theme_bw() +
        theme_custom()

    print(p)
    dev.off()
}


## plot the random effects
if (!file.exists(here::here("images", "MRA-random-effects.png"))) {
    png(file = here::here("images", "MRA-random-effects.png"), width = 16, height = 9, units = "in", res = res)

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
        theme_bw() +
        theme_custom()

    plot(sim_alphas)
    dev.off()
}


## plot the simulated data
if (!file.exists(here::here("images", "MRA-simulated-data.png"))) {
    png(file = here::here("images", "MRA-simulated-data.png"), width = 16, height = 9, units = "in", res = res)

    zlims <- c(
        min(c(X %*% beta), c(y), c(W %*% alpha)),
        max(c(X %*% beta), c(y), c(W %*% alpha))
    )


    g_fixed <- data.frame(
        temp = c(X %*% beta),
        x    = locs[, 1],
        y    = locs[, 2]
    ) %>%
        ggplot(aes(x = x, y = y, fill = temp)) +
        geom_raster() +
        xlab("x") +
        ylab("y") +
        ggtitle("Fixed effects") +
        colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = zlims) +
        theme_bw() +
        theme_custom(0.70)

    dat_map <- data.frame(
        temp = c(y),
        re   = c(W %*% alpha),
        x    = locs[, 1],
        y    = locs[, 2]
    )

    g_full <- dat_map %>%
        ggplot(aes(x = x, y = y, fill = temp)) +
        geom_raster() +
        xlab("x") +
        ylab("y") +
        ggtitle("Observed surface") +
        colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = zlims) +
        theme_bw() +
        theme_custom(0.70)

    g_re <- dat_map %>%
        ggplot(aes(x = x, y = y, fill = re)) +
        geom_raster() +
        xlab("x") +
        ylab("y") +
        ggtitle("Spatial process") +
        colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = zlims) +
        theme_bw() +
        theme(legend.position = "none") +
        theme_custom(0.70)

    print(g_fixed + g_re + g_full + plot_layout(guides = "collect"))
    dev.off()
}



# plot the MRA spatial process
if (!file.exists(here::here("images", "MRA-resolutions.png"))) {
    png(file = here::here("images", "MRA-resolutions.png"), width = 16, height = 9, units = "in", res = res)

    W_alpha_res = unlist(sapply(1:M, function(m) W[, dims_idx == m] %*% alpha[dims_idx == m], simplify = "matrix"))
    dimnames(W_alpha_res) <- list(
        site = 1:N,
        res  = paste("resolution", 1:M)
    )
    dat_locs <- data.frame(
        site = 1:N,
        x    = locs[, 1],
        y    = locs[, 2]
    )

    sim_MRA <- rbind(
        merge(dat_locs, as.data.frame.table(W_alpha_res, responseName = "W_alpha")),
        data.frame(
            site    = 1:N,
            x       = locs[, 1],
            y       = locs[, 2],
            W_alpha = W %*% alpha,
            res     = "full process"
        )
    ) %>%
        ggplot(aes(x = x, y = y, fill = W_alpha)) +
        geom_raster() +
        facet_wrap( ~ res, ncol = 2) +
        xlab("x") +
        ylab("y") +
        ggtitle("Multi-resolution process") +
        colorspace::scale_fill_continuous_diverging("Blue-Red 3") +
        theme_bw() +
        theme_custom()
    print(sim_MRA)
    dev.off()
}


# Sample the simulated data --------------------------------------------------

set.seed(111)
n_sites   <- 500
s         <- sample(N, n_sites)
y_s       <- y[s, ]
y_oos     <- y[ - s, ]
X_s       <- X[s, ]
X_oos     <- X[ - s, ]
locs_s    <- locs[s, ]
locs_oos  <- locs[ - s, ]


if (!file.exists(here::here("images", "sample-data.png"))) {
    png(file = here::here("images", "sample-data.png"), width = 16, height = 9, units = "in", res = res)
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
        theme_custom() +
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
    Sigma_beta   = 5 * diag(ncol(X_s)))


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
        n_cores       = 1L,
        verbose       = FALSE
    )
    end     <- Sys.time()
    runtime <- end - start
    save(out, runtime, file = here::here("results", "fit-mra-sim.RData"))
}

# The fitting of `mcmc_mra()` took `r format(runtime, nsmall = 2, digits = 2)`.

## Check the output to make sure the sum to 0 constraints are satisfied
iter_idx <- 224
sum(out$MRA$W %*% out$alpha[iter_idx, ])
sum(out$MRA$W[, out$MRA$dims_idx == 1] %*% out$alpha[iter_idx, out$MRA$dims_idx == 1])
sum(out$MRA$W[, out$MRA$dims_idx == 2] %*% out$alpha[iter_idx, out$MRA$dims_idx == 2])
sum(out$MRA$W[, out$MRA$dims_idx == 3] %*% out$alpha[iter_idx, out$MRA$dims_idx == 3])


# Examining the MRA MCMC output -----------------------------------------------


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

png(file = here::here("images", "trace-plots.png"), width = 16, height = 9, units = "in", res = 100)
(p_tau2 + p_lambda) /  (p_sigma2 + p_beta)
dev.off()


png(file = here::here("images", "trace-plots-alpha.png"), width = 16, height = 9, units = "in", res = 100)
## posterior plots for alpha
dat_alpha <- data.frame(
    alpha = c(out$alpha),
    iteration = rep(1:nrow(out$alpha), each = ncol(out$alpha)),
    knot      = factor(unlist(sapply(1:M, function(m) 1:n_dims[m]))),
    resolution = factor(paste("resolution", dims_idx))
)

plot_idx <- sapply(1:M, function(m) sample(n_dims[m], 10))

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
dev.off()



png(file = here::here("images", "estimated-vs-predicted.png"), width = 16, height = 9, units = "in", res = res)
## predicted vs. estimated beta
dat_plot <- data.frame(
    truth = beta,
    mean  = apply(out$beta, 2, mean),
    lower = apply(out$beta, 2, quantile, prob = 0.025),
    upper = apply(out$beta, 2, quantile, prob = 0.975),
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
p_sigma2 <- data.frame(sigma2 = out$sigma2) %>%
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
dat_plot <- data.frame(
    truth = tau2,
    mean  = apply(out$tau2, 2, mean),
    lower = apply(out$tau2, 2, quantile, prob = 0.025),
    upper = apply(out$tau2, 2, quantile, prob = 0.975),
    parameter = factor(1:length(tau2))
)


p_tau2 <- dat_plot %>%
    ggplot(aes(x = truth, y = mean, color = parameter)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    ylab("estimate") +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. predicted tau2") +
    theme(legend.position = "none")

# plot alpha estimates
dat_plot <- data.frame(
    truth = alpha,
    mean = apply(out$alpha, 2, mean),
    lower = apply(out$alpha, 2, quantile, prob = 0.025),
    upper = apply(out$alpha, 2, quantile, prob = 0.975),
    parameter = factor(1:length(alpha)),
    res       = factor(dims_idx)
)

p_alpha <- dat_plot %>%
    ggplot(aes(x = truth, y = mean, color = res)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. predicted alpha /n non-identifialbe but not really a concern") +
    ylab("estimate")

print((p_beta + p_sigma2) / (p_tau2 + p_alpha))
dev.off()



png(file = here::here("images", "fitted-spatial-random-effects.png"), width = 16, height = 9, units = "in", res = res)
fitted_alphas <- data.frame(
    x = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 1])),
    y = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 2])),
    res = factor(unlist(sapply(1:M, function(i) rep(i, each = nrow(MRA$locs_grid[[i]]))))),
    alpha = unlist(sapply(1:M, function(i) apply(out$alpha, 2, mean)[dims_idx == i]))
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
    theme_custom()
print(fitted_alphas)
dev.off()


png(file = here::here("images", "fitted-process.png"), width = 16, height = 9, units = "in", res = res)
## Estimated MRA spatial proces
W_alpha_res <- unlist(sapply(1:M, function(m)
    W[, dims_idx == m] %*% apply(out$alpha, 2, mean)[dims_idx == m],
    simplify = "matrix"))
dimnames(W_alpha_res) <- list(
    site = 1:N,
    res  = paste("resolution", 1:M)
)
dat_locs <- data.frame(
    site = 1:N,
    lat  = locs[, 2],
    lon  = locs[, 1]
)

estimated_MRA <- rbind(
    merge(dat_locs, as.data.frame.table(W_alpha_res, responseName = "W_alpha")),
    data.frame(
        site = 1:N,
        lat  = locs[, 2],
        lon  = locs[, 1],
        W_alpha = apply(W %*% t(out$alpha), 1, mean),
        res = "full process"
    )
) %>%
    ggplot(aes(x = lon, y = lat, fill = W_alpha)) +
    geom_raster() +
    facet_wrap( ~ res, ncol = 2) +
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle("Estimated Multi-resolution process") +
    colorspace::scale_fill_continuous_diverging("Blue-Red 3") +
    theme_bw() +
    theme_custom()

print(estimated_MRA)
dev.off()



png(file = here::here("images", "sim-vs-fitted-process.png"), width = 16, height = 9, units = "in", res = res)

## plot estimated mean response
Xbeta_post <- t(X_s %*% t(out$beta))
Walpha_post <- t(out$MRA$W %*% t(out$alpha))
mu_post <- Xbeta_post + Walpha_post

dat_plot <- data.frame(
    mean_Xbeta   = apply(Xbeta_post, 2, mean),
    lower_Xbeta  = apply(Xbeta_post, 2, quantile, prob = 0.025),
    upper_Xbeta  = apply(Xbeta_post, 2, quantile, prob = 0.975),
    truth_Xbeta  = X_s %*% beta,

    mean_Walpha  = apply(Walpha_post, 2, mean),
    lower_Walpha = apply(Walpha_post, 2, quantile, prob = 0.025),
    upper_Walpha = apply(Walpha_post, 2, quantile, prob = 0.975),
    truth_Walpha = (W %*% alpha)[s],

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

plot_Walpha <- dat_plot %>%
    ggplot(aes(x = truth_Walpha, y = mean_Walpha)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_errorbar(aes(ymin = lower_Walpha, ymax = upper_Walpha)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated spatial process")

plot_mu <- dat_plot %>%
    ggplot(aes(x = truth_mu, y = mean_mu)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_errorbar(aes(ymin = lower_mu, ymax = upper_mu)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated mean response")

print(plot_Xbeta + plot_Walpha + plot_mu)
dev.off()



# spBayes fit to subsample data -------------------------



set.seed(999)

knots.pp <- as.matrix(
        expand.grid(
            seq(0.01, 0.99, length.out = 4),
            seq(0.01, 0.99, length.out = 4)
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
}

if (file.exists(here::here("results", "fit-pp-recover.RData"))) {
    load(here::here("results", "fit-pp-recover.RData"))
} else {
    fit_out <- spRecover(fit, start = 2500, thin = 2, verbose = FALSE)
    save(fit_out, file = here::here("results", "fit-pp-recover.RData"))
}

if (file.exists(here::here("results", "predict-pp-recover.RData"))) {
    load(here::here("results", "predict-pp-recover.RData"))
} else {
    ## predict using the tuned Matern fit
    preds <- spPredict(
        fit_out,
        start = 2500,
        thin = 10,
        pred.coords = locs,
        pred.covars = X
    )
    save(preds, file = here::here("results", "predict-pp-recover.RData"))
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
        W    = apply(preds$p.y.predictive.samples, 1, mean),
        site = 1:N,
        lat  = locs[, 2],
        lon  = locs[, 1]
    ) %>%
        ggplot(aes(x = lon, y = lat, fill = W)) +
        geom_raster() +
        xlab("Longitude") +
        ylab("Latitude") +
        ggtitle("Estimated Predictive process") +
        colorspace::scale_fill_continuous_diverging("Blue-Red 3") +
        theme_bw() +
        theme_custom()
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
