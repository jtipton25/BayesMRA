##
## Need to speed up the model initialization
##

library(BayesMRA)
library(tidyverse)
library(patchwork)
set.seed(404)

N <- 100^2
locs <- as.matrix(expand_grid(lon = seq(0, 1, length.out = sqrt(N)),
                    lat  = seq(0, 1, length.out = sqrt(N))))

# explore simulations from the MRA-vecchia process
M <- 3
grid_id <- partition_observations(locs, M)
dat <- data.frame(lon = locs[, 1], lat = locs[, 2], grid_id = grid_id)

ggplot(dat, aes(x = lon, y = lat, fill = grid_id)) +
    geom_raster() +
    scale_fill_viridis_c(option = "plasma") +
    ggtitle("Grid IDs")

# Simulate data from the Vecchia MRA process

N <- 50^2
locs <- as.matrix(expand_grid(lon = seq(0, 1, length.out = sqrt(N)),
                              lat  = seq(0, 1, length.out = sqrt(N))))

# setup the process
M <- 1
W_list <- vector(mode = "list", length = M)
tau2 <- 100 / (2^(1:M - 1))
Walpha_full <- rep(0, N)
for (m in 1:M) {
    W_list[[m]] <- vector(mode = "list", length = 4^(m-1))

    # this is the slow function
    if (m == 1) {
        grid_idx <- partition_observations(locs, m)
        MRA <- mra_wendland_2d(locs, M = 1)
        Q_alpha <- list(make_Q_alpha_2d(sqrt(MRA$n_dims), phi = 0.99))
        class(Q_alpha) <- "Q_alpha"
        class(Q_alpha[[1]]) <- "spam"
        Q <- make_Q_alpha_tau2(Q_alpha, tau2[m])
        constraints <- make_constraint(MRA, constraint = "overall", joint = TRUE)
        A_constraint <- constraints$A_constraint
        a_constraint <- constraints$a_constraint
        alpha <- drop(spam::rmvnorm.prec.const(1, rep(0, MRA$n_dims), Q, A = A_constraint, a = a_constraint))
        Walpha = MRA$W %*% alpha
        W_list[[m]] <- list(MRA = MRA, Q = Q, grid_idx = grid_idx, alpha = alpha, Walpha = Walpha,
                            A_constraint = A_constraint, a_constraint = a_constraint)
        Walpha_full <- Walpha_full + Walpha
    } else {
        grid_idx <- partition_observations(locs, m)
        for (j in 1:(4^(m-1))) {
            MRA <- mra_wendland_2d(locs[grid_idx == j, ], M = 1)
            Q_alpha <- list(make_Q_alpha_2d(sqrt(MRA$n_dims), phi = 0.99))
            class(Q_alpha) <- "Q_alpha"
            class(Q_alpha[[1]]) <- "spam"
            Q <- make_Q_alpha_tau2(Q_alpha, tau2[m])
            constraints <- make_constraint(MRA, constraint = "overall", joint = TRUE)
            A_constraint <- constraints$A_constraint
            a_constraint <- constraints$a_constraint
            alpha <- drop(spam::rmvnorm.prec.const(1, rep(0, MRA$n_dims), Q, A = A_constraint, a = a_constraint))
            Walpha = MRA$W %*% alpha
            W_list[[m]][[j]] <- list(MRA = MRA, Q = Q, grid_idx = grid_idx, alpha = alpha, Walpha = Walpha,
                                     A_constraint = A_constraint, a_constraint = a_constraint)
            Walpha_full[grid_idx == j] <- Walpha_full[grid_idx == j] + Walpha
        }
    }
}

# simulate from the process
Walpha_full <- rep(0, N)
for (m in 1:M) {

    if (m == 1) {
        alpha <- drop(spam::rmvnorm.prec.const(1, rep(0, W_list[[m]]$MRA$n_dims), W_list[[m]]$Q,
                                         A = W_list[[m]]$A_constraint, a = W_list[[m]]$a_constraint))
        Walpha = W_list[[m]]$MRA$W %*% alpha
        Walpha_full <- Walpha_full + Walpha
    } else {
        for (j in 1:(4^(m-1))) {
            alpha <- drop(spam::rmvnorm.prec.const(1, rep(0, W_list[[m]][[j]]$MRA$n_dims), W_list[[m]][[j]]$Q,
                                                   A = W_list[[m]][[j]]$A_constraint, a = W_list[[m]][[j]]$a_constraint))
            Walpha = W_list[[m]][[j]]$MRA$W %*% alpha
            Walpha_full[W_list[[m]][[j]]$grid_idx == j] <- Walpha_full[W_list[[m]][[j]]$grid_idx == j] + Walpha
        }
    }
}

sigma2 <- 0.05^2
y <- rnorm(N, Walpha_full, sqrt(sigma2))

dat <- data.frame(lon = locs[, 1], lat = locs[, 2], Walpha = Walpha_full, y = y)

ggplot(dat, aes(x = lon, y = lat, fill = Walpha)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Simulated latent process")

ggplot(dat, aes(x = lon, y = lat, fill = y)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Simulated observations")

# Fit the MCMC model
M <- 3
tau2 <- 100 / (2^(1:M - 1))

params <- list(n_adapt = 100, n_mcmc = 400, n_thin = 1, n_message = 50)

if (file.exists(here::here("results", "mra-vecchia-normalize.RData"))) {
    load(here::here("results", "mra-vecchia-normalize.RData"))
} else {
    start <- Sys.time()
    out <- mcmc_mra_vecchia(y, locs, params, M = 3, normalize = TRUE)
    time_vecchia <- Sys.time() - start
    save(out, time_vecchia, file = here::here("results", "mra-vecchia-normalize.RData"))
}


if (file.exists(here::here("results", "mra-vecchia-non-normalize.RData"))) {
    load(here::here("results", "mra-vecchia-non-normalize.RData"))
} else {
    start <- Sys.time()
    out_non <- mcmc_mra_vecchia(y, locs, params, M = 3, normalize = FALSE)
    time_vecchia_non <- Sys.time() - start
    save(out_non, time_vecchia_non, file = here::here("results", "mra-vecchia-non-normalize.RData"))
}


dat <- data.frame(lon = locs[, 1], lat = locs[, 2],
                  mu = apply(Walpha_full, 1, sum),
                  mu_post= apply(out$mu, 2, mean),
                  mu_sd= apply(out$mu, 2, sd))
dat_non <- data.frame(lon = locs[, 1], lat = locs[, 2],
                  mu = apply(Walpha_full, 1, sum),
                  mu_post= apply(out_non$mu, 2, mean),
                  mu_sd= apply(out_non$mu, 2, sd))

layout(matrix(1:4, 2, 2))
plot(out$sigma2, type = 'l')
abline(h = sigma2, col = "red")
matplot(out$tau2, type = 'l')
matplot(out$mu[, sample(1:ncol(out$mu), 20)], type = 'l')

layout(matrix(1:4, 2, 2))
plot(out_non$sigma2, type = 'l')
abline(h = sigma2, col = "red")
matplot(out_non$tau2, type = 'l')
matplot(out_non$mu[, sample(1:ncol(out_non$mu), 20)], type = 'l')

p_predict <- ggplot(dat, aes(x = mu, y = mu_post)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red")


p_mu <- ggplot(dat, aes(x = lon, y = lat, fill = mu)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Simulated mu")

p_mu_post <- ggplot(dat, aes(x = lon, y = lat, fill = mu_post)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Estimated mu")

p_mu_sd <- ggplot(dat, aes(x = lon, y = lat, fill = mu_sd)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("sd(mu)")

(p_mu + p_mu_post) / (p_predict + p_mu_sd)

p_predict_non <- ggplot(dat_non, aes(x = mu, y = mu_post)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red")


p_mu_non <- ggplot(dat_non, aes(x = lon, y = lat, fill = mu)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Simulated mu")

p_mu_post_non <- ggplot(dat_non, aes(x = lon, y = lat, fill = mu_post)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Estimated mu")

p_mu_sd_non <- ggplot(dat_non, aes(x = lon, y = lat, fill = mu_sd)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("sd(mu)")

(p_mu_non + p_mu_post_non) / (p_predict_non + p_mu_sd_non)

p_mu_post + p_mu_post_non
p_mu_sd + p_mu_sd_non

plot(dat$mu_post, dat_non$mu_post)
abline(0, 1, col = "red")
dat %>%
    mutate(normalize = "yes") %>%
    left_join(dat_non %>%
                  mutate(normalize = "no")) %>%
    pivot_wider(cols = mu, names_from = normalize)
    ggplot(aes(x = yes, y = no))



# ggsave(filename = "~/Desktop/mcmc_mra_vecchia.png",
#        p_mu + p_mu_post,
#        device = "png",
#        height = 9, width = 16, units = "in")

if (file.exists(here::here("results", "mra-vecchia-mra.RData"))) {
    load(here::here("results", "mra-vecchia-mra.RData"))
} else {
    start_mra <- Sys.time()
    out_mra <- mcmc_mra(y, matrix(1, N, 1), locs, params, joint = FALSE, constraint = "resolution")
    time_mra <- Sys.time() - start_mra
    save(out_mra, time_mra, file = here::here("results", "mra-vecchia-mra.RData"))
}



mu_mra <- apply(out_mra$data$X %*% t(out_mra$beta) + out_mra$MRA$W %*% t(out_mra$alpha), 1, mean) * out_mra$data$sd_y +  out_mra$data$mu_y
mu_sd_mra <- apply(out_mra$data$X %*% t(out_mra$beta) + out_mra$MRA$W %*% t(out_mra$alpha), 1, sd) * out_mra$data$sd_y

dat <- data.frame(lon = locs[, 1], lat = locs[, 2],
                  mu = apply(Walpha_full, 1, sum),
                  mu_post = mu_mra,
                  mu_sd = mu_sd_mra)

plot(out_mra$sigma2, type = 'l')
abline(h = sigma2, col = "red")
matplot(out_mra$tau2, type = 'l')

matplot(t(out_mra$MRA$W %*% t(out_mra$alpha))[, sample(1:ncol(out$mu), 20)], type = 'l')

p_predict_mra <- ggplot(dat, aes(x = mu, y = mu_post)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red")


p_mu_mra <- ggplot(dat, aes(x = lon, y = lat, fill = mu)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Simulated mu")

p_mu_post_mra <- ggplot(dat, aes(x = lon, y = lat, fill = mu_post)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("Estimated mu")

p_mu_sd_mra <- ggplot(dat, aes(x = lon, y = lat, fill = mu_sd)) +
    geom_raster() +
    scale_fill_viridis_c() +
    ggtitle("sd(mu)")

(p_mu_mra + p_mu_post_mra) / (p_predict_mra + p_mu_sd_mra)


(p_mu_post + p_mu_post_mra) / (p_mu_sd + p_mu_sd_mra)

message("The Normalized Vecchia MRA model took ", round(time_vecchia, digits = 2), " ", units(time_vecchia))
message("The Vecchia MRA model took ", round(time_vecchia_non, digits = 2), " ", units(time_vecchia_non))
message("The MRA model took ", round(time_mra, digits = 2), " ", units(time_mra))
