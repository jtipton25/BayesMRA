#' Plot the observed and "true" (holdout) spatial data
#'
#' Plot the fitted model from the output of `mcmc_mra()` or `mcmc_mra_integrated()`. Also add in the observed (or true data if using simulated/hold out data) to compare model fit to the data.
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param data A dataframe with four columns: Lat, Lon, Observed, Fitted
#' @param base_size The base size for the plot
#' @param ...
#'
#' @return
#' @import ggplot2
#' @import patchwork
#' @export
#'
#'
plot_fitted_process <- function(out, data, base_size = 12, ...) {

    X <- out$data$X
    for(i in 2:ncol(X)) {
        X[, i] <- (X[, i] - out$data$mu_X[i-1]) / out$data$sd_X[i-1]
    }
    Xbeta_post <- t(X %*% t(out$beta)) * out$data$sd_y + out$data$mu_y
    Walpha_post <- t(out$MRA$W %*% t(out$alpha)) * out$data$sd_y
    mu_post <- Xbeta_post + Walpha_post

    dat_plot <- data.frame(
        Lat          = out$data$locs[, 1],
        Lon          = out$data$locs[, 2],
        mean_Xbeta   = apply(Xbeta_post, 2, mean),
        mean_Walpha  = apply(Walpha_post, 2, mean),
        mean_mu      = apply(mu_post, 2, mean))
    zlims = range(
        range(dat_plot$mean_mu),
        range(data$Observed, na.rm = TRUE))

    plot_Xbeta <- dat_plot %>%
        ggplot(aes(x = Lat, y = Lon, fill = mean_Xbeta)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B")

    plot_Walpha <- dat_plot %>%
        ggplot(aes(x = Lat, y = Lon, fill = mean_Walpha)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B")

    plot_mu <- dat_plot %>%
        ggplot(aes(x = Lat, y = Lon, fill = mean_mu)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B", limits = zlims)

    plot_obs <- data %>%
        filter(!is.na(Observed)) %>%
        ggplot(aes(x = Lat, y = Lon, fill = Observed)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B", limits = zlims)

    (plot_Xbeta + plot_Walpha) / (plot_mu + plot_obs)

    return((plot_Xbeta + plot_Walpha) / (plot_mu + plot_obs))
}
