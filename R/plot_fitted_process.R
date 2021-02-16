#' Plot the observed and "true" (holdout) spatial data
#'
#' Plot the fitted model from the output of `mcmc_mra()` or `mcmc_mra_integrated()`. Also add in the observed (or true data if using simulated/hold out data) to compare model fit to the data.
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param data A dataframe with four columns: Lat, Lon, Observed, Fitted
#' @param base_size The base size for the plot
#' @param file If `file = NULL`, the ggplot object is returned. If `file` is not NULL, an image is saved to the file path specified by `file`
#' @param width If a file path is specified, `width` determines the width of the saved image (in inches)
#' @param height If a file path is specified, `height` determines the height of the saved image (in inches)
#'
#' @return Either a ggplot object (if `file = NULL`) or a saved image file with no return (`file` is not NULL)
#' @import ggplot2
#' @import patchwork
#' @export
#'
#'
plot_fitted_process <- function(out, data, base_size = 12, file = NULL, width = 16, height = 9) {

    if (!(class(out) %in% c("mcmc_mra", "mcmc_mra_integrated")))
        stop('out must be of class "mcmc_mra" or "mcmc_mra_integrated"')
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")


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
        ggplot(aes(x = .data$Lat, y = .data$Lon, fill = .data$mean_Xbeta)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B")

    plot_Walpha <- dat_plot %>%
        ggplot(aes(x = .data$Lat, y = .data$Lon, fill = .data$mean_Walpha)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B")

    plot_mu <- dat_plot %>%
        ggplot(aes(x = .data$Lat, y = .data$Lon, fill = .data$mean_mu)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B", limits = zlims)

    plot_obs <- data %>%
        filter(!is.na(.data$Observed)) %>%
        ggplot(aes(x = .data$Lat, y = .data$Lon, fill = .data$Observed)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B", limits = zlims)

    (plot_Xbeta + plot_Walpha) / (plot_mu + plot_obs)

    if (is.null(file)) {
        return((plot_Xbeta + plot_Walpha) / (plot_mu + plot_obs))
    } else {
        ggsave(filename = file,
               plot = (plot_Xbeta + plot_Walpha) / (plot_mu + plot_obs),
               device = "png",
               width = width,
               height = height,
               units = "in")
    }
}
