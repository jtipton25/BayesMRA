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
plot_fitted_alphas <- function(out, base_size = 12, file = NULL, width = 16, height = 9, ...) {

    if (!(class(out) %in% c("mcmc_mra", "mcmc_mra_integrated")))
        stop('out must be of class "mcmc_mra" or "mcmc_mra_integrated"')
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")



    fitted_alphas <- data.frame(
        x = unlist(sapply(1:out$MRA$M, function(i) out$MRA$locs_grid[[i]][, 1])),
        y = unlist(sapply(1:out$MRA$M, function(i) out$MRA$locs_grid[[i]][, 2])),
        res = factor(unlist(sapply(1:out$MRA$M, function(i) rep(i, each = nrow(out$MRA$locs_grid[[i]]))))),
        alpha = unlist(sapply(1:out$MRA$M, function(i) apply(out$alpha, 2, mean)[out$MRA$dims_idx == i]))
    ) %>%
        ggplot(aes(x = x, y = y, fill = alpha)) +
        geom_raster() +
        scale_fill_viridis_c() +
        ggtitle("Posterior mean spatial random effects") +
        facet_wrap( ~ res, ncol = 2) +
        geom_rect(
            data = NULL,
            aes(xmin = min(out$data$locs[, 1]),
                xmax = max(out$data$locs[, 1]),
                ymin = min(out$data$locs[, 2]),
                ymax = max(out$data$locs[, 2])
            ),
            fill  = NA,
            color = "black",
            alpha = 0.5
        ) +
        theme_bw(base_size)

    if (is.null(file)) {
        return(fitted_alphas)
    } else {
        ggsave(filename = file,
               plot = fitted_alphas,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }
}
