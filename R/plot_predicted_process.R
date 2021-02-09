#' Plot the observed (or true) and predicted spatial data
#'
#' Plot the fitted model from the output of `mcmc_mra()` or `mcmc_mra_integrated()`. Also add in the observed (or true data if using simulated/hold out data) to compare model fit to the data.
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param data A dataframe with four columns: Lat, Lon, Observed, Truth
#' @param preds A prediction object that is the output of `predict_mra()`
#' @param base_size The base size for the plot
#' @param ...
#'
#' @return
#' @import ggplot2
#' @import patchwork
#' @export
#'
#'
plot_predicted_process <- function(out, data, preds, base_size = 12, file = NULL, width = 16, height = 9, ...) {

    ## add in titles, grouping of legends, change sd color scheme

    # add in error checking

    if(class(preds) != "mcmc_mra_pred")
        stop('preds must be of class "mcmc_mra_pred" which is the output of the function predict_mra()')
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")


    Wpred_mean <- apply(preds$Walpha_pred, 2, mean)
    y_pred_mean <- apply(preds$y_pred, 2, mean)
    y_pred_sd   <- apply(preds$y_pred, 2, sd)

    dat_pred <- data.frame(
        y = y_pred_mean,
        sd = y_pred_sd,
        Lat = preds$locs_pred[, 1],
        Lon = preds$locs_pred[, 2])

    zlims = range(
        range(y_pred_mean),
        range(data$Observed, na.rm = TRUE))

    p_obs <- data %>%
        ggplot(aes(x = Lat, y = Lon, fill = Observed)) +
        geom_raster() +
        scale_fill_viridis_c(option = "B", limits = zlims)
    p_fit <- dat_pred %>%
        ggplot(aes(x = Lat, y = Lon, fill = y)) +
        geom_raster() +
        scale_fill_viridis_c(option = "B", limits = zlims)
    p_sd <- dat_pred %>%
        ggplot(aes(x = Lat, y = Lon, fill = sd)) +
        geom_raster() +
        scale_fill_viridis_c(option = "B")

    if (is.null(file)) {
        return(p_obs + p_fit + p_sd)
    } else {
        ggsave(filename = file,
               plot = p_obs + p_fit + p_sd,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }
}

