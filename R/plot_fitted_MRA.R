#' Plot the fitted MRA process
#'
#' Plot the fitted spatial process from the output of `mcmc_mra()` or `mcmc_mra_integrated()` for each resolution.
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
plot_fitted_MRA <- function(out, preds, base_size = 12) {


    W_alpha_res = unlist(sapply(1:out$MRA$M, function(m) preds$MRA_pred$W_pred[, out$MRA$dims_idx == m] %*% apply(out$alpha, 2, mean)[out$MRA$dims_idx == m], simplify = "matrix"))
    dimnames(W_alpha_res) <- list(
        site = 1:nrow(W_alpha_res),
        res  = paste("resolution", 1:out$MRA$M)
    )
    dat_locs <- data.frame(
        site    = 1:nrow(W_alpha_res),
        Lat     = preds$locs_pred[, 1],
        Lon     = preds$locs_pred[, 2]
    )
    dat_full <- data.frame(
        site    = 1:nrow(W_alpha_res),
        Lat     = preds$locs_pred[, 1],
        Lon     = preds$locs_pred[, 2],
        W_alpha = apply(W_alpha_res, 1, sum),
        res     = "full process"
    )

    estimated_MRA <- merge(
        dat_locs,
        as.data.frame.table(W_alpha_res, responseName = "W_alpha")
    ) %>%
        rbind(dat_full) %>%
        ggplot(aes(x = Lon, y = Lat, fill = W_alpha)) +
        geom_raster() +
        facet_wrap( ~ res, ncol = 2) +
        xlab("Longitude") +
        ylab("Latitude") +
        ggtitle("Estimated Multiresolution process") +
        scale_fill_viridis_c(option = "B") +
        theme_bw(base_size = base_size)

    return(estimated_MRA)
}


