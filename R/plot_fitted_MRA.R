#' Plot the fitted MRA process
#'
#' Plot the fitted spatial process from the output of `mcmc_mra()` or `mcmc_mra_integrated()` for each resolution.
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param preds A prediction object that is the output of `predict_mra()`
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
plot_fitted_MRA <- function(out, preds, base_size = 12, file = NULL, width = 16, height = 9) {

    if (!(class(out) %in% c("mcmc_mra", "mcmc_mra_integrated")))
        stop('out must be of class "mcmc_mra" or "mcmc_mra_integrated"')
    if(class(preds) != "mcmc_mra_pred")
        stop('preds must be of class "mcmc_mra_pred" which is the output of the function predict_mra()')
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is_positive_numeric(base_size, 1))
        stop("base_size must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")



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
        ggplot(aes(x = .data$Lon, y = .data$Lat, fill = .data$W_alpha)) +
        geom_raster() +
        facet_wrap( ~ res, ncol = 2) +
        xlab("Longitude") +
        ylab("Latitude") +
        ggtitle("Estimated Multiresolution process") +
        scale_fill_viridis_c(option = "B") +
        theme_bw(base_size = base_size)

    if (is.null(file)) {
        return(estimated_MRA)
    } else {
        ggsave(filename = file,
               plot = estimated_MRA,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }

}


