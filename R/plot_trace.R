#' Trace plots from the output of `mcmc_mra()` or `mcmc_mra_integrated()`
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param base_size The base size for the plot
#' @param file If `file = NULL`, the ggplot object is returned. If `file` is not NULL, an image is saved to the file path specified by `file`
#' @param width If a file path is specified, `width` determines the width of the saved image (in inches)
#' @param height If a file path is specified, `height` determines the height of the saved image (in inches)
#'
#' @return Either a ggplot object of the model trace plots (not including the spatial coefficients alpha) (if `file = NULL`) or a saved image file with no return (`file` is not NULL)
#' @export
#'
#' @import patchwork
#' @import tidyverse
#' @import latex2exp
#'

plot_trace <- function(out, base_size = 12, file = NULL, width = 16, height = 9) {

    if (!inherits(out, c("mcmc_mra", "mcmc_mra_integrated")))
        stop('out must be of class "mcmc_mra" or "mcmc_mra_integrated"')
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is_positive_numeric(base_size, 1))
        stop("base_size must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")

    M <- out$MRA$M

    p_tau2 <- data.frame(
        tau2      = c(out$tau2),
        iteration = rep(1:nrow(out$tau2), times = M),
        resolution = factor(rep(1:M, each = nrow(out$tau2)))
    ) %>%
        ggplot(aes(x = .data$iteration, y = .data$tau2, group = .data$resolution, color = .data$resolution)) +
        geom_line(alpha = 0.75) +
        scale_color_viridis_d(end = 0.8) +
        ggtitle(TeX("Trace plots for $\\tau2$")) +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\tau^2$"))

    p_sigma2 <- data.frame(sigma2 = c(out$sigma2), iteration = 1:length(out$sigma2)) %>%
        ggplot(aes(x = .data$iteration, y = .data$sigma2)) +
        geom_line(alpha = 0.75) +
        scale_color_viridis_d(end = 0.8) +
        ggtitle(TeX("Trace plots for $\\sigma2$")) +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\sigma^2$"))

    p_beta  <- data.frame(
        beta      = c(out$beta),
        iteration = rep(1:nrow(out$beta), times = ncol(out$beta)),
        covariate = factor(rep(1:ncol(out$beta), each = nrow(out$beta)))
    ) %>%
        ggplot(aes(x = .data$iteration, y = .data$beta, group = .data$covariate, color = .data$covariate)) +
        geom_line(alpha = 0.75) +
        scale_color_viridis_d(end = 0.8) +
        ggtitle(TeX("Trace plots for $\\beta$")) +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\beta$"))


    if (is.null(file)) {
        return(p_tau2 / p_sigma2 / p_beta)
    } else {
        ggsave(filename = file,
               plot = p_tau2 / p_sigma2 / p_beta,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }
}
