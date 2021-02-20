#' Posterior distribution plots from the output of `mcmc_mra()` or `mcmc_mra_integrated()`
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param base_size The base size for the plot
#' @param alpha The plot transparency
#' @param file If `file = NULL`, the ggplot object is returned. If `file` is not NULL, an image is saved to the file path specified by `file`
#' @param width If a file path is specified, `width` determines the width of the saved image (in inches)
#' @param height If a file path is specified, `height` determines the height of the saved image (in inches)
#'
#' @return Either a ggplot object of the posterior distribution of some model parameters (if `file = NULL`) or a saved image file with no return (`file` is not NULL)
#'
##' @export
#'
#' @import patchwork
#' @import tidyverse
#' @import latex2exp
#'

plot_posterior_params <- function(out, base_size = 12, alpha = 0.1, file = NULL, width = 16, height = 9) {

    if (!(class(out) %in% c("mcmc_mra", "mcmc_mra_integrated")))
        stop('out must be of class "mcmc_mra" or "mcmc_mra_integrated"')
    if (!is_numeric(alpha, 1))
        stop("alpha must be a number between 0 and 1")
    if (alpha < 0 | alpha > 1)
        stop("alpha must be a number between 0 and 1")
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is_positive_numeric(base_size, 1))
        stop("base_size must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")



    M <- out$MRA$M
    p_beta <- data.frame(
        beta = c(out$beta),
        parameter = factor(rep(1:ncol(out$beta), each = nrow(out$beta)))
        ) %>%
        ggplot(aes(x = .data$parameter, y = .data$beta)) +
        geom_violin() +
        geom_point(position = "jitter", alpha = alpha) +
        geom_hline(yintercept = 0, col = "red", lty = 2) +
        ggtitle(TeX("Posterior Distribution for $\\beta$")) +
        theme(legend.position = "none") +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\beta$"))


    p_sigma2 <- data.frame(sigma2 = out$sigma2) %>%
        ggplot(aes(y = .data$sigma2, x = "")) +
        geom_violin() +
        geom_point(position = "jitter", alpha = alpha) +
        xlab("") +
        ggtitle(TeX("Posterior Distribution for $\\sigma^2$")) +
        theme_bw(base_size = base_size) +
        ylab(TeX("$\\sigma^2$"))


    p_tau2 <- dat_plot <- data.frame(
        tau2 = c(out$tau2),
        resolution = factor(rep(1:M, each = nrow(out$tau2)))) %>%
        ggplot(aes(x = .data$resolution, y = .data$tau2)) +
        geom_violin() +
        geom_point(position = "jitter", alpha = alpha) +
        geom_hline(yintercept = 0, col = "red", lty = 2) +
        ggtitle(TeX("Posterior Distribution for $\\tau^2$")) +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\tau^2$"))


    if (is.null(file)) {
        return(p_beta / p_sigma2 / p_tau2)
    } else {
        ggsave(filename = file,
               plot = p_beta / p_sigma2 / p_tau2,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }
}
