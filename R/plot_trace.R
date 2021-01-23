#' Trace plots from the output of `mcmc_mra()` or `mcmc_mra_integrated()`
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param base_size The base size for the plot
#'
#' @return Trace plots of the parameters
#' @export
#'
#' @import patchwork
#' @import tidyverse
#' @import latex2exp
#'

plot_trace <- function(out, base_size = 12) {
    if (class(out) != "mcmc_mra")
        stop('out must be of class "mcmc_mra"')

    M <- out$MRA$M

    p_tau2 <- data.frame(
        tau2      = c(out$tau2),
        iteration = rep(1:nrow(out$tau2), times = M),
        resolution = factor(rep(1:M, each = nrow(out$tau2)))
    ) %>%
        ggplot(aes(x = iteration, y = tau2, group = resolution, color = resolution)) +
        geom_line(alpha = 0.75) +
        scale_color_viridis_d(end = 0.8) +
        ggtitle(TeX("Trace plots for $\\tau2$")) +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\tau^2$"))

    p_sigma2 <- data.frame(sigma2 = c(out$sigma2), iteration = 1:length(out$sigma2)) %>%
        ggplot(aes(x = iteration, y = sigma2)) +
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
        ggplot(aes(x = iteration, y = beta, group = covariate, color = covariate)) +
        geom_line(alpha = 0.75) +
        scale_color_viridis_d(end = 0.8) +
        ggtitle(TeX("Trace plots for $\\beta$")) +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\beta$"))

    return(p_tau2 / p_sigma2 / p_beta)
}
