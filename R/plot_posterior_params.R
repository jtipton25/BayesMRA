#' Posterior distribution plots from the output of `mcmc_mra()` or `mcmc_mra_integrated()`
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param base_size The base size for the plot
#' @param alpha The plot transparency
#'
#' @return Plots of the posterior distribution of model parameters
#' @export
#'
#' @import patchwork
#' @import tidyverse
#' @import latex2exp
#'

plot_posterior_params <- function(out, base_size = 12, alpha = 0.1, ...) {
    if (class(out) != "mcmc_mra")
        stop('out must be of class "mcmc_mra"')

    M <- out$MRA$M
    p_beta <- data.frame(
        beta = c(out$beta),
        parameter = factor(rep(1:ncol(out$beta), each = nrow(out$beta)))
        ) %>%
        ggplot(aes(x = parameter, y = beta)) +
        geom_violin() +
        geom_point(position = "jitter", alpha = alpha) +
        geom_hline(yintercept = 0, col = "red", lty = 2) +
        ggtitle(TeX("Posterior Distribution for $\\beta$")) +
        theme(legend.position = "none") +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\beta$"))


    p_sigma2 <- data.frame(sigma2 = out$sigma2) %>%
        ggplot(aes(y = sigma2, x = "")) +
        geom_violin() +
        geom_point(position = "jitter", alpha = alpha) +
        xlab("") +
        ggtitle(TeX("Posterior Distribution for $\\sigma^2$")) +
        theme_bw(base_size = base_size) +
        ylab(TeX("$\\sigma^2$"))


    p_tau2 <- dat_plot <- data.frame(
        tau2 = c(out$tau2),
        resolution = factor(rep(1:M, each = nrow(out$tau2)))) %>%
        ggplot(aes(x = resolution, y = tau2)) +
        geom_violin() +
        geom_point(position = "jitter", alpha = alpha) +
        geom_hline(yintercept = 0, col = "red", lty = 2) +
        ggtitle(TeX("Posterior Distribution for $\\tau^2$")) +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\tau^2$"))


    return(p_beta / p_sigma2 / p_tau2)
}
