#' Plot the fitted model parameters vs. simulated model parameters
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param X The simulated regression covariates
#' @param beta The simulated regression coefficients
#' @param W The simulated spatial basis functions
#' @param alpha The simulated spatial basis coefficients
#' @param tile_size The size of the ggplot2 tile
#' @param base_size The base size for the plot
#' @param title The title for the plot
#'
#' @return
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import spam
#' @export
#'

plot_sim_vs_fitted <- function(out, X, beta, W, alpha, base_size = 12, file = NULL, width = 16, height = 9, ...) {

    if (!(class(out) %in% c("mcmc_mra", "mcmc_mra_integrated")))
        stop('out must be of class "mcmc_mra" or "mcmc_mra_integrated"')
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")
    if(!is_numeric_matrix(X, nrow(X), ncol(X)))
        stop("X must be a numeric matrix")
    if(!is_numeric_vector(beta, length(beta)))
        stop("beta must be a numeric vector")
    if (ncol(X) != length(beta))
        stop("X must have the same number of columns as the length of beta")

    if (nrow(W) != nrow(out$MRA$W))
        stop("The basis function matrix W must have the same number of rows as the fitted basis function matrix out$MRA$W")
    if (nrow(X) != nrow(X))
        stop("X must have the same number of rows as the basis function matrix W")

    beta_post <- unscale_beta(out)
    Xbeta_post <- t(X %*% t(beta_post))
    Walpha_post <- t(out$MRA$W %*% t(out$alpha) * out$data$sd_y)
    mu_post <- Xbeta_post + Walpha_post

    dat_plot <- data.frame(
        mean_Xbeta   = apply(Xbeta_post, 2, mean),
        lower_Xbeta  = apply(Xbeta_post, 2, quantile, prob = 0.025),
        upper_Xbeta  = apply(Xbeta_post, 2, quantile, prob = 0.975),
        truth_Xbeta  = X %*% beta,

        mean_Walpha  = apply(Walpha_post, 2, mean),
        lower_Walpha = apply(Walpha_post, 2, quantile, prob = 0.025),
        upper_Walpha = apply(Walpha_post, 2, quantile, prob = 0.975),
        truth_Walpha = (W %*% alpha),

        mean_mu      = apply(mu_post, 2, mean),
        lower_mu     = apply(mu_post, 2, quantile, prob = 0.025),
        upper_mu     = apply(mu_post, 2, quantile, prob = 0.975),
        truth_mu     = X_s %*% beta + (W %*% alpha)[s]
    )

    plot_Xbeta <- dat_plot %>%
        ggplot(aes(x = truth_Xbeta, y = mean_Xbeta)) +
        scale_color_viridis_d(begin = 0, end = 0.8) +
        geom_point(alpha = 0.5) +
        geom_errorbar(aes(ymin = lower_Xbeta, ymax = upper_Xbeta)) +
        geom_abline(intercept = 0, slope = 1, col = "red") +
        ggtitle("Estimated vs. simulated fixed effects") +
        xlab("Simulated fixed effects") +
        ylab("Estimated fixed effects") +
        theme_bw(base_size = base_size)

    plot_Walpha <- dat_plot %>%
        ggplot(aes(x = truth_Walpha, y = mean_Walpha)) +
        scale_color_viridis_d(begin = 0, end = 0.8) +
        geom_point(alpha = 0.5) +
        geom_errorbar(aes(ymin = lower_Walpha, ymax = upper_Walpha)) +
        geom_abline(intercept = 0, slope = 1, col = "red") +
        ggtitle("Estimated vs. simulated spatial process") +
        xlab("Simulated spatial process") +
        ylab("Estimated spatial process") +
        theme_bw(base_size = base_size)

    # plot_mu <- dat_plot %>%
    #     ggplot(aes(x = truth_mu, y = mean_mu)) +
    #     scale_color_viridis_d(begin = 0, end = 0.8) +
    #     geom_point(alpha = 0.5) +
    #     geom_errorbar(aes(ymin = lower_mu, ymax = upper_mu)) +
    #     geom_abline(intercept = 0, slope = 1, col = "red") +
    #     ggtitle("Estimated vs. simulated mean response")
    if (is.null(file)){
        return(plot_Xbeta + plot_Walpha)
    } else {
        ggsave(filename = file,
               plot = plot_Xbeta + plot_Walpha,
               device = "png",
               width = width,
               height = height,
               units = "in")

    }
}
