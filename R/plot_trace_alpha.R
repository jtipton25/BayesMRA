#' Trace plots from the output of `mcmc_mra()` or `mcmc_mra_integrated()`
#'
#' @param out The output from `mcmc_mra()` or `mcmc_mra_integrated()`
#' @param base_size The base size for the plot
#' @param alpha The plot transparency
#' @param n_per_res The number of alpha parameters to sample per resolutoin
#'
#' @return Trace plots of the alpha parameter
#' @export
#'
#' @import patchwork
#' @import tidyverse
#' @import latex2exp
#'

plot_trace_alpha <- function(out, base_size = 12, alpha = 0.1, n_per_res = 10, ...) {
    if (class(out) != "mcmc_mra")
        stop('out must be of class "mcmc_mra"')

    M        <- out$MRA$M
    n_dims   <- out$MRA$n_dims
    dims_idx <- out$MRA$dims_idx

    alpha_post <- out$alpha
    dimnames(alpha_post) <- list(
        iteration = 1:nrow(alpha_post),
        par_id    = 1:ncol(alpha_post))
    dat_alpha <- as.data.frame.table(alpha_post, responseName = "alpha")
    dat_id <- data.frame(
        par_id     = factor(1:sum(n_dims)),
        knot       = factor(unlist(sapply(1:M, function(m) 1:n_dims[m]))),
        resolution = factor(paste("resolution", dims_idx))
    )

    # programmatically subsample over the resolutions
    plot_idx <- sapply(1:M, function(m) sample(n_dims[m], n_per_res))
    for (m in 1:M) {
        assign(paste("filter", m, sep = "_"), parse(text = paste0('(resolution == "resolution ', m, '" & knot %in% plot_idx[, ', m, '])')))
    }

    p_alpha <- dat_alpha %>%
        left_join(dat_id, by = "par_id") %>%
        filter(
            eval(parse(text = paste(paste0("eval(filter_", 1:M, ")"), collapse = " | ")))
        ) %>%
        ggplot(aes(x = iteration, y = alpha, group = knot, color = knot)) +
        geom_line(alpha = alpha) +
        facet_wrap( ~ resolution, ncol = 1) +
        scale_color_viridis_d(end = 0.8) +
        ggtitle(TeX("Trace plots for $\\alpha$")) +
        theme_bw(base_size = base_size) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\alpha$"))

    return(p_alpha)
}
