#' Plot the observed and "true" (holdout) spatial data
#'
#' @param data A dataframe with four columns: Lat, Lon, Observed, Fitted
#' @param base_size The base size for the plot
#' @param ... Other ggplot2 options
#'
#' @return
#' @import ggplot2
#' @import patchwork
#' @export
#'
#'
plot_test_data <- function(data, base_size = 12, file = NULL, width = 16, height = 9, ...) {
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")

    plot_obs <- data %>%
        filter(!is.na(Observed)) %>%
        ggplot(aes(x = Lat, y = Lon, fill = Observed)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B")
    plot_full <- data %>%
        filter(!is.na(Truth)) %>%
        ggplot(aes(x = Lat, y = Lon, fill = Truth)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B")

    if (is.null(file)) {
        return(plot_obs + plot_full)
    } else {
        ggsave(filename = file,
               plot = plot_obs + plot_full,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }
}
