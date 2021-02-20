#' Plot the observed and "true" (holdout) spatial data
#'
#' @param data A dataframe with four columns: Lat, Lon, Observed, Fitted
#' @param base_size The base size for the plot
#' @param file If `file = NULL`, the ggplot object is returned. If `file` is not NULL, an image is saved to the file path specified by `file`
#' @param width If a file path is specified, `width` determines the width of the saved image (in inches)
#' @param height If a file path is specified, `height` determines the height of the saved image (in inches)
#'
#' @return Either a ggplot object of the test data (if `file = NULL`) or a saved image file with no return (`file` is not NULL)
#' @import ggplot2
#' @import patchwork
#' @export
#'
#'
plot_test_data <- function(data, base_size = 12, file = NULL, width = 16, height = 9) {
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is_positive_numeric(base_size, 1))
        stop("base_size must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")

    ## write checks for the data object
    if (!is.data.frame(data))
        stop("data must be a data.frame with variables Observed, Lat, and Lon")
    if (!is_numeric_with_na(data$Observed, nrow(data)))
        stop("The data.frame data must contain a numeric vector named Observed of the observed values")
    if (!is_numeric_vector(data$Truth, nrow(data)))
        stop("The data.frame data must contain a numeric vector named Truth of the full values")
    if (!is_numeric_vector(data$Lat, nrow(data)))
        stop("The data.frame data must contain a numeric vector named Lat of the latitude locations")
    if (!is_numeric_vector(data$Lon, nrow(data)))
        stop("The data.frame data must contain a numeric vector named Lon of the longitude locations")


    plot_obs <- data %>%
        filter(!is.na(.data$Observed)) %>%
        ggplot(aes(x = .data$Lat, y = .data$Lon, fill = .data$Observed)) +
        geom_raster() +
        theme_bw(base_size = base_size) +
        scale_fill_viridis_c(option = "B")
    plot_full <- data %>%
        filter(!is.na(.data$Truth)) %>%
        ggplot(aes(x = .data$Lat, y = .data$Lon, fill = .data$Truth)) +
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
