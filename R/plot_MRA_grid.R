#' Plot the MRA spatial grid
#'
#' Plot the MRA  spatial grid from the output of `mra_wendland_2d()`
#'
#' @param MRA The output from `mra_wendland_2d()`
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
plot_MRA_grid <- function(MRA, base_size = 12, file = NULL, width = 16, height = 9) {

    if (!inherits(MRA, c("mra_wendland_2d", "mra_wendland_2d_pred")))
        stop('MRA must be of class "mra_wendland_2d"')
    if (!is_positive_numeric(width, 1))
        stop("width must be a positive number")
    if (!is_positive_numeric(height, 1))
        stop("height must be a positive number")
    if (!is_positive_numeric(base_size, 1))
        stop("base_size must be a positive number")
    if (!is.null(file) & !is.character(file))
        stop("file must be a character string")



    M <- MRA$M


    p_plot <- data.frame(
        x = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 1])),
        y = unlist(sapply(1:M, function(i) MRA$locs_grid[[i]][, 2])),
        res = factor(unlist(sapply(1:M, function(i) rep(paste("resolution", i), each = nrow(MRA$locs_grid[[i]])))))
    ) %>%
        ggplot(aes(x = .data$x, y = .data$y, color = .data$res)) +
        geom_point(alpha = 0.75, size = 1) +
        coord_cartesian(
            xlim = c(
                min(MRA$locs_grid[[1]][, 1]),
                max(MRA$locs_grid[[1]][, 1])
            ),
            ylim = c(
                min(MRA$locs_grid[[1]][, 2]),
                max(MRA$locs_grid[[1]][, 2])
            )
        ) +
        geom_rect(data = NULL,
                  aes(xmin = min(MRA$locs[, 1]), xmax = max(MRA$locs[, 1]), ymin = min(MRA$locs[, 2]), ymax = max(MRA$locs[, 2])),
                  fill = NA, color = "black", alpha = 0.5
        ) +
        ggtitle("MRA grid points") +
        scale_color_viridis_d(begin = 0, end = 0.8) +
        xlab("Longitude") +
        ylab("Latitude") +
        theme_bw(base_size = base_size) +
        labs(color = "Resolution") +
        guides(alpha = "none")

    if (is.null(file)) {
        return(p_plot)
    } else {
        ggsave(filename = file,
               plot = p_plot,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }

}


