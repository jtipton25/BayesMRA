#' Plot the sparsity structure of a spam matrix object using ggplot2
#'
#' Plot the sparsity of a spam matrix object using ggplot2
#'
#' @param x The sparse matrix which must be of class "spam"
#' @param tile_size The size of the ggplot2 tile
#' @param base_size The base size for the plot
#' @param title The title for the plot
#'
#' @return A ggplot of the sparse matrix structure of a matrix of class "spam"
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import spam
#' @export
#'
#'
plot_sparse_matrix <- function(x, tile_size = 1, base_size = 12, title = "") {
    if (!inherits(x, "spam"))
        stop('x must be of class "spam"')

    # extract the sparsity structure
    nrow <- x@dimension[1]
    ncol <- x@dimension[2]

    data.frame(
        column   = x@colindices,
        row      = rep.int((1:nrow), diff(x@rowpointers)),
        sparsity = 1) %>%
        drop_na() %>%
        # plot the sparse matrix
        ggplot(aes(y = .data$row, x = .data$column, fill = .data$sparsity)) +
        geom_tile(aes(width = tile_size, height = tile_size)) +
        scale_fill_gradient(low = "black", high = "black") +
        scale_y_reverse() +
        theme_bw(base_size = base_size) +
        theme(legend.position = "none") +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5))
}
