#' Partition observations into MRA grid
#'
#' Partition the observation locations into MRA grids for use in \code{\link{mcmc_mra_vecchia()}}
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param M The number of resolutions.
#'
#' @return A data.frame of partition ids for each location
#' @export
#'
partition_observations <- function(locs, M) {
    # add in error checks for locs
    if (M < 0)
        stop ("M must be a positive integer")

    N <- nrow(locs)
    range_x <- range(locs[, 1])
    range_y <- range(locs[, 2])
    n_grids <- 4^(M-1)
    cuts_x <- seq(range_x[1], range_x[2], length.out = 2^(M-1)+1)
    cuts_y <- seq(range_y[1], range_y[2], length.out = 2^(M-1)+1)

    # programmatically define the grids
    grid_labels <- matrix(1:n_grids, 2^(M-1), 2^(M-1), byrow = TRUE)
    dimnames(grid_labels) <- list(row_id = 1:(2^(M-1)), col_id = 1:(2^(M-1)))
    dat_grid_labels <- as.data.frame.table(grid_labels, responseName = "grid_id") %>%
        mutate(row_id = as.numeric(row_id), col_id = as.numeric(col_id))

    grid_idx <- NULL

    if (M == 1) {
        grid_idx <- rep(1, N)
        return(grid_idx)
    }

    row_id <- as.numeric(cut(locs[, 1], breaks = cuts_x, include.lowest = TRUE))
    col_id <- as.numeric(cut(locs[, 2], breaks = cuts_x, include.lowest = TRUE))
    grid_idx <- data.frame(row_id = row_id, col_id = col_id) %>%
        left_join(dat_grid_labels, by = c("row_id", "col_id")) %>%
        pull(grid_id)

    return(grid_idx)
}
