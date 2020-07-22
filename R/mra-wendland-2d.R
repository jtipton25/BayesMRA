#' Code to construct the mutli-resolution sparse basis function representation for fitting spatial processes
#'
#' @param locs The location variables in 2 dimensions over which to construct the basis function representation
#' @param M The number of resolutions
#' @param n_neighbors 
#' @param padding 
#' @param use_spam 
#'
#' @importFrom fields rdist
#' @importFrom spam spam
#' @importFrom Matrix Matrix
#' @return
#' @export
mra_wendland_2d <- function(
    locs,
    M           = 4,
    n_neighbors = 68,
    n_coarse_grid = 10,
    # n_max_fine_grid = 2^12,
    # radius      = 25,
    n_padding     = 5L,
    use_spam    = TRUE
) {
    ## helper function
    wendland_basis <- function(d, radius) {
        if (any(d < 0)) {
            stop("d must be nonnegative")
        }
        if (radius <= 0) {
            stop("radius must be positive")
        }
        d_rad <- d / radius
        return(((1 - d_rad)^6 * (35 * d_rad^2 + 18 * d_rad + 3)) / 3 * (d_rad < 1))
    }

    if (!is_positive_integer(n_padding, 1)) {
        stop("the padding must be a positive integer")
    }

    N <- nrow(locs)
    ## Assign as many gridpoints (approximately) as data
    # n_grid <- ceiling(sqrt(N / 2^(M:1 - 1)))  
    
    ## finest grid is the smaller of n_max_fine_grid
    ## or the largest power of 2 larger than N
    # n_grid <- ceiling(min(n_max_fine_grid, 2^ceiling(log(N, base = 2)))^(0.5) / 2^(M:1 - 1))
    n_grid <- n_coarse_grid * 2^(1:M - 1)
    if (min(n_grid) < 4) {
        stop("There are too many resolutions to form a reliable grid. Reduce M and try again.")
    }

    ## define radius so that each basis function covers approximately 5 neighbors
    area_domain      <- diff(range(locs[, 1])) * diff(range(locs[, 2]))
    density_domain   <- max(n_grid^2) / area_domain
    radius_fine_grid <- sqrt(n_neighbors / (density_domain * base::pi))
    # radius           <- radius_fine_grid * (2^(1:M - 1))^2
    radius           <- radius_fine_grid * (2^(M:1 - 1))
    # radius <- radius / 2^(1:M - 1)

    ## generate a set of grid locations for the Wendland basis
    locs_grid    <- vector(mode = "list", length = M)
    out          <- vector(mode = "list", length = M)
    W            <- vector(mode = "list", length = M)

    for (m in 1:M) {
        ## right now assuming a 2D grid -- can generalize later
        seq_x <- seq(
            range(locs[, 1])[1],
            range(locs[, 1])[2],
            length.out = n_grid[m]
        )
        seq_y <- seq(
            range(locs[, 2])[1],
            range(locs[, 2])[2],
            length.out = n_grid[m]
        )
        delta_x <- seq_x[2] - seq_x[1]
        delta_y <- seq_y[2] - seq_y[1]
        seq_x <- c(
            min(seq_x) - delta_x * (n_padding:1), 
            seq_x,
            max(seq_x) + delta_x * (1:n_padding)
        )
        seq_y <- c(
            min(seq_y) - delta_y * (n_padding:1), 
            seq_y,
            max(seq_y) + delta_y * (1:n_padding)
        )
        
        locs_grid[[m]] <- expand.grid(seq_x, seq_y)
        
        D <- rdist(locs, locs_grid[[m]])
        if (use_spam) {
            ## use the spam sparse matrix package
            W[[m]] <- spam(c(wendland_basis(D, radius[m])), nrow = nrow(D), ncol = ncol(D))
        } else {
            ## use the Matrix sparse matrix package
            W[[m]] <- Matrix(wendland_basis(D, radius[m]), sparse = TRUE)
        }
    }
    return(
        list(
            locs_grid = locs_grid,
            W         = W,
            radius    = radius
        )
    )
}
