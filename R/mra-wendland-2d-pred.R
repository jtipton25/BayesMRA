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
mra_wendland_2d_pred <- function(
    locs,
    locs_pred,
    MRA,
    # M           = 4,
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

    # if (padding < 0 | padding > 1) {
    #     stop("the padding must be a number between 0 and 1")
    # }

    N      <- nrow(locs)
    N_pred <- nrow(locs_pred)
    M      <- length(MRA$radius)
    
    # ## Assign as many gridpoints (approximately) as data
    # n_grid <- ceiling(sqrt(N / 2^(M:1 - 1)))
    # if (min(n_grid) < 4) {
    #     stop("There are too many resolutions to form a reliable grid. Reduce M and try again.")
    # }
    # 
    # ## define radius so that each basis function covers approximately 5 neighbors
    # area_domain      <- diff(range(locs[, 1])) * diff(range(locs[, 2]))
    # density_domain   <- N / area_domain
    # radius_fine_grid <- sqrt(n_neighbors / (density_domain * base::pi))
    # # radius           <- radius_fine_grid * (2^(1:M - 1))^2
    # radius           <- radius_fine_grid * (2^(M:1 - 1))
    # # radius <- radius / 2^(1:M - 1)
    # 
    # ## generate a set of grid locations for the Wendland basis
    # locs_grid    <- vector(mode = "list", length = M)
    # out          <- vector(mode = "list", length = M)
    W_pred       <- vector(mode = "list", length = M)

    for (m in 1:M) {
        ## right now assuming a 2D grid -- can generalize later
        # locs_grid[[m]] <- expand.grid(
        #     seq(
        #         range(locs[, 1])[1] - padding * diff(range(locs[, 1])),
        #         range(locs[, 1])[2] + padding * diff(range(locs[, 1])),
        #         length.out = n_grid[m]
        #     ),
        #     seq(
        #         range(locs[, 2])[1] - padding * diff(range(locs[, 2])),
        #         range(locs[, 2])[2] + padding * diff(range(locs[, 2])),
        #         length.out = n_grid[m]
        #     )
        # )
        D_pred <- rdist(locs_pred, MRA$locs_grid[[m]])
        if (use_spam) {
            ## use the spam sparse matrix package
            W_pred[[m]] <- spam(c(wendland_basis(D_pred, MRA$radius[m])), nrow = nrow(D_pred), ncol = ncol(D_pred))
        } else {
            ## use the Matrix sparse matrix package
            W_pred[[m]] <- Matrix(wendland_basis(D_pred, MRA$radius[m]), sparse = TRUE)
        }
    }
    return(
        list(
            locs_grid = MRA$locs_grid,
            W_pred    = W_pred,
            radius    = MRA$radius
        )
    )
}
