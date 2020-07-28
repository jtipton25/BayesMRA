#' Code to construct the mutli-resolution sparse basis function representation for fitting spatial processes
#'
#' @param locs The location variables in 2 dimensions over which to construct the basis function representation in the fitting stage.
#' @param locs_pred The location variables in 2 dimensions over which to construct the basis function representation in the prediction stage.
#' @param MRA The multiresolution basis expansion at the observed locations. This object is the output of \code{mra_wendland-2d()}.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE} (see spam and Matrix packages for details).
#'
#' @importFrom fields rdist
#' @importFrom spam spam
#' @importFrom Matrix Matrix
#' @return A list of objects including the MRA knots locations \code{locs_grid},
#' the Wendland basis representation matrix \code{W_pred} at the prediction locations, and the basis radius \code{radius}
#' @export
mra_wendland_2d_pred <- function(
    locs,
    locs_pred,
    MRA,
    use_spam = TRUE
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

    N      <- nrow(locs)
    N_pred <- nrow(locs_pred)
    M      <- length(MRA$radius)


    W_pred       <- vector(mode = "list", length = M)

    for (m in 1:M) {

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
