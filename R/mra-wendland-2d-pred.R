#' Code to construct the mutli-resolution sparse basis function representation for fitting spatial processes
#'
#' @param locs_pred The location variables in 2 dimensions over which to construct the basis function representation in the prediction stage.
#' @param MRA The multi-resolution basis expansion at the observed locations. This object is the output of `mra_wendland-2d()` and is of class "mra_wendland_2d".
#' @param max_points The expected number of pairs less than or equal to the radius. Default is `nrow(locs)` * `num_neighbors`.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (`use_spam = TRUE`) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" `use_spam = FALSE` (see spam and Matrix packages for details).
#'
#' @importFrom fields fields.rdist.near
#' @importFrom Matrix Matrix
#' @import spam
#' @import spam64
#' @return A list of objects including the observation locations `locs`,
#' the MRA knots locations `locs_grid`,
#' the Wendland basis representation matrix `W_pred` at the prediction locations,
#' the basis radius `radius`,
#' the numbers of resolutions `M`,
#' the basis function dimensions `n_dims`,
#' the basis function resolution indices `dims_idx`,
#' the number of expected neighbors in the interior of each grid `n_neighbors`,
#' the number of interior basis functions in one direction `n_coarse_grid`,
#' the number of additional padding basis functions given by `n_padding`,
#' and the setting `use_spam` which determines whether the MRA output uses the `spam` format.
#'
#' @examples
#' set.seed(111)
#' locs <- matrix(runif(20), 10, 2)
#' locs_pred <- matrix(runif(20), 10, 2)
#' MRA <- mra_wendland_2d(locs, M = 2, n_coarse_grid = 4)
#' MRA_pred <- mra_wendland_2d_pred(locs_pred, MRA)
#'
#' ## plot the MRA prediction grid at different resolutions
#' plot_MRA_grid(MRA_pred)
#'
#' @export
mra_wendland_2d_pred <- function(
    locs_pred,
    MRA,
    max_points    = NULL,
    use_spam = TRUE
) {
    n_pred <- nrow(locs_pred)

    # if (class(MRA) != "mra_wendland_2d")
    if (!inherits(MRA, "mra_wendland_2d"))
        stop('MRA must be of class "mra_wendland_2d"')

    if (!is_numeric_matrix(locs_pred, n_pred, 2)) {
        stop("locs_pred must be a numeric matrix with n_pred rows and 2 columns")
    }

    if (!(is.null(max_points) | is_positive_integer(max_points, 1))) {
        stop("max_points must be either NULL or a positive numeric integer")
    }
    if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
        stop("use_spam must be either TRUE or FALSE")
    }


    M      <- MRA$M

    ## Define max_points parameter
    if (is.null(max_points)) {
        max_points <- n_pred * MRA$n_neighbors
    }

    W_pred       <- vector(mode = "list", length = M)

    for (m in 1:M) {

        # D_pred <- rdist(locs_pred, MRA$locs_grid[[m]])
        D_pred <- fields.rdist.near(locs_pred, MRA$locs_grid[[m]],
                                    delta = MRA$radius[m], max.points = max_points)

        D_pred$ra <- wendland_basis(D_pred$ra, MRA$radius[m])
        if (use_spam) {
            ## use the spam sparse matrix package
            # W_pred[[m]] <- spam(c(wendland_basis(D_pred, MRA$radius[m])), nrow = nrow(D_pred), ncol = ncol(D_pred))
            W_pred[[m]] <- spam(D_pred[c("ind", "ra")], nrow = D_pred$da[1], ncol = D_pred$da[2])
        } else {
            stop("The Matrix package is not currently supported")
            ## use the Matrix sparse matrix package
            W_pred[[m]] <- Matrix(wendland_basis(D_pred, MRA$radius[m]), sparse = TRUE)
        }
    }

    W_pred <- do.call(cbind, W_pred)
    out <- list(
        locs          = MRA$locs,
        locs_grid     = MRA$locs_grid,
        W_pred        = W_pred,
        radius        = MRA$radius,
        M             = MRA$M,
        n_dims        = MRA$n_dims,
        dims_idx      = MRA$dims_idx,
        n_neighbors   = MRA$n_neighbors,
        n_coarse_grid = MRA$n_coarse_grid,
        n_padding     = MRA$n_padding,
        use_spam      = MRA$use_spam


    )

    class(out) <- "mra_wendland_2d_pred"

    return(out)
}


