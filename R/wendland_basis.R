#' calculate the Wendland basis function
#'
#' @param d The distance over which to calculate the Wendland basis
#' @param radius The effective radius over which the Wendland basis is defined
#'
#' @return The output of the Wendland basis applied to the distance `d` for a given radius `radius`.
#' @export
#'
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
