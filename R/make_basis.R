#' Generate the basis function expansion
#'
#' @param d A vector of distances
#' @param radius A radius
#' @param basis_type The basis function type. Currently, only "wendland" basis functions are allowed
#' @param use_spam is a boolean flag to determine whether the output is a list of `spam::spam` matrix objects (`use_spam = TRUE`) or a an \eqn{n \times n}{n x n} sparse Matrix of class `Matrix::dgCMatrix` `use_spam = FALSE` (see spam and Matrix packages for details).
#'
#' @return The output of the basis function
#' @export
#'
make_basis <- function(d, radius, basis_type = "wendland", use_spam = TRUE) {
    if (!(basis_type %in% c("wendland")))
        stop('The only currently supported basis functions are "wendland"')
    if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
        stop("use_spam must be either TRUE or FALSE")
    }

    if (basis_type == "wendland") {
        return(wendland_basis(d, radius))
    }


    # if (basis_type == "") {
    #
    # }
}
