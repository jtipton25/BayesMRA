#' Title
#'
#' @param Q_alpha a list of length M composed of matrices that are the correlation structure of the CAR prior on beta
#' @param tau2 a vector of length M that contains the CAR prior precisions
#' @param use_spam a boolean that determines if the output matrix is of class "spam" (\code{use_spam = TRUE}) or of class "dgCMatrix" (\code{use_spam = FALSE}; see Matrix package for details)
#'
#' @importFrom spam bdiag.spam
#' @importFrom Matrix bdiag
#'
#' @return
#' @export
make_Q_alpha_tau2 <- function(Q_alpha, tau2, use_spam = TRUE) {

    if (length(Q_alpha) != length(tau2))
        stop("Q_alpha must be a list of length M and tau2 must be a vector of length M")
    M <- length(Q_alpha)
    Q_alpha_tau2 <- vector(mode = "list", length = M)
    for (m in 1:M) {
        Q_alpha_tau2[[m]] <- Q_alpha[[m]] * tau2[m]
    }
    if (use_spam) {
        ## use the spam package
        Q_alpha_tau2 <- do.call(bdiag.spam, Q_alpha_tau2)
    } else {
        ## use the Matrix package
        Q_alpha_tau2 <- do.call(bdiag, Q_alpha_tau2)
    }
    return(Q_alpha_tau2)
}
