#' Construct the constraint matrices for the spatial random effects alpha
#'
#' @param MRA An MRA object. Usually contained in the output from `mcmc_mra()` or `mcmc_mra_integrated()` or for developers, the MRA object can be created directly from the function `mra_wendland_2d()`
#' @param constraint What constraint should be applied to the spatial process? Options include no constraint (`constraint = "unconstrained"`), a constraint so the entire MRA process sums to 0 (`constraint = "overall"`), a constraint so that each of the M levels of the MRA process sum to 0 (`constraint = "resolution"`), or whether the predicte process must sum to 0 (`constraint = "predicted"`). Note: `constraint = "predicted"` is NOT currently supported.
#'
#' @return
#' @export
#'
make_constraint <- function(MRA, constraint = "unconstrained") {

    ## check the constraints on alpha
    if (!(constraint %in% c("unconstrained", "overall", "resolution", "predicted"))) {
        stop('constraint must be either "unconstrained", "overall", "resolution", or "predicted"')
    }
    if (constraint == "predicted") {
        stop('constraint = "predicted" is not currently supported -- developer note: add W_pred to function call to enable this in future results')
    }
    A_constraint <- NULL
    a_constraint <- NULL
    if (constraint == "overall") {
        ## overall constraint on random effect
        A_constraint <- rep(1, nrow(MRA$W)) %*% MRA$W
        a_constraint <- 0
    } else if (constraint == "resolution") {
        ## resolution level constraint on random effect W_m %*% alpha_m
        A_constraint_list <- vector(mode = "list", length = MRA$M)
        a_constraint <- NULL
        A_constraint_tmp <- sapply(1:MRA$M, function(i){
            tmp <- rep(1, nrow(MRA$W)) %*% MRA$W[, MRA$dims_idx == i]
            return(tmp)
        })
        A_constraint <- matrix(0, MRA$M, sum(MRA$n_dims))
        for (i in 1:MRA$M) {
            A_constraint[i, MRA$dims_idx == i] <- A_constraint_tmp[[i]]
        }

        a_constraint <- rep(0, MRA$M)
    }
    return(list(A_constraint = A_constraint, a_constraint = a_constraint))
}
