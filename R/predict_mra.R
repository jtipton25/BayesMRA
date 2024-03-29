#' Predict for new data using the MRA model
#'
#' @param out The fitted mra object from `mcmc_mra()`
#' @param new_data A list of new data at predicted locations with two inputs: `locs_pred` and `X_pred` representing the prediction locations and covariates, respectively
#'
#' @return A list of prediction objects: `MRA_pred` is the MRA grid information for the predicted locations. `Xbeta_pred` is the predicted fixed effects, `Walpha_pred` is the predicted spatial process, `y_pred` is the predicted response
#' @export
#'
predict_mra <- function(out, new_data) {
    # check new_data input for locs_pred
    # to do
    if (is.null(new_data$locs_pred)) {
        stop("new_data must be a list that contains a n_pred by 2 matrix of locations in the variable locs_pred at which to make predictions")
    }
    if (is.null(new_data$X_pred)) {
        stop("new_data must be a list that contains a n_pred by p matrix of covariates in the variable X_pred which are used to make predictions")
    }
    if (ncol(out$data$X) != ncol(new_data$X_pred))
        stop("The columns of X_pred in the list new_data must be the same number as X used to fit the model")
    if (nrow(new_data$locs_pred) != nrow(new_data$X_pred))
        stop("The number of rows in locs_pred in the list new_data must be equal to the number of rows in X_pred in the list new_data")
    # check new_data input for X_pred
    # to do

    # check if out is class MRA
    if (!(class(out) %in% c("mcmc_mra", "mcmc_mra_integrated")))
        stop('out must be of class "mcmc_mra" or "mcmc_mra_integrated"')

    # generate the prediction MRA grid
    MRA_pred <- mra_wendland_2d_pred(new_data$locs_pred, out$MRA)

    # center and scale the predictors using the same scaling as the model fit
    X_pred <- new_data$X_pred
    if (ncol(X_pred) >= 2) {
        for(i in 2:ncol(X_pred)) {
            X_pred[, i] <- (X_pred[, i] - out$data$mu_X[i-1]) / out$data$sd_X[i-1]
        }
    }

    Xbeta_pred  <- t(X_pred %*% t(out$beta)) * out$data$sd_y + out$data$mu_y
    Walpha_pred <- t(MRA_pred$W_pred %*% t(out$alpha)) * out$data$sd_y
    y_pred      <- Xbeta_pred + Walpha_pred

    out <- list(
        MRA_pred    = MRA_pred,
        Xbeta_pred  = Xbeta_pred,
        Walpha_pred = Walpha_pred,
        y_pred      = y_pred,
        X_pred      = new_data$X_pred,
        locs_pred   = new_data$locs_pred)

    class(out) <- "mcmc_mra_pred"
    return(out)
}
