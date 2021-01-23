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
        stop("new_data must contain a n_pred by 2 matrix of locations at which to make predictions")
    }

    # check new_data input for X_pred
    # to do

    # check if out is class MRA
    if (class(out)!= "mcmc_mra") {
        stop("out must be of class mcmc_mra which is the output of mcmc_mra()")
    }


    # generate the prediction MRA grid
    MRA_pred <- mra_wendland_2d_pred(new_data$locs_pred, out$MRA)

    # center and scale the predictors using the same scaling as the model fit
    X_pred <- new_data$X_pred
    for(i in 2:ncol(X)) {
        X_pred[, i] <- (X_pred[, i] - out$mu_X[i-1]) / out$sd_X[i-1]
    }

    Xbeta_pred  <- t(X_pred %*% t(out$beta)) * out$sd_y + out$mu_y
    Walpha_pred <- t(MRA_pred$W_pred %*% t(out$alpha)) * out$sd_y
    y_pred      <- Xbeta_pred + Walpha_pred

    for(i in 2:ncol(X)) {
        X_pred[, i] <- X_pred[, i] * out$sd_X[i-1] + out$mu_X[i-1]
    }

    out <- list(
        MRA_pred    = MRA_pred,
        Xbeta_pred  = Xbeta_pred,
        Walpha_pred = Walpha_pred,
        y_pred      = y_pred,
        X_pred = X_pred,
        locs_pred = locs_pred)

    class(out) <- "mcmc_mra_pred"
    return(out)
}
