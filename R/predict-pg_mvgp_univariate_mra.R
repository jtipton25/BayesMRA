#' Bayesian Polya-gamma regression prediction
#' 
#' this function generates predictions from the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param out is a list of MCMC outputs from pgSPLM
#' @param X is a \eqn{n \times p}{n x p} matrix of covariates at the observed locations.
#' @param X_pred is a \eqn{n_{pred} \times p}{n_{pred} x p} matrix of covariates at the locations where predictions are to be made. 
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of locations where observations were taken.
#' @param locs_pred is a \eqn{n_pred \times 2}{n_pred x 2} matrix of locations where predictions are to be made.
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param progress is a logicial input that determines whether to print a progress bar.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' @importFrom stats toeplitz
#' @importFrom fields rdist
#' 
#' @export 

predict_pg_mvgp_univariate_mra <- function(
    out,
    X,
    X_pred,
    locs,
    locs_pred,
    n_cores = 1L,
    progress = TRUE, 
    verbose = FALSE
    # posterior_mean_only = TRUE
) {
    
    # if (posterior_mean_only) {
    #     message("For now, this function generates a posterior predictive draws for the posterior mean only, not the full posterior predictive distribution")
    # }
    

    if (class(out) != "pg_mvgp_univariate_mra")
        stop("THe MCMC object out must be of class pg_mvgp_univariate_mra which is the output of the pg_mvgp_univariate_mra() function.")

    
    ## 
    ## extract the parameters 
    ##
    
    ## polya-gamma variables -- not needed
    beta      <- out$beta
    sigma2    <- out$sigma2
    eta       <- out$eta
    
    gamma     <- out$gamma
    alpha     <- out$alpha
    tau2      <- out$tau2
    Z         <- out$Z
    rho       <- out$rho
    n_samples <- nrow(beta)  
    N         <- nrow(X)
    M         <- ncol(tau2)
    n_time    <- dim(Z)[3]
    n_pred    <- nrow(X_pred)
    J         <- dim(beta)[3] + 1
    MRA       <- out$MRA
    
    if (n_pred > 15000) {
        stop("Number of prediction points must be less than 15000")
    }
    
    ## generate the MRA spatial basis
    W_pred      <- mra_wendland_2d_pred(locs, locs_pred, MRA = MRA)$W
    W_pred      <- do.call(cbind, W_pred)
    
    ## intialize the covariance matric
    Z_pred <- array(0, dim = c(n_samples, n_pred, n_time))
    
    if (progress) {
        message("Beginning Kriging estimates")
        progressBar <- utils::txtProgressBar(style = 3)
    }
    percentage_points <- round((1:100 / 100) * n_samples)   
    
    ## parallelize this later
    
    ## the comments below are to verify that the faster calculations are equivalent
    ## to the slower but simpler mathematical representations
    ## \begin{pmatrix} \boldsymbol{\eta}_o \\ \boldsymbol{\eta}_{oos} \end{pmatrix} & \sim \operatorname{N} \left( \begin{pmatrix} \mathbf{X}_o \\ \mathbf{X}_{oos}  \end{pmatrix} \boldsymbol{\beta}, \boldsymbol{\Sigma}_{time} \otimes \begin{pmatrix} \boldsymbol{\Sigma}_o & \boldsymbol{\Sigma}_{o, oos} \\ \boldsymbol{\Sigma}_{oos, o} & \boldsymbol{\Sigma}_{oos} \end{pmatrix} \right)
    
    for (k in 1:n_samples) {

        X_gamma_pred <- X_pred %*% gamma[k, ]
        W_alpha_pred <- W_pred %*% alpha[k, , ]
        
        Z_pred[k, , ] <- matrix(X_gamma_pred, n_pred, n_time) + W_alpha_pred
        
        if (k %in% percentage_points && progress) {
            utils::setTxtProgressBar(progressBar, k / n_samples)
        }
    }
    
    if (progress) {
        close(progressBar)
    }
    
    return(list(Z = Z_pred))
}

