#' Code for the elliptical slice sampler
#'
#' @param current The current state of the parameter of interest in the 
#' MCMC algorithm.
#' @param prior A sample from the mean zero, Gaussian prior of the 
#' parameter of interest.
#' @param prior_mean The mean, if nonzero, of the Gaussian prior 
#' distribution for the parameter of interest.
#' @param pars The set of parameters needed to evaluation the 
#' log-likelihood function for the parameter of interest. At a minimum, 
#' the pars object must contain the variables \code{parameter} which is 
#' the parameter to be updated, \code{par_name} which is a string name
#' of the variable to update, \code{verbose} which is a boolean option for verbose output,
#' and the ess performance tracking variables \code{num_calls_ess}, 
#' \code{num_max_contract}, and \code{num_unusual_contract}.
#' @param log_like_fun The log-likelihood function for the parameter of 
#' interest.
#' @param update_pars_fun The function to update the parameters 
#' \code{pars}.
#' @param max_iter The maximum number of iterations to perform for the 
#' elliptical slice sampler.
#'
#' @return
#' @export
#'

ess <- function(current, prior, prior_mean, pars, log_like_fun, update_pars_fun, max_iter = 500) {
    
    pars$num_calls_ess <- pars$num_calls_ess + 1
    
    ## calculate log likelihood of current value
    
    current_log_like <- log_like_fun(pars)
    
    hh <- log(runif(1, 0.0, 1.0)) + current_log_like
    
    ## Setup a bracket and pick a first proposal
    ## Bracket whole ellipse with both edges at first proposed point
    phi_angle <- runif(1, 0.0, 1.0) * 2.0 * base::pi
    phi_angle_min <- phi_angle - 2.0 * base::pi
    phi_angle_max <- phi_angle
    
    updated_pars <- pars
    updated_pars$parameter <- current
    
    test <- TRUE
    count <- 0
    
    ## Slice sampling loop
    
    while (test) {
        
        if (count >= max_iter) {
            ## ESS failed to generate an acceptable proposal within the 
            ## maximum number of iterations. Return the initial values
            if (pars$verbose)
                message("maximum number of contractions in ess reached, if this happens often, try increasing max_iter")
            pars$num_max_contract <- pars$num_max_contract + 1
            pars$parameter        <- current
            return(pars)
        }
        
        ## compute proposal for angle difference and check to see if it is on the slice
        ## adjust X so that it is mean 0 for the ellipse to be valid
        
        proposal     <- (current - prior_mean) * cos(phi_angle) + prior * sin(phi_angle) + prior_mean
        updated_pars <- update_pars_fun(proposal, pars)
        
        ## calculate log likelihood of proposed value
        
        proposal_log_like <- log_like_fun(updated_pars)
        
        if (!is.finite(proposal_log_like)) {
            ## The log-likelihood for the proposal is invalid 
            ## this happens sometimes when the proposal is far from the posterior
            ## and if it happens occasionally it isn't a problem
            if (pars$verbose)
                message("An unusual log-likelihood was evaluated for ", pars$par_name, " using the elliptical slice sampler. If this warning is rare, it should be safe to ignore")
            pars$num_unusual_ll         <- pars$num_unusual_ll + 1
            updated_pars$num_unusual_ll <- updated_pars$num_unusual_ll + 1
            ## Skip this iteration and propose new angle difference
            phi_angle <- runif(1, 0.0, 1.0) * (phi_angle_max - phi_angle_min) + phi_angle_min
            count     <- count + 1
        } else if (proposal_log_like > hh) {
            ## proposal is on the slice
            updated_pars$parameter <- proposal
            test <- FALSE
        } else if (phi_angle > 0.0) {
            ## proposal is not on the slice and adjust the proposal window
            phi_angle_max <- phi_angle
        } else if (phi_angle < 0.0) {
            ## proposal is not on the slice and adjust the proposal window
            phi_angle_min <- phi_angle
        } else {
            warning("Bug - ESS for eta_star shrunk to current position \n")
            test <- FALSE
        }
        ## Propose new angle difference
        phi_angle <- runif(1, 0.0, 1.0) * (phi_angle_max - phi_angle_min) + phi_angle_min
        count <- count + 1
        
    }
    return(updated_pars)
}