% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize-basis.R
\name{calc_W_star}
\alias{calc_W_star}
\title{Calculate normalized basis function}
\usage{
calc_W_star(W, omega, tau2)
}
\arguments{
\item{W}{The sparse basis function matrix}

\item{omega}{The normalizing weights for the basis functions which is the output of \code{\link[=calc_omega_norm]{calc_omega_norm()}}}

\item{tau2}{The basis function resolution precision}
}
\value{
A normalized basis function matrix
}
\description{
Calculate the normalized basis function matrix so that all basis functions have the same marginal variance
}
