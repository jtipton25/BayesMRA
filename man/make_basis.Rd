% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_basis.R
\name{make_basis}
\alias{make_basis}
\title{Generate the basis function expansion}
\usage{
make_basis(d, radius, basis_type = "wendland", use_spam = TRUE)
}
\arguments{
\item{d}{A vector of distances}

\item{radius}{A radius}

\item{basis_type}{The basis function type. Currently, only "wendland" basis functions are allowed}

\item{use_spam}{is a boolean flag to determine whether the output is a list of \code{spam::spam} matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class \code{Matrix::dgCMatrix} \code{use_spam = FALSE} (see spam and Matrix packages for details).}
}
\value{
The output of the basis function
}
\description{
Generate the basis function expansion
}
