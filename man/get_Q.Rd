% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Q_functions.R
\name{get_Q}
\alias{get_Q}
\title{Get movement transition matrix from fitted moveMMP object}
\usage{
get_Q(fit, sparse = TRUE)
}
\arguments{
\item{fit}{A fitted model object from \code{\link[walk]{fit_ctmc}}.}

\item{sparse}{Logical. Should the matrix be returned in a sparse format from the \code{Matrix}
package. Defaults to \code{sparse = TRUE}.}
}
\description{
Get movement transition matrix from fitted moveMMP object
}
\author{
Devin S. Johnson
}
