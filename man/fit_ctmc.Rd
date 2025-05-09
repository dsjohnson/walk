% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_ctmc.R
\name{fit_ctmc}
\alias{fit_ctmc}
\title{Fit CTMC movement model to telemetry data}
\usage{
fit_ctmc(
  walk_data,
  model_parameters = ctmc_control(),
  pen_fun = NULL,
  hessian = TRUE,
  reals = FALSE,
  start = NULL,
  method = "nlminb",
  fit = TRUE,
  eq_prec = 1e-08,
  check_rho = NULL,
  debug = 0,
  ...
)
}
\arguments{
\item{walk_data}{A design data list produced by the function \code{\link{make_walk_data}}.}

\item{model_parameters}{A named list giving specifications for the model. See \code{\link[walk]{ctmc_control}}.}

\item{pen_fun}{An optional penalty function. Should be on the scale of a log-prior distribution.}

\item{hessian}{Logical. Should the Hessian matrix be calculated to obtain the parameter
variance-covariance matrix.}

\item{reals}{Calculate real values for expected residency, cell transition probabilities, and
outlier proportion for observations.}

\item{start}{Optional starting values for the parameter must be a list of the
form \code{list(beta_l=c(), beta_q_r=c(), beta_q_r=c())}.}

\item{method}{Optimization method. See \code{\link[optimx]{optimr}}}

\item{fit}{Logical. Should the likelihood be optimized?}

\item{eq_prec}{Error rate of matrix exponential calculation. Defaults to \code{1.0e-8}. This is
a generous value. If the model is running slow, you can try reducing it to, say, \code{1.0e-4}.}

\item{check_rho}{Check if rho is too big for uniformitazation calculation of exp\{Qt\}. Value is the size for which rho is too big, causing likelihood calculation error.}

\item{debug}{Integer from 1-4. Opens browser() at various points in the function call. Mostly for
package developers.}

\item{...}{Additional arguments passed to the optimization function
\code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}.}
}
\description{
Fit CTMC movement model to telemetry data
}
\author{
Devin S. Johnson
}
