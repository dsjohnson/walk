% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_walk_data.R
\name{make_walk_data}
\alias{make_walk_data}
\title{Produce design data for use in fitting MMPP movement models}
\usage{
make_walk_data(
  proc_data,
  cell_data,
  grad = NULL,
  rast_mask = NULL,
  directions = "rook",
  debug = 0,
  ...
)
}
\arguments{
\item{proc_data}{A sparse matrix with rows corresponding to time steps and columns corresponding to cells.
The entries are probabilities that the animal is located in the corresponding cell and time.}

\item{cell_data}{A \code{SpatRaster} object from the \code{terra} package.}

\item{grad}{A character vector of names of \code{cell_data} layers for which gradient covariates will be constructed.}

\item{rast_mask}{Raster mask for inaccessible cells when \code{cell_data} is of type \code{SpatRaster} from the \code{terra} package. This is ignored
if \code{cell_data} is an \code{POLYGON} data frame from the \code{sf} package.}

\item{directions}{Neighborhood structure, one of \code{"rook"} or \code{"queen"}.}

\item{debug}{Debugging level: 1-3 mainly for package developers.}

\item{...}{Ignored arguments.}
}
\description{
Produce design data for use in fitting MMPP movement models
}
