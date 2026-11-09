#' @title Convert `simLangevin`  or `dataLangevin' object to a `SpatRaster` stack from the `terra` package.
#' @param data A `simLangevin` or `dataLangevin' object from the `langevinSSM` package, or a list
#' of them (e.g. one per individual).
#' @param cell_data A `\link[terra]{SpatRaster}` stack of covariates that will be used in CTMC 
#' movement modeling. Cells with `NA` values will be considered areas to the animal cannot travel, i.e., 
#' likelihood surfaces will be `0` for those cells.
#' @param return_type Type of object returned. One if `"data.frame"`, `"sparse"` (sparse matrix), 
#' `"matrix_df"` (matrix form of `"data.frame"`), or `"dense"` (dense matrix).
#' @param max_err The maximum error in meters. If unspecified it will be set to
#' 4 times the maximum error standard deviation as determined by the UERE and HDOP
#' of the telemetry data.
#' @param trunc The smallest probability value that is considered to be > 0. Defaults to 1.0e-8.
#' @param  return_inadmissible Logical. If there are locations outside the admissible areas of movement return the observations and raster cell numbers.
#' @param nthreads Number of OpenMP threads used for the cell probability calculations. Defaults to 1.
#' Has no effect if `walk` was built without OpenMP (a warning is given).
#' @param neighborhood How the error neighborhood of each location is found. `"disc"` (default,
#' fastest) uses every raster cell that touches an exact circle of radius `max_err`, computed
#' arithmetically from the raster geometry. `"buffer"` uses `terra::buffer()` +
#' `terra::cells(touches=TRUE)` as in earlier versions. The polygonised buffer sits slightly
#' inside the true circle, so it can miss a few edge cells; otherwise the results are identical.
#' Cell probabilities are computed in C++ either way.
#' @param id Name of the column in `data` identifying individuals (default `"id"`). If the column
#' is absent, list elements are named by their list names (or position), and a single object is
#' individual `"1"`. Time differences are computed within individual and every individual
#' starts a new, independent Markov chain in the likelihood. Rows of one individual must be
#' contiguous and in time order. Set to `NULL` to ignore any `id` column.
#' @param time_unit Only used when the `date` column is POSIXct (or Date): unit of the time
#' differences `dt` (any unit accepted by \code{\link[base]{difftime}}). Defaults to the
#' `"time.unit"` attribute of `data` (set by `langevinSSM::formatData()`), or `"hours"` if there
#' is none. Numeric dates are used as they are.
#' @param max_gap Optional. Gaps longer than `max_gap` (in the units of `dt`) start a new,
#' independent segment (the chain is restarted from `delta`).
#' @details This function takes the error information in `data` to 
#' produce a likelihood surface over the `SpatRaster` defined by 
#' `cell_data` for each location. This can then be passed to 
#' CTMC HMM fitting functions. The `times` data frame returned with the likelihood has
#' columns `obs`, `id`, `seg`, `timestamp` and `dt`; `dt` is `0` at the first observation of
#' every segment (individual or `max_gap` burst). Several outputs can be combined with
#' \code{\link{bind_proc_data}}.
#' @author Devin S. Johnson
#' @importFrom Matrix sparseMatrix
#' @importFrom terra vect buffer cells res values
#' @importFrom methods as
#' @export
proc_Langevin <- function(data, cell_data, return_type="sparse", max_err=NULL, trunc=1.0e-8,
                          return_inadmissible=FALSE, nthreads=1, neighborhood=c("disc","buffer"),
                          id="id", time_unit=NULL, max_gap=NULL){
  neighborhood <- match.arg(neighborhood)
  proc_check_threads(nthreads)
  # Check arguments
  dl <- proc_as_list(data, c("simLangevin", "dataLangevin"),
                     "a 'simLangevin' or 'dataLangevin' object from the langevinSSM package")
  if(!inherits(cell_data, "SpatRaster")) stop("'cell_data' must be a terra::SpatRaster object!")
  xl <- dl$x
  
  # Individual / segment structure
  idl <- proc_ids(xl, id, dl$is_list)
  if(any(vapply(xl, function(d) is.null(d$date), TRUE))) stop("'data' has no 'date' column.")
  ts <- do.call(c, lapply(xl, function(d) d$date))
  is_time <- inherits(ts, c("POSIXt", "Date"))
  if(is_time && is.null(time_unit)){
    time_unit <- attr(xl[[1]], "time.unit", exact=TRUE)
    if(is.null(time_unit)) time_unit <- "hours"
  }
  times <- proc_times(unlist(idl, use.names=FALSE), ts, time_unit, max_gap)
  if(is_time) attr(times, "time_unit") <- time_unit
  
  # Stack only the columns used below (list elements may carry different extra columns)
  cols <- c("x", "y", "smaj", "smin", "eor", "x.err", "y.err")
  d <- do.call(rbind, lapply(xl, function(z){
    z <- as.data.frame(z)
    for(cc in setdiff(cols, names(z))) z[[cc]] <- NA_real_
    z[, cols]
  }))
  n <- nrow(d)
  
  # Evaluate error ellipse covariance
  cov_df <- argos_diag_to_cov(d$smaj, d$smin, 180*d$eor/pi)
  
  # x and y standard deviation errors (assumes no correlation, i.e., rho_o = 0)
  sderrInd <- which(!is.na(d$x.err) & !is.na(d$y.err))
  if(length(sderrInd)){
    rho_o <- 0
    v1 <- d$x.err[sderrInd]^2
    v2 <- d$y.err[sderrInd]^2
    c12 <- d$x.err[sderrInd] * d$y.err[sderrInd] * rho_o
    cov_df$cov.x.x[sderrInd] <- v1
    cov_df$cov.y.y[sderrInd] <- v2
    cov_df$cov.x.y[sderrInd] <- c12
    cov_df$cov.check[sderrInd] <-  (v1*v2-c12^2) > 0
  }
  
  if(is.null(max_err)) {
    max_err <- 4 * sqrt(pmax(cov_df$cov.x.x, cov_df$cov.y.y))
  } else if(length(max_err) == 1) {
    max_err <- rep(max_err, n)
  }
  
  proc_lik(as.numeric(d$x), as.numeric(d$y), cov_df$cov.x.x, cov_df$cov.y.y, cov_df$cov.x.y,
           max_err, cell_data, times, return_type, trunc, return_inadmissible, nthreads,
           neighborhood)
}

#' @title Combine processed telemetry from several individuals or data sources
#' @description Stack the output of several calls to \code{\link{proc_telem}} and/or
#' \code{\link{proc_Langevin}} (e.g. one per individual) into a single object that
#' \code{\link{make_walk_data}} and \code{\link{fit_ctmc}} treat as a joint,
#' multi-individual data set.
#' @param ... Objects returned by `proc_telem(..., return_type="sparse")` or
#' `proc_Langevin(..., return_type="sparse")`, in any combination, or a single list of them.
#' @param ids Optional character vector of individual names, one per input. Overrides the
#' `id` column of each input's `times` (each input is then treated as one individual).
#' @details All inputs must have been built on the same `cell_data` and use the same time
#' units. Segments are renumbered so they stay distinct, and an individual may not appear in
#' more than one input (use `max_gap`, not separate inputs, to split one animal's track).
#' @return A list with `L` (sparse) and `times` (columns `obs`, `id`, `seg`, `timestamp`, `dt`),
#' the same form as the output of \code{\link{proc_telem}} and \code{\link{proc_Langevin}}.
#' @export
bind_proc_data <- function(..., ids=NULL){
  x <- list(...)
  if(length(x)==1 && is.list(x[[1]]) && is.null(x[[1]]$L)) x <- x[[1]]
  if(length(x)==0) stop("No inputs.")
  if(!all(vapply(x, function(z) is.list(z) && (inherits(z$L, "Matrix") || is.matrix(z$L)) && is.data.frame(z$times), TRUE)))
    stop("All inputs must be output of proc_telem() or proc_Langevin() with return_type = \"sparse\" or \"dense\".")
  if(!is.null(ids) && length(ids) != length(x)) stop("'ids' must have one entry per input.")
  nc <- vapply(x, function(z) ncol(z$L), 1)
  if(length(unique(nc))!=1) stop("All inputs must have been processed with the same 'cell_data' (different numbers of cells).")
  tu <- unique(vapply(x, function(z){ u <- attr(z$times, "time_unit"); if(is.null(u)) NA_character_ else u }, ""))
  if(length(tu) != 1) stop("All inputs must have the same 'time_unit' attribute.")
  seg_off <- 0L
  tl <- vector("list", length(x))
  for(k in seq_along(x)){
    tk <- x[[k]]$times
    if(nrow(tk) != nrow(x[[k]]$L)) stop("Input ", k, ": 'times' and 'L' have different numbers of rows.")
    if(is.null(tk$seg)) tk$seg <- 1L
    if(is.null(tk$id)) tk$id <- as.character(k)
    if(!is.null(ids)) tk$id <- as.character(ids[k])
    tk$id <- as.character(tk$id)
    tk$dt[1] <- 0
    tk$seg <- as.integer(factor(tk$seg, levels=unique(tk$seg))) + seg_off
    seg_off <- max(tk$seg)
    tl[[k]] <- tk[, c("obs", "id", "seg", "timestamp", "dt")]
  }
  times <- do.call(rbind, tl)
  rownames(times) <- NULL
  id_in <- unlist(lapply(tl, function(tk) unique(tk$id)), use.names=FALSE)
  if(any(duplicated(id_in))) stop("The same individual id appears in more than one input: ",
                                  paste(unique(id_in[duplicated(id_in)]), collapse=", "),
                                  ". Use 'ids' to rename, or process each animal in one call (use 'max_gap' to split a track).")
  times$obs <- seq_len(nrow(times))
  if(!is.na(tu)) attr(times, "time_unit") <- tu
  L <- do.call(rbind, lapply(x, function(z) as(z$L, "CsparseMatrix")))
  list(L=as(as(L, "generalMatrix"), "dgCMatrix"), times=times)
}
