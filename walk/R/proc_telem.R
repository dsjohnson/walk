
#' @title Convert `ctmm` telemetry data object to a `SpatRaster` stack from the `terra` package.
#' @param data A `telemetry` object from the `ctmm` package, or a list of them (one per
#' individual), as returned by `ctmm::as.telemetry()` when the data contain several animals.
#' @param cell_data A `\link[terra]{SpatRaster}` stack of covariates that will be used in CTMC 
#' movement modeling. Cells with `NA` values will be considered areas to the animal cannot travel, i.e., 
#' likelihood surfaces will be `0` for those cells.
#' @param time_unit Unit of measurement for numeric conversion of the POSIX timestamps 
#' in the telemetry data. Defaults to `"hours"`. Can specify, e.g., 
#'  `"secs"` or `"days"` (any unit accepted by \code{\link[base]{difftime}}).
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
#' @param id Name of a column in `data` identifying individuals (default `"id"`). A ctmm
#' `telemetry` object normally has no such column: its identity is stored in
#' `attr(data, "info")$identity`, and that is used when the column is absent. If neither is
#' available, list elements are named by their list names (or position), and a single object
#' is individual `"1"`. Time differences are computed within individual and every individual
#' starts a new, independent Markov chain in the likelihood. Set to `NULL` to ignore any `id` column.
#' @param max_gap Optional. Gaps longer than `max_gap` (in `time_unit`) start a new, independent
#' segment (the chain is restarted from `delta`).
#' @details This function takes the HDOP information in the `telemetry` object to 
#' produce a `SpatRaster` likelihood surface over the `SpatRaster` defined by 
#' the `raster` argument for each location. This can then be passed to 
#' CTMC HMM fitting functions. The `times` data frame returned with the likelihood has
#' columns `obs`, `id`, `seg`, `timestamp` and `dt`; `dt` is `0` at the first observation of
#' every segment (individual or `max_gap` burst). Several outputs can be combined with
#' \code{\link{bind_proc_data}}.
#' @author Devin S. Johnson
#' @importFrom Matrix sparseMatrix
#' @importFrom terra vect buffer cells crds is.lonlat crs project res values
#' @importFrom ctmm uere
#' @importFrom methods as
#' @export
proc_telem <- function(data, cell_data, time_unit="hours", return_type="sparse", max_err=NULL,
                       trunc=1.0e-8, return_inadmissible=FALSE, nthreads=1,
                       neighborhood=c("disc","buffer"), id="id", max_gap=NULL){
  neighborhood <- match.arg(neighborhood)
  proc_check_threads(nthreads)
  # Check arguments
  dl <- proc_as_list(data, "telemetry", "a ctmm 'telemetry' object")
  if(!inherits(cell_data, "SpatRaster")) stop("'cell_data' must be a terra::SpatRaster object!")
  if(is.lonlat(cell_data)) stop("'cell_data' must be projected with units in meters.")
  xl <- dl$x
  
  # Individual / segment structure
  idl <- proc_ids(xl, id, dl$is_list, fallback=function(d){
    v <- tryCatch(attr(d, "info")$identity, error=function(e) NULL)
    if(is.null(v)) NULL else as.character(v)[1]
  })
  ts <- do.call(c, lapply(xl, function(d) d$timestamp))
  times <- proc_times(unlist(idl, use.names=FALSE), ts, time_unit, max_gap)
  attr(times, "time_unit") <- time_unit
  
  # Evaluate error covariance (per telemetry object: the UERE is per object)
  err_nms <- c("COV.x.x", "COV.x.y", "COV.y.y")
  cv <- lapply(xl, function(d){
    if(all(err_nms %in% colnames(d))){
      cbind(d$COV.x.x, d$COV.y.y, d$COV.x.y)
    } else {
      if(! "HDOP" %in% colnames(d)){
        warning("'HDOP' not detected in data it will be set to 1.")
        d$HDOP <- 1
      }
      u <- summary(uere(d))[,,"horizontal"]
      v <- (d$HDOP*u['est'])^2/2
      cbind(v, v, rep(0, length(v)))
    }
  })
  cv <- do.call(rbind, cv)
  cov_xx <- as.numeric(cv[,1]); cov_yy <- as.numeric(cv[,2]); cov_xy <- as.numeric(cv[,3])
  n <- nrow(times)
  
  # Project locations to the raster CRS
  lon <- unlist(lapply(xl, function(d) d$longitude), use.names=FALSE)
  lat <- unlist(lapply(xl, function(d) d$latitude), use.names=FALSE)
  telem_pts <- terra::vect(cbind(lon, lat), crs="epsg:4326") |>
    terra::project(crs(cell_data, proj=TRUE))
  xy <- terra::crds(telem_pts)
  
  if(is.null(max_err)) max_err <- 4*sqrt(pmax(cov_xx,cov_yy))
  max_err <- rep_len(as.numeric(max_err), n)
  
  proc_lik(xy[,1], xy[,2], cov_xx, cov_yy, cov_xy, max_err, cell_data, times, return_type, trunc,
           return_inadmissible, nthreads, neighborhood)
}
