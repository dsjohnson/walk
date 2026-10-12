
#' @title Convert `simLangevin`  object to a `SpatRaster` stack from the `terra` package.
#' @param data A `simLangevin` object from the `langevinSSM` package
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
#' @details This function takes the HDOP information in the `telemetry` object to 
#' produce a `SpatRaster` likelihood surface over the `SpatRaster` defined by 
#' the `raster` argument for each location. This can then be passed to 
#' CTMC HMM fitting functions. 
#' @author Devin S. Johnson
#' @importFrom mvtnorm pmvnorm
#' @importFrom Matrix sparseMatrix
#' @importFrom terra vect buffer cells crds is.lonlat xyFromCell crs project
#' @importFrom methods as
#' @export
proc_simLangevin<- function(data, cell_data, return_type="sparse", max_err=NULL, trunc=1.0e-8, return_inadmissible=FALSE){
  # Check arguments
  if(!inherits(data, "simLangevin")) stop("'data' must be a langevinSSM::simLangevin object!")
  if(!inherits(cell_data, "SpatRaster")) stop("'cell_data' must be a terra::SpatRaster object!")
  
  # Evaluate error covariance
  cov_df <- argos_diag_to_cov(data$smaj, data$smin, 180*data$eor/pi)
  
  # Extract error neighborhood and get cells 
  telem_pts <- terra::vect(cbind(data$x, data$y)) 
  if(is.null(max_err)) max_err <- 4*sqrt(pmax(cov_df$cov.x.x,cov_df$cov.y.y))
  t_buf <- terra::buffer(telem_pts, max_err)
  t_err <- terra::cells(cell_data, t_buf, touches=TRUE) 
  t_err <- split(t_err[,'cell'], t_err[,'ID'])
  xy <- terra::crds(telem_pts) 
  
  lik_list <- vector("list", nrow(xy))
  out <- NULL
  out_inad <- NULL
  zero_ind <- FALSE
  cell_data_val <- values(cell_data)
  
  for(i in 1:nrow(xy)){
    # i <- 1
    if(any(!is.finite(t_err[[i]]))) stop("There are error buffered locations completely outside of raster area.")
    sigma <- matrix(c(cov_df$cov.x.x[i], cov_df$cov.x.y[i], cov_df$cov.x.y[i], cov_df$cov.y.y[i]), 2, 2)
    mean <- as.vector(xy[i,])
    lower <- get_corner(t_err[[i]], cell_data, "ll")
    upper <- get_corner(t_err[[i]], cell_data, "ur")
    dfi <- cbind(obs=i, cell=t_err[[i]], lik=NA)
    dfi[,'lik'] <- sapply(1:length(t_err[[i]]), 
                          \(j) pmvnorm(lower=lower[j,], upper=upper[j,], mean=mean, sigma=sigma, keepAttr=FALSE)
    )
    m <- ifelse(is.na(cell_data_val[t_err[[i]]]), 0, 1)
    dfi[,3] <- dfi[,3]*m
    if(sum(dfi[,3]) <=sqrt(.Machine$double.eps)){
      zero_ind <- TRUE
      out_inad <- rbind(out_inad, dfi)
      next
    }
    dfi[,3] <- dfi[,3]/sum(dfi[,3])
    dfi[,3] <- ifelse(dfi[,3]<trunc, 0, dfi[,3])
    dfi <- dfi[dfi[,3]>0,,drop=FALSE]
    out <- rbind(out, dfi)
  }
  
  if(return_inadmissible) return(out_inad)
  if(zero_ind) warning("Some observations have likelihood values of 0 in all 'cell_data' cells!")
  
  times <- data.frame(obs=1:nrow(data), timestamp = data$date)
  times <- cbind(times, dt=c(0,diff(times$timestamp)))
  
  if(return_type=="data.frame"){
    out <- as.data.frame(out)
    out <- merge(out, times, by="obs")
    return(out)
  } else if(return_type=="sparse"){
    M <- Matrix::sparseMatrix(i = out[,'obs'], j = out[,'cell'], x = out[,'lik'], dims = c(nrow(xy), prod(dim(cell_data)[1:2])))
    M <- as(M, "dgCMatrix")
    return(list(L=M, times=times))
  } else if(return_type=="matrix_df"){
    return(out)
  } else if(return_type=="dense"){
    M <- Matrix::sparseMatrix(i = out[,'obs'], j = out[,'cell'], x = out[,'lik'], dims = c(nrow(xy), prod(dim(cell_data)[1:2])))
    M <- as.matrix(M)
    return(list(L=M, times=times))
  } else{
    warning("Unknown 'return_type' returning 'matrix_df'")
    return(out)
  }
  
  # 
  # Possible extension for doing this in C++ directly:
  # https://stackoverflow.com/questions/51290014/rcpp-implementation-of-mvtnormpmvnorm-slower-than-original-r-function
  #  source code: https://rdrr.io/cran/mvtnorm/src/R/mvt.R
  #
  
}
