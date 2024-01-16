
#' @title Convert {ctmm} telemetry data object to a {terra} `SpatRaster` stack
#' @param data A `telemetry` object from the {ctmm} package
#' @param ras A `\link[terra]{SpatRaster}` stack of covariates that will be used in CTMC 
#' movement modeling. Cells with `NA` values will be considered areas to the animal cannot travel, i.e., 
#' likelihood surfaces will be `0` for those cells.
#' @param return_type Type of object returned. One if `"data.frame"`, `"sparse"` (sparse matrix), 
#' `"matrix_df"` (matrix form of `"data.frame"`), or `"dense"` (dense matrix).
#' @param max_err The maximum error in meters. If unspecified it will be set to
#' 4 times the maximum error standard deviation as determined by the UERE and HDOP
#' of the telemetry data.
#' @param trunc The smallest probability value that is considered to be > 0. Defaults to 1.0e-8.
#' @details This function takes the HDOP information in the `telemetry` object to 
#' produce a `SpatRaster` likelihood surface over the `SpatRaster` defined by 
#' the `raster` argument for each location. This can then be passed to 
#' CTMC HMM fitting functions. 
#' @author Devin S. Johnson
#' @importFrom mvtnorm pmvnorm
#' @importFrom Matrix sparseMatrix
#' @importFrom terra vect buffer cells crds is.lonlat xyFromCell crs project
#' @importFrom ctmm uere
#' @importFrom methods as
#' @export
telem_to_ras <- function(data, ras, return_type="sparse", max_err=NULL, trunc=1.0e-8){
  # Check arguments
  if(!inherits(data, "telemetry")) stop("'data' must be a ctmm::telemery object!")
  if(!inherits(ras, "SpatRaster")) stop("'ras' must be a terra::SpatRaster object!")
  if(is.lonlat(ras)) stop("'ras' must be projected with units in meters.")
  
  # Evaluate error covariance
  err_nms <- c("COV.x.x", "COV.x.y", "COV.y.y")
  if(all(err_nms %in% colnames(data))){
    cov_xx <- data$COV.x.x
    cov_yy <- data$COV.y.y
    cov_xy <- data$COV.x.y
  } else {
    if(! "HDOP" %in% colnames(data)){
      warning("'HDOP' not detected in data it will be set to 1.")
      data$HDOP <- 1
    }
    u <- summary(uere(data))[,,"horizontal"]
    cov_xx <- (data$HDOP*u['est'])^2/2
    cov_yy <- (data$HDOP*u['est'])^2/2
    cov_xy <- rep(0, length(cov_xx))
  }
  
  # Extract error neighborhood and get cells 
  ras_prj <- crs(ras, proj=TRUE)
  telem_pts <- terra::vect(cbind(data$longitude, data$latitude), crs="epsg:4326") |>
    terra::project(ras_prj)
  if(is.null(max_err)) max_err <- 4*sqrt(pmax(cov_xx,cov_yy))
  t_buf <- terra::buffer(telem_pts, max_err)
  t_err <- terra::cells(ras, t_buf, touches=TRUE) 
  t_err <- split(t_err[,'cell'], t_err[,'ID'])
  xy <- terra::crds(telem_pts) 
  
  lik_list <- vector("list", nrow(xy))
  out <- NULL
  zero_ind <- FALSE
  ras_val <- values(ras)
  
  for(i in 1:nrow(xy)){
    # i <- 1
    if(any(!is.finite(t_err[[i]]))) stop("There are error buffered locations completely outside of raster area.")
    sigma <- matrix(c(cov_xx[i], cov_xy[i], cov_xy[i], cov_yy[i]), 2, 2)
    mean <- as.vector(xy[i,])
    lower <- get_corner(t_err[[i]], ras, "ll")
    upper <- get_corner(t_err[[i]], ras, "ur")
    dfi <- cbind(obs=i, cell=t_err[[i]], lik=NA)
    dfi[,'lik'] <- sapply(1:length(t_err[[i]]), 
                          \(j) pmvnorm(lower=lower[j,], upper=upper[j,], mean=mean, sigma=sigma, keepAttr=FALSE)
    )
    m <- ifelse(is.na(ras_val[t_err[[i]]]), 0, 1)
    dfi[,3] <- dfi[,3]*m
    dfi[,3] <- dfi[,3]/sum(dfi[,3])
    dfi[,3] <- ifelse(dfi[,3]<trunc, 0, dfi[,3])
    dfi <- dfi[dfi[,3]>0,,drop=FALSE]
    dfi[,3] <- dfi[,3]/sum(dfi[,3])
    if(sum(dfi[,3]) <=sqrt(.Machine$double.eps)) zero_ind <- TRUE
    out <- rbind(out, dfi)
  }
  
  if(zero_ind) warning("Some observations have likelihood values of 0 in all ras cells!")
  
  if(return_type=="data.frame"){
    return(as.data.frame(out))
  } else if(return_type=="sparse"){
    M <- Matrix::sparseMatrix(i = out[,'obs'], j = out[,'cell'], x = out[,'lik'], dims = c(nrow(xy), prod(dim(ras)[1:2])))
    M <- as(M, "dgCMatrix")
    return(M)
  } else if(return_type=="matrix_df"){
    return(out)
  } else if(return_type=="dense"){
    M <- Matrix::sparseMatrix(i = out[,'obs'], j = out[,'cell'], x = out[,'lik'], dims = c(nrow(xy), prod(dim(ras)[1:2])))
    M <- as.matrix(M)
    return(M)
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


#' @importFrom terra res
get_corner <- function(cells, r, which){
  sz <- res(r)
  xy <- xyFromCell(r, cells)
  if(which=="ur"){
    xy[,1] <- xy[,1] + sz[1]/2
    xy[,2] <- xy[,2] + sz[2]/2
    colnames(xy) <- c("x_ur","y_ur")
  } else if(which=="ll"){
    xy[,1] <- xy[,1] - sz[1]/2
    xy[,2] <- xy[,2] - sz[2]/2
    colnames(xy) <- c("x_ll","y_ll")
  }
  return(xy)
}
