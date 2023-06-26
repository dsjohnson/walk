
#' @title Convert {ctmm} telemtry data object to a {terra} `SpatRaster` stack
#' @param data A `telemetry` object from the {ctmm} package
#' @param cov A `\link[terra]{SpatRaster}` stack of covariates that will be used in CTMC 
#' movement modeling. Cells with `NA` values will be considered areas to the animal cannot travel, i.e., 
#' likelihood surfaces will be `0` for those cells.
#' @param return_type Type of object returned. One if `"data.frame"`, `"sparse"` (sparse matrix), 
#' `"matrix_df"` (matrix form of `"data.frame"`), or `"dense"` (dense matrix).
#' @param max_err The maximum error in meters. If unspecified it will be set to
#' 4 times the maximum error standard deviation as determined by the UERE and HDOP
#' of the telemetry data.
#' @param trunc The smallest probability value that is considered to be > 0. Deafults to 1.0e-8.
#' @details This function takes the HDOP information in the `telemetry` object to 
#' produce a `SpatRaster` likelihood surface over the `SpatRaster` defined by 
#' the `raster` argument for each location. This can then be passed to 
#' CTMC HMM fitting functions. 
#' @author Devin S. Johnson
#' @importFrom mvtnorm pmvnorm
#' @importFrom Matrix sparseMatrix
#' @importFrom terra vect buffer cells crds is.lonlat xyFromCell crs
#' @importFrom ctmm uere
#' @importFrom methods as
#' @export
telem_to_sparse <- function(data, cov, return_type="sparse", max_err=NULL, trunc=1.0e-8){
  # Check arguments
  if(!inherits(data, "telemetry")) stop("'data' must be a ctmm::telemery object!")
  if(!inherits(cov, "SpatRaster")) stop("'cov' must be a terra::SpatRaster object!")
  if(is.lonlat(cov)) stop("'cov' must be projected with units in meters.")
  
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
  cov_prj <- crs(cov, proj=TRUE)
  telem_pts <- terra::vect(cbind(data$longitude, data$latitude), crs="epsg:4326") |>
    terra::project(cov_prj)
  if(is.null(max_err)) max_err <- 4*sqrt(pmax(cov_xx,cov_yy))
  t_buf <- terra::buffer(telem_pts, max_err)
  t_err <- terra::cells(cov, t_buf, touches=TRUE) 
  t_err <- split(t_err[,'cell'], t_err[,'ID'])
  xy <- terra::crds(telem_pts) 
  
  lik_list <- vector("list", nrow(xy))
  
  t_err <- t_err |> bind_cols(get_corner(t_err$cell, cov, "ur"))
  t_err <- t_err |> bind_cols(get_corner(t_err$cell, cov, "ll"))
  
  out <- NULL
  for(i in 1:nrow(xy)){
    # i <- 1
    sigma <- matrix(c(cov_xx[i], cov_xy[i], cov_xy[i], cov_yy[i]), 2, 2)
    mean <- as.vector(xy[i,])
    lower <- get_corner(t_err[[i]], cov, "ll")
    upper <- get_corner(t_err[[i]], cov, "ur")
    dfi <- cbind(obs=i, cell=t_err[[i]], lik=NA)
    dfi[,'lik'] <- sapply(1:length(t_err[[i]]), 
                          \(j) pmvnorm(lower=lower[j,], upper=upper[j,], mean=mean, sigma=sigma)[1]
    )
    dfi[,3] <- dfi[,3]/sum(dfi[,3])
    dfi[,3] <- ifelse(dfi[,3]<trunc, 0, dfi[,3])
    dfi <- dfi[dfi[,3]>0,,drop=FALSE]
    dfi[,3] <- dfi[,3]/sum(dfi[,3])
    out <- rbind(out, dfi)
    
    if(return_type=="data.frame"){
      return(as.data.frame(out))
    } else if(return_type=="sparse"){
      M <- Matrix::sparseMatrix(i = out[,'obs'], j = out[,'cell'], x = out[,'lik'], dims = c(nrow(xy), prod(dim(cov)[1:2])))
      M <- as(M, "dgCMatrix")
      return(M)
    } else if(return_type=="matrix_df"){
      return(out)
    } else if(return_type=="dense"){
      M <- Matrix::sparseMatrix(i = out[,'obs'], j = out[,'cell'], x = out[,'lik'], dims = c(nrow(xy), prod(dim(cov)[1:2])))
      M <- as.matrix(M)
      return(M)
    } else{
      warning("Unknown 'return_type' returning 'matrix_df'")
      return(out)
    }
    
  }
  
  # 
  # Possible extention:
  # https://stackoverflow.com/questions/51290014/rcpp-implementation-of-mvtnormpmvnorm-slower-than-original-r-function
  #  source code: https://rdrr.io/cran/mvtnorm/src/R/mvt.R
  
  
}


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
