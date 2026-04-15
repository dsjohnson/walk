#' @title Z transform a raster
#' @param r A spatRaster object from the terra package
#' @importFrom terra global
#' @export
raster_z <- function(r){
  return(
  (r - global(r, "mean", na.rm=TRUE)[1,1]) / global(r, "sd", na.rm=TRUE)[1,1]
  )
}

#' @title Bhattacharyya distance for 2 raster probability densities
#' @param r SpatRaster object
#' @param s SpatRaster object
#' @details
#' Computes Bhattacharyya distance for 2 raster UDs. 
#' @importFrom terra global
#' @export
bc_raster <- function(r,s){
  r <- r/(global(r, "sum", na.rm=T)[1,1])
  s <- s/(global(s, "sum", na.rm=T)[1,1])
  xx <- sqrt(r*s)
  bc <- global(xx,"sum",na.rm=TRUE)[1,1]
  return(bc)
}
