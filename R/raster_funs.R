#' @title Z transform a raster
#' @param r A spatRaster object from the terra package
#' @importFrom terra global
#' @export
raster_z <- function(r){
  return(
  (r - global(r, "mean", na.rm=TRUE)[1,1]) / global(r, "sd", na.rm=TRUE)[1,1]
  )
}
