#' @title Denote cells where the animal cannot travel
#' @param x a polygon object that denotes places where the animal cannot travel
#' @param y A `SpatRaster` where cells will be set to `NA` that are inside the polygon object.
#' @param overlap The amount of overlap betwen the raster cell and polygon such that the cell is considered 
#' "inside" the polygon. Defaults to `overlap = 0.95`. If `overlap = 1` then a cell has to be completely inside the 
#' polygon to be maked.
#' @export
#' @importFrom sf st_geometry_type
#' @importFrom terra vect geomtype rasterize
set_barrier_mask <- function(x, y, overlap=0.95){
  if(!inherits(x, "SpatRaster")) stop("'x' must be of class 'SpatRaster'")
  if(inherits(y, "sf")){
    if(!all(st_geometry_type(y) %in% c("POLYGON","MULTIPOLYGON"))){stop("y must be polygons!")}
    y <- terra::vect(y)
  } else if(inherits(y, "SpatVector")){
    if(terra::geomtype(y)!="polygons") stop("y must be polygons!")
  } else{
    stop("y must be of class 'sf' or 'SpatVector")
  }
  b <- terra::rasterize(y, x, touches=TRUE, cover=TRUE)
  x[b > overlap] <- NA
  return(x)
}