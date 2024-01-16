#' @title Produce design data for use in fitting MMPP movement models
#' @param cell_data ---.
#' @param grad ---.
#' @param rast_mask Raster mask for inaccessible cells when \code{cell_data} is of type \code{SpatRaster} from the \code{terra} package. This is ignored
#' if \code{cell_data} is an \code{POLYGON} data frame from the \code{sf} package.
#' @param directions ---.
#' @param debug Debugging level: 1-3 mainly for package developers.
#' @param ... Ignored arguments.
#' @import dplyr 
#' @importFrom units set_units
#' @export
make_design_data <- function(cell_data, grad=NULL, rast_mask=NULL, directions="rook", debug=0,...){
  
  cells <- obs <- timestamp <- quad <- NULL
  
  if(debug==1) browser()
  cell_name <- NULL
  quad_data <- NULL
  
  if(debug==2) browser()
  
  if(inherits(cell_data, "SpatRaster")){
    q_list <- make_q_data_rast(cell_data, grad=grad, rast_mask=rast_mask, directions=directions)
  } else if(inherits(cell_data, "sf")){
    q_list <- make_q_data_sf(cell_data, cell_name)
  }
  
  class(q_list) <- c(class(q_list), "Qdf")
  if(debug==3) browser()
  return(q_list)
  
}




