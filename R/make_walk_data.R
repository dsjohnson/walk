#' @title Produce design data for use in fitting MMPP movement models
#' @param proc_data ---.
#' @param cell_data ---.
#' @param grad ---.
#' @param rast_mask Raster mask for inaccessible cells when \code{cell_data} is of type \code{SpatRaster} from the \code{terra} package. This is ignored
#' if \code{cell_data} is an \code{POLYGON} data frame from the \code{sf} package.
#' @param directions ---.
#' @param debug Debugging level: 1-3 mainly for package developers.
#' @param ... Ignored arguments.
#' @import dplyr 
#' @importFrom Matrix rowSums
#' @export
make_walk_data <- function(proc_data, cell_data, grad=NULL, rast_mask=NULL, 
                           directions="rook", debug=0,...){
  
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
  out <- list()
  out$L <- proc_data$L[,q_list$q_r$cell]
  if(min(Matrix::rowSums(out$L))==0) warning("Some observations have likelihood values of 0 in all ras cells after 'rast_mask' cells are removed!")
  out$times <- proc_data$times
  out$q_r <- q_list$q_r
  out$q_m <- q_list$q_m
  class(out) <- c(class(out), "walk_ddl")
  if(debug==3) browser()
  return(out)
}




