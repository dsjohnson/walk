#' @title Produce design data for use in fitting MMPP movement models
#' @param proc_data A sparse matrix with rows corresponding to time steps and columns corresponding to cells. 
#' The entries are probabilities that the animal is located in the corresponding cell and time.
#' @param cell_data A `SpatRaster` object from the `terra` package.
#' @param grad A character vector of names of `cell_data` layers for which gradient covariates will be constructed.
#' @param rast_mask Raster mask for inaccessible cells when `cell_data` is of type `SpatRaster` from the `terra` package. This is ignored
#' if \code{cell_data} is an `POLYGON` data frame from the `sf` package.
#' @param directions Neighborhood structure, one of `"rook"` or `"queen"`.
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
    attr(q_list$q_m, "directions") <- directions
    attr(q_list$q_m, "cov_class") <- "SpatRaster"
  } else if(inherits(cell_data, "sf")){
    stop("In the current version cell_data must be a 'SpatRaster' from the terra package. In future versions we hope to incorporate sf polygons data sets.")
    # q_list <- make_q_data_sf(cell_data, cell_name)
    # attr(q_list$q_m, "directions") <- NULL
    # attr(q_list$q_m, "cov_class") <- "sf"
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




