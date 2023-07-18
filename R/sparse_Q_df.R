#' @title Produce design data for use in fitting MMPP movement models
#' @param ras A `\link[terra]{SpatRaster}` stack of covariate rasters that will be used in CTMC 
#' movement modeling. Cells with `NA` values will be considered areas to the animal cannot travel, i.e., 
#' movement into those cells is not allowed. See the `\link[terra]{mask}` function to accomplish masking of
#' areas where the animal cannot go.
#' @param debug Debugging level: 1-3 mainly for package developers.
#' @param ... Ignored arguments.
#' @importFrom terra adjacent values xyFromCell terrain
#' @export
sparse_Q_df <- function(ras, grad=NULL, directions="rook", debug=0,...){
  
  cell <- cellx <- timestamp <- quad <- period <- fix <- from_cellx <- NULL
  neighborhood <- to_cellx <- boundary <- NULL
  
  # if(debug==1) browser()
  # 
  # quad_data <- select(sighting_data, timestamp, quad, period) %>% 
  #   filter(quad==1) %>% distinct() %>% arrange(timestamp) %>% select(-quad)
  
  
  ### Create Q matrix data frame
  if(debug==1) browser()
  
  nb <- terra::adjacent(ras, 1:ncell(ras), directions)
  nonzeros <- c(1:ncell(ras))[is.na(values(ras[[1]]))]
  
  q_data <- data.frame(from_cell=rep(1:ncell(ras), each=ncol(nb)), to_cell=as.vector(t(nb))) |>
    subset(!is.na(to_cell))
  
  from_df <- as.data.frame(terra::xyFromCell(ras, q_data$from_cell))
  from_df <- cbind(from_df, values(ras, data.frame=TRUE)[q_data$from_cell,])
  colnames(from_df) <- paste0("from_", colnames(from_df))
  
  to_df <- as.data.frame(terra::xyFromCell(ras, q_data$to_cell))
  to_df <- cbind(to_df, values(ras, data.frame=TRUE)[q_data$to_cell,])
  colnames(to_df) <- paste0("to_", colnames(to_df))
  
  q_data <- cbind(q_data, from_df, to_df)
  
  q_data$dist <- with(q_data,
                      sqrt((from_x-to_x)^2 + (from_y-to_y)^2)
  )
  q_data$w_x <- (q_data$to_x - q_data$from_x)/q_data$dist
  q_data$w_y <- (q_data$to_y - q_data$from_y)/q_data$dist
  n <- ncol(q_data)
  
  if(!is.null(grad)){
    for(i in 1:length(grad)){
      cov_tmp <- ras[[grad[i]]]
      grad_tmp <- get_grad(cov_tmp) |> terra::values()
      grad_tmp <- grad_tmp[q_data$from_cell,]
      grad_tmp <- grad_tmp[,"dx"]*q_data$w_x + grad_tmp[,"dy"]*q_data$w_y
      q_data <- cbind(q_data, grad_tmp)
      
    }
  }
  
  # if(dynamic_movement){
  #   # add code here to expand q_data by quad_data if time indexed movement is 
  #   # desired.
  # }
  
  
  
  if(debug==2) browser()
  
  return(out)
  
}


get_grad <- function(r){
  sa <- terra::terrain(r, v=c('slope', 'aspect'), unit="radians")
  sa[["dx"]] <- sa[['slope']]*sin(sa[['aspect']])
  sa[["dy"]] <- sa[['slope']]*cos(sa[['aspect']])
  sa <- sa[[-c(1:2)]]
  return(sa)
}
