#' @title Produce design data for use in fitting MMPP movement models
#' @param ras A `\link[terra]{SpatRaster}` stack of covariate rasters that will be used in CTMC 
#' movement modeling. Cells with `NA` values will be considered areas to the animal cannot travel, i.e., 
#' movement into those cells is not allowed. See the `\link[terra]{mask}` function to accomplish masking of
#' areas where the animal cannot go. If `ras` has multiple layers, only the first is used to determine masking.
#' @param debug Debugging level: 1-3 mainly for package developers.
#' @param ... Ignored arguments.
#' @importFrom terra adjacent values xyFromCell terrain
#' @export
make_Q_data <- function(ras, grad=NULL, directions="rook", debug=0,...){
  
  cell <- cellx <- timestamp <- quad <- period <- fix <- r_cellx <- NULL
  neighborhood <- m_cellx <- boundary <- NULL
  
  # if(debug==1) browser()
  # 
  # quad_data <- select(sighting_data, timestamp, quad, period) %>% 
  #   filter(quad==1) %>% distinct() %>% arrange(timestamp) %>% select(-quad)
  
  
  ### Create Q matrix data frame
  if(debug==1) browser()
  
  nb <- terra::adjacent(ras, 1:ncell(ras), directions)
  nonzeros <- c(1:ncell(ras))[!is.na(values(ras[[1]]))] |> sort()
  nb <- nb[nonzeros,]
  
  q_data <- data.frame(r_cell=rep(nonzeros, each=ncol(nb)), m_cell=as.vector(t(nb))) |>
    subset(!is.na(m_cell))
  
  r_df <- as.data.frame(terra::xyFromCell(ras, q_data$r_cell))
  r_df <- cbind(r_df, values(ras, data.frame=TRUE)[q_data$r_cell,])
  colnames(r_df) <- paste0("r_", colnames(r_df))
  
  m_df <- as.data.frame(terra::xyFromCell(ras, q_data$m_cell))
  m_df <- cbind(m_df, values(ras, data.frame=TRUE)[q_data$m_cell,])
  colnames(m_df) <- paste0("m_", colnames(m_df))
  
  q_data <- cbind(q_data, r_df, m_df)
  
  q_data$dist <- with(q_data,
                      sqrt((r_x-m_x)^2 + (r_y-m_y)^2)
  )
  q_data$w_x <- (q_data$m_x - q_data$r_x)/q_data$dist
  q_data$w_y <- (q_data$m_y - q_data$r_y)/q_data$dist
  n <- ncol(q_data)
  
  if(!is.null(grad)){
    grad_df <- NULL
    for(i in 1:length(grad)){
      cov_tmp <- ras[[grad[i]]]
      grad_tmp <- get_grad(cov_tmp) |> terra::values()
      grad_tmp <- grad_tmp[q_data$r_cell,]
      grad_tmp <- grad_tmp[,"dx"]*q_data$w_x + grad_tmp[,"dy"]*q_data$w_y
      grad_tmp <- ifelse(!is.finite(grad_tmp), 0, grad_tmp)
      grad_df <- cbind(grad_df, grad_tmp)
    }
    colnames(grad_df) <- paste0("m_", grad, "_grad")
    q_data <- cbind(q_data, grad_df)
  }
  
  # if(dynamic_movement){
  #   # add code here to expand q_data by quad_data if time indexed movement is 
  #   # desired.
  # }
  
  if(debug==2) browser()
  
  return(q_data)
  
}


get_grad <- function(r){
  sa <- terra::terrain(r, v=c('slope', 'aspect'), unit="radians")
  sa[["dx"]] <- sa[['slope']]*sin(sa[['aspect']])
  sa[["dy"]] <- sa[['slope']]*cos(sa[['aspect']])
  sa <- sa[[-c(1:2)]]
  return(sa)
}
