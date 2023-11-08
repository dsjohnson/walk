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
  cell_conv <- data.frame(site_idx=1:length(nonzeros), cell=nonzeros)
  
  q_r_data <- cell_conv
  q_r_data <- cbind(q_r_data, as.data.frame(terra::xyFromCell(ras, q_r_data$cell)))
  q_r_data <- cbind(q_r_data, values(ras, data.frame=TRUE)[q_r_data$cell,])
  
  q_m_data <- data.frame(
    from_cell = rep(nonzeros, each=ncol(nb)),
    from_site_idx = rep(cell_conv$site_idx, each=ncol(nb)),
    cell=as.vector(t(nb))
    ) |> subset(!is.na(cell))
  q_m_data <- merge(q_m_data, cell_conv, sort=FALSE)
  q_m_data <- q_m_data[,c("from_site_idx","site_idx","from_cell","cell")]
  q_m_data <- q_m_data[order(q_m_data$from_site_idx, q_m_data$site_idx),]
  q_m_data <- cbind(q_m_data, as.data.frame(terra::xyFromCell(ras, q_m_data$cell)))
  q_m_data <- cbind(q_m_data, values(ras, data.frame=TRUE)[q_m_data$cell,])
  from_df <- as.data.frame(terra::xyFromCell(ras, q_m_data$from_cell))
  from_df <- cbind(from_df, values(ras, data.frame=TRUE)[q_m_data$from_cell,])
  colnames(from_df) <- paste0("from_", colnames(from_df))
  q_m_data <- cbind(q_m_data, from_df)
  q_m_data$dist <- with(q_m_data,
                      sqrt((from_x-x)^2 + (from_y-y)^2)
  )
  q_m_data$w_x <- (q_m_data$x - q_m_data$from_x)/q_m_data$dist
  q_m_data$w_y <- (q_m_data$y - q_m_data$from_y)/q_m_data$dist
  # n <- ncol(q_data)
  
  if(!is.null(grad)){
    grad_df <- NULL
    for(i in 1:length(grad)){
      cov_tmp <- ras[[grad[i]]]
      grad_tmp <- get_grad(cov_tmp) |> terra::values()
      grad_tmp <- grad_tmp[q_m_data$from_cell,]
      grad_tmp <- grad_tmp[,"dx"]*q_m_data$w_x + grad_tmp[,"dy"]*q_m_data$w_y
      grad_tmp <- ifelse(!is.finite(grad_tmp), 0, grad_tmp)
      grad_df <- cbind(grad_df, grad_tmp)
    }
    colnames(grad_df) <- paste0(grad, "_grad")
    q_m_data <- cbind(q_m_data, grad_df)
  }
  
  # if(dynamic_movement){
  #   # add code here to expand q_data by quad_data if time indexed movement is 
  #   # desired.
  # }
  
  if(debug==2) browser()
  
  q_data <- list(q_r = q_r_data, q_m = q_m_data)
  class(q_data) <- c(class(q_data), "Qdf")
  
  return(q_data)
  
}


get_grad <- function(r){
  sa <- terra::terrain(r, v=c('slope', 'aspect'), unit="radians")
  sa[["dx"]] <- sa[['slope']]*sin(sa[['aspect']])
  sa[["dy"]] <- sa[['slope']]*cos(sa[['aspect']])
  sa <- sa[[-c(1:2)]]
  return(sa)
}
