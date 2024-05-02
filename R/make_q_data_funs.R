
###-----------------------------------------------------------------------------
### Q data for {terra} rasters
###-----------------------------------------------------------------------------
#' @importFrom terra compareGeom xFromCell yFromCell adjacent ncell values
#' @import dplyr
make_q_data_rast <- function(cell_data, grad=NULL, rast_mask=NULL, directions="rook", ...){
  
  mask <- cell <- from_cell <- NULL
  
  ### data checks
  if(!is.null(rast_mask)){
    if(!terra::compareGeom(cell_data, rast_mask)) stop("'cell_data' and 'rast_mask' must have the same geometry!")
  }
  
  N <- terra::ncell(cell_data)
  
  ### Make Q_r data
  q_r_data <- dplyr::tibble(
    cell = 1:N,
    x = terra::xFromCell(cell_data, 1:N),
    y = terra::yFromCell(cell_data, 1:N),
  ) 
  q_r_data <- q_r_data %>% bind_cols(terra::values(cell_data, dataframe=TRUE))
  
  if(!is.null(rast_mask)){
    q_r_data$mask <- terra::values(rast_mask, dataframe=TRUE)[,1]
    q_r_data <- q_r_data %>% dplyr::filter(!is.na(mask)) %>% dplyr::select(-mask)
  } 
  
  q_r_data <- q_r_data %>% mutate(
    cellx =  as.numeric(factor(cell))
  )
  
  ### Make Q_m data
  adj_mat <- terra::adjacent(cell_data, 1:N, directions=directions)
  q_m_data <- tibble(
    from_cell = rep(1:N, each=ncol(adj_mat)), 
    cell=as.vector(t(adj_mat))
  ) %>% filter(!is.na(cell)) %>% filter(from_cell %in% q_r_data$cell) %>% 
    filter(cell %in% q_r_data$cell)
  
  covs <- q_r_data 
  q_m_data <- q_m_data %>% left_join(covs, by="cell")
  colnames(covs) <- paste0("from_", colnames(covs))
  q_m_data <- q_m_data %>% left_join(covs, by="from_cell")
  rm(covs)
  # browser()
  q_m_data$dist <- with(q_m_data,sqrt((from_x-x)^2 + (from_y-y)^2))
  q_m_data$w_x <- (q_m_data$x - q_m_data$from_x)/q_m_data$dist
  q_m_data$w_y <- (q_m_data$y - q_m_data$from_y)/q_m_data$dist
  
  if(!is.null(grad)){
    grad_df <- NULL
    for(i in 1:length(grad)){
      cov_tmp <- cell_data[[grad[i]]]
      grad_tmp <- get_grad(cov_tmp) |> terra::values()
      grad_tmp <- grad_tmp[q_m_data$from_cell,]
      grad_tmp <- grad_tmp[,"dx"]*q_m_data$w_x + grad_tmp[,"dy"]*q_m_data$w_y
      grad_tmp <- ifelse(!is.finite(grad_tmp), 0, grad_tmp)
      grad_df <- cbind(grad_df, grad_tmp)
    }
    colnames(grad_df) <- paste0(grad, "_grad")
    q_m_data <- cbind(q_m_data, grad_df)
  }
  
  return(list(q_r = q_r_data, q_m = q_m_data))
  
}


###-----------------------------------------------------------------------------
### Q data for {sf} polygons
###-----------------------------------------------------------------------------

#' @importFrom sf st_drop_geometry st_geometry_type
#' @importFrom spdep poly2nb
#' @import dplyr
#' @importFrom tidyr unnest
make_q_data_sf <- function(cell_data, cell_name,...){
  
  cellx <- from_cellx <- neighborhood <- NULL
  
  ### checks
  if(!all(st_geometry_type(cell_data) %in% c("POLYGON","MULTIPOLYGON"))){stop("'cell_data' all geometry types must be POLYGON or MULTIPOLYGON!")}
  if(is.null(cell_name)){
    cell_name <- "cell"
    cell_data$cell <- 1:nrow(cell_data)
  }
  cell_data$cellx <- 1:nrow(cell_data) #as.integer(factor(cell_data[[cell_name]]))
  cell_data <- cell_data %>% rename(cell=.data[[cell_name]])
  ### Make Q_r data
  q_r_data <- cell_data %>% st_drop_geometry()
  
  ### Make Q_m data
  nb <- spdep::poly2nb(cell_data, queen=FALSE)
  q_m_data <- cell_data %>% select(cellx) %>% st_drop_geometry() %>% rename(from_cellx = cellx) %>%
    rowwise() %>%
    mutate(
      neighborhood = list(
        tibble(cellx = nb[[from_cellx]], num_neigh=length(nb[[from_cellx]]))
      )
    ) %>% ungroup() %>% arrange(from_cellx) %>% unnest(cols=neighborhood) 
  q_m_data <- q_m_data %>% arrange(from_cellx, cellx)
  covs <- q_r_data 
  q_m_data <- q_m_data %>% left_join(covs, by="cellx")
  colnames(covs) <- paste0("from_", colnames(covs))
  q_m_data <- q_m_data %>% left_join(covs, by="from_cellx")
  rm(covs)
  
  return(list(q_r = q_r_data, q_m = q_m_data))
}


get_grad <- function(r){
  h <- res(r)
  sa <- terra::terrain(r, v=c('slope', 'aspect'), unit="radians")
  sa[["dx"]] <- h[1]*sa[['slope']]*sin(sa[['aspect']])
  sa[["dy"]] <- h[2]*sa[['slope']]*cos(sa[['aspect']])
  sa <- sa[[-c(1:2)]]
  
  return(sa)
}
