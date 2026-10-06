#' @title Create neighborhood indices for subsetting cells between movement locations to improve 
#' computational speed at the cost of some likelihood approximation. 
#' @param walk_data A `walk_ddl` object created by `walk::make_walk_data`
#' @param vmax A single numerical value describing a rough estimate of the maximum velocity of the animal. This will be used for 
#' for calculating the buffer around movement subsets for each time step, `buffer[i] = vmax * dt[i`, for `i = 1...,N`.
#' @param subset_buffer A value > 0 such that the bounding box of cells created 
#' by union of the location uncertainty layers (`walk_data$L`) are buffered by `subset_buffer*walk_data$times$dt[i]` to reduce the likelihood 
#' approximation. Once `subset_buffer` is large enough such that all cells are involved 
#' at each time step, there is no longer any approximation. 
#' @import dplyr
#' @export 
make_subset_neigh <- function(walk_data, vmax=NULL, subset_buffer=0){
  if(!inherits(walk_data, "walk_ddl")) stop("'walk_data' must be a 'walk_ddl' object!")
  # Pre-computation in R before calling C++ likelihood
  xy_df <- walk_data$q_r |> select(cell, cellx, x, y)
  dt <- walk_data$times$dt
  N <- length(dt)
  L <- walk_data$L
  active_indices_list <- vector("list", N)
  if(is.numeric(vmax)){
    if(length(vmax)>1) stop("The 'vmax' argument must be a single value.")
    subset_buffer <- vmax*sqrt(dt)
  } else{
    if(length(subset_buffer)==1) subset_buffer <- rep(subset_buffer, N)
  }
  if(length(subset_buffer)!=N) stop("Something is wrong. Resulting 'subset_buffer' is not the correct length.")
  
  
  for(i in 2:N) {
    # 1. Identify non-zero likelihood cells for L(i-1) and L(i)
    cells_prev <- which(L[i-1, ] > 0)
    cells_curr <- which(L[i, ] > 0)
    
    # 2. Get min/max x, y bounds from both observations
    coords_prev <- filter(xy_df, cellx%in%cells_prev)
    coords_curr <- filter(xy_df, cellx%in%cells_curr)
    
    x_min <- min(coords_prev$x, coords_curr$x) - subset_buffer[i-1]*dt[i]
    x_max <- max(coords_prev$x, coords_curr$x) + subset_buffer[i-1]*dt[i]
    y_min <- min(coords_prev$y, coords_curr$y) - subset_buffer[i-1]*dt[i]
    y_max <- max(coords_prev$y, coords_curr$y) + subset_buffer[i-1]*dt[i]
    
    # 3. Generate contiguous 1D grid cell indices for this bounding box
    active_indices_list[[i]] <- filter(xy_df, x>=x_min, x<=x_max, y>=y_min, y<=y_max)$cellx-1
  }
  
  # Step 1 active set is based solely on L[1, ]
  active_indices_list[[1]] <- 0
  
  # Flatten into contiguous vectors for C++
  active_indices <- unlist(active_indices_list)
  active_lengths <- sapply(active_indices_list, length)
  active_offsets <- c(0, cumsum(active_lengths)) # C++ offsets
  return(list(active_offsets=active_offsets, active_indices=active_indices))
}