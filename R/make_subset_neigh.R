#' @title Create neighborhood indices for subsetting cells between movement locations to improve 
#' computational speed at the cost of some likelihood approximation. 
#' @param walk_data A `walk_ddl` object created by `walk::make_walk_data`
#' @param vmax A single numerical value describing a rough estimate of the maximum velocity of the animal. This will be used for 
#' for calculating the buffer around movement subsets for each time step, `buffer[i] = vmax * sqrt(dt[i])`, for `i = 1...,N`.
#' The `vmax` argument is assumed to be in units of prjected coordinate units per time units in the original telemetry data, e.g., m/h usually. 
#' @param min_buffer A value >= 0 such that the buffer is guaranteed to be at least `min_buffer` to prevent numerical issues
#' for small time steps and high precision locations. If not specified, calculates an approximate 2 cell minimum buffer.  
#' @import dplyr
#' @export 
make_subset_neigh <- function(walk_data, vmax, min_buffer=NULL){
  w_x <- w_y <- h <- NULL
  if(!inherits(walk_data, "walk_ddl")) stop("'walk_data' must be a 'walk_ddl' object!")
  # Pre-computation in R before calling C++ likelihood
  cell <- cellx <- x <- y <- NULL
  xy_df <- walk_data$q_r |> select(cell, cellx, x, y)
  dt <- walk_data$times$dt
  N <- length(dt)
  L <- walk_data$L
  res <- c(
    dx = walk_data$q_m |> filter(w_x==1, w_y==0) |> pull(h) |> max(),
    dy = walk_data$q_m |> filter(w_y==1, w_x==0) |> pull(h) |> max()
  )
  if(is.null(min_buffer)) min_buffer <- 2*max(res) + 0.1

  if(is.numeric(vmax)){
    if(length(vmax)>1) stop("The 'vmax' argument must be a single value.")
    buffer <- vmax*sqrt(dt)
    buffer <- pmax(buffer, min_buffer)
  }
  active_indices_list <- vector("list", N)
  
  for(i in 2:N) {
    # 1. Identify non-zero likelihood cells for L(i-1) and L(i)
    cells_prev <- which(L[i-1, ] > 0)
    cells_curr <- which(L[i, ] > 0)
    # 2. Get min/max x, y bounds from both observations
    coords_prev <- filter(xy_df, cellx%in%cells_prev)
    coords_curr <- filter(xy_df, cellx%in%cells_curr)
    x_min <- min(coords_prev$x, coords_curr$x) - buffer[i-1]
    x_max <- max(coords_prev$x, coords_curr$x) + buffer[i-1]
    y_min <- min(coords_prev$y, coords_curr$y) - buffer[i-1]
    y_max <- max(coords_prev$y, coords_curr$y) + buffer[i-1]
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