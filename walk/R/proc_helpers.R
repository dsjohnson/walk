# Internal helpers shared by proc_telem() and proc_Langevin().
# The two functions only differ in how they read locations / error covariances
# from their input class; individual/segment bookkeeping, the C++ cell
# likelihoods and output assembly live here so the two behave identically.

# Accept a single object or a list of objects (one per individual).
proc_as_list <- function(data, cls, what){
  if(inherits(data, cls)) return(list(x=list(data), is_list=FALSE))
  if(is.list(data) && !is.data.frame(data) && length(data) > 0 &&
     all(vapply(data, function(z) inherits(z, cls), TRUE))) return(list(x=data, is_list=TRUE))
  stop("'data' must be ", what, " or a list of them (one per individual).")
}

# Individual ids for each row of each list element.
# Priority: 'id' column -> 'fallback' (e.g. ctmm identity) -> list names -> list index.
# A single (non-list) object without an 'id' column and without a fallback is individual "1".
proc_ids <- function(xl, id, is_list, fallback=NULL){
  nms <- names(xl)
  out <- vector("list", length(xl))
  for(k in seq_along(xl)){
    d <- xl[[k]]
    if(!is.null(id) && id %in% names(d)){
      v <- as.character(d[[id]])
    } else {
      v <- if(!is.null(fallback)) fallback(d) else NULL
      if(is.null(v) || length(v) != 1 || is.na(v) || !nzchar(v)){
        v <- if(is_list){
          if(!is.null(nms) && !is.na(nms[k]) && nzchar(nms[k])) nms[k] else as.character(k)
        } else "1"
      }
      v <- rep(as.character(v), nrow(d))
    }
    out[[k]] <- v
  }
  if(is_list){
    u <- unlist(lapply(out, unique), use.names=FALSE)
    if(any(duplicated(u)))
      stop("The same individual id (", paste(unique(u[duplicated(u)]), collapse=", "),
           ") appears in more than one list element of 'data'. Use 'id = NULL' to name ",
           "individuals by the list names instead.", call.=FALSE)
  }
  out
}

# Build the 'times' data frame (obs, id, seg, timestamp, dt).
# dt is computed within individual (0 at the first row of each individual and at
# every 'max_gap' break, where a new independent segment starts).
# POSIXct/Date time stamps are converted to 'time_unit'; numeric ones are used as they are.
proc_times <- function(idv, timestamp, time_unit, max_gap){
  n <- length(idv)
  if(any(is.na(idv))) stop("Missing values in the individual ids.")
  if(any(duplicated(rle(idv)$values)))
    stop("Rows for each individual must be contiguous. Sort 'data' by id, then time ",
         "(or give each individual once in a list).")
  if(length(timestamp) != n) stop("Time stamps are missing.")
  if(any(is.na(timestamp))) stop("Missing values in the time stamps.")
  first <- !duplicated(idv)
  if(inherits(timestamp, c("POSIXt", "Date"))){
    timestamp <- as.POSIXct(timestamp)
    dt <- c(0, as.numeric(difftime(timestamp[-1], timestamp[-n], units=time_unit)))
  } else {
    dt <- c(0, diff(as.numeric(timestamp)))
  }
  dt[first] <- 0
  if(any(dt < 0)) stop("Observations are not in time order within individual.")
  new_seg <- first
  if(!is.null(max_gap)){
    if(!is.numeric(max_gap) || length(max_gap) != 1 || is.na(max_gap) || max_gap <= 0)
      stop("'max_gap' must be a single positive number (in the units of 'dt').")
    new_seg <- new_seg | (dt > max_gap)
    dt[new_seg] <- 0
  }
  times <- data.frame(obs=seq_len(n), id=idv, seg=as.integer(cumsum(new_seg)),
                      timestamp=timestamp, dt=dt, stringsAsFactors=FALSE)
  rownames(times) <- NULL
  times
}

# Warn about nthreads when OpenMP is unavailable / limited.
proc_check_threads <- function(nthreads){
  if(nthreads > 1){
    oi <- omp_info()
    if(!oi$enabled){
      warning("'walk' was compiled without OpenMP support, so 'nthreads' is ignored and 1 thread is used.",
              call.=FALSE)
    } else if(oi$thread_limit < nthreads){
      warning("OpenMP thread limit (OMP_THREAD_LIMIT) is ", oi$thread_limit,
              "; using at most that many threads.", call.=FALSE)
    }
  }
  invisible(NULL)
}

# Cell likelihoods + output assembly shared by proc_telem() and proc_Langevin().
# x, y: projected locations (units of cell_data); sxx, syy, sxy: error covariance;
# max_err: per-location neighbourhood radius (NA / <= 0 = no usable error information:
# the location's own cell gets likelihood 1).
proc_lik <- function(x, y, sxx, syy, sxy, max_err, cell_data, times, return_type, trunc,
                     return_inadmissible, nthreads, neighborhood){
  n <- length(x)
  valid <- !is.na(max_err) & max_err > 0
  me  <- ifelse(valid, max_err, 0)
  sxx <- ifelse(valid, sxx, 1)
  syy <- ifelse(valid, syy, 1)
  sxy <- ifelse(valid, sxy, 0)

  # Raster geometry and inadmissible-cell mask (NA in first layer)
  e <- as.vector(terra::ext(cell_data))
  geom <- c(e[["xmin"]], e[["ymax"]], terra::res(cell_data),
            terra::nrow(cell_data), terra::ncol(cell_data))
  mask <- is.na(terra::values(cell_data[[1]], mat=FALSE))

  # Cell likelihoods (computed, masked, normalised and truncated in C++)
  if(neighborhood == "disc"){
    res <- bvn_lik_disc(x, y, sxx, syy, sxy, me, geom, mask, trunc, as.integer(nthreads))
  } else {
    vi <- which(valid)
    obs <- integer(0); cell <- numeric(0)
    if(length(vi)){
      t_buf <- terra::buffer(terra::vect(cbind(x[vi], y[vi]), crs=terra::crs(cell_data)), max_err[vi])
      extr <- terra::cells(cell_data, t_buf, touches = TRUE)
      obs  <- as.integer(vi[extr[, 'ID']])
      cell <- as.numeric(extr[, 'cell'])
      o <- order(obs); obs <- obs[o]; cell <- cell[o]
    }
    res <- bvn_lik_list(obs, cell, x, y, sxx, syy, sxy, me, geom, mask, trunc, as.integer(nthreads))
  }

  zero_ind <- length(res$inad_obs) > 0
  if(return_inadmissible){
    if(!zero_ind) return(NULL)
    return(cbind(obs=res$inad_obs, cell=res$inad_cell, lik=res$inad_lik))
  }
  if(zero_ind) warning("Some observations have likelihood values of 0 in all 'cell_data' cells!")

  ncells <- prod(dim(cell_data)[1:2])
  if(return_type %in% c("sparse", "dense")){
    M <- Matrix::sparseMatrix(i = res$obs, j = res$cell, x = res$lik, dims = c(n, ncells))
    if(return_type == "sparse") M <- as(M, "dgCMatrix") else M <- as.matrix(M)
    return(list(L=M, times=times))
  }
  out <- cbind(obs=res$obs, cell=res$cell, lik=res$lik)
  rm(res)
  if(return_type == "data.frame"){
    out <- as.data.frame(out)
    out <- merge(out, times, by="obs")
    return(out)
  } else if(return_type != "matrix_df"){
    warning("Unknown 'return_type' returning 'matrix_df'")
  }
  out
}
