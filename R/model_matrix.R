#' @title Create sparse representations of design matrices for 
#' movement and resight models
#' @param formula An \code{R} formula object for the paramter vector.
#' @param ddl Design data list
#' @details This function is not designed for end-users but is exported for posterity.
#' @name dm_matrix
NULL

#' @rdname dm_matrix
#' @importFrom stats model.matrix sd
#' @export
dm_q_m <- function(formula, ddl){
  X <- model.matrix(formula, ddl$q_m)
  # if(is.null(par_list$offset)){
  #   par_list$offset <- ~0 + log(1/num_neigh)
  #   offset = model.matrix(par_list$offset, ddl$q)
  #   offset = as.vector(offset)
  # } else if(par_list$offset==0){
  #   offset=rep(0,nrow(X))
  # }else{
  #   offset = model.matrix(par_list$offset, ddl$q)
  #   if(ncol(offset)>1) stop("q offset formula must result in model.matrix of only 1 column.")
  #   offset = as.vector(offset)
  #   offset = ifelse(is.na(offset), 0, offset)
  # }
  keep_col <- !(apply(X, 2, sd)==0)
  X <- X[,keep_col,drop=FALSE]
  #if(ncol(X)==0) X <- NA
  # uX <- unique(X)
  # uX <- uX[rowSums(uX)!=0,,drop=FALSE]
  # dX <- data.frame(cbind(from_cellx=ddl$q$from_cellx, to_cellx=ddl$q$to_cellx, X))
  # duX <- data.frame(cbind(idx_q=1:nrow(uX), uX))
  # mX <- merge(dX, duX)
  # mX <- with(mX, mX[order(from_cellx, to_cellx),])
  # lookup <- mX[,c('from_cellx','to_cellx','idx_q')]
  # return(list(X_q = uX, idx_q=lookup, off_q=offset))
  return(list(X_q_m = X))
}

#' @rdname dm_matrix
#' @export
dm_q_r <- function(formula, ddl){
  X <- model.matrix(formula, ddl$q_r)
  keep_col <- !((colMeans(X)!=1) & (apply(X, 2, sd)==0))
  X <- X[,keep_col,drop=FALSE]
  
  #if(ncol(X)==0) X <- NA
  # uX <- unique(X)
  # uX <- uX[rowSums(uX)!=0,,drop=FALSE]
  # dX <- data.frame(cbind(from_cellx=ddl$q$from_cellx, to_cellx=ddl$q$to_cellx, X))
  # duX <- data.frame(cbind(idx_q=1:nrow(uX), uX))
  # mX <- merge(dX, duX)
  # mX <- with(mX, mX[order(from_cellx, to_cellx),])
  # lookup <- mX[,c('from_cellx','to_cellx','idx_q')]
  # return(list(X_q = uX, idx_q=lookup, off_q=offset))
  return(list(X_q_r = X))
}