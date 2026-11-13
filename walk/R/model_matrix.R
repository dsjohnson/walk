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
  return(X)
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
  return(X)
}

#' @rdname dm_matrix
#' @details `dm_q_cell` (used by `form = "logistic"`) evaluates the movement formula on the cell data
#' (`ddl$q_r`), giving the cell-level potential covariates `x_i`; the intercept and constant columns
#' are removed (they cancel in the differences `x_j - x_i`).
#' @export
dm_q_cell <- function(formula, ddl){
  X <- tryCatch(model.matrix(formula, ddl$q_r), error = function(e)
    stop("For form = 'logistic' the 'q_m' formula must use cell covariates (columns of walk_data$q_r), ",
         "not edge covariates such as gradients: ", conditionMessage(e), call. = FALSE))
  if(nrow(X) != nrow(ddl$q_r)) stop("Missing values in the cell covariates of the 'q_m' formula.")
  keep_col <- apply(X, 2, function(x) diff(range(x)) > 0)
  X[, keep_col, drop=FALSE]
}
