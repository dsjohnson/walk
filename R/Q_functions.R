#' @title Get movement transition matrix from fitted moveMMP object
#' @param fit A moveMMPP fitted model object from \code{\link[moveMMPP]{fit_mmpp_dir}}.
#' @param sparse Logical. Should the matrix be returned in a sparse format from the \code{Matrix}
#' package. Defaults to \code{sparse = TRUE}.
#' @author Devin S. Johnson
#' @export
get_Q <- function(fit, sparse=TRUE){
  dl <- fit$data_list
  Xb_q_r <- dl$X_q_r %*% fit$results$beta$q_r$est
  Xb_q_m <- dl$X_q_m %*% fit$results$beta$q_m$est
  from_to <- t(cbind(dl$from, dl$to))
  Q <- load_Q(from_to, Xb_q_r, Xb_q_m, dl$ns, norm = TRUE)
  if(!sparse) Q <- as.matrix(Q)
  return(Q)
}

#' @title Get the limiting utilization distribution of the CTMC movement process
#' @param fit A moveMMPP fitted model object from \code{\link[moveMMPP]{fit_mmpp_dir}}.
#' @param Q A movement rate matrix. If provided, will ignore the \code{fit} argument.
#' @author Devin S. Johnson
#' @export
#' @importFrom Matrix t
get_lim_ud <- function(fit=NULL, Q=NULL){
  if(is.null(Q)) Q <- get_Q(fit)
  eigen_list <- eigen(Matrix::t(Q))
  idx <- which.min(abs(eigen_list$values))
  u <- eigen_list$vectors[,idx]
  u <- u/sum(u)
  return(u)
}