#' @title Get movement transition matrix from fitted moveMMP object
#' @param fit A ctmc fitted model object from \code{\link[ctmc_arma]{ctmc_arma}}.
#' @param sparse Logical. Should the matrix be returned in a sparse format from the \code{Matrix}
#' package. Defaults to \code{sparse = TRUE}.
#' @author Devin S. Johnson
#' @export
get_Q <- function(fit, sparse=TRUE){
  dl <- fit$data_list
  beta_q <- fit$results$beta$q$est
  X_q <- dl$X_q
  Xb_q <- X_q%*%beta_q
  from_to_q <- t(cbind(dl$from_q, dl$to_q))
  idx_q <- dl$idx_q
  off_q <- dl$off_q
  ns <- dl$ns
  Q <- load_Q(from_to_q, idx_q, Xb_q, off_q, ns)
  if(!sparse) Q <- as.matrix(Q)
  return(Q)
}

#' @title Get the limiting utilization distribution of the CTMC movement process
#' @param fit A ctmc fitted model object from \code{\link[ctmc_arma]{ctmc_arma}}.
#' @author Devin S. Johnson
#' @export
#' @importFrom Matrix t
get_lim_ud <- function(fit){
  Q <- get_Q(fit)
  eigen_list <- eigen(Matrix::t(Q))
  idx <- which.min(abs(eigen_list$values))
  u <- eigen_list$vectors[,idx]
  u <- u/sum(u)
  return(u)
}