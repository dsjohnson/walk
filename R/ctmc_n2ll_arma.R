#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param ... Extra wiggle room for ignored arguments.
#' @author Devin S. Johnson
#' @export
ctmc_n2ll_arma <- function(par, data_list,...){
from_to_q <- t(cbind(data_list$from, data_list$to))
beta_q <- par[(ncol(data_list$X_l)+1):(ncol(data_list$X_l)+ncol(data_list$X_q))]

ctmc_arma(
  data_list$id, 
  data_list$period, 
  data_list$dt, 
  data_list$cell, 
  data_list$ns, 
  data_list$np, 
  from_to_q, 
  data_list$X_q, 
  data_list$off_q,
  data_list$idx_q, 
  beta_q
  )$n2ll
}



