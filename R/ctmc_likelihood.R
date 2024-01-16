#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param debug For developers only, leave in the default setting.
#' @param ... Extra wiggle room for ignored arguments.
#' @importFrom stats aggregate
#' @author Devin S. Johnson
#' @export
ctmc_n2ll <- function(par, data_list, debug=0, ...){
  if(debug>0) browser()
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  
  Xb_q_r <- data_list$X_q_r %*% beta_q_r
  Xb_q_m <- data_list$X_q_m %*% beta_q_m
  
  mx <- (aggregate(Xb_q_m, list(data_list$from), max)[,2])[data_list$from+1]
  Xb_q_m <- Xb_q_m-mx
    
  ctmc_n2ll_arma(
    data_list$dt, 
    data_list$cell, 
    data_list$ns, 
    from_to, 
    Xb_q_r, 
    Xb_q_m
  )$n2ll
}



