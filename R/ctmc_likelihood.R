#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param debug For developers only, leave in the default setting.
#' @param ... Extra wiggle room for ignored arguments.
#' @importFrom stats aggregate plogis
#' @author Devin S. Johnson
#' @export
ctmc_n2ll <- function(par, data_list, debug=0, ...){
  if(debug==1) browser()
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  logit_p <- par[par_map$logit_p]
  
  Xb_q_r <- data_list$X_q_r %*% beta_q_r
  Xb_q_m <- data_list$X_q_m %*% beta_q_m
  
  # mx <- (aggregate(Xb_q_m, list(data_list$from), max)[,2])[data_list$from+1]
  # Xb_q_m <- Xb_q_m-mx
  
  # Q <- load_Q(from_to, Xb_q_r, Xb_q_m, ns=data_list$ns, norm = TRUE)
  
  if(!is.null(par_map$logit_p)){
    p = plogis(logit_p)
  } else{
    p = 0
  }
  
  if(data_list$delta=="stationary"){
    delta <- get_lim_ud(list(par = par, data_list = data_list))
    delta <- delta$ud
    delta <- delta/sum(delta)
  }
  
  
  if(debug==2) browser()
  
  ctmc_n2ll_arma(
    L = data_list$L, 
    dt = data_list$dt, 
    ns = data_list$ns, 
    from_to = from_to, 
    Xb_q_r = Xb_q_r, 
    Xb_q_m = Xb_q_m,
    p = p,
    delta = matrix(delta, nrow=1),
    eq_prec = data_list$eq_prec,
    link_r = which(data_list$link_r==c("soft_plus", "log")),
    link_m = which(data_list$link_m==c("soft_plus", "log")),
    form = which(data_list$form==c("mult", "add", "sde")),
    a_r = data_list$a_r,
    a_m = data_list$a_m,
    k = data_list$k,
    norm = data_list$norm
  )$n2ll
  
}



