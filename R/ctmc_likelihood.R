#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param debug For developers only, leave in the default setting.
#' @param check_rho Check if rho is too big for uniformitazation calculation of exp{Q}. Value is the size of rho to check. 
#' @param ... Extra wiggle room for ignored arguments.
#' @importFrom stats aggregate plogis
#' @importFrom Matrix diag
#' @author Devin S. Johnson
#' @export
ctmc_n2ll <- function(par, data_list, check_rho=NULL, debug=0, ...){
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
  
  if(is.character(data_list$delta) && data_list$delta=="stationary"){
    delta <- get_lim_ud(list(par = par, data_list = data_list))
    delta <- delta$ud
    delta <- delta/sum(delta)
  } else{
    delta <- data_list$delta
  }
  
  if(!is.null(check_rho)){
    if(is.logical(check_rho)) stop("'check_rho' needs to be a numeric value!")
    Q <- get_Q(list(par = par, data_list = data_list))
    qii <- abs(Matrix::diag(Q))
    qii_cell <- c(1:data_list$ns)[which(qii==max(qii))]
    qii <- max(qii)
    rho <- ceiling(qpois(data_list$eq_prec, qii*data_list$dt))
    if(any(rho>=check_rho)){
      # cat("rho is too large for these times: ", c(1:length(data_list$dt))[which(rho>check_rho)], "\n")
      cat("Large rho, cell(s): ", qii_cell, "\n")
      cat("max Qii: ", qii, "\n")
      cat("par: ", par, "\n")
      return(NA)
    }
    
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
    link_r = which(data_list$link_r==c("soft_plus", "log", "logit")),
    a_r = data_list$a_r,
    l_r = data_list$l_r,
    u_r = data_list$u_r,
    link_m = which(data_list$link_m==c("soft_plus", "log")),
    a_m = data_list$a_m,
    form = which(data_list$form==c("mult", "add", "sde")),
    k = data_list$k,
    norm = data_list$norm,
    clip = data_list$clip
  )$n2ll
  
}



