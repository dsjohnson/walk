#' @title Predict Movement Path From Fitted CTMC Model
#' @param fit A fitted model object produced by \code{\link[walk]{fit_ctmc}}
#' @param walk_data A design data list produced by the function \code{\link{make_walk_data}}.
#' @param aux_timestamp Additional POSIX times for which location prediction is desired.
#' @param debug Developer debugging. 
#' @param ... Aditional arguments passed to internal function \code{ctmc_predict_arma}. 
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois qlogis
#' @export
predict_ctmc <- function(fit, walk_data, aux_timestamp=NULL, debug=0, ...){
  if(debug==1) browser()
  data_list <- fit$data_list
  times <- walk_data$times
  time_unit <- attr(times, "time_unit")
  times$type <- "o"
  # st <- lubridate::ceiling_date(times$timestamp[1],"day")
  # end <- lubridate::floor_date(tail(times$timestamp,1),"day")
  # aug_timestamp <- seq.POSIXt(st, end, by="day")
  if(!is.null(aux_timestamp)){
    aux_ts <- data.frame(timestamp=aux_timestamp, type="p")
    times <- merge(times, aux_ts, by=c("timestamp","type"), all=TRUE)
    dt <- diff(times$timestamp)
    units(dt) <- time_unit
    times$dt <- c(0,dt)
  }
  times <- times[order(times$timestamp, times$type),]
  times <- times[!duplicated(times[,c("timestamp","type")]),]
  
  obs_time <- (times$type=="o")
  Lpred <- Matrix::Matrix(0, nrow=nrow(times), ncol = data_list$ns)
  Lpred[obs_time,] <- data_list$L
  
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  par <- fit$par
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  logit_p <- par[par_map$logit_p]
  
  Xb_q_r <- data_list$X_q_r %*% beta_q_r
  Xb_q_m <- data_list$X_q_m %*% beta_q_m
  
  # mx <- (aggregate(Xb_q_m, list(data_list$from), max)[,2])[data_list$from+1]
  # Xb_q_m <- Xb_q_m-mx
  
  if(!is.null(par_map$logit_p)){
    p <- plogis(logit_p)
  } else{
    p <- 0
  }
  
  if(is.character(data_list$delta) && data_list$delta=="stationary"){
      delta <- get_lim_ud(fit)$ud
      delta <- delta/sum(delta)
  } else{
    delta <- data_list$delta
  }
  
  link <- ifelse(fit$data_list$link=="soft_plus", 1, 0)
  if(debug==2) browser()
  #(L, dt, ns, from_to, Xb_q_r, Xb_q_m, p, delta, eq_prec = 1.0e-8, trunc_tol = 1.0e-8, link = 1L, row_sweep = TRUE)
  out <- ctmc_predict_arma(
    L=Lpred, 
    obs=1.0*(rowSums(Lpred)>0), 
    dt=times$dt, 
    ns=data_list$ns, 
    from_to=from_to, 
    Xb_q_r=Xb_q_r, Xb_q_m=Xb_q_m, p=p, 
    delta = matrix(delta, nrow=1),
    eq_prec = data_list$eq_prec, 
    link_r = which(data_list$link_r==c("soft_plus", "log")),
    link_m = which(data_list$link_m==c("soft_plus", "log")),
    form = which(data_list$form==c("mult", "add", "sde")),
    a_r = data_list$a_r,
    a_m = data_list$a_m,
    k = data_list$k,
    norm = data_list$norm
  )
  
  times$row = 1:nrow(Lpred)
  
  if(debug==3) browser()
  
  return(list(local_state_prob=out$local_state_prob, times=times))
  
}
