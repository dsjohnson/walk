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
    dt <- diff(times$timestamp) |> `units<-`(time_unit)
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
  
  mx <- (aggregate(Xb_q_m, list(data_list$from), max)[,2])[data_list$from+1]
  Xb_q_m <- Xb_q_m-mx
  
  if(!is.null(par_map$logit_p)){
    p = plogis(logit_p)
  } else{
    p = 0
  }
  if(debug==2) browser()
  out <- ctmc_predict_arma(
    L=Lpred, dt=times$dt, ns=data_list$ns, from_to=from_to, 
    Xb_q_r=Xb_q_r, Xb_q_m=Xb_q_m, p=p, 
    delta = matrix(data_list$delta, nrow=1),
    eq_prec = data_list$eq_prec, ...
  )
  
  if(debug==3) browser()
  
  return(list(local_state_prob=out$local_state_prob, times=times))
  
}
