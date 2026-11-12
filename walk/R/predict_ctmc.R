#' @title Predict Movement Path From Fitted CTMC Model
#' @param fit A fitted model object produced by \code{\link[walk]{fit_ctmc}}
#' @param walk_data A design data list produced by the function \code{\link{make_walk_data}}.
#' @param aux_timestamp Additional POSIX times for which location prediction is desired.
#' For data with several segments (individuals or `max_gap` bursts) each segment gets the
#' `aux_timestamp` times that fall within its time span.
#' @param trunc_tol Level for truncating small location prediction probabilities. 
#' @param debug Developer debugging. 
#' @param ... Ignored.
#' @details Each independent segment (`walk_data$times$seg`) is smoothed separately with the
#' shared fitted rate matrix, and the results are stacked in the rows of `local_state_prob`.
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois qlogis
#' @export
predict_ctmc <- function(fit, walk_data, aux_timestamp=NULL, trunc_tol=1.0e-8, debug=0, ...){
  if(debug==1) browser()
  data_list <- fit$data_list
  times_all <- walk_data$times
  time_unit <- attr(times_all, "time_unit")
  if(is.null(times_all$seg)) times_all$seg <- 1L
  if(is.null(times_all$id)) times_all$id <- "1"
  times_all$row_obs <- seq_len(nrow(times_all))
  segs <- unique(times_all$seg)
  
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  par <- fit$par
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  logit_p <- par[par_map$logit_p]
  
  Xb_q_r <- data_list$X_q_r %*% beta_q_r
  Xb_q_m <- data_list$X_q_m %*% beta_q_m
  
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
  
  # time differences in the units of the fit
  get_dt <- function(ts){
    if(length(ts) < 2) return(rep(0, length(ts)))
    if(inherits(ts, "POSIXt") && !is.null(time_unit)){
      dt <- diff(ts)
      units(dt) <- time_unit
      c(0, as.numeric(dt))
    } else c(0, diff(as.numeric(ts)))
  }
  
  if(debug==2) browser()
  res_P <- vector("list", length(segs)); res_T <- vector("list", length(segs))
  for(k in seq_along(segs)){
    tk <- times_all[times_all$seg == segs[k], , drop=FALSE]
    tk$type <- "o"
    if(!is.null(aux_timestamp)){
      ak <- aux_timestamp
      if(length(segs) > 1) ak <- ak[ak >= min(tk$timestamp) & ak <= max(tk$timestamp)]
      if(length(ak)){
        aux_ts <- data.frame(timestamp=ak, type="p")
        tk <- merge(tk, aux_ts, by=c("timestamp","type"), all=TRUE)
        tk$seg <- segs[k]; tk$id <- tk$id[!is.na(tk$id)][1]
      }
    }
    tk <- tk[order(tk$timestamp, tk$type),]
    tk <- tk[!duplicated(tk[,c("timestamp","type")]),]
    tk$dt <- get_dt(tk$timestamp)
    obs_time <- (tk$type=="o")
    Lpred <- Matrix::Matrix(0, nrow=nrow(tk), ncol = data_list$ns, sparse=TRUE)
    Lpred[obs_time,] <- data_list$L[tk$row_obs[obs_time], , drop=FALSE]
    res_P[[k]] <- ctmc_predict_arma(
      L=as(Matrix::t(Lpred), "dgCMatrix"), # arma code needs the observations as columns
      dt=tk$dt, 
      ns=data_list$ns, 
      from_to=from_to, 
      Xb_q_r=Xb_q_r, Xb_q_m=Xb_q_m, p=p, 
      delta = matrix(delta, nrow=1),
      hij = data_list$hij,
      eq_prec = data_list$eq_prec, 
      trunc_tol = trunc_tol,
      link_r = which(data_list$link_r==c("soft_plus", "log", "logit")),
      a_r = data_list$a_r,
      l_r = data_list$l_r,
      u_r = data_list$u_r,
      link_m = which(data_list$link_m==c("soft_plus", "log")),
      a_m = data_list$a_m,
      form = which(data_list$form==c("mult", "sde")),
      k = data_list$k,
      norm = data_list$norm,
      clip = data_list$clip
    )
    tk$row_obs <- NULL
    res_T[[k]] <- tk
  }
  out <- do.call(rbind, res_P)
  times <- do.call(rbind, res_T)
  rownames(times) <- NULL
  times$row <- seq_len(nrow(times))
  attr(times, "time_unit") <- time_unit
  
  if(debug==3) browser()
  
  return(list(local_state_prob=out, times=times))
  
}
