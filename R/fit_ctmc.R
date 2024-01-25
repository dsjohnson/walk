#' @title Fit CTMC movement model to telemetry data
#' @param walk_data A design data list produced by the function \code{\link{make_walk_data}}.
#' @param model_parameters Model formula for the detection and movement portions
#' of the MMPP model. 
#' @param pen_fun An optional penalty function. Should be on the scale of a log-prior distribution.
#' @param hessian Logical. Should the Hessian matrix be calculated to obtain the parameter
#' variance-covariance matrix.
#' @param get_reals Calculate real values for expected residency, cell transition probabilities, and 
#' outlier proportion for observations.  
#' @param start Optional starting values for the parameter must be a list of the 
#' form \code{list(beta_l=c(), beta_q_r=c(), beta_q_r=c())}.
#' @param method Optimization method. See \code{\link[optimx]{optimr}}
#' @param fit Logical. Should the likelihood be optimized?
#' @param eq_prec Error rate of matrix exponential calculation. Defaults to \code{1.0e-8}. This is 
#' a generous value. If the model is running slow, you can try reducing it to, say, \code{1.0e-4}.
#' @param debug Integer from 1-4. Opens browser() at various points in the function call. Mostly for 
#' package developers. 
#' @param ... Additional arguments passed to the optimization function 
#' \code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}.
#' @details
#' Two model forms are available \code{list(lambda=list(form=~1, offset=NULL), q=list(form=~1, offset=~log(1.num_neigh)))}. For the 
#' \code{q} model you must use \code{offset=0} to not have one. If it is left off, \code{offset=~log(1.num_neigh)))} 
#' will be used. To use the movement rate form of Hewitt et al. (2023), one must use, e.g., 
#' `q = list(res_form = ~from_var, mov_form=~to_var, offset=...)`, where `from_var` is a variable
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois qlogis
#' @export
fit_ctmc <- function(walk_data, 
                     model_parameters = list(q_r = ~1, q_m = ~1, p = FALSE, delta=NULL), 
                     pen_fun = NULL, hessian=TRUE, get_reals=FALSE, start=NULL, method="nlminb", 
                     fit=TRUE, eq_prec = 1.0e-8, debug=0, ...){
  
  if(debug==1) browser()
  # cell_idx_df <- select(walk_data$q_r, cell, cellx) %>% distinct()
  # data <- data %>% left_join(cell_idx_df, by="cell")
  
  X_q_r <- dm_q_r(model_parameters$q_r, walk_data)
  X_q_m <- dm_q_m(model_parameters$q_m, walk_data)
  
  par_map <- list(beta_q_r = 1:ncol(X_q_r))
  if(ncol(X_q_m)!=0) par_map$beta_q_m <- c(1:ncol(X_q_m)) + ncol(X_q_r)
  
  if(model_parameters$p){
    par_map$logit_p <- ncol(X_q_m) + ncol(X_q_r) + 1
  }

  if(is.null(model_parameters$delta)){
    delta <- walk_data$L[1,]
  }
  delta <- delta/sum(delta)
  
  data_list <- list(
    N = nrow(walk_data$L),
    ns = nrow(walk_data$q_r),
    dt = walk_data$times$dt,
    L = walk_data$L,
    delta = delta,
    ### Q
    from = as.integer(walk_data$q_m$from_cellx-1),
    to = as.integer(walk_data$q_m$cellx-1),
    X_q_r = X_q_r,
    X_q_m = X_q_m,
    par_map = par_map,
    eq_prec = eq_prec,
    cell_map = walk_data$q_r[,c("cell","cellx")]
  )
  
  
  if(is.null(start$beta_q_r)) start$beta_q_r <- rep(0, ncol(X_q_r))
  if(is.null(start$beta_q_m)) start$beta_q_m <- rep(0, ncol(X_q_m))
  if(is.null(start$logit_p) & model_parameters$p) start$logit_p <- qlogis(0.05)
  
  par_start <- c(start$beta_q_r, start$beta_q_m, start$logit_p)
  
  # ctmc_n2ll(par, data_list)
  
  if(is.null(pen_fun)){
    obj_fun <- function(par, data_list, ...){ctmc_n2ll(par, data_list, debug=0, ...)}
  } else{
    obj_fun <- function(par, data_list, ...){ctmc_n2ll(par, data_list, debug=0, ...) - 2*pen_fun(par)}
  }
  
  if(debug==2) browser()
  
  # opt <- optimx::optimr(par=par_start, fn=obj_fun, method=method, data_list=data_list, control = list(trace=10, rel.tol=1.0e-3))
  
  
  if(fit){
    message('Optimizing likelihood...')  
    if(debug==2) browser()
    opt <- optimx::optimr(par=par_start, fn=obj_fun, method=method, data_list=data_list, ...)
    
    if(opt$convergence!=0){
      message("There was a problem with optimization... See output 'optimx' object.")
      # return(list(opt=opt, data_list=data_list))
      # hessian <- FALSE
      # V <- NULL
    }
    if(hessian){
      message('Calculating Hessian and variance-covariance matrices...')  
      H <- numDeriv::hessian(ctmc_n2ll, opt$par, data_list=data_list)
      V <- 2*solve(H)
    } else{
      V <- NULL
      H <- NULL 
    }
  } else{
    hessian <- FALSE
    V <- NULL
    opt <- list(par=start, objective=ctmc_n2ll(start, data_list))
  }
  
  if(debug==3) browser()
  
  par <- as.vector(opt$par)
  beta <- get_betas(par, V, data_list)
  if(get_reals){
    reals <- get_reals(par, V, data_list, walk_data, model_parameters)
  } else{
    reals <- NULL
  }
  
  
  if(!hessian) V <- NULL
  
  out <- list(
    par = par,
    vcov = V,
    log_lik = -0.5*opt$value,
    aic = opt$value + 2*length(par),
    results = list(beta = beta, real = reals),
    opt = opt,
    start = start,
    data_list = data_list,
    pen_fun = pen_fun,
    model_parameters = model_parameters
  )
  
  if(debug==4) browser()
  
  return(out)
  
}