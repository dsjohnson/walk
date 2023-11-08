#' @title Fit CTMC movement model
#' @param Q_dd Design data list produced by the function \code{\link{make_Q_data}}.
#' @param Lik_mat Sparse matrix where each row is the likelihood surface for the corresponding location observation.x
#' @param model_parameters Model formula for the residency and movement portions of the model, e.g., 
#'  \code{list(Q = list(q_r=~1, q_m=~1, separable=TRUE), L = NULL)}.
#' @param hessian Logical. Should the Hessian matrix be calculated to obtain the parameter
#' variance-covariance matrix.
#' @param start Optional starting values for the parameter must be a list of the 
#' form \code{list(beta_l=c(), beta_q=c())}.
#' @param method Optimization method. See \code{\link[optimx]{optimr}}
#' @param fit Logical. Should the likelihood be optimized?
#' @param debug Integer from 1-4. Opens browser() at various points in the function call. Mostly for 
#' package developers. 
#' @param ... Additional arguments passed to the optimization function 
#' \code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}.
#' @details
#' The \code{separable = TRUE} element of the \code{Q} list for \code{model.parameters} indicates that the 
#' separable parameteriztion of Hewitt et al. (2023) will be used. If set to \code{separable = FALSE}
#' the traditional loglinear formulation of Johnson et al. (2021) and Hanks et al. (2015) will be used. If  \code{separable = FALSE} the 
#' \code{q_r} formula is ignored. 
#' @references Hanks, E. M., Hooten, M. B., and Alldredge, M. W. (2015) Continuous-time discrete-space models for animal movement. Annals of Applied Statistics. 9:145-165.
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @references Johnson, D. S., Pelland, N. A., and Sterling, J. T. (2021) A Continuous-Time Semi-Markov Model for Animal Movement in a Dynamic Environment. The Annals of Applied Statistics, 15:797-812.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois
#' @export
fit_ctmc <- function(Q_dd, Lik_mat, 
                     model_parameters = list(
                       Q = list(q_r=~1, q_m=~1, separable=TRUE),
                       L = NULL
                     ), 
                     hessian=TRUE, start=NULL, method="nlminb", fit=TRUE, 
                     debug=0, ...){
  
  if(debug==1) browser()
  
  if(!inherits(Q_dd, "Qdf")) stop(" 'Q_dd' is not a Q design data list. See '?walk::make_Q_data'")
  
  separable <- model_parameters$Q$separable
  if(is.null(separable)) separable <- TRUE
  # Design matrices
  Xr <- model.matrix(model_parameters$Q$r_form, Q_dd)
  Xr <- Xr[,rcheck_cols(Xr),drop=FALSE]
  Xm <- model.matrix(model_parameters$Q$m_form, Q_dd)
  Xm <- Xm[,mcheck_cols(Xm),drop=FALSE]
  if(separable){
    Xr <- cbind(cell=Q_dd$r_cell, Xr)
    Xr <- unique(Xr)
    Xr <- Xr[,-1]
  }
  
  data_list <- list(
    Xr = Xr, 
    Xm = Xm,
    L = L,
    nb_idx = as.matrix(Q_dd[,c("r_site_idx","m_site_idx")]),
    sep = separable,
    r_idx <- 1:ncol(Xr),
    m_idx <- c(1:ncol(Xm)) + ncol(Xr)
    )
  
  if(is.null(start)){
    par_list <- list(
      beta_r=rep(0, ncol(Xr)), 
      beta_m=rep(0, ncol(Xm))
    )
  } else{
    par_list=start
  }
  start <- c(par_list$beta_r, par_list$beta_m)
  
  # message('Building model...')
  # foo <- MakeADFun(
  #   data=append(list(model="mmpp"), data_list),
  #   parameters=par_list,
  #   #random=c(),
  #   DLL="moveMMPP_TMBExports"
  # )
  
  if(debug==2) browser()
  
  if(fit){
    message('Optimizing likelihood...')  
    if(debug==2) browser()
    # opt <- nlminb(start=start, objective=mmpp_ll, data_list=data_list, ...)
    opt <- optimx::optimr(par=start, fn=ctmc_n2ll, method=method, data_list=data_list, ...)
    
    if(opt$convergence!=0){
      message("There was a problem with optimization... See output 'optimx' object.")
      # return(list(opt=opt, data_list=data_list))
      hessian <- FALSE
      V <- NULL
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
  
  ### Get real lambda values
  par <- opt$par
  # real <- get_reals(par, V, data_list, ddl, model_parameters)
  # beta <- get_betas(par, V, data_list)
  
  if(!hessian) V <- NULL
  
  out <- list(
    # par = c(beta_l,beta_q),
    par = par,
    vcov = V,
    log_lik = -0.5*opt$value,
    aic = opt$value + 2*length(par),
    # results = list(
    #   beta = beta,
    #   real = real
    # ),
    opt = opt,
    start=start,
    data_list=data_list
  )
  
  if(debug==4) browser()
  
  return(out)
  
}
