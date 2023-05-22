#' @title Fit CTMC movement model
#' @param data A processed data frame produced by the function \code{\link{process_data}}
#' @param ddl A design data list produced by the function \code{\link{make_design_data}}.
#' @param model_parameters Model formula for the detection and movement portions
#' of the MMPP model. Must be of the form, e.g., 
#'  \code{list(lambda=list(form=~1, offset=NULL), q=list(form=~1, offset=~log(1.num_neigh)))}. For the 
#'  \code{q} model you must use \code{offset=0} to not have one. If it is left off, \code{offset=~log(1.num_neigh)))}
#'  will be used. 
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
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois
#' @export
ctmc_arma <- function(data, ddl, 
                     model_parameters = list(
                       lambda = list(form=~1, offset=~NULL),
                       q = list(form=~1, offset=~log(1/num_neigh)-1)
                     ), 
                     hessian=TRUE, start=NULL, method="nlminb", fit=TRUE, 
                     debug=0, ...){
  
  cell <- cellx <- fix <- NULL
  
  if(debug==1) browser()
  
  cell_idx_df <- select(ddl$lambda, cell, cellx) %>% distinct()
  data <- data %>% left_join(cell_idx_df, by="cell")
  
  dml_list <- dm_lambda(model_parameters$lambda, ddl)
  dmq_list <- dm_q(model_parameters$q, ddl)
  
  data$period <- ifelse(is.na(data$cell), data$period-1, data$period)
  
  data_list <- list(
    N = as.integer(nrow(data)),
    ns = as.integer(length(unique(ddl$lambda$cell))),
    np = as.integer(max(ddl$quad_pts$period)),
    # detection
    id = data$idx-1,
    period = as.integer(data$period-1),
    dt = data$delta,
    cell = as.integer(data$cellx-1),

    # Q
    from_q = as.integer(dmq_list$idx_q$from_cellx-1),
    to_q = as.integer(dmq_list$idx_q$to_cellx-1),
    X_q = dmq_list$X_q,
    off_q = dmq_list$off_q,
    idx_q = as.integer(dmq_list$idx_q$idx_q-1)
  )
  
  if(is.null(start)){
    par_list <- list(
      beta_l=rep(0,ncol(dml_list$X_l)), 
      beta_q=rep(0, ncol(dmq_list$X_q))
    )
  } else{
    par_list=start
  }
  
  start <- c(par_list$beta_l, par_list$beta_q)
  
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
    opt <- optimx::optimr(par=start, fn=mmpp_ll, method=method, data_list=data_list, ...)
    
    if(opt$convergence!=0){
      message("There was a problem with optimization... See output 'optimx' object.")
      # return(list(opt=opt, data_list=data_list))
      hessian <- FALSE
      V <- NULL
    }
    if(hessian){
      message('Calculating Hessian and variance-covariance matrices...')  
      H <- numDeriv::hessian(mmpp_ll, opt$par, data_list=data_list)
      V <- 2*solve(H)
    } else{
      V <- NULL
      H <- NULL 
    }
  } else{
    hessian <- FALSE
    V <- NULL
    opt <- list(par=start, objective=mmpp_ll(start, data_list))
  }
  
  if(debug==3) browser()
  
  ### Get real lambda values
  par <- opt$par
  real <- get_reals(par, V, data_list, ddl, model_parameters)
  # X_l <- data_list$X_l
  # lll <- 1:ncol(X_l)
  # beta_l <- opt$par[lll]
  # if(hessian) V_l <- V[lll,lll]
  # l_vals <- exp(X_l %*% beta_l)
  # # L <- load_L(data_list$period_l, data_list$cell_l, data_list$idx_l, 
  # #                        data_list$fix_l, l_vals, data_list$ns, data_list$np)
  # vars_l <- unique(
  #   unique(c(
  #     c('cell', 'cellx', 'period', 'fix'), 
  #     all.vars(model_parameters$lambda$form)
  #   ))
  # )
  # df_l <- ddl$lambda[,vars_l]
  # df_l$real <- l_vals[data_list$idx_l+1]
  # df_l$real <- ifelse(is.na(df_l$real), df_l$fix, df_l$real)
  # if(hessian){
  #   xbVbx_l <- X_l %*% V_l %*% t(X_l)
  #   se_real_l <- as.vector(l_vals) * sqrt(diag(xbVbx_l)) 
  #   ci_lower_l <- exp(X_l %*% beta_l - 1.96*se_real_l)
  #   ci_upper_l <- exp(X_l %*% beta_l + 1.96*se_real_l)
  #   df_l$se_real <- se_real_l[data_list$idx_l+1]
  #   df_l$se_real <- ifelse(is.na(df_l$se_real) & !is.na(df_l$fix), 0, df_l$se_real)
  #   df_l$ci_lower <- ci_lower_l[data_list$idx_l+1]
  #   df_l$ci_upper <- ci_upper_l[data_list$idx_l+1]
  # }
  # df_l$prob_det <- ppois(0, df_l$real, lower.tail=FALSE)
  # if(hessian){
  #   df_l$ci_det_prob_lower <- ppois(0, df_l$ci_lower, lower.tail=FALSE)
  #   df_l$ci_det_prob_upper <- ppois(0, df_l$ci_upper, lower.tail=FALSE)
  # }
  # df_l <- dplyr::distinct(df_l)
  
  ### Get beta lambda values
  beta <- get_betas(par, V, data_list)
  # l_nms <- colnames(X_l)
  # if(hessian){
  #   df_beta_l <- data.frame(parameter = l_nms, est=beta_l, se_beta=diag(V_l))
  # } else{
  #   df_beta_l <- data.frame(parameter = l_nms, est=beta_l)
  # }
  # 
  # ### Get real Q  values
  # X_q <- data_list$X_q
  # qqq <- (ncol(X_l)+1):(ncol(X_l)+ncol(X_q))
  # beta_q <- opt$par[qqq]
  # if(hessian) V_q <- V[qqq,qqq]
  # from_to_q <- t(cbind(data_list$from, data_list$to))
  # q_vals <- exp(X_q %*% beta_q)
  # q_nms <- colnames(X_q)
  # # Q <- load_Q(from_to_q, data_list$idx_q, q_vals, data_list$ns)
  # vars_q <- all.vars(model_parameters$q$form)
  # df_q <- ddl$q[,vars_q]
  # df_q$real <- q_vals[data_list$idx_q+1]
  # if(hessian){
  #   xbVbx_q <- X_q %*% V_q %*% t(X_q)
  #   se_real_q <- as.vector(q_vals) * sqrt(diag(xbVbx_q)) 
  #   df_q$se_real <- se_real_q[data_list$idx_q+1]
  #   ci_lower_q <- exp(X_q %*% beta_q - 1.96*se_real_q)
  #   ci_upper_q <- exp(X_q %*% beta_q + 1.96*se_real_q)
  #   df_q$ci_lower <- ci_lower_q[data_list$idx_q+1]
  #   df_q$ci_upper <- ci_upper_q[data_list$idx_q+1]
  # }
  # df_q <- dplyr::distinct(df_q)
  # 
  # ### Get beta Q values
  # q_nms <- colnames(X_q)
  # if(hessian){
  #   df_beta_q <- data.frame(parameter = q_nms, est=beta_q, se_beta=diag(V_q))
  # } else{
  #   df_beta_q <- data.frame(parameter = q_nms, est=beta_q)
  # }
  
  ####
  # statd <- eigen(t(Q))$vectors[,78] /sum(eigen(t(Q))$vectors[,78])
  # zones$ppp <- statd
  # mapview::mapview(zones, zcol='ppp')
  
  if(!hessian) V <- NULL
  
  out <- list(
    # par = c(beta_l,beta_q),
    par = par,
    vcov = V,
    log_lik = -0.5*opt$value,
    aic = opt$value + 2*length(par),
    results = list(
      # beta = list(
      #   lambda = df_beta_l,
      #   q = df_beta_q
      # ),
      beta = beta,
      real = real
    #   real = list(
    #     lambda = df_l,
    #     q = df_q
    #   )
    ),
    opt = opt,
    start=start,
    data_list=data_list
  )
  
  if(debug==4) browser()
  
  return(out)
  
}