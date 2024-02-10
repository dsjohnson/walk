#' @importFrom stats aggregate 
normvec <- function(v, ind){
  v <- as.vector(v)
  v <- v/aggregate(v, list(ind), sum)$x[ind]
}


#' @title Derive real parameter values
#' @description Takes fitted model and data object and produces estimates of the
#' real parameter values, i.e., movement and detection rates
#' @param par Model estimated parameters
#' @param V Optional variance-covariance matrix for the parameters. If none is given
#' no standard errors will be provided for the real parameter estimates
#' @param data_list The data list from a fitted model object
#' @param walk_data The design data list fed into a `fit_mmpp` function.
#' @param model_parameters The model formula for the fitted model.
#' @author Devin S. Johnson
#' @importFrom dplyr distinct
#' @importFrom mvnfast rmvn
#' @importFrom stats quantile plogis rnorm
#' @export

get_reals <- function(par, V=NULL, data_list, walk_data, model_parameters){
  have_V <- !is.null(V)
  if(have_V && any(eigen(V)$values<=0)){
    warning("Parameter variance covariance is not proper!")
    have_V <- FALSE
  }
  
  link <- model_parameters$link
  if(is.null(link)) link <- "soft_plus"
  if(link=="soft_plus"){
    a <- model_parameters$a
    if(is.null(a)) a <- 1
    lfun <- function(x){soft_plus(x, a = a)}
    gr_lfun <- function(x){plogis(a*x)}
  } else{
    lfun <- function(x) exp(x)
    gr_lfun <- function(x) exp(x)
  }
  
  ### Q reals
  # q_r
  X_q_r <- data_list$X_q_r
  beta_q_r <- par[data_list$par_map$beta_q_r]
  if(have_V) V_q_r <- V[data_list$par_map$beta_q_r, data_list$par_map$beta_q_r]
  q_r_vals <- as.vector(lfun(X_q_r %*% beta_q_r))
  
  vars_q_r <- unique(c("cell", "cellx", all.vars(model_parameters$q_r)))
  df_q_r <- walk_data$q_r[,vars_q_r]
  df_q_r$real <- 1/q_r_vals
  if(have_V){
    xbVbx_q_r <- sqrt(diag(X_q_r %*% V_q_r %*% t(X_q_r)))
    df_q_r$se_real <- as.vector(gr_lfun(X_q_r%*%beta_q_r) * xbVbx_q_r * df_q_r$real^2)
    df_q_r$ci_lower <- as.vector(1/lfun(X_q_r %*% beta_q_r + 2*xbVbx_q_r))
    df_q_r$ci_upper <- as.vector(1/lfun(X_q_r %*% beta_q_r - 2*xbVbx_q_r))
  }
  df_q_r <- dplyr::distinct(df_q_r)
  
  # browser()
  
  #q_m
  if(!is.null(data_list$par_map$beta_q_m)){
    X_q_m <- data_list$X_q_m
    beta_q_m <- par[data_list$par_map$beta_q_m]
    if(have_V) V_q_m <- V[data_list$par_map$beta_q_m, data_list$par_map$beta_q_m]
    q_m_vals <- as.vector(lfun(X_q_m %*% beta_q_m))
    
    vars_q_m <- unique(c("from_cell","cell", "from_cellx", "cellx", all.vars(model_parameters$q_m)))
    df_q_m <- walk_data$q_m[,vars_q_m]
    df_q_m$real <- normvec(q_m_vals, df_q_m$from_cellx)
    if(have_V){
      bsamp <- mvnfast::rmvn(2000, beta_q_m, V_q_m)
      q_m_samp <- lfun(X_q_m %*% t(bsamp)) |> matrix(nrow=nrow(X_q_m), ncol=2000)
      q_m_samp <- apply(q_m_samp, 2, normvec, ind=df_q_m$from_cellx)
      df_q_m$se_real <- apply(q_m_samp, 1, sd)
      df_q_m$ci_lower <- apply(q_m_samp, 1, quantile, prob=0.025)
      df_q_m$ci_upper <- apply(q_m_samp, 1, quantile, prob=0.975)
    }
    df_q_m <- dplyr::distinct(df_q_m)
  } else{
    df_q_m <- NULL
  }
  
  # p
  if(!is.null(data_list$par_map$logit_p)){
    qqq <- data_list$par_map$logit_p
    logit_p <- par[qqq]
    if(have_V){
      V_p <- as.matrix(V[qqq,qqq])
      df_p <- data.frame(parameter = "p", 
                         real = plogis(logit_p), 
                         se_real = sqrt(plogis(logit_p)*(1-plogis(logit_p)))*sqrt(V_p) , 
                         ci_lower = plogis(logit_p - 2*sqrt(V_p)), 
                         ci_upper = plogis(logit_p + 2*sqrt(V_p))
                         )
      rownames(df_p) <- NULL
    } else{
      df_p <- data.frame(parameter = "p", est=plogis(logit_p))
    }
  } else{
    df_p <- NULL
  }
  
  
  
  return(list(residency=df_q_r, movement=df_q_m, outlier=df_p))
  
}




#' @title Summarize beta parameter values
#' @description Takes fitted model and produces a summary of the beta estimates. 
#' @param par Model estimated parameters
#' @param V Optional variance-covariance matrix for the parameters. If none is given
#' no standard errors will be provided for the real parameter estimates
#' @param data_list The data list from a fitted model object
#' @author Devin S. Johnson
#' @export

get_betas <- function(par, V=NULL, data_list){
  have_V <- !is.null(V)
  
  ### Get beta Q_r values
  qqq <- data_list$par_map$beta_q_r
  beta_q_r <- par[qqq]
  q_nms <- colnames(data_list$X_q_r)
  if(have_V){
    V_q_r <- as.matrix(V[qqq,qqq])
    df_beta_q_r <- data.frame(parameter = q_nms, est=beta_q_r, se_beta=diag(V_q_r))
  } else{
    df_beta_q_r <- data.frame(parameter = q_nms, est=beta_q_r)
  }
  
  ### Get beta Q_m values
  if(!is.null(data_list$par_map$beta_q_m)){
    qqq <- data_list$par_map$beta_q_m
    beta_q_m <- par[qqq]
    q_nms <- colnames(data_list$X_q_m)
    if(have_V){
      V_q_m <- as.matrix(V[qqq,qqq])
      df_beta_q_m <- data.frame(parameter = q_nms, est=beta_q_m, se_beta=diag(V_q_m))
    } else{
      df_beta_q_m <- data.frame(parameter = q_nms, est=beta_q_m)
    }
  } else{
    df_beta_q_m <- NULL
  }
  ### Get ZI parameter  
  if(!is.null(data_list$par_map$logit_p)){
    qqq <- data_list$par_map$logit_p
    logit_p <- par[qqq]
    if(have_V){
      V_p <- as.matrix(V[qqq,qqq])
      df_logit_p <- data.frame(parameter = "logit_p", est=logit_p, se_beta=diag(V_p))
    } else{
      df_logit_p <- data.frame(parameter = "logit_p", est=logit_p)
    }
  } else{
    df_logit_p <- NULL
  }
  
  return(list(q_r=df_beta_q_r, q_m=df_beta_q_m, p=df_logit_p))
  
}