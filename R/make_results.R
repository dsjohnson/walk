#' @importFrom stats aggregate 
normvec <- function(v, ind){
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
  
  ### Q reals
  # q_r
  X_q_r <- data_list$X_q_r
  beta_q_r <- -1*par[data_list$par_map$beta_q_r]
  if(have_V) V_q_r <- V[data_list$par_map$beta_q_r, data_list$par_map$beta_q_r]
  q_r_vals <- as.vector(exp(X_q_r %*% beta_q_r))
  
  vars_q_r <- unique(c("cell", "cellx", all.vars(model_parameters$q_r)))
  df_q_r <- walk_data$q_r[,vars_q_r]
  df_q_r$real <- round(q_r_vals,3)
  if(have_V){
    xbVbx_q_r <- X_q_r %*% V_q_r %*% t(X_q_r)
    se_real_q_r <- as.vector(q_r_vals) * sqrt(diag(xbVbx_q_r))
    df_q_r$se_real <- round(se_real_q_r,3)
    df_q_r$ci_lower <- round(as.vector(exp(X_q_r %*% beta_q_r - 1.96*se_real_q_r)),3)
    df_q_r$ci_upper <- round(as.vector(exp(X_q_r %*% beta_q_r + 1.96*se_real_q_r)),3)
  }
  df_q_r <- dplyr::distinct(df_q_r)
  
  #q_m
  if(!is.null(data_list$par_map$beta_q_m)){
    X_q_m <- data_list$X_q_m
    beta_q_m <- par[data_list$par_map$beta_q_m]
    if(have_V) V_q_m <- V[data_list$par_map$beta_q_m, data_list$par_map$beta_q_m]
    q_m_vals <- as.vector(exp(X_q_m %*% beta_q_m))
    
    vars_q_m <- unique(c("from_cell","cell", "from_cellx", "cellx", all.vars(model_parameters$q_m)))
    df_q_m <- walk_data$q_m[,vars_q_m]
    df_q_m$real <- round(normvec(q_m_vals, df_q_m$from_cellx),3)
    if(have_V){
      bsamp <- mvnfast::rmvn(2000, beta_q_m, V_q_m)
      q_m_samp <- exp(X_q_m %*% t(bsamp))
      q_m_samp <- apply(q_m_samp, 2, normvec, ind=df_q_m$from_cellx)
      df_q_m$ci_lower <- round(apply(q_m_samp, 1, quantile, prob=0.025),3)
      df_q_m$ci_upper <- round(apply(q_m_samp, 1, quantile, prob=0.975),3)
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
      samp <- plogis(rnorm(5000, logit_p, sqrt(V_p)))
      df_p <- data.frame(parameter = "p", 
                         real=round(plogis(logit_p),3), 
                         se_real=round(sd(samp),3), 
                         ci_lower=round(quantile(samp, 0.025),3), 
                         ci_upper=round(quantile(samp, 0.975),3)
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