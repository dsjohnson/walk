#' @title Fit CTMC movement model to telemetry data
#' @param walk_data A design data list produced by the function \code{\link{make_walk_data}}.
#' @param model_parameters A named list giving specifications for the model. See \code{\link[walk]{ctmc_control}}. 
#' @param pen_fun An optional penalty function. Should be on the scale of a log-prior distribution.
#' @param hessian Logical. Should the Hessian matrix be calculated to obtain the parameter
#' variance-covariance matrix.
#' @param start Optional starting values for the parameter must be a list of the 
#' form \code{list(beta_l=c(), beta_q_r=c(), beta_q_r=c())}.
#' @param method Optimization method. See \code{\link[optimx]{optimr}}
#' @param fit Logical. Should the likelihood be optimized?
#' @param eq_prec Error rate of matrix exponential calculation. Defaults to \code{1.0e-8}. This is 
#' a generous value. If the model is running slow, you can try reducing it to, say, \code{1.0e-4}.
#' @param nthreads Number of OpenMP threads used to evaluate independent segments
#' (individuals / `max_gap` bursts) of the joint likelihood in parallel. Defaults to 1.
#' @param prune Optional threshold for dropping filtered state probabilities smaller than
#' `prune * max(phi)` during the forward recursion. `0` (default) is exact.
#' @param drop Probability mass that may be dropped per time step while the state distribution
#' is propagated between observations: after each uniformization step, cells holding less than
#' `drop / (number of steps * number of occupied cells)` are set to zero, so at most `drop` of the
#' probability is lost per time step. Defaults to `eq_prec`, the truncation error that the
#' matrix exponential already has. `drop = 0` keeps the exact support, which is much slower for
#' long time steps on large grids.
#' @param gradient Logical. Use the exact gradient of the log-likelihood (computed in the same
#' forward pass as the likelihood, see \code{\link{ctmc_n2ll_grad}}) for the optimization and
#' compute the Hessian by central differences of it? This needs far fewer likelihood evaluations
#' than finite differences. Used when the model allows it (`form = "mult"`, `clip = 0`, `delta`
#' not `"stationary"`) and no `gr` is given in `...`; otherwise finite differences are used.
#' @param debug Integer from 1-4. Opens browser() at various points in the function call. Mostly for 
#' package developers. 
#' @param ... Additional arguments passed to the optimization function 
#' \code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}.
#' @details Data from several individuals (see the `id` argument of \code{\link{proc_Langevin}}
#' and \code{\link{proc_telem}}, and \code{\link{bind_proc_data}}) are fitted as a single joint
#' model: all individuals share the movement parameters, each segment (`walk_data$times$seg`;
#' an individual, or a burst created with `max_gap`) is an independent chain started from
#' `delta`, and the objective is the sum of the per-segment -2 log-likelihoods. A data set
#' without a `seg` column is a single segment. `seg_info` gives the segment number, individual
#' `id`, first row and number of observations of each segment.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois qlogis
#' @export
fit_ctmc <- function(walk_data, 
                     model_parameters = ctmc_control(), 
                     pen_fun = NULL, hessian=TRUE, start=NULL, method="nlminb", 
                     fit=TRUE, eq_prec = 1.0e-8, nthreads = 1, prune = 0, drop = eq_prec,
                     gradient = TRUE, debug=0, ...){
  
  if(debug==1) browser()
  # cell_idx_df <- select(walk_data$q_r, cell, cellx) %>% distinct()
  # data <- data %>% left_join(cell_idx_df, by="cell")
  
  if(is.null(model_parameters$q_r$form) | is.null(model_parameters$q_m$form)) stop("Model formulas must be specified for both residency and movement.")
  X_q_r <- dm_q_r(model_parameters$q_r$form, walk_data)
  X_q_m <- dm_q_m(model_parameters$q_m$form, walk_data)
  
  par_map <- list(beta_q_r = 1:ncol(X_q_r))
  if(ncol(X_q_m)!=0) par_map$beta_q_m <- c(1:ncol(X_q_m)) + ncol(X_q_r)
  
  if(is.null(model_parameters$p)) model_parameters$p <- FALSE
  if(model_parameters$p){
    par_map$logit_p <- ncol(X_q_m) + ncol(X_q_r) + 1
  }
  
  delta <- model_parameters$delta
  if(is.numeric(delta)){
    if(length(delta) != nrow(walk_data$q_r)) stop("The length of 'delta' vector is not equal to the number of cells.")
    delta <- delta/sum(delta)
  }
  if(delta=="uniform"){
    delta <- rep(1,nrow(walk_data$q_r))
    delta <- delta/sum(delta)
  }
  
  if(is.null(model_parameters$q_r$link)){
    link_r <- "soft_plus"
  } else{
    link_r <- model_parameters$q_r$link
  }
  if(is.null(model_parameters$q_m$link)){
    link_m <- "soft_plus"
  } else{
    link_m <- model_parameters$q_m$link
  }
  if(!link_r%in%c("soft_plus","log","logit")) stop("The 'link_r' objects in must be either 'soft_plus', 'log', or 'logit'")
  if(!link_m%in%c("soft_plus","log")) stop("The 'link_m' objects in must be either 'soft_plus' or 'log'.")
  
  
  form <- model_parameters$form
  if(!is.null(form)){
    if(!form%in%c("mult","sde")) stop("The 'form' object in must be either 'mult',or 'sde'.")
  } else{
    form<- "mult"
  }
  
  
  if(is.null(model_parameters$q_r$a)){
    a_r <- 1
  } else{
    a_r <- model_parameters$q_r$a
  }
  
  if(is.null(model_parameters$q_r$U)){
    u_r <- 0
  } else{
    u_r <- model_parameters$q_r$U
  }
  if(is.null(model_parameters$q_r$L)){
    l_r <- 0
  } else{
    l_r <- model_parameters$q_r$L
  }
  if(link_r=="logit"){
    if(l_r < 0) stop("L for the logit link must be >= 0")
    if(u_r <= 0) stop("U for the logit link must be > 0")
    if(u_r <= l_r) stop("U must be > L for the logit link")
  }
  if(is.null(model_parameters$q_m$a)){
    a_m <- 1
  } else{
    a_m <- model_parameters$q_m$a
  }
  
  if(a_r<1 | a_m<1) stop("The 'a' parameter for the 'soft_plus' and 'sde' link functions must be >1.")
  
  if(form %in% c("mult","add")){
    k <- 0
  } else {
    k <- if(attr(walk_data$q_m, "directions")=="rook"){
      k <- 2
    } else if(attr(walk_data$q_m, "directions")=="queen"){
      k <- 4
    } else {
      stop("'sde' model currently only functions with 'SpatRaster' class habitat data.")
    }
  }
  
  norm <- model_parameters$norm
  if(is.null(norm) | !is.logical(norm)) norm <- TRUE
  if(model_parameters$form=="sde") norm <- FALSE
  
  clip <- model_parameters$clip
  if(is.null(clip) || clip<0) clip <- 0
  
  # independent segments (individuals / max_gap bursts); 0-based first row of each
  times <- walk_data$times
  seg <- times$seg
  if(is.null(seg)) seg <- rep(1L, nrow(times))
  seg_start <- which(c(TRUE, seg[-1] != seg[-length(seg)])) - 1L
  if(any(duplicated(seg[seg_start + 1L]))) stop("Segments in 'walk_data$times' are not contiguous.")
  dt <- times$dt
  dt[seg_start + 1L] <- 0
  if(any(dt < 0)) stop("Negative time differences within a segment. Re-run proc_telem() or proc_Langevin() with the 'id' argument.")
  
  data_list <- list(
    N = nrow(walk_data$L),
    seg_start = as.integer(seg_start),
    Lt = as(Matrix::t(walk_data$L), "dgCMatrix"),
    nthreads = as.integer(nthreads),
    prune = prune,
    drop = drop,
    ns = nrow(walk_data$q_r),
    dt = dt,
    L = walk_data$L,
    delta = delta,
    hij=walk_data$q_m$h,
    ### Q
    from = as.integer(walk_data$q_m$from_cellx-1),
    to = as.integer(walk_data$q_m$cellx-1),
    X_q_r = X_q_r,
    X_q_m = X_q_m,
    par_map = par_map,
    eq_prec = eq_prec,
    link_r=link_r,
    a_r = a_r,
    l_r = l_r,
    u_r = u_r,
    link_m=link_m,
    a_m = a_m,
    form=form,
    k = k,
    norm = norm,
    clip = clip,
    cell_map = walk_data$q_r[,c("cell","cellx")]
  )
  
  
  if(is.null(start$beta_q_r)) start$beta_q_r <- rep(0, ncol(X_q_r))
  if(is.null(start$beta_q_m)) start$beta_q_m <- rep(0, ncol(X_q_m))
  if(is.null(start$logit_p) & model_parameters$p) start$logit_p <- qlogis(0.05)
  
  par_start <- c(start$beta_q_r, start$beta_q_m, start$logit_p)
  
  # ctmc_n2ll(par, data_list)
  
  pen <- if(is.null(pen_fun)) function(par) 0 else pen_fun
  # The value at the last gradient point is cached: optimisers usually ask for f and gr at the
  # same parameters, and ctmc_n2ll_grad() returns both.
  cache <- new.env(parent = emptyenv())
  use_grad <- isTRUE(gradient) && grad_supported(data_list) && !("gr" %in% names(list(...)))
  # Rate guard: the cost of one evaluation grows with the uniformization rate (the fastest total
  # leaving rate), and an optimiser trial step can jump to a rate many times higher than anything
  # seen so far, where a single evaluation can take hours on a large grid. Trial points whose
  # rate exceeds 'rate_guard' times the largest rate at an evaluated point get the value Inf
  # without computing the likelihood, so the optimiser shortens the step.
  rate_guard <- 10
  guard_hit <- make_rate_guard(data_list, rate_guard)
  gr_fun <- function(par, data_list, ...){
    if(!is.null(cache$par) && identical(par, cache$par) && !is.null(cache$g)) return(cache$g)
    if(guard_hit(par)){ cache$par <- par; cache$f <- Inf; cache$g <- rep(NA_real_, length(par)); return(cache$g) }
    v <- ctmc_n2ll_grad(par, data_list)
    g <- attr(v, "gradient")
    if(!is.null(pen_fun)) g <- g - 2*numDeriv::grad(pen_fun, par)
    cache$par <- par; cache$f <- as.numeric(v) - 2*pen(par)
    if(!is.finite(cache$f)) g[] <- NA_real_     # zero likelihood: let the optimiser back off
    cache$g <- g
    g
  }
  # With the analytic gradient the objective value comes from the same computation as the
  # gradient, so the optimiser sees one consistent function.
  obj_fun <- function(par, data_list, ...){
    if(!is.null(cache$par) && identical(par, cache$par)) return(cache$f)
    if(use_grad){ gr_fun(par, data_list); return(cache$f) }
    if(guard_hit(par)) return(Inf)
    as.numeric(ctmc_n2ll_joint(par, data_list, debug=0)) - 2*pen(par)
  }
  
  if(debug==2) browser()
  
  # opt <- optimx::optimr(par=par_start, fn=obj_fun, method=method, data_list=data_list, control = list(trace=10, rel.tol=1.0e-3))
  # ctmc_n2ll(par_start, data_list, debug=1)
  
  if(fit){
    message('Optimizing likelihood...')  
    if(debug==2) browser()
    opt <- if(use_grad){
      optimx::optimr(par=par_start, fn=obj_fun, gr=gr_fun, method=method, data_list=data_list, ...)
    } else {
      optimx::optimr(par=par_start, fn=obj_fun, method=method, data_list=data_list, ...)
    }
    
    if(opt$convergence!=0){
      warning("There might be a problem with optimization... See output 'optimx' object.")
      # return(list(opt=opt, data_list=data_list))
      # hessian <- FALSE
      # V <- NULL
    }
    if(hessian){
      message('Calculating Hessian and variance-covariance matrices...')  
      H <- if(use_grad) hess_from_grad(function(b) gr_fun(b, data_list), opt$par) else
        numDeriv::hessian(obj_fun, opt$par, data_list=data_list)
      V <- 2*solve(H)
    } else{
      V <- NULL
      H <- NULL 
    }
  } else{
    hessian <- FALSE
    V <- NULL
    opt <- list(par=par_start, objective=obj_fun(par_start, data_list))
  }
  
  if(debug==3) browser()
  
  par <- as.vector(opt$par)
  beta <- get_betas(par, V, data_list)
  
  if(!hessian){
    V <- NULL
    H <- NULL
  }
  
  seg_info <- data.frame(seg = seg[seg_start + 1L],
                         id = if(is.null(times$id)) "1" else as.character(times$id[seg_start + 1L]),
                         start = seg_start + 1L, n = diff(c(seg_start, length(seg))),
                         stringsAsFactors = FALSE)
  
  out <- list(
    par = par,
    n_seg = length(seg_start),
    seg_info = seg_info,
    vcov = V,
    log_lik = -0.5*opt$value,
    aic = opt$value + 2*length(par),
    hessian = 0.5*H,
    results = beta, #list(beta = beta, real = reals),
    opt = opt,
    start = start,
    data_list = data_list,
    pen_fun = pen_fun,
    model_parameters = model_parameters
  )
  
  if(debug==4) browser()
  
  return(out)
  
}

# Hessian by central differences of an exact gradient (symmetrised); step 1e-4 * max(|x|, 1)
#' @noRd
hess_from_grad <- function(gfun, x, eps = 1e-4){
  n <- length(x); H <- matrix(0, n, n)
  for(k in seq_len(n)){
    h <- eps * max(abs(x[k]), 1); e <- numeric(n); e[k] <- h
    H[, k] <- (gfun(x + e) - gfun(x - e)) / (2 * h)
  }
  (H + t(H)) / 2
}

# Rate guard used by fit_ctmc(): returns a function g(par) that is TRUE (reject) if the
# uniformization rate at 'par' exceeds 'factor' times the largest rate at a point accepted so
# far (the first point is always accepted).
#' @noRd
make_rate_guard <- function(data_list, factor = 10){
  amax <- NA_real_
  function(par){
    a <- unif_rate(par, data_list)
    if(is.na(amax)) { amax <<- a; return(FALSE) }
    if(is.finite(a) && a <= factor * amax) { amax <<- max(amax, a); return(FALSE) }
    TRUE
  }
}

# Uniformization rate alpha_0 = max_i |Q_ii| at 'par' (the cost of one likelihood evaluation is
# proportional to it).
#' @noRd
unif_rate <- function(par, data_list){
  a <- tryCatch(max(abs(Matrix::diag(get_Q(list(par = par, data_list = data_list))))),
                error = function(e) Inf)
  if(is.na(a)) Inf else a
}
