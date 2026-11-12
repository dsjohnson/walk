#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param debug For developers only, leave in the default setting.
#' @param ... Extra wiggle room for ignored arguments.
#' @importFrom stats aggregate plogis
#' @importFrom Matrix diag
#' @author Devin S. Johnson
#' @export
ctmc_n2ll <- function(par, data_list, debug=0, ...){
  if(debug==1) browser()
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  logit_p <- par[par_map$logit_p]
  
  Xb_q_r <- data_list$X_q_r %*% beta_q_r
  Xb_q_m <- data_list$X_q_m %*% beta_q_m
  
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
  
  if(debug==2) browser()
  
  ctmc_n2ll_arma(
    L = t(data_list$L), #arma code needs the obs on the columns
    dt = data_list$dt, 
    ns = data_list$ns, 
    from_to = from_to, 
    Xb_q_r = Xb_q_r, 
    Xb_q_m = Xb_q_m,
    p = p,
    delta = matrix(delta, nrow=1),
    hij = data_list$h,
    eq_prec = data_list$eq_prec,
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
  
}

#' @title Joint (multi-individual) movement MMPP log-likelihood
#' @description Evaluates the CTMC likelihood for data containing several
#' independent chains (individuals, or segments of individuals split at long gaps)
#' that share one set of movement parameters. The total is the sum of the
#' per-segment log-likelihoods, each segment starting from `delta`.
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood. Must contain
#' `seg_start` (0-based first observation of each segment). If absent the data are
#' treated as one segment.
#' @param by_segment Logical. Return the vector of per-segment -2 log-likelihoods
#' (with attribute `zero_obs`) instead of their sum.
#' @param debug For developers only, leave in the default setting.
#' @param ... Ignored.
#' @details The forward recursion only touches cells that can carry probability
#' (the support of the observation likelihood, expanded by the uniformization
#' steps), which is exact, and segments are evaluated in parallel with OpenMP using
#' `data_list$nthreads` threads. A segment containing an observation whose error support
#' does not overlap the predicted distribution has -2 log-likelihood `Inf`; the first such
#' observation of each segment is returned in `attr(, "zero_obs")` (rows of
#' `walk_data$times`).
#'
#' The uniformization series for \eqn{\exp(Q\,dt)} is truncated so that at most `eq_prec` of
#' the probability mass is left out. For an observation far in the tail of its prediction (e.g. a
#' fast displacement) the omitted terms are exactly the ones that reach it, so the step is then
#' recomputed with more terms until the omitted mass is below `1e-4 / 20` of the observation's
#' predictive probability (down to 1e-300). An observation is therefore only impossible (`Inf`)
#' when it cannot be reached at all (up to underflow), not because the series was cut short.
#' @author Devin S. Johnson, Brett T. McClintock
#' @export
ctmc_n2ll_joint <- function(par, data_list, by_segment=FALSE, debug=0, ...){
  if(debug==1) browser()
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  logit_p <- par[par_map$logit_p]
  
  Xb_q_r <- as.vector(data_list$X_q_r %*% beta_q_r)
  Xb_q_m <- if(length(beta_q_m)) as.vector(data_list$X_q_m %*% beta_q_m) else rep(0, length(data_list$from))
  p <- if(!is.null(par_map$logit_p)) plogis(logit_p) else 0
  
  if(is.character(data_list$delta) && data_list$delta=="stationary"){
    delta <- get_lim_ud(list(par = par, data_list = data_list),)$ud
    delta <- delta/sum(delta)
  } else{
    delta <- data_list$delta
  }
  seg_start <- data_list$seg_start
  if(is.null(seg_start)) seg_start <- 0L
  Lt <- data_list$Lt
  if(is.null(Lt)) Lt <- Matrix::t(data_list$L)
  nthreads <- data_list$nthreads; if(is.null(nthreads)) nthreads <- 1L
  prune <- data_list$prune; if(is.null(prune)) prune <- 0
  drop <- data_list$drop; if(is.null(drop)) drop <- 0
  if(debug==2) browser()
  r <- ctmc_n2ll_joint_arma(
    L = Lt, dt = data_list$dt, seg_start = as.integer(seg_start),
    ns = data_list$ns, from_to = from_to,
    Xb_q_r = Xb_q_r, Xb_q_m = Xb_q_m, p = p,
    delta = matrix(delta, nrow=1), hij = data_list$hij,
    eq_prec = data_list$eq_prec,
    link_r = which(data_list$link_r==c("soft_plus", "log", "logit")),
    a_r = data_list$a_r, l_r = data_list$l_r, u_r = data_list$u_r,
    link_m = which(data_list$link_m==c("soft_plus", "log")),
    a_m = data_list$a_m,
    form = which(data_list$form==c("mult", "sde")),
    k = data_list$k, norm = data_list$norm, clip = data_list$clip,
    nthreads = as.integer(nthreads), prune = prune, drop = drop
  )
  out <- if(by_segment) r$n2ll else sum(r$n2ll)
  if(length(r$zero_obs)) attr(out, "zero_obs") <- r$zero_obs
  out
}

#' @title Joint CTMC -2 log-likelihood with its analytic gradient
#' @description Same value as \code{\link{ctmc_n2ll_joint}}, plus the exact gradient with respect
#' to `par`, computed in the same forward pass (forward-mode differentiation of the
#' uniformization series). One call costs much less than the `length(par) + 1` likelihood
#' evaluations of a finite-difference gradient. Available for `form = "mult"`, `clip = 0`, and a
#' fixed `delta` (`"uniform"` or numeric); see `grad_supported()`.
#' @param par Parameter vector.
#' @param data_list `data_list` element of a \code{\link{fit_ctmc}} object (or of `fit_ctmc(..., fit = FALSE)`).
#' @param by_segment Logical. Return per-segment values (and an `n_par x n_seg` gradient matrix)?
#' @return The -2 log-likelihood with attribute `gradient` (and `zero_obs` if some observation has
#' zero predictive probability; the gradient of such a segment is set to 0).
#' @export
ctmc_n2ll_grad <- function(par, data_list, by_segment = FALSE){
  if(!grad_supported(data_list)) stop("The analytic gradient needs form = \"mult\", clip = 0 and a fixed 'delta' (not \"stationary\").")
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_q_r <- par[par_map$beta_q_r]; beta_q_m <- par[par_map$beta_q_m]
  Xb_q_r <- as.vector(data_list$X_q_r %*% beta_q_r)
  Xm <- data_list$X_q_m
  if(is.null(Xm) || !length(beta_q_m)) Xm <- matrix(0, length(data_list$from), 0)
  Xb_q_m <- if(length(beta_q_m)) as.vector(Xm %*% beta_q_m) else rep(0, length(data_list$from))
  has_p <- !is.null(par_map$logit_p)
  p <- if(has_p) plogis(par[par_map$logit_p]) else 0
  seg_start <- data_list$seg_start; if(is.null(seg_start)) seg_start <- 0L
  Lt <- data_list$Lt; if(is.null(Lt)) Lt <- Matrix::t(data_list$L)
  nthreads <- data_list$nthreads; if(is.null(nthreads)) nthreads <- 1L
  prune <- data_list$prune; if(is.null(prune)) prune <- 0
  drop <- data_list$drop; if(is.null(drop)) drop <- 0
  r <- ctmc_n2ll_grad_arma(
    L = Lt, dt = data_list$dt, seg_start = as.integer(seg_start), ns = data_list$ns,
    from_to = from_to, Xb_q_r = Xb_q_r, Xb_q_m = Xb_q_m,
    X_q_r = as.matrix(data_list$X_q_r), X_q_m = as.matrix(Xm), p = p, has_p = has_p,
    delta = matrix(data_list$delta, nrow = 1), eq_prec = data_list$eq_prec,
    link_r = which(data_list$link_r == c("soft_plus", "log", "logit")),
    a_r = data_list$a_r, l_r = data_list$l_r, u_r = data_list$u_r,
    link_m = which(data_list$link_m == c("soft_plus", "log")), a_m = data_list$a_m,
    norm = data_list$norm, nthreads = as.integer(nthreads), prune = prune, drop = drop)
  # parameter order: beta_q_r, beta_q_m, logit_p (as in par_map)
  ord <- c(par_map$beta_q_r, par_map$beta_q_m, par_map$logit_p)
  G <- r$grad; Gp <- matrix(0, length(par), ncol(G)); Gp[ord, ] <- G
  if(by_segment){ out <- r$n2ll; attr(out, "gradient") <- Gp }
  else { out <- sum(r$n2ll); attr(out, "gradient") <- rowSums(Gp) }
  if(length(r$zero_obs)) attr(out, "zero_obs") <- r$zero_obs
  out
}

#' @rdname ctmc_n2ll_grad
#' @export
grad_supported <- function(data_list){
  identical(data_list$form, "mult") && (is.null(data_list$clip) || data_list$clip == 0) &&
    !is.character(data_list$delta)
}
