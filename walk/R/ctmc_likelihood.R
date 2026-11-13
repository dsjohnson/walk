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
    link_r = link_r_code(data_list),
    a_r = data_list$a_r,
    l_r = data_list$l_r,
    u_r = data_list$u_r,
    link_m = link_m_code(data_list),
    a_m = data_list$a_m,
    form = form_code(data_list),
    k = data_list$k,
    norm = data_list$norm,
    clip = data_list$clip,
    sde_a = sde_a_of(data_list)
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
    link_r = link_r_code(data_list),
    a_r = data_list$a_r, l_r = data_list$l_r, u_r = data_list$u_r,
    link_m = link_m_code(data_list),
    a_m = data_list$a_m,
    form = form_code(data_list),
    k = data_list$k, norm = data_list$norm, clip = data_list$clip,
    nthreads = as.integer(nthreads), prune = prune, drop = drop, sde_a = sde_a_of(data_list)
  )
  out <- if(by_segment) r$n2ll else sum(r$n2ll)
  if(length(r$zero_obs)) attr(out, "zero_obs") <- r$zero_obs
  out
}

#' @title Joint CTMC -2 log-likelihood with its analytic gradient
#' @description Same value as \code{\link{ctmc_n2ll_joint}}, plus the exact gradient with respect
#' to `par`, computed in the same forward pass (forward-mode differentiation of the
#' uniformization series). One call costs much less than the `length(par) + 1` likelihood
#' evaluations of a finite-difference gradient. Available for all model forms (`"mult"`,
#' `"sde"` and `"logistic"`), with or without `clip`, and for every `delta` (`"uniform"`, numeric, or
#' `"stationary"`); see `grad_supported()`.
#' @details
#' Rate matrix: for `form = "mult"`, q_ij = r_i m_ij (or r_i m_ij / sum_k m_ik with `norm = TRUE`);
#' for `form = "sde"`, q_ij = exp(Xb_r,i) / (k h_ij^2) + max(0, Xb_m,ij) / ((k/2) h_ij); for
#' `form = "logistic"`, q_ij = lambda_i logistic((x_j - x_i)'b). With
#' `clip > 0`, rows whose total rate q_i exceeds `clip` are scaled by `clip / q_i`.
#'
#' The log-likelihood is continuous but only piecewise differentiable in two cases: the
#' `max(0, x)` of the `"sde"` drift (kink where an edge's drift changes sign) and the `clip`
#' threshold (kink where a row's total rate crosses `clip`). The gradient returned is that of the
#' smooth piece containing `par`; exactly at a kink the drift derivative uses 1/2 (the average of
#' the one-sided derivatives) and a row with q_i = clip is treated as unclipped.
#'
#' With `delta = "stationary"` the initial distribution pi(theta) solves pi Q = 0, sum(pi) = 1, and
#' its derivative solves dpi Q = -pi dQ, sum(dpi) = 0 (one sparse solve with
#' `length(beta_q_r) + length(beta_q_m)` right-hand sides per call). For `form = "logistic"` the closed
#' form pi_i propto exp(x_i'b)/lambda_i and its derivative are used instead (see \code{\link{lim_ud_logistic}}).
#' @param par Parameter vector.
#' @param data_list `data_list` element of a \code{\link{fit_ctmc}} object (or of `fit_ctmc(..., fit = FALSE)`).
#' @param by_segment Logical. Return per-segment values (and an `n_par x n_seg` gradient matrix)?
#' @return The -2 log-likelihood with attribute `gradient` (and `zero_obs` if some observation has
#' zero predictive probability; the gradient of such a segment is set to 0).
#' @export
ctmc_n2ll_grad <- function(par, data_list, by_segment = FALSE){
  if(!grad_supported(data_list)) stop("The analytic gradient is not available for this model.")
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_q_r <- par[par_map$beta_q_r]; beta_q_m <- par[par_map$beta_q_m]
  Xb_q_r <- as.vector(data_list$X_q_r %*% beta_q_r)
  Xm <- data_list$X_q_m
  if(is.null(Xm) || !length(beta_q_m)) Xm <- matrix(0, length(data_list$from), 0)
  Xb_q_m <- if(length(beta_q_m)) as.vector(Xm %*% beta_q_m) else rep(0, length(data_list$from))
  has_p <- !is.null(par_map$logit_p)
  p <- if(has_p) plogis(par[par_map$logit_p]) else 0
  form_i <- form_code(data_list)
  link_r <- link_r_code(data_list)
  link_m <- link_m_code(data_list)
  hij <- if(form_i == 2) as.vector(data_list$hij) else numeric(0)
  k_sde <- if(is.null(data_list$k)) 0 else data_list$k
  clip <- if(is.null(data_list$clip)) 0 else data_list$clip
  X_q_r <- as.matrix(data_list$X_q_r); X_q_m <- as.matrix(Xm)
  if(is.character(data_list$delta) && data_list$delta == "stationary"){
    # same pi as ctmc_n2ll_joint(); d pi from  d pi Q = - pi dQ,  sum(d pi) = 0
    if(ud_closed_ok(data_list)){
      # form = "logistic": closed-form pi and d pi (no sparse solve)
      cf <- lim_ud_logistic(par, data_list, grad = TRUE)
      delta <- cf$ud; ddelta <- cf$dud
    } else {
    delta <- get_lim_ud(list(par = par, data_list = data_list))$ud
    delta <- delta/sum(delta)
    qg <- ctmc_Q_grad_arma(from_to = from_to, Xb_q_r = Xb_q_r, Xb_q_m = Xb_q_m, X_q_r = X_q_r,
                           X_q_m = X_q_m, hij = hij, ns = data_list$ns, link_r = link_r,
                           a_r = data_list$a_r, l_r = data_list$l_r, u_r = data_list$u_r,
                           link_m = link_m, a_m = data_list$a_m, form = form_i, k_sde = k_sde,
                           norm = data_list$norm, clip = clip, pi = delta,
                           sde_a = sde_a_of(data_list))
    ddelta <- stat_dist_grad(qg$Q, qg$B)
    }
  } else {
    delta <- data_list$delta
    ddelta <- matrix(0, data_list$ns, 0)
  }
  seg_start <- data_list$seg_start; if(is.null(seg_start)) seg_start <- 0L
  Lt <- data_list$Lt; if(is.null(Lt)) Lt <- Matrix::t(data_list$L)
  nthreads <- data_list$nthreads; if(is.null(nthreads)) nthreads <- 1L
  prune <- data_list$prune; if(is.null(prune)) prune <- 0
  drop <- data_list$drop; if(is.null(drop)) drop <- 0
  r <- ctmc_n2ll_grad_arma(
    L = Lt, dt = data_list$dt, seg_start = as.integer(seg_start), ns = data_list$ns,
    from_to = from_to, Xb_q_r = Xb_q_r, Xb_q_m = Xb_q_m,
    X_q_r = X_q_r, X_q_m = X_q_m, p = p, has_p = has_p,
    delta = matrix(delta, nrow = 1), eq_prec = data_list$eq_prec,
    link_r = link_r, a_r = data_list$a_r, l_r = data_list$l_r, u_r = data_list$u_r,
    link_m = link_m, a_m = data_list$a_m, norm = data_list$norm,
    hij = hij, form = form_i, k_sde = k_sde, clip = clip, ddelta = ddelta,
    nthreads = as.integer(nthreads), prune = prune, drop = drop, sde_a = sde_a_of(data_list))
  # parameter order: beta_q_r, beta_q_m, logit_p (as in par_map)
  ord <- c(par_map$beta_q_r, par_map$beta_q_m, par_map$logit_p)
  G <- r$grad; Gp <- matrix(0, length(par), ncol(G)); Gp[ord, ] <- G
  if(by_segment){ out <- r$n2ll; attr(out, "gradient") <- Gp }
  else { out <- sum(r$n2ll); attr(out, "gradient") <- rowSums(Gp) }
  if(length(r$zero_obs)) attr(out, "zero_obs") <- r$zero_obs
  out
}

# Derivative of the stationary distribution: solves dpi Q = -B^T (B[, k] = (pi dQ_k)^T) subject to
# sum(dpi) = 0, i.e. t(Q) dpi^T = -B with the last equation replaced by the constraint
# (the same replacement as stat_dist()).
#' @noRd
stat_dist_grad <- function(Q, B){
  if(!ncol(B)) return(B)
  A <- Matrix::t(Q)
  n <- nrow(A)
  A[n, ] <- 1
  rhs <- -B; rhs[n, ] <- 0
  as.matrix(Matrix::solve(methods::as(A, "CsparseMatrix"), rhs))
}

# drift transform for form = "sde" stored in a data_list (Inf = hard plus; also for older objects
# without the element)
#' @noRd
sde_a_of <- function(data_list){
  a <- data_list$sde_a
  if(is.null(a)) Inf else a
}

#' @rdname ctmc_n2ll_grad
#' @export
grad_supported <- function(data_list){
  isTRUE(data_list$form %in% c("mult", "sde", "logistic"))
}

# Integer codes passed to the C++ code. The "logistic" form is built by the "mult" code (form 1)
# with the logistic movement link (link_m 3) and norm = FALSE.
#' @noRd
form_code <- function(data_list) if(identical(data_list$form, "sde")) 2L else 1L
#' @noRd
link_r_code <- function(data_list) match(data_list$link_r, c("soft_plus", "log", "logit"))
#' @noRd
link_m_code <- function(data_list) match(data_list$link_m, c("soft_plus", "log", "logistic"))

# residency rate lambda = g_r(eta) and its derivative, as link_fn() in ctmc_grad_joint.cpp
#' @noRd
link_r_eval <- function(eta, data_list){
  switch(data_list$link_r,
    soft_plus = { a <- data_list$a_r
      list(f = pmax(0, eta) + log1p(exp(-abs(a*eta)))/a, df = stats::plogis(a*eta)) },
    log = { f <- exp(pmin(eta, log(.Machine$double.xmax))); list(f = f, df = f) },
    logit = { s <- stats::plogis(eta); d <- data_list$u_r - data_list$l_r
      list(f = data_list$l_r + d*s, df = d*s*(1 - s)) })
}

# Can the closed-form limiting distribution of form = "logistic" be used? Needs every edge i -> j
# to have its reverse j -> i (checked once in fit_ctmc(), stored as data_list$ud_closed).
#' @noRd
ud_closed_ok <- function(data_list){
  identical(data_list$form, "logistic") && isTRUE(data_list$ud_closed) && !is.null(data_list$X_q_m_cell)
}

#' @title Closed-form limiting distribution of the "logistic" CTMC form
#' @description For `form = "logistic"` (see \code{\link{ctmc_control}}) the CTMC is reversible and its
#' limiting (utilization) distribution is
#' `u_i = exp(x_i'b) / lambda_i / sum_j exp(x_j'b) / lambda_j`, where `x_i'b` is the movement linear
#' predictor of cell `i` (the potential) and `lambda_i` its residency rate. With `clip`, a clipped
#' row is scaled by `clip / q_i`, i.e. its `lambda_i` is replaced by `clip / sum_j logistic((x_j - x_i)'b)`.
#' @param par Parameter vector.
#' @param data_list `data_list` element of a \code{\link{fit_ctmc}} object with `form = "logistic"`.
#' @param grad Logical. Also return the derivative of `u` with respect to `c(beta_q_r, beta_q_m)`?
#' @param clip Logical. Use the `clip` of the model (if any)?
#' @return A list with `ud` (vector over cells, in `cellx` order), `log_potential` (`x_i'b`),
#' `lambda` (effective residency rates) and, if `grad = TRUE`, `dud` (`ns x (n_beta_q_r + n_beta_q_m)`).
#' @export
lim_ud_logistic <- function(par, data_list, grad = FALSE, clip = TRUE){
  dl <- data_list
  if(!identical(dl$form, "logistic")) stop("'data_list' must use form = \"logistic\".")
  pm <- dl$par_map
  X_r <- as.matrix(dl$X_q_r); Xc <- dl$X_q_m_cell
  b <- par[pm$beta_q_m]
  pot <- if(length(b)) as.vector(Xc %*% b) else rep(0, dl$ns)
  lr <- link_r_eval(as.vector(X_r %*% par[pm$beta_q_r]), dl)
  lam <- lr$f
  clip_v <- if(clip && !is.null(dl$clip)) dl$clip else 0
  clipped <- rep(FALSE, dl$ns)
  from <- dl$from + 1L
  if(clip_v > 0 && dl$link_r != "logit"){       # as load_Q(): no clip with the logit residency link
    d <- if(length(b)) as.vector(dl$X_q_m %*% b) else rep(0, length(from))
    sg <- stats::plogis(d)
    S <- as.vector(tapply(sg, factor(from, levels = 1:dl$ns), sum)); S[is.na(S)] <- 0
    clipped <- lam*S > clip_v
    lam[clipped] <- clip_v / S[clipped]
  }
  lu <- pot - log(lam)
  u <- exp(lu - max(lu)); u <- u/sum(u)
  out <- list(ud = u, log_potential = pot, lambda = lam)
  if(grad){
    kr <- ncol(X_r); km <- length(b)
    G <- matrix(0, dl$ns, kr + km)
    G[, seq_len(kr)] <- -(lr$df/lr$f) * X_r * (!clipped)
    if(km){
      G[, kr + seq_len(km)] <- Xc
      if(any(clipped)){
        # clipped: log lambda_eff = log(clip) - log S_i, dS_i = sum_j s(1-s) X_q_m
        dS <- rowsum(sg*(1 - sg) * as.matrix(dl$X_q_m), factor(from, levels = 1:dl$ns), reorder = TRUE)
        dS <- as.matrix(dS)[as.character(1:dl$ns), , drop = FALSE]; dS[is.na(dS)] <- 0
        G[clipped, kr + seq_len(km)] <- G[clipped, kr + seq_len(km)] + dS[clipped, , drop = FALSE]/S[clipped]
      }
    }
    out$dud <- u * sweep(G, 2, colSums(u * G))
  }
  out
}
