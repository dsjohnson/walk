#' @title Extra diffusion added by the soft-plus drift of the SDE form
#' @description For a fit (or `fit_ctmc(..., fit = FALSE)` object) with `form = "sde"`, computes on every
#' edge the drift `mu` (movement linear predictor) and the extra symmetric rate introduced by the
#' soft-plus drift transform (`ctmc_control(sde_a = )`), relative to the diffusion of the model:
#' `h (f(mu) + f(-mu) - |mu|) / exp(beta_q_r)`. For the default hard plus (`sde_a = Inf`) it is 0.
#' @param fit A fitted model object from \code{\link{fit_ctmc}}.
#' @param sde_a Drift transform scale to evaluate. Defaults to the one used in `fit`; supply other values
#' to see how large `sde_a` must be at the fitted parameters.
#' @param par Parameter vector at which to evaluate. Defaults to `fit$par`.
#' @return A data frame with one row per value of `sde_a`: quantiles of the relative extra diffusion over
#' edges, and the proportion of edges with `sde_a |mu| < 1` (drift in the smoothed region).
#' @export
sde_drift_check <- function(fit, sde_a = NULL, par = fit$par){
  dl <- fit$data_list
  if(!identical(dl$form, "sde")) stop("'fit' must use form = \"sde\".")
  if(is.null(sde_a)) sde_a <- sde_a_of(dl)
  pm <- dl$par_map
  s2 <- exp(as.vector(dl$X_q_r %*% par[pm$beta_q_r]))[dl$from + 1]
  mu <- if(length(pm$beta_q_m)) as.vector(dl$X_q_m %*% par[pm$beta_q_m]) else rep(0, length(dl$from))
  h <- as.vector(dl$hij)
  do.call(rbind, lapply(sde_a, function(a){
    ex <- if(is.finite(a)) h * (sde_drift(mu, a) + sde_drift(-mu, a) - abs(mu)) / s2 else 0 * mu
    q <- stats::quantile(ex, c(.5, .9, 1), names = FALSE)
    data.frame(sde_a = a, rel_extra_median = q[1], rel_extra_q90 = q[2], rel_extra_max = q[3],
               prop_smoothed = if(is.finite(a)) mean(a * abs(mu) < 1) else 0)
  }))
}
