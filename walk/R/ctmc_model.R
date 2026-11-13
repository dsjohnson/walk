#' @title Define specifics for model structure and fitting control
#' @name arg_funs
#' @aliases ctmc_model
#' @aliases ctmc_control
#' @param link A link funtion for the residency or movement model. One of `"soft_plus"` or `"log"`.
#' @param a Scale parameter for the `"soft_plus"` link function. Ignored for `link="log"`. The `a` parameter determines the approximation to 
#' a hard plus function, i.e., as `a` becomes large the soft plus function converges to `g^{-1}(x) = max(0,x)`. For this specification, `a` must be greater than or equal to 1.
#' @param L Lower bound for general logit link
#' @param U Upper bound for general logit link
#' @param q_r A named list for the residency model. Must contain elements: `form`, `link`, and `a`. 
#' The easiest way to construct this is the `ctmc_model` function.
#' @param q_m A named list for the movement portion of the model. See `q_r`.
#' @param p Logical. Should a zero-inflation parameter be used. Necessary when there might be outliers when the location likelihood does not
#' overlap likely locations. 
#' @param delta Vector of probabilities for the initial location. Can be one of `"uniform"` for a uniform distribution over all cells, `"stationary"` to use the 
#' stationary distribution implied by the parameters, or a numeric vector the same length as the number of cells. 
#' @param form For `ctmc_model`, the model formula for the residency or movement component. For `ctmc_control`, the form of the rate matrix entries: `"mult"` (default) for residency times movement, 
#' i.e., q_\{ij\} = q_r(i) * pi_m(i,j), `"sde"` for a Langevin diffusion approximation (requires `q_r` form `~1`; sets `norm = FALSE`), 
#' or `"logistic"` for q_\{ij\} = lambda_i * logistic((x_j - x_i)'b), which has a closed-form limiting distribution (see Details; sets `norm = FALSE`
#' and the movement link to logistic). The exact likelihood gradient used by \code{\link{fit_ctmc}} is available for all forms.
#' @param norm Should the movement portion be adjusted to sum to 1. This is the parameterization suggested by Hewitt et al. (2023). 
#' @param clip Positive value for the maximum movement rate for all cells. All cells with modeled total movement rates above this value will be truncated to `clip`.
#' @param sde_a Drift transform for `form = "sde"` (ignored otherwise). `Inf` (default) uses the hard plus
#' `max(0, x)` of the original model. A finite `sde_a > 0` uses the soft plus
#' `f(x) = max(0, x) + log(1 + exp(-sde_a |x|))/sde_a`, which makes the likelihood smooth in `beta_q_m`.
#' See Details.
#' @details
#' \strong{Soft-plus drift for the SDE form.} With `form = "sde"` the rate from cell i to neighbour j is
#' `q_ij = exp(beta_q_r)/(k h^2) + f(mu_ij)/((k/2) h)`, where `mu_ij` is the linear predictor of the
#' movement model (the drift along the edge) and `h` the edge length. Both transforms satisfy
#' `f(x) - f(-x) = x`, so the net drift between two cells, and hence the meaning of `beta_q_m`
#' (coefficients of the gradient of the log potential), is the same. The soft plus adds
#' `f(mu) + f(-mu) - |mu|` (between 0 and `2 log(2)/sde_a`, largest where `mu` is near 0) of symmetric rate
#' to each pair of edges, i.e. extra diffusion that is not part of the Langevin model. The diffusion
#' parameter `exp(beta_q_r)` is then biased downward to compensate. Relative to the diffusion of the
#' model, the extra diffusion on an edge is `h (f(mu) + f(-mu) - |mu|) / exp(beta_q_r)`, at most
#' `2 log(2) h / (sde_a exp(beta_q_r))`. `sde_a` has units of 1/(units of `mu`), so a value that is large
#' enough for one landscape need not be for another; check a fit with \code{\link{sde_drift_check}}
#' (e.g. require the median relative extra diffusion to be < 0.01). The hard plus (default) gives a likelihood that is continuous but has
#' a kink wherever an edge drift changes sign; the exact gradient of \code{\link{fit_ctmc}} handles this,
#' so the soft plus is mainly useful when optimisation stalls at a kink.
#'
#' \strong{The "logistic" form.} With `form = "logistic"` the rate from cell i to neighbour j is
#' `q_ij = lambda_i * logistic((x_j - x_i)'b)`, where `lambda_i` is the residency model (`q_r`, with
#' coefficients `beta_q_r`) and `x_i'b` is a cell-level potential: the `q_m` formula is evaluated on the
#' cell data (`walk_data$q_r`, so it uses cell covariates, not gradient/edge covariates; the intercept drops
#' out) and the design of edge i -> j is `x_j - x_i`. Because `logistic(d)/logistic(-d) = exp(d)`, the chain
#' is reversible (when every neighbour relation is symmetric, which is checked) and its limiting
#' (utilization) distribution is
#' \deqn{u_i \propto \exp(x_i'b) / \lambda_i .}
#' With the bounded residency link `q_r = ctmc_model(~ z, link = "logit", U = U)` (L = 0), i.e.
#' `lambda_i = U logistic(z_i'a)`, this is `u_i propto exp(x_i'b) (1 + exp(-z_i'a))`: `U` is a hard upper
#' bound on the rate of leaving any cell (each edge rate is below `lambda_i`, so the total leaving rate is below
#' `U` times the number of neighbours), `b` has the usual resource-selection interpretation
#' (`exp(x'b)` when `lambda` is constant, e.g. `q_r` form `~1`), and a positive coefficient in `a` means faster
#' movement out of (shorter residency in), and hence lower use of, cells with large `z`. Compared with
#' `clip`, which caps rates by a kink, the bound here is smooth. The closed form is used by
#' \code{\link{get_lim_ud}} and for `delta = "stationary"` (with its exact derivative), see
#' \code{\link{lim_ud_logistic}}. With `clip > 0` (not applied with the logit residency link, as for
#' `"mult"`) a clipped row keeps the same structure with `lambda_i` replaced by `clip / sum_j logistic(.)`,
#' so the closed form still holds.
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @export
ctmc_model <- function(form=~1, link="soft_plus", a=1, L=0, U=0){
  if(link=="logit" & U==0) stop("U>0 and U>L must be specified when using the logit link.") 
  out <- list(
    form = form, link = link, a = a, L=L, U=U
  )
  class(out) <- c("ctmc_mod","list")
  return(out)
}

#' @name arg_funs
#' @importFrom stats terms
#' @export
ctmc_control <- function(q_r = ctmc_model(), q_m=ctmc_model(), 
                         p=FALSE, delta="uniform", form="mult", norm=TRUE, clip=0, sde_a=Inf){
  if(q_m$link=="logit") stop("'logit' link not implimented for movement portion of the CTMC model. Please select either 'log' or 'soft_plus'")
  form <- match.arg(form, c("mult", "sde", "logistic"))
  if(clip<0) stop("'clip' must be >= 0!")
  if(!is.numeric(sde_a) || length(sde_a) != 1 || is.na(sde_a) || sde_a <= 0) stop("'sde_a' must be a single number > 0 (Inf for the hard plus).")
  if(length(attr(terms(q_r$form),"term.labels"))>0 & form=="sde") stop("The SDE cannot handle q_r forms other than ~1 at this time.")
  if(form=="sde") norm <- FALSE
  if(form=="logistic"){ norm <- FALSE; q_m$link <- "logistic" }
  out <- list(
    q_r = q_r, q_m=q_m, 
    p=p, delta=delta, form=form, norm=norm, clip=clip, sde_a=sde_a
  )
  class(out) <- c("ctmc_cont","list")
  return(out)
}
