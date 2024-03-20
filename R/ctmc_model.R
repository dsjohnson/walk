#' @title Define specifics for model structure and fitting control
#' @name arg_funs
#' @aliases ctmc_model
#' @aliases ctmc_control
#' @param form Model formula for CTMC model components, residency or movement
#' @param link A link funtion for the residency or movement model. One of `"soft_plus"` or `"log"`.
#' @param a Scale parameter for the `"soft_plus"` link function. Ignored for `link="log"`.
#' @export
ctmc_model <- function(form=~1, link="soft_plus", a=1){
    list(
      form = form,
      link = link,
      a = a
    )
}

#' @name arg_funs
#' @export
ctmc_control <- function(q_r = ctmc_model(), q_m=ctmc_model(), 
                          p=FALSE, delta=NULL, form="mult", norm=TRUE){
  list(
    q_r = q_r, q_m=q_m, 
    p=p, delta=delta, form=form, norm=norm
  )
}
