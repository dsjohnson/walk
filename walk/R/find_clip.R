#' @title Find Clip Value
#' @param f Fraction of time unit in processed telemetry data
#' @param prob Cumulative probability of residency time
#' @param upper The maximum value for the clip value
#' @description
#' This function allows the user to calculate the clipping value for the CTMC rate matrix Q,
#' such that Prob(residency in a cell < (time unit)/f) = prob at the maximum movement rate. This allows 
#' a definition, or sorts, for non-habitat. That is, if an animal would move so fast through a cell at to almost certainly
#' reside in a cell less than (time unit)/f length of time. Then it is essentially not habitat for the
#' individual.
#' @importFrom stats optimize pexp
#' @export
 find_clip <- function(f, prob=0.9, upper=1000){
   foo <- function(lam){
     2*abs(pexp(lam*f)-prob)
   }
   opt <- optimize(foo, c(0,1000))
   return(opt$minimum)
 }
 