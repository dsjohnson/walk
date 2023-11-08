#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param ... Extra wiggle room for ignored arguments.
#' @author Devin S. Johnson
#' @export
ctmc_n2ll <- function(par, data_list, ...){
  
  Q <- load_rate_mat()
  
  n2ll <- ctmc_arma(
    Q = , 
    delta = , 
    L = , 
    dt = ,
  )$n2ll
  
  return(n2ll)
}



