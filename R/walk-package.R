#' @rawNamespace useDynLib(walk, .registration=TRUE); useDynLib(walk_TMBExports)
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL


###
### Misc functions
###

#' @title Soft-plus function
#' @param x numeric value
#' @param a scale parameter, `a` must be >1.
#' @export
softPlus <- function(x, a = 1.0) {
  soft_plus(x, a)
}

#' @title General inverse logit function
#' @param x numeric value
#' @param L Lower bound
#' @param U Upper bound
#' @export
gen_invLogit <- function(x, L = 0.0, U = 1.0){
  logit(x, L = L, U =U)
}



#' @title Penlization specification
#' @param type Character vector describing the type of penalty, one of `"lasso"` or `"ridge"`
#' @param group_dm A design matrix for grouping shrinkage 

# mcheck_cols <- function(Xm){
#   ind1 <- apply(Xm, 2, sd)!=0
#   v <- colnames(Xm)
#   ind2 <- unlist(gregexpr('r_', v))
#   ind2 <- (ind2!=1)
#   return(ind1 & ind2)
# }
# 
# rcheck_cols <- function(Xr){
#   v <- colnames(Xr)
#   ind2 <- unlist(gregexpr('m_', v))
#   ind2 <- (ind2!=1)
#   return(ind2)
# }