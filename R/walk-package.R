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
mcheck_cols <- function(Xm){
  ind1 <- apply(Xm, 2, sd)!=0
  v <- colnames(Xm)
  ind2 <- unlist(gregexpr('r_', v))
  ind2 <- (ind2!=1)
  return(ind1 & ind2)
}

rcheck_cols <- function(Xr){
  v <- colnames(Xr)
  ind2 <- unlist(gregexpr('m_', v))
  ind2 <- (ind2!=1)
  return(ind2)
}