#' @rawNamespace useDynLib(walk, .registration=TRUE);
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
    paste(package, version, paste("(",date, ")", sep=""))
  )
}



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


#' @import Matrix
create_Q <- function(n, density = 0.2, rate_range = c(0.1, 1.0)) {
  if (n < 2) stop("n must be at least 2")
  
  # Number of possible off-diagonal elements
  num_off_diag <- n * (n - 1)
  num_to_fill <- ceiling(num_off_diag * density)
  
  # Generate random (i, j) pairs where i ≠ j
  possible_indices <- which(matrix(1, n, n) - diag(n) == 1, arr.ind = TRUE)
  sampled_indices <- possible_indices[sample(nrow(possible_indices), num_to_fill), , drop = FALSE]
  
  # Generate random rates
  rates <- runif(num_to_fill, rate_range[1], rate_range[2])
  
  # Create sparse matrix with off-diagonal values
  Q <- sparseMatrix(
    i = sampled_indices[, 1],
    j = sampled_indices[, 2],
    x = rates,
    dims = c(n, n)
  )
  
  # Set diagonals so each row sums to zero
  row_sums <- rowSums(Q)
  Q <- Q - Diagonal(n, row_sums)
  
  return(Q)
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