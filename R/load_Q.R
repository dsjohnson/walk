#' @title Load the rate matrix for the CTMC movelment model
#' @param Xr The design matrix for the residency only covariates
#' @param Xm The design matrix for the movement covariates
#' @param beta_r Coefficient vector for the residency only covariates
#' @param beta_m Coefficent vector for the movement covariates
#' @param nb_idx A 2-column vector of indecies for nonzero elements of the rate matrix. 
#' @param ns Number of sites in the study area.
#' @param separable Logical. Should the separable parameterization of Hewitt et al. (2023)
#' be used. Defaults to `separable = TRUE`
#' @references Hewitt et al. (2023) xxx. Annals of Applied Statistics. 
#' @import Matrix

load_Q <- function(Xr, Xm, beta_r, beta_m, nb_idx, separable){
  if(is.null(Xm)){
    pi_vals = rep(1, nrow(nb_idx))
  } else{
    pi_vals <- exp(Xm %*% beta_m) |> as.vector()
  }
  
  if(separable){
    q_vals <- Matrix::Diagonal(x = exp(Xr %*% beta_r))
    Q <- Matrix::sparseMatrix(i=nb_idx[,1], j=nb_idx[,2], x=pi_vals, repr = "C") 
    Qrs <- Matrix::Diagonal(x=1/Matrix::rowSums(Q) )
    Q <- q_vals %*% (Qrs %*% Q)
    return(Q)
  } else{
    q_vals <- exp(Xr%*%beta_r + Xm%*%beta_m) |> as.vector()
    Q <- Matrix::sparseMatrix(i=nb_idx[,1], j=nb_idx[,2], x=q_vals, repr = "C") 
    return(Q)
  }
}