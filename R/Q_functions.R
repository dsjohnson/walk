#' @title Get movement transition matrix from fitted moveMMP object
#' @param fit A moveMMPP fitted model object from \code{\link[moveMMPP]{fit_mmpp_dir}}.
#' @param sparse Logical. Should the matrix be returned in a sparse format from the \code{Matrix}
#' package. Defaults to \code{sparse = TRUE}.
#' @author Devin S. Johnson
#' @export
get_Q <- function(fit, sparse=TRUE){
  par <- fit$par
  dl <- fit$data_list
  par_map <- dl$par_map
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  Xb_q_r <- dl$X_q_r %*% beta_q_r
  Xb_q_m <- dl$X_q_m %*% beta_q_m
  from_to <- t(cbind(dl$from, dl$to))
  link <- which(dl$link==c("soft_plus", "log"))
  if(dl$form=="mult"){
    Q <- load_Q_mult(from_to, Xb_q_r, Xb_q_m, dl$ns, 
                     link_r = which(dl$link_r==c("soft_plus", "log")), 
                     link_m = which(dl$link_m==c("soft_plus", "log")), 
                     a_r = dl$a_r, a_m=dl$a_m,
                     dl$norm)
  } else if(dl$form=="add"){
    Q <- load_Q_add(from_to, Xb_q_r, Xb_q_m, dl$ns,  
                    link_r = which(dl$link_r==c("soft_plus", "log")), 
                    link_m = which(dl$link_m==c("soft_plus", "log")), 
                    a_r = dl$a_r, a_m=dl$a_m) 
  } else{
    Q <- load_Q_sde(from_to, Xb_q_r, Xb_q_m, dl$ns, dl$k, dl$a_r) 
  }
  if(!sparse) Q <- as.matrix(Q)
  return(Q)
}

#' @title Get the limiting utilization distribution of the CTMC movement process
#' @param fit A moveMMPP fitted model object from \code{\link[moveMMPP]{fit_mmpp_dir}}.
#' @param hpd A vector of probabilities. Will return columns with highest probability area for each specified probability. E.g., 
#' \code{hpd=c(0.5, 0.95)} will return 2 extra columns with 50 and 95% HPD densities. 
#' @param method Method used for eigen decomposition. One of \code{"lu"} or \code{"arpack"}.
#' @param ... Extra arguments to pass to \code{\link[rARPACK]{eigs}}
#' @author Devin S. Johnson
#' @export
#' @importFrom Matrix t lu expand Matrix
#' @importFrom stats median
#' @importFrom rARPACK eigs
get_lim_ud <- function(fit=NULL, hpd=NULL, method="lu",...){
  
  tQ <- t(get_Q(fit))
  
  if(method=="lu"){
    # browser()
    n <- nrow(tQ)
    LU.decomp=Matrix::lu(tQ)
    LU=Matrix::expand(LU.decomp)
    P=LU$P
    Q=LU$Q
    L=LU$L
    U=LU$U
    ## so t(G)=P'LUQ
    ## max(abs(t(G)-t(P)%*%L%*%U%*%Q)/max(G))
    
    ## unit vector
    e.n=Matrix(0,n,1,sparse=T)
    e.n[n,1]=1
    
    ## solve (UQ)pi=e.n
    U[n,n]=1
    Qpi=Matrix::solve(U,e.n)
    ud=t(Q)%*%Qpi |> as.vector()
    ud <-Re(ud) * sign(median(Re(ud)))
    # ud <- pmax(0,ud)
  } else if(method=="arpack"){
    ud <- rARPACK::eigs(tQ, 3, "LM", sigma=0, ...)$vectors 
    ud <-Re(ud) * sign(median(Re(ud)))
    # ud <- pmax(0,ud)
  } else {
    stop("Unknown method for calculation! Should be either 'lu' or 'arpack'")
  }
  ud <- ud/sum(ud)
  ud <- pmax(0,ud)
  ud <- ud/sum(ud)
  ud <- cbind(fit$data_list$cell_map, ud)
  if(!is.null(hpd)){
    hpd <- round(hpd*100,0)
    ud <- ud[order(ud$ud, decreasing=TRUE),]
    val <- cumsum(ud$ud)
    hpd_df <- NULL
    for(i in 1:length(hpd)){
      hpd_df <- cbind(hpd_df, ud$ud)
      hpd_df[val>hpd[i]/100,i] <- NA
    }
    colnames(hpd_df) <- paste0("hpd_",hpd)
    ud <- cbind(ud, hpd_df)
    ud <- ud[order(ud$cellx),]
  }
  return(ud)
}