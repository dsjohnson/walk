#' @importFrom stats nlminb 
#' @importFrom TMB sdreport
fit_norm <- function(x){
  data_list <- list(model="NormalNLL", x=x)
  obj <- TMB::MakeADFun(data = data_list,
                         parameters = list(mu=0, sigma=1),
                         DLL = "walk_TMBExports") # package's DLL
  opt <- nlminb(obj$par,obj$fn,obj$gr)
  sdr <- sdreport(obj,opt$par)
  return(sdr)
}