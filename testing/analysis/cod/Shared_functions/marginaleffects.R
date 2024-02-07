# Make objects
get_coef.custom_tmb = function(model, param, ...){
  model$parhat[[param]]
}
get_vcov.custom_tmb = function(model, param, ...){
  whichrows = which( names(model$opt$par) == param )
  array(model$opt$SD$cov.fixed[whichrows,whichrows],dim=rep(length(whichrows),2))
}
set_coef.custom_tmb = function(model, newpar, param, ...){
  model$parhat[[param]] <- newpar
  return(model)
}
get_predict.custom_tmb = function(model, newdata, param, center=FALSE, ...){
  # build original model.frame
  frame0 = model.frame( formula=model$formula, data=model$data )
  terms0 = terms( frame0 )
  xlevels = .getXlevels( terms0, frame0 )
  # get new design matrix
  terms1 = delete.response( terms0 )
  frame1 = model.frame( terms1, newdata, xlev=xlevels )
  X_ik = model.matrix( terms1, frame1 )
  gamma_k = get_coef.custom_tmb(model, param)
  yhat_i = X_ik %*% gamma_k
  if(center==TRUE) yhat_i = yhat_i - mean(yhat_i)
  out = data.frame( rowid=seq_along(yhat_i[,1]), estimate=yhat_i )
  return(out)
}
options("marginaleffects_model_classes" = "custom_tmb")

