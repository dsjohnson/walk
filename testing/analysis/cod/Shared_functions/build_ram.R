build_ram <-
function( model,
          vars ){

  # Extract info
  vars = sapply( vars, FUN=function(char){gsub("-", "", gsub(" ", "", char))} )
  n.paths = nrow(model)
  par.names = model[, 2]
  startvalues = model[,3]

  # EXCERPT FROM `getAnywhere("sem.semmod")`
  heads = from = to = rep(0, n.paths)
  for (p in 1:n.paths) {
    path = sem:::parse.path(model[p, 1])
    heads[p] = abs(path$direction)
    to[p] = path$second
    from[p] = path$first
    if (path$direction == -1) {
      to[p] = path$first
      from[p] = path$second
    }
  }
  missing_vars = setdiff( c(from,to), vars )
  if( length(missing_vars) > 0 ) stop( "Check `build_ram`:", paste0(missing_vars,sep=", ") )

  # Constrct RAM
  ram = data.frame(matrix(0, nrow=p, ncol=5))
  pars = na.omit(unique(par.names))
  ram[, 1] = heads
  ram[, 2] = apply(outer(vars, to, "=="), 2, which)
  ram[, 3] = apply(outer(vars, from, "=="), 2, which)
  par.nos = apply(outer(pars, par.names, "=="), 2, which)
  if(length(par.nos) > 0){
    ram[, 4] = unlist(lapply(par.nos, function(x) if (length(x)==0){0}else{x}))
  }
  ram[, 5] = startvalues
  colnames(ram) = c("heads", "to", "from", "parameter", "start")
  return(ram)
}
