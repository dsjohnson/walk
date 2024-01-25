###~~~~~~~~~~~~~~~~~~~~~~~
### Make some testing data
###~~~~~~~~~~~~~~~~~~~~~~~
source("/Users/devin.johnson/research/projects/r_packages/walk/testing/functions/make_test_data.R")
# source("/Users/devin.johnson/research/projects/r_packages/walk/testing/functions/make_Q_data_test.R")

library(ctmm)
data(pelican)
telem_data <- ctmm::tbind(pelican)
telem_data$COV.x.x <- ifelse(telem_data$COV.x.x==Inf, 20^2, telem_data$COV.x.x)
telem_data$COV.y.y <- ifelse(telem_data$COV.y.y==Inf, 20^2, telem_data$COV.y.y)

r2 <- ras
values(r2) <- values(r2) + rnorm(ncell(ras)) + 5*(1:ncell(ras))/ncell(ras)
cell_data <- c(ras, r2)
names(cell_data) <- c("cov1","cov2")


pt_data <- walk::proc_telem(telem_data, cell_data)
walk_data <- make_walk_data(pt_data, cell_data, grad="cov2", rast_mask = ras[[1]])

fit <- fit_ctmc(walk_data,
         model_parameters = list(q_r = ~1, q_m = ~w_x + w_y, p = TRUE, delta=NULL), 
         pen_fun = NULL, hessian=TRUE, get_reals=TRUE, 
         start=list(beta_q_r = -1.3, beta_q_m=c(0.17,-0.19), logit_p=-3.9),
         method="nlminb", 
         fit=TRUE, eq_prec = 1.0e-8, debug=0, control=list(rel.tol=1.0e-6, trace=10) )

Q <- get_Q(fit)
ud <- get_lim_ud(fit)
cell_data$ud <- 0
cell_data$ud[ud$cell] <- ud$ud

cell_data <- wrap(cell_data)

save(list=ls(), file="~/research/projects/r_packages/walk/testing/functions/.RData")

st <- lubridate::ceiling_date(walk_data$times$timestamp[1],"12 hours")
end <- lubridate::floor_date(tail(walk_data$times$timestamp,1),"12 hours")
aux_timestamp <- seq.POSIXt(st, end, by="12 hours")

pred <- predict_ctmc(fit, walk_data, aux_timestamp)



