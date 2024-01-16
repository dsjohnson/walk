###~~~~~~~~~~~~~~~~~~~~~~~
### Make some testing data
###~~~~~~~~~~~~~~~~~~~~~~~
source("/Users/devin.johnson/research/projects/r_packages/walk/testing/functions/make_test_data.R")
# source("/Users/devin.johnson/research/projects/r_packages/walk/testing/functions/make_Q_data_test.R")

library(ctmm)
data(pelican)
data <- ctmm::tbind(pelican)
data$COV.x.x <- ifelse(data$COV.x.x==Inf, 20^2, data$COV.x.x)
data$COV.y.y <- ifelse(data$COV.y.y==Inf, 20^2, data$COV.y.y)
L <- walk::telem_to_ras(data, ras)

r2 <- ras
values(r2) <- values(r2) + rnorm(ncell(ras)) + 5*(1:ncell(ras))/ncell(ras)
ras <- c(ras, r2)
names(ras) <- c("cov1","cov2")

ddl <- make_design_data(ras, grad="cov2", rast_mask = ras[[1]])



