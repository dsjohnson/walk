
###~~~~~~~~~~~~~~~~~~~~~~~
### Make some testing data
###~~~~~~~~~~~~~~~~~~~~~~~
source("/Users/devin.johnson/research/projects/r_packages/walk/testing/functions/make_test_data.R")


r2 <- ras
values(r2) <- values(r2) + rnorm(ncell(ras)) + 5*(1:ncell(ras))/ncell(ras)

ras <- c(ras, r2)
names(ras) <- c("cov1","cov2")

ddl <- make_design_data(ras, grad="cov2", rast_mask = ras[[1]])
inherits(ddl, "Qdf")


