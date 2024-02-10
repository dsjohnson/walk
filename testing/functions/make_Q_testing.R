devtools::load_all("~/research/projects/r_packages/walk")
library(terra)

n <- 10

R <- rast(matrix(1:n^2, n, n))

adj_mat <- terra::adjacent(R, 1:ncell(R))
adj_df <- cbind(from_cell = rep(1:ncell(R), each=4), to_cell=as.vector(t(adj_mat)))
from_to <- t(adj_df[is.finite(adj_df[,2]),]) -1 

Xb_q_r <- log(4)*rep(1,ncell(R))
Xb_q_m <- rep(0, ncol(from_to))
Q <- load_Q(from_to, Xb_q_r, Xb_q_m, ns=ncell(R), norm = TRUE)


Xb_q_r <- 4*rep(1,ncell(R))
Xb_q_m <- rep(1, ncol(from_to))
Qsp <-  load_Q_sp(from_to, Xb_q_r, Xb_q_m, ns=ncell(R), a=10, norm = TRUE)
Qsp[55,]


D_val <- 4*rep(1,ncell(R))
Z_val<- rep(1, ncol(from_to))
Q_dz <- load_Q_DZ(from_to, D_val, Z_val, ns=ncell(R), a =10.0)
Q_dz[55,]
