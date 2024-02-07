library(expQ2)
library(Matrix)
library(magrittr)
# library(walk)
devtools::load_all("~/research/projects/r_packages/walk")
library(terra)

v <- c(90,1) %>% {t(./sum(.))}
M <- matrix(c(2,1,2,1),2,2); diag(M) <- -diag(M)
M <- Matrix(M) %>% as("dgCMatrix")

walk:::my_test(v, M*0.01, prec=1.0e-10)


R <- rast(nrows=71, ncol=71)
v <- matrix(1/ncell(R), 1, ncell(R))
adj <- adjacent(R, cells=1:ncell(R))
from_to <- cbind(from=rep(1:ncell(R), each=4), to=as.vector(t(adj)))
from_to <- from_to[!is.na(from_to[,2]),] |> t()
Xb_q_m <- rnorm(ncol(from_to), 0, 1)
Xb_q_r <- rnorm(ncell(R), 10, 3)
Qsp <- load_Q_sp(from_to-1, Xb_q_r, Xb_q_m, ns=ncell(R))
Q <- load_Q(from_to-1, Xb_q_r, Xb_q_m, ns=ncell(R))
-min(diag(Qsp))

system.time({walk:::my_test(v,Qsp,1e-8, checks=F)})



