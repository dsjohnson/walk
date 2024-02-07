library(expQ2)
library(Matrix)
library(magrittr)
library(walk)
library(terra)
v <- c(90,1) %>% {t(./sum(.))}
M <- matrix(c(2,1,2,1),2,2); diag(M) <- -diag(M)
M <- Matrix(M) %>% as("dgCMatrix")

walk:::my_test(v, M*0.01, prec=1.0e-10)


R <- rast(nrows=20, ncol=20)
adj <- adjacent(R, cells=1:ncell(R))
from_to <- cbind(from=rep(1:ncell(R), each=4), to=as.vector(t(adj)))
from_to <- from_to[!is.na(from_to[,2]),] |> t()
Xb_q_m <- rnorm(ncol(from_to), 0, 1)
Xb_q_r <- rnorm(ncell(R), 2, 0.00001)

Q <- load_Q(from_to, Xb_q_r, Xb_q_m, )


# Qs <- zapsmall(Qs)
# diag(Qs)=0-1*Qsum
diag(Qs) <- -1
rate <- 2
Qs <- rate*Qs
out <- as.vector(walk:::my_test(v,Qs,1e-8))

### Walk



M <- as.matrix(1/dist(1:10000))
M[M<0.25] <- 0
# M <- M/4
diag(M) <- -rowSums(M)
M <- Matrix(M) %>% as("dgCMatrix")
v <- matrix(c(1,rep(0,999)), 1, 10000)

walk:::my_test(v,M*500,1e-15) %>% hist()
expQ2::v_exp_Q(v,M, 1.0e-8)


