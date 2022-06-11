library(expQ2)
library(Matrix)
library(magrittr)
library(walk)
v <- c(90,1) %>% {t(./sum(.))}
M <- matrix(c(2,1,2,1),2,2); diag(M) <- -diag(M)
M <- Matrix(M)

walk:::my_test(v, M*0.01, prec=1.0e-10)


d=10000; f=0.03; ones=rep(1,d); v=runif(d)%>% {t(./sum(.))}
Qs=abs(rsparsematrix(d,d,f))
diag(Qs)=0; Qs=round(Qs,1); Qsum=rowSums(Qs); diag(Qs)=0-1*Qsum
out <- as.vector(walk:::my_test(v,Qs,1e-15))

### Walk



M <- as.matrix(1/dist(1:10000))
M[M<0.25] <- 0
# M <- M/4
diag(M) <- -rowSums(M)
M <- Matrix(M) %>% as("dgCMatrix")
v <- matrix(c(1,rep(0,999)), 1, 10000)

walk:::my_test(v,M*500,1e-15) %>% hist()
expQ2::v_exp_Q(v,M, 1.0e-8)


