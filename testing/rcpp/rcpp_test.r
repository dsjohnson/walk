library(expQ2)
library(Matrix)
library(magrittr)
library(walk)

v <- matrix(rep(c(rep(1.0e-7,9),1), 10), nrow=1, ncol=100)
walk:::dense_to_sparse(v)

### Walk



