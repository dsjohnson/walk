###~~~~~~~~~~~~~~~~~~~~~~~
### Make some testing data
###~~~~~~~~~~~~~~~~~~~~~~~
source("/Users/devin.johnson/research/projects/r_packages/walk/testing/functions/make_test_data.R")
source("/Users/devin.johnson/research/projects/r_packages/walk/testing/functions/make_Q_data_test.R")

library(ctmm)
data(pelican)
data <- ctmm::tbind(pelican)
data$COV.x.x <- ifelse(data$COV.x.x==Inf, 20^2, data$COV.x.x)
data$COV.y.y <- ifelse(data$COV.y.y==Inf, 20^2, data$COV.y.y)
L <- walk::telem_to_ras(data, ras)

### Function bits
model_parameters = list(
  Q = list(
    r_form=~r_cov2, 
    m_form=~m_x+m_y+m_cov2_grad, separable=TRUE
    ),
  L = NULL
)

Xr <- model.matrix(model_parameters$Q$r_form, Q_dd)
Xm <- model.matrix(model_parameters$Q$m_form, Q_dd)
Xm <- Xm[,mcheck_cols(Xm),drop=FALSE]
separable=TRUE

