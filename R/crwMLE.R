#' @title Fit continuous-time correlated random walk models to animal
#' telemetry dat
#' @export

crwMLE <- function(
  data,
  cov_data,
  mov.model = ~ 1,
  err.model = NULL,
  activity = NULL,
  time.scale = NULL,
  theta,
  fixPar,
  optim_method="optim"
){

}
