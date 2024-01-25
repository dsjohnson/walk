// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


//[[Rcpp::export]]
arma::sp_mat dense_to_sparse(const arma::mat& M, const double& tol=1.0e-8){
  arma::sp_mat out(M);
  out.clean(tol);
  return(out);
}



// [[Rcpp::export]]
arma::mat my_test(const arma::mat v, SEXP Q, double prec, bool renorm=true, bool t2=true, bool checks=true) {
  arma::mat out;
  out = v_exp_Q(v, Q, prec, true, true, true);
  return(out);
}