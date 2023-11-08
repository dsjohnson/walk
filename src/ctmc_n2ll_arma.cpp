// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// function prototypes
arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8);
arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const bool& row_sweep=true);


// Calculate likelihood ///////////////
// [[Rcpp::export]]
Rcpp::List ctmc_n2ll_arma(const arma::vec& id, 
                          const arma::vec& dt, const arma::sp_mat& L, 
                          const int& ns, 
                          const arma::umat& from_to, 
                          const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
                          const arma::vec& delta, 
                          const bool& row_sweep=true)
{
  int N = dt.size();
  double u = 0.0;
  arma::vec log_lik_v(N);
  arma::sp_mat Qt(ns,ns);
  
  arma::sp_mat Q = load_Q(from_to, Xb_q_r, Xb_q_m, ns, row_sweep);

  // Start forward loop
  arma::rowvec v(ns);
  arma::rowvec phi = delta;
  arma::sp_mat P(ns,ns);
  
  // Start Foward alg loop (index = i)
  for(int i=0; i<N; i++){
      Qt = Q*dt(i);
      v = phi_exp_lnG(phi, Qt);
      v = v % L.row(i);
      u = accu(v);
      log_lik_v(i) = log(u);
      phi = v/u;
  } // end i
  
  double n2ll = -2*accu(log_lik_v);
  
  return Rcpp::List::create(
    // Rcpp::Named("log_lik_v") = log_lik_v,
    Rcpp::Named("n2ll") = n2ll
  );
  
}