// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;



// [[Rcpp::export]]
arma::mat v_exp_M(const arma::mat& v, const arma::sp_mat&  M, 
                    const double& prec=1.0e-8) {
  arma::mat out = expQ2::sv_exp_Q(v, M, prec, false, true);
  return out;
}


// // [[Rcpp::export]]
// arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& q_vals,
//                     const int& ns){
//   arma::sp_mat Q(from_to, q_vals, ns, ns);
//   arma::colvec row_sums = Q * ones(ns);
//   Q.diag() = -1.0*row_sums;
//   return Q;
// }

// // [[Rcpp::export]]
// arma::sp_mat load_Q_hp(const arma::umat& from_to, const arma::vec& pi_vals, 
//                        const arma::vec& q_vals, const int& ns){
//   int n = from_to.n_cols;
//   arma::sp_mat Qr(ns,ns);
//   arma::vec qvals_m(n);
//   for(int i=0; i<ns; i++){Qr(i,i) = exp(off_q(i) + Xb_q_r(idx_q_r(i)));}
//   for(int i=0; i<n; i++){qvals_m(i) = exp(Xb_q_m(idx_q_m(i)));}
//   arma::sp_mat Qm(from_to, qvals_m, ns, ns);
//   Qm = normalise(Qm, 1, 1);
//   arma::sp_mat Q = Qr * Qm;
//   Q.diag() = -1.0*Qr.diag();
//   return Q;
// }



// Calculate likelihood ///////////////
// [[Rcpp::export]]
Rcpp::List ctmc_n2ll_arma(const arma::sp_mat& Q, const arma::vec& delta, 
                     const arma::sp_mat& L, const arma::vec dt)
{
  int N = dt.size();
  int ns = Q.n_cols;
  double u = 0.0;
  arma::vec log_lik_v(N);
  arma::sp_mat Qt(ns,ns);

  // Start forward loop
  arma::rowvec v(ns);
  arma::rowvec phi = delta;
  arma::sp_mat P(ns,ns);
  
  // Start Foward alg loop (index = i)
  for(int i=0; i<N; i++){
      Qt = Q*dt(i);
      v = v_exp_M(phi, Qt);
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