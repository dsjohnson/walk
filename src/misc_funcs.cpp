#define arma_64bit_word 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[rcpp::plugins(cpp11)]] 
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;



// [[Rcpp::export]]
arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8) {
  arma::mat out = expQ2::sv_exp_Q(phi, lnG, prec, false, true);
  return out;
}


// [[Rcpp::export]]
arma::sp_mat load_Q(const arma::umat& from_to,
                    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
                    const int& ns, const bool& norm=true) {
  arma::sp_mat Qr(ns,ns);
  Qr.diag() = trunc_exp(Xb_q_r);
  arma::sp_mat Qm(from_to, trunc_exp(Xb_q_m), ns, ns);
  if(norm){
    Qm = normalise(Qm, 1, 1);
  }
  Qm.diag().ones();
  Qm.diag() = -1*Qm.diag();
  arma::sp_mat Q = Qr * Qm;
  return Q;
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

