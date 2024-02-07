#define arma_64bit_word 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[rcpp::plugins(cpp11)]] 
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// [[Rcpp::export]]
double sp1(const double& x, const double& a){
  double y = std::max(0.0, x) + log1p(exp(-abs(a*x)))/a;
  return y;
} 

//[[Rcpp::export]]
arma::vec soft_plus(const arma::vec& x, const double& a = 1){
  arma::vec out(x);
  for(int i=0; i<x.size(); i++){
    out(i) = sp1(x(i),a);
  }
  return out;
} 


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

// [[Rcpp::export]]
arma::sp_mat load_Q_sp(const arma::umat& from_to,
                    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
                    const int& ns, const bool& norm=true) {
  arma::sp_mat Qr(ns,ns);
  Qr.diag() = soft_plus(Xb_q_r);
  arma::sp_mat Qm(from_to, soft_plus(Xb_q_m), ns, ns);
  if(norm){
    Qm = normalise(Qm, 1, 1);
  }
  Qm.diag().ones();
  Qm.diag() = -1*Qm.diag();
  arma::sp_mat Q = Qr * Qm;
  return Q;
}

