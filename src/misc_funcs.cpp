#define arma_64bit_word 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[rcpp::plugins(cpp11)]] 
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;


double sp1(const double& x, const double& a=1.0){
  double y = std::max(0.0, x) + log1p(exp(-abs(a*x)))/a;
  return y;
} 

//[[Rcpp::export]]
arma::vec soft_plus(const arma::vec& x, const double& a = 1.0){
  arma::vec out(x);
  for(int i=0; i<x.size(); i++){
    out(i) = sp1(x(i),a);
  }
  return out;
} 

//[[Rcpp::export]]
arma::vec hard_plus(const arma::vec& x){
  arma::vec out(x);
  for(int i=0; i<x.size(); i++){
    out(i) = std::max(0.0, x(i));
  }
  return out;
} 


//[[Rcpp::export]]
arma::vec trimmed_logistic(const arma::vec& x, const double& a=0.0, const double& b=1.0) {
  arma::vec out(x);
  for(int i=0; i<x.size(); i++){
    out(i) = a + (b-a)/(1+trunc_exp(-x(i)));
  }
  return(out);
}



// [[Rcpp::export]]
arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8) {
  arma::mat out = expQ2::sv_exp_Q(phi, lnG, prec, false, true);
  return out;
}

// [[Rcpp::export]]
arma::sp_mat load_Q_mult(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const int& link_r=1, const int& link_m=1, const double& a_r=1.0, const double& a_m=1.0, const bool& norm=true) {
  arma::sp_mat Qr(ns,ns);
  arma::sp_mat Q(ns, ns);
  arma::vec Qm_vals;
  if(link_r==1){
    Qr.diag() = soft_plus(Xb_q_r, a_r);
  } else{
    Qr.diag() = trunc_exp(Xb_q_r);
  }
  if(link_m==1){
    Qm_vals = soft_plus(Xb_q_m, a_m);
  } else{
    Qm_vals = trunc_exp(Xb_q_m);
  }
  
  arma::sp_mat Qm(from_to, Qm_vals, ns, ns);
  if(norm){
    Qm = normalise(Qm, 1, 1);
    Qm.diag().ones();
    Qm.diag() *= -1;
    Q = Qr * Qm;
  } else{
    Q = Qr * Qm;
    arma::sp_mat qii = sum(Q,1);
    Q.diag() -= 1*qii;
  }

  return Q;
}

// [[Rcpp::export]]
arma::sp_mat load_Q_add(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const int& link_r=1, const int& link_m=1, const double& a_r=1.0, const double& a_m=1.0) {
  arma::vec aij(from_to.n_cols, fill::ones);
  arma::sp_mat A(from_to, aij, ns, ns);  
  arma::sp_mat Dfill(ns,ns);
  arma::vec Z_vals;
  if(link_r==1){
    Dfill.diag() = soft_plus(Xb_q_r, a_r);
  } else{
    Dfill.diag() = trunc_exp(Xb_q_r);
  }
  if(link_m==1){
    Z_vals = soft_plus(Xb_q_m, a_m);
  } else{
    Z_vals = trunc_exp(Xb_q_m);
  }
  arma::sp_mat D = Dfill * A;
  arma::sp_mat Z(from_to, Z_vals, ns, ns);
  arma::sp_mat Q = D + Z;
  Q.diag() -= sum(Q,1);
  return Q;
}


// [[Rcpp::export]]
arma::sp_mat load_Q_sde(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const double& k,const double& a_r=1.0) {
  
  // make D
  arma::vec aij(from_to.n_cols, fill::ones);
  arma::sp_mat A(from_to, aij, ns, ns);  
  arma::sp_mat Dfill(ns,ns);
  Dfill.diag() = soft_plus(Xb_q_r, a_r);
  arma::sp_mat D = Dfill * A;
  D.diag() -= 1*sum(D,1);
  D /= k;
  // make Z
  arma::sp_mat Z(from_to, hard_plus(Xb_q_m), ns, ns);
  Z.diag() -= 1*sum(Z,1);
  Z /= (k/2.0);
  
  arma::sp_mat Q = D + Z;
  return Q;
}
