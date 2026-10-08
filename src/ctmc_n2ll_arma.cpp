// [[Rcpp::depends(RcppArmadillo)]]

#include "walk_types.h"

using namespace Rcpp;
using namespace arma;

// Calculate likelihood ///////////////

// [[Rcpp::export]]
double ctmc_n2ll_arma(
    const arma::sp_mat& L, 
    const arma::vec& dt, 
    int ns,
    const arma::umat& from_to, 
    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
    double p,
    const arma::rowvec& delta, 
    const arma::vec& hij,
    double eq_prec,
    int link_r,
    double a_r, 
    double l_r,
    double u_r,
    int link_m,
    double a_m, 
    int form,
    double k,
    bool norm,
    double clip)
{
  int N = dt.size();
  double u = 0.0;
  arma::vec log_lik_v(N);
  arma::sp_mat Q;
  if(form==1){
    Q = load_Q(from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, l_r, u_r, link_m, a_m, norm, clip);
    // } else if(form==2){
    //   Q = load_Q_add(from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, link_m, a_m, clip);
  } else {
    Q = load_Q_sde(from_to, Xb_q_r, Xb_q_m, hij, ns, k, clip);
  }
  
  // Start forward loop
  arma::sp_mat Li = L.col(0).t();
  arma::rowvec v = delta % ((1-p)*Li) + (p/ns)*delta;
  u = accu(v);
  log_lik_v(0) = log(u);
  arma::rowvec phi = v/u;
  // Start Forward alg loop (index = i)
  for(int i=1; i<N; i++){
    Li = L.col(i).t();
    v = v_exp_Q_t(phi, Q, dt(i), eq_prec);
    if(accu(Li)>0) v = v % ((1-p)*Li) + (p/ns)*v;
    u = accu(v);
    log_lik_v(i) = log(u);
    phi = v/u;
  } // end i
  
  double n2ll = -2*accu(log_lik_v);
  return n2ll;
  
}