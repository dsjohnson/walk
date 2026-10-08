// [[Rcpp::depends(RcppArmadillo)]]

#include "walk_types.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat ctmc_predict_arma(
    const arma::sp_mat& L, 
    const arma::vec& dt, 
    int ns, 
    const arma::umat& from_to, 
    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
    double p,
    const arma::rowvec& delta, 
    const arma::vec& hij,
    double eq_prec,
    double trunc_tol,
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
  arma::sp_mat Q;
  if(form == 1){
    Q = load_Q(from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, l_r, u_r, link_m, a_m, norm, clip);
  } else {
    Q = load_Q_sde(from_to, Xb_q_r, Xb_q_m, hij, ns, k, clip);
  }
  
  arma::mat A(N, ns, fill::zeros);
  arma::mat G(N, ns, fill::zeros);
  arma::vec scale_factors(N, fill::zeros);
  
  // Forward Pass (Alpha) ----------------------------------------------------
  arma::sp_mat Li = L.col(0).t();
  arma::rowvec v = delta;
  if(accu(Li) > 0) {
    v = v % ((1.0 - p) * Li) + (p / ns) * v;
  }
  
  double u = accu(v);
  scale_factors(0) = u;
  A.row(0) = v / u;
  
  for(int i = 1; i < N; i++) {
    Li = L.col(i).t();
    // Propagate forward through continuous time transition matrix
    v = v_exp_Q_t(A.row(i - 1), Q, dt(i), eq_prec);
    
    // Apply emission / observation model
    if(accu(Li) > 0) {
      v = v % ((1.0 - p) * Li) + (p / ns) * v;
    }
    
    u = accu(v);
    scale_factors(i) = u;
    A.row(i) = v / u;
  }
  
  // Backward Pass (Beta) ---------------------------------------------------
  // Initialize backward vector at terminal time N-1
  arma::rowvec beta_curr(ns, fill::ones);
  arma::rowvec  b_emission(ns, fill::zeros);
  
  // Compute smoothed probability for the last state
  arma::rowvec ab = A.row(N - 1) % beta_curr;
  G.row(N - 1) = ab / accu(ab);
  
  for(int i = N - 1; i > 0; i--) {
    Li = L.col(i).t();
    
    // 1. Incorporate emission probability at step i
    b_emission = beta_curr;
    if(accu(Li) > 0) {
      b_emission = b_emission % ((1.0 - p) * Li) + (p / ns) * b_emission;
    }
    
    // 2. Backpropagate vector across continuous interval dt(i)
    // Note: v_exp_Q_t(b, Q, dt) computes b * exp(Q * dt).
    // Transposition relation: exp(Q * dt) * b^T = (b * exp(Q^T * dt))^T
    beta_curr = v_exp_Q_t(b_emission, Q.t(), dt(i), eq_prec) / scale_factors(i);
    
    // 3. Compute smoothed state probabilities Gamma = Forward % Backward
    ab = A.row(i - 1) % beta_curr;
    ab = ab / accu(ab);
    
    // Apply threshold truncation clean-up
    ab.elem(find(ab < trunc_tol)).zeros();
    G.row(i - 1) = ab / accu(ab);
  }
  
  return G;
}