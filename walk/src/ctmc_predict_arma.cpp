// [[Rcpp::depends(RcppArmadillo)]]

#include "ctmc_common.h"

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
  arma::sp_mat Q = build_Q(from_to, Xb_q_r, Xb_q_m, hij, ns, link_r, a_r, l_r, u_r,
                           link_m, a_m, form, k, norm, clip);
  
  // P = I + Q/alpha_0 for forward (Q) and backward (Q^T) passes, built once.
  arma::sp_mat Pf, Pb;
  double alpha_0 = uniformization_P(Q, Pf);
  arma::sp_mat Qt = Q.t();
  uniformization_P(Qt, Pb);

  arma::mat A(N, ns, fill::zeros);
  arma::mat G(N, ns, fill::zeros);
  arma::vec scale_factors(N, fill::zeros);
  
  // Forward Pass (Alpha) ----------------------------------------------------
  arma::rowvec v = delta;
  emission_update(v, L.col(0).t(), p, ns);
  
  double u = accu(v);
  scale_factors(0) = u;
  A.row(0) = v / u;
  
  // per-step series length: as in the likelihood (ctmc_n2ll_joint.cpp), a step whose observation
  // is in the far tail of the prediction is recomputed with more terms (kTruncRelTol); the
  // backward pass uses the same number of terms for that step
  arma::vec ltol_step(N, fill::zeros);
  for(int i = 1; i < N; i++) {
    arma::sp_mat Li = L.col(i).t();
    double lmax = Li.n_nonzero ? arma::max(arma::nonzeros(Li)) : 0.0;
    double wmax = (accu(Li) > 0) ? ((p > 0) ? (1.0 - p) * lmax + p / ns : lmax) : 1.0;
    double ltol = 0.0;
    for(;;) {
      // Propagate forward through continuous time transition matrix
      v = v_exp_Q_t_pre(A.row(i - 1), Pf, alpha_0, dt(i), eq_prec, ltol);
      // Apply emission / observation model
      emission_update(v, Li, p, ns);
      u = accu(v);
      double cur_tol = (ltol < 0.0) ? std::exp(ltol) : eq_prec;
      if (dt(i) > 0.0 && alpha_0 > 0.0 && ltol > kTruncMinLogTol && !(2.0 * cur_tol * wmax <= kTruncRelTol * u)) {
        double target = (u > 0.0) ? std::log(kTruncRelTol * u / (20.0 * wmax)) : std::log(cur_tol) - 46.0;
        ltol = std::max(kTruncMinLogTol, std::min(target, std::log(cur_tol) - 2.3));
        continue;
      }
      break;
    }
    ltol_step(i) = ltol;
    
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
    // 1. Incorporate emission probability at step i
    b_emission = beta_curr;
    emission_update(b_emission, L.col(i).t(), p, ns);
    
    // 2. Backpropagate vector across continuous interval dt(i)
    // Note: v_exp_Q_t_pre(b, P, alpha_0, dt) computes b * exp(Q * dt).
    // Transposition relation: exp(Q * dt) * b^T = (b * exp(Q^T * dt))^T
    beta_curr = v_exp_Q_t_pre(b_emission, Pb, alpha_0, dt(i), eq_prec, ltol_step(i)) / scale_factors(i);
    
    // 3. Compute smoothed state probabilities Gamma = Forward % Backward
    ab = A.row(i - 1) % beta_curr;
    ab = ab / accu(ab);
    
    // Apply threshold truncation clean-up
    ab.elem(find(ab < trunc_tol)).zeros();
    G.row(i - 1) = ab / accu(ab);
  }
  
  return G;
}
