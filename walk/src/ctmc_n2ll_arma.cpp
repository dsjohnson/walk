// [[Rcpp::depends(RcppArmadillo)]]

#include "ctmc_common.h"

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
  // Reference single-chain forward algorithm (dense vectors). fit_ctmc() uses
  // ctmc_n2ll_joint_arma(); this one is kept as an independent check of it.
  // All observations are treated as ONE chain (segments are ignored).
  int N = dt.size();
  arma::vec log_lik_v(N);
  arma::sp_mat Q = build_Q(from_to, Xb_q_r, Xb_q_m, hij, ns, link_r, a_r, l_r, u_r,
                           link_m, a_m, form, k, norm, clip);
  arma::sp_mat P;
  double alpha_0 = uniformization_P(Q, P);   // built once per Q

  arma::rowvec v = delta;
  emission_update(v, L.col(0).t(), p, ns);
  double u = accu(v);
  log_lik_v(0) = log(u);
  arma::rowvec phi = v/u;
  for(int i=1; i<N; i++){
    v = v_exp_Q_t_pre(phi, P, alpha_0, dt(i), eq_prec);
    emission_update(v, L.col(i).t(), p, ns);
    u = accu(v);
    log_lik_v(i) = log(u);
    phi = v/u;
  }
  return -2*accu(log_lik_v);
}
