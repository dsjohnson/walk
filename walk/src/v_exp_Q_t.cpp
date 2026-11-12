
// [[Rcpp::depends(RcppArmadillo)]]

#include "ctmc_common.h"

// ---------------------------------------------------------------------------
// Build the uniformization matrix P = I + Q / alpha_0 once, and return alpha_0.
//
// P depends only on Q, NOT on t, so it is built once per Q and reused for every
// time step (rebuilding it per step would allocate an ns x ns sparse matrix N
// times per likelihood evaluation).
// ---------------------------------------------------------------------------
double uniformization_P(const arma::sp_mat& Q, arma::sp_mat& P) {
  arma::spdiagview<double> diag_Q = Q.diag();
  double alpha_0 = 0.0;
  for (arma::uword i = 0; i < diag_Q.n_elem; ++i) {
    double abs_val = std::abs(diag_Q(i));
    if (abs_val > alpha_0) alpha_0 = abs_val;
  }
  if (alpha_0 <= 0.0) {           // no movement at all: exp(Q t) = I
    P = arma::sp_mat(Q.n_rows, Q.n_cols);
    return 0.0;
  }
  arma::sp_mat I;
  I.eye(Q.n_rows, Q.n_cols);
  P = I + (Q / alpha_0);
  return alpha_0;
}

// ---------------------------------------------------------------------------
// v * exp(Q t) given a pre-built P and alpha_0.
// ---------------------------------------------------------------------------
arma::mat v_exp_Q_t_pre(const arma::rowvec& v, const arma::sp_mat& P,
                        double alpha_0, double t, double tolerance, double log_tol_ext) {
  int N = P.n_rows;

  // Handle the trivial edge case instantly
  if (t <= 0.0 || alpha_0 <= 0.0) {
    return arma::conv_to<arma::mat>::from(v);
  }

  // 1. Calculate the time-dependent uniformization rate
  double alpha = alpha_0 * t;

  // 2.-3. Poisson truncation window [m_lower, m] (shared with ctmc_n2ll_joint.cpp)
  int m_lower, m;
  // log_tol_ext < 0: all terms from n = 0, upper tail below exp(log_tol_ext) (see kTruncRelTol)
  if (log_tol_ext < 0.0) { m_lower = 0; m = poisson_m_upper(alpha, log_tol_ext); }
  else poisson_window(alpha, tolerance, m_lower, m);

  // 4. Forward accumulate the Poisson series
  arma::rowvec v_curr = v;
  arma::rowvec result(N, arma::fill::zeros);

  for (int n = 0; n <= m; ++n) {
    if (n >= m_lower) {
      double log_p_n = n * std::log(alpha) - alpha - std::lgamma(n + 1);
      result += std::exp(log_p_n) * v_curr;
    }
    if (n < m) {
      // NB: `v_curr * P` (dense rowvec x sparse) is ~14x slower than
      // `P.t() * v_curr.t()` in Armadillo, for bit-identical results.
      v_curr = (P.t() * v_curr.t()).t();
    }
  }

  return arma::conv_to<arma::mat>::from(result);
}
