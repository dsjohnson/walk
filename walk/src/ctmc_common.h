#ifndef WALK_CTMC_COMMON_H
#define WALK_CTMC_COMMON_H

// Pieces shared by the likelihood (ctmc_n2ll_arma.cpp, ctmc_n2ll_joint.cpp),
// prediction (ctmc_predict_arma.cpp) and matrix-exponential (v_exp_Q_t.cpp) code.

#include "walk_types.h"

// Rate matrix for either model form (1 = "mult", 2 = "sde").
inline arma::sp_mat build_Q(const arma::umat& from_to, const arma::vec& Xb_q_r,
                            const arma::vec& Xb_q_m, const arma::vec& hij, int ns,
                            int link_r, double a_r, double l_r, double u_r,
                            int link_m, double a_m, int form, double k, bool norm,
                            double clip) {
  if (form == 1)
    return load_Q(from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, l_r, u_r, link_m, a_m, norm, clip);
  return load_Q_sde(from_to, Xb_q_r, Xb_q_m, hij, ns, k, clip);
}

// Observation update with the outlier mixture:  v <- v % ((1-p) L_i) + (p/ns) v.
// An observation with no information (all-zero L_i) leaves v unchanged.
inline void emission_update(arma::rowvec& v, const arma::sp_mat& Li, double p, int ns) {
  if (arma::accu(Li) > 0) v = v % ((1.0 - p) * Li) + (p / ns) * v;
}

// Poisson truncation window [m_lower, m] for uniformization rate alpha at
// tolerance tol: the terms n < m_lower carry < tol of the mass on the left,
// and the terms n <= m carry >= 1 - tol in total.
inline void poisson_window(double alpha, double tol, int& m_lower, int& m) {
  m = 0;
  double rs = 0.0;
  while (rs < (1.0 - tol)) {
    rs += std::exp(m * std::log(alpha) - alpha - std::lgamma(m + 1));
    m++;
    if (m > 1000000) break;        // safety guard
  }
  m_lower = 0;
  double ls = 0.0;
  for (int n = 0; n <= m; ++n) {
    ls += std::exp(n * std::log(alpha) - alpha - std::lgamma(n + 1));
    if (ls < tol) m_lower = n; else break;
  }
}

// Upper truncation point for a RELATIVE accuracy requirement: the smallest m >= floor(alpha)
// whose Poisson upper tail sum_{n>m} p_n is below exp(log_tol), using the bound
// sum_{n>m} p_n <= p_{m+1} / (1 - alpha/(m+2)) (valid for m + 2 > alpha). Works for any log_tol
// (e.g. log(1e-250)), unlike poisson_window(), whose 1 - tol test needs tol >~ 1e-15.
inline int poisson_m_upper(double alpha, double log_tol) {
  if (alpha <= 0.0) return 0;
  const double la = std::log(alpha);
  int m = (int)std::floor(alpha);
  for (;; ++m) {
    double lp = (m + 1) * la - alpha - std::lgamma(m + 2.0);
    double lb = lp - std::log1p(-alpha / (m + 2.0));
    if (lb < log_tol || m > 10000000) return m;
  }
}

// The truncated series misses at most 2 * tol of probability mass, which can change an
// observation's predictive probability u by at most 2 * tol * max_j L_ij. If that exceeds
// kTruncRelTol * u (an observation far in the tail of the prediction, e.g. a fast displacement,
// for which the omitted terms are the ones that reach it), the step is recomputed with all terms
// from n = 0 and a tail bound of kTruncRelTol * u / (20 max_j L_ij), so a far observation never
// gets probability exactly 0 (or a badly underestimated one) just because the series was cut
// short. Ordinary observations (u >= 2 tol max L / kTruncRelTol, i.e. 2e-4 max L with the default
// eq_prec = 1e-8) are not affected; at the switch, -2 log L can change by at most 2 kTruncRelTol.
const double kTruncRelTol = 1e-4;
const double kTruncMinLogTol = -690.0;   // ~ log(1e-300): give up beyond this (u underflows anyway)

// Validate 0-based segment starts: seg_start[0] == 0, strictly increasing, < N.
inline void check_seg_start(const arma::uvec& seg_start, arma::uword N) {
  if (seg_start.n_elem == 0 || seg_start(0) != 0) Rcpp::stop("seg_start must start at 0");
  for (arma::uword s = 1; s < seg_start.n_elem; ++s)
    if (seg_start(s) <= seg_start(s - 1) || seg_start(s) >= N)
      Rcpp::stop("seg_start must be strictly increasing and < N");
}

#endif // WALK_CTMC_COMMON_H
