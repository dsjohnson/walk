// [[Rcpp::depends(RcppArmadillo)]]
#include "ctmc_common.h"
#include <vector>
#include <deque>
#ifdef _OPENMP
#include <omp.h>
#endif

// ===========================================================================
// Joint CTMC -2 log-likelihood AND its exact gradient in one forward pass.
//
// Same model and the same forward recursion as ctmc_n2ll_joint_arma() (support
// tracking, uniformization, mass dropping), but every filtered vector carries
// its derivatives with respect to the parameters ("forward-mode" / tangent
// propagation).  For uniformization with a FIXED rate alpha,
//      exp(Q t) = sum_n Pois(n; alpha t) P^n,   P = I + Q / alpha,
// holds for every Q, so d/dtheta exp(Q t) is obtained by differentiating the
// mat-vecs only:  d(v P) = dv P + v dP,  dP = dQ / alpha.  The index work of a
// mat-vec (the expensive part) is shared by the value and all derivatives, so
// one call costs much less than the 1 + n_par evaluations of a forward-
// difference gradient.
//
// Supported: form = "mult", links soft_plus / log / logit (q_r), soft_plus /
// log (q_m), norm TRUE / FALSE, outlier mixture p, fixed delta, prune, drop.
// Not supported (the R side falls back to numerical derivatives): form "sde",
// clip > 0, delta = "stationary".
// ===========================================================================

namespace {

// must equal kDropRelTol in ctmc_n2ll_joint.cpp so value and gradient describe the same function
const double kDropRelTolG = 1e-6;
  // carried-error check, as in segment_n2ll() (ctmc_n2ll_joint.cpp): kLookback / kCarryRelTol
  const int kLookbackG = 4;
  const double kCarryRelTolG = 1e-4;
  struct CkptG {
    arma::uword step;
    std::vector<arma::uword> idx;
    std::vector<double> val;          // value + derivatives, stride S
    double ll;
    std::vector<double> g;            // gradient accumulated before that step
    double e, amp;
  };
  
  void chk_int_fn_g(void*) { R_CheckUserInterrupt(); }
  inline bool user_interrupt_g() { return R_ToplevelExec(chk_int_fn_g, NULL) == FALSE; }
  
  struct CtlG {
    volatile int* abort;
    bool main;
    inline bool stop() { if (main && !*abort && user_interrupt_g()) *abort = 1; return *abort != 0; }
  };
  
  // P = I + Q/alpha in CSR form (row i = cells reachable from i, diagonal first),
  // with dP/dtheta_k (k < npq) stored row-major per entry.
  struct PMat {
    std::vector<arma::uword> rp, ci;
    std::vector<double> v, dv;   // dv[e * npq + k]
    int npq = 0;
  };
  
  // stride S = 1 + np: component 0 = value, 1..np = derivatives
  struct WorkG {
    int S;
    std::vector<double> vc, vn, res;
    std::vector<unsigned int> mark_n, mark_r;
    std::vector<arma::uword> act, nact, ract;
    std::vector<arma::uword> sav_i;   // copy of phi (value + derivatives) before a dropping propagate
    std::vector<double> sav_v;
    unsigned int stamp_n = 0, stamp_r = 0;
    WorkG(arma::uword ns, int S_) : S(S_), vc(ns * S_, 0.0), vn(ns * S_, 0.0), res(ns * S_, 0.0),
    mark_n(ns, 0), mark_r(ns, 0) { act.reserve(1024); nact.reserve(1024); ract.reserve(1024); }
  };
  
  // a cell matters if its value OR any derivative is non-zero (a cell can have value 0 and a
  // non-zero derivative, e.g. parity zeros of P^n when P has zero diagonal entries)
  inline bool any_nz(const double* x, int S) { for (int c = 0; c < S; ++c) if (x[c] != 0.0) return true; return false; }
  // drop only cells whose value AND derivatives are all below thr
  // thr_d[k] = thr * (L1 norm of derivative k of the incoming vector): at most 'drop' of the
  // mass and 'drop' of each derivative's total |.| is removed per call (scale invariant)
  inline bool negligible(const double* x, int S, double thr, const double* thr_d) {
    if (x[0] >= thr) return false;
    for (int c = 1; c < S; ++c) if (std::abs(x[c]) >= thr_d[c - 1]) return false;
    return true;
  }
  inline void zero_cell(std::vector<double>& x, arma::uword j, int S) {
    double* p = &x[j * S]; for (int c = 0; c < S; ++c) p[c] = 0.0;
  }
  
  template <int NP, int NPQ>
  double propagate_t(WorkG& ws, const PMat& P, int np_rt, double alpha_0, double t,
                     double tol, double drop, CtlG& ctl, double log_tol_ext) {
    const int np = (NP >= 0) ? NP : np_rt;
    const int npq = (NPQ >= 0) ? NPQ : P.npq;
    const int S = np + 1;
    double dropped = 0.0;
    ws.ract.clear();
    ++ws.stamp_r;
    if (t <= 0.0 || alpha_0 <= 0.0) {
      for (arma::uword i : ws.act) {
        for (int c = 0; c < S; ++c) { ws.res[i * S + c] = ws.vc[i * S + c]; ws.vc[i * S + c] = 0.0; }
        ws.mark_r[i] = ws.stamp_r; ws.ract.push_back(i);
      }
      ws.act.clear();
      return 0.0;
    }
    double alpha = alpha_0 * t;
    int m_lower, m;
    if (log_tol_ext < 0.0) { m_lower = 0; m = poisson_m_upper(alpha, log_tol_ext); }   // see kTruncRelTol
    else poisson_window(alpha, tol, m_lower, m);
    double la = std::log(alpha);
    std::vector<double> Dk(np, 0.0), thr_d(np, 0.0);
    if (drop > 0.0) for (arma::uword i : ws.act) for (int k = 0; k < np; ++k) Dk[k] += std::abs(ws.vc[i * S + 1 + k]);
    for (int n = 0; n <= m; ++n) {
      if (n >= m_lower) {
        double pn = std::exp(n * la - alpha - std::lgamma(n + 1));
        for (arma::uword i : ws.act) {
          double* r = &ws.res[i * S];
          if (ws.mark_r[i] != ws.stamp_r) {
            ws.mark_r[i] = ws.stamp_r; ws.ract.push_back(i);
            for (int c = 0; c < S; ++c) r[c] = 0.0;
          }
          const double* x = &ws.vc[i * S];
          for (int c = 0; c < S; ++c) r[c] += pn * x[c];
        }
      }
      if (n < m) {
        ++ws.stamp_n;
        ws.nact.clear();
        for (arma::uword i : ws.act) {
          const double* x = &ws.vc[i * S];
          const double x0 = x[0];
          for (arma::uword e = P.rp[i]; e < P.rp[i + 1]; ++e) {
            arma::uword j = P.ci[e];
            double* y = &ws.vn[j * S];
            if (ws.mark_n[j] != ws.stamp_n) {
              ws.mark_n[j] = ws.stamp_n; ws.nact.push_back(j);
              for (int c = 0; c < S; ++c) y[c] = 0.0;
            }
            const double pv = P.v[e];
            const double* dp = npq ? &P.dv[e * npq] : nullptr;
            y[0] += x0 * pv;
            for (int k = 0; k < npq; ++k) y[1 + k] += x[1 + k] * pv + x0 * dp[k];
            for (int k = npq; k < np; ++k) y[1 + k] += x[1 + k] * pv;
          }
        }
        for (arma::uword i : ws.act) zero_cell(ws.vc, i, S);
        std::swap(ws.vc, ws.vn);
        std::swap(ws.act, ws.nact);
        if (drop > 0.0 && !ws.act.empty()) {
          const double thr = drop / ((double)std::max(m, 1) * (double)ws.act.size());
          for (int k = 0; k < np; ++k) thr_d[k] = thr * Dk[k];
          std::size_t w = 0;
          for (std::size_t q = 0; q < ws.act.size(); ++q) {
            arma::uword j = ws.act[q];
            if (negligible(&ws.vc[j * S], S, thr, thr_d.data())) { dropped += ws.vc[j * S]; zero_cell(ws.vc, j, S); } else ws.act[w++] = j;
          }
          ws.act.resize(w);
        }
        if ((n & 63) == 63 && ctl.stop()) break;
      }
    }
    for (arma::uword i : ws.act) zero_cell(ws.vc, i, S);
    ws.act.clear();
    return dropped;
  }
  
  // dispatch to a version compiled for the parameter count (common sizes), else generic
  double propagate_g(WorkG& ws, const PMat& P, int np, double alpha_0, double t,
                     double tol, double drop, CtlG& ctl, double log_tol_ext) {
    const int q = P.npq;
#define WALK_PG(A, B) if (np == A && q == B) return propagate_t<A, B>(ws, P, np, alpha_0, t, tol, drop, ctl, log_tol_ext);
    WALK_PG(1,1) WALK_PG(2,2) WALK_PG(3,3) WALK_PG(4,4) WALK_PG(5,5) WALK_PG(6,6) WALK_PG(7,7) WALK_PG(8,8)
      WALK_PG(2,1) WALK_PG(3,2) WALK_PG(4,3) WALK_PG(5,4) WALK_PG(6,5) WALK_PG(7,6) WALK_PG(8,7) WALK_PG(9,8)
#undef WALK_PG
      return propagate_t<-1, -1>(ws, P, np, alpha_0, t, tol, drop, ctl, log_tol_ext);
  }
  
  // normalise the filtered vector on ws.act (value + derivatives) by its mass u;
  // adds log(u) and dlog(u) to ll / g
  inline void normalise_g(WorkG& ws, int np, double u, const std::vector<double>& du,
                          double& ll, double* g) {
    const int S = ws.S;
    ll += std::log(u);
    for (int k = 0; k < np; ++k) g[k] += du[k] / u;
    for (arma::uword j : ws.act) {
      double* x = &ws.vc[j * S];
      x[0] /= u;
      for (int k = 0; k < np; ++k) x[1 + k] = x[1 + k] / u - x[0] * du[k] / u;
    }
  }
  
  double segment_g(arma::uword i0, arma::uword i1, const arma::uword* Lcp, const arma::uword* Lri,
                   const double* Lv, const arma::vec& dt, const arma::rowvec& delta, double p,
                   double dp_dlogit, int ip, arma::uword ns, const PMat& P, int np, double alpha_0,
                   double tol, double prune, double drop, CtlG& ctl, std::vector<int>& zobs,
                   WorkG& ws, double* g) {
    const int S = ws.S;
    double ll = 0.0;
    bool full = (p > 0.0);
    std::vector<double> du(np, 0.0);
    {
      ws.act.clear();
      double u = 0.0;
      std::fill(du.begin(), du.end(), 0.0);
      arma::uword k0 = Lcp[i0], k1 = Lcp[i0 + 1];
      if (k1 > k0 && !full) {
        for (arma::uword k = k0; k < k1; ++k) {
          arma::uword j = Lri[k];
          double val = delta[j] * Lv[k];
          if (val > 0.0) { zero_cell(ws.vc, j, S); ws.vc[j * S] = val; ws.act.push_back(j); u += val; }  // delta, L do not depend on par
        }
      } else {
        std::vector<double> Ld;
        if (k1 > k0) { Ld.assign(ns, 0.0); for (arma::uword k = k0; k < k1; ++k) Ld[Lri[k]] = Lv[k]; }
        for (arma::uword j = 0; j < ns; ++j) {
          double val = delta[j], dval = 0.0;
          if (k1 > k0) {
            dval = delta[j] * (-Ld[j] + 1.0 / ns) * dp_dlogit;
            val = val * ((1.0 - p) * Ld[j]) + (p / ns) * val;
          }
          if (val > 0.0) {
            zero_cell(ws.vc, j, S); ws.vc[j * S] = val; ws.act.push_back(j); u += val;
            if (ip >= 0) { ws.vc[j * S + 1 + ip] = dval; du[ip] += dval; }
          }
        }
      }
      if (!(u > 0.0)) { for (arma::uword j : ws.act) zero_cell(ws.vc, j, S); ws.act.clear(); return R_PosInf; }
      normalise_g(ws, np, u, du, ll, g);
    }
    std::deque<CkptG> hist;
    arma::uword exact_until = 0;
    for (arma::uword i = i0 + 1; i < i1; ++i) {
      if (*ctl.abort || (((i - i0) & 7) == 0 && ctl.stop())) break;
      const double ll_before = ll;
      std::vector<double> g_before;
      if (drop > 0.0) g_before.assign(g, g + np);
      // same safeguard as segment_n2ll(): if the dropped mass could change this observation's
      // predictive probability u by more than kDropRelTolG * u, redo the step with drop = 0
      double dcur = (i <= exact_until) ? 0.0 : drop;
      double ltol = 0.0;                  // < 0: extended series (kTruncRelTol), as in segment_n2ll()
      double last_dropped = 0.0, last_wmax = 1.0;
      ws.sav_i.assign(ws.act.begin(), ws.act.end()); ws.sav_v.clear();
      for (arma::uword j : ws.act) ws.sav_v.insert(ws.sav_v.end(), &ws.vc[j * S], &ws.vc[j * S] + S);
      arma::uword k0 = Lcp[i], k1 = Lcp[i + 1];
      double u = 0.0;
      for (;;) {
        const double dropped = propagate_g(ws, P, np, alpha_0, dt(i), tol, dcur, ctl, ltol);
        double wmax = 1.0;
        last_dropped = dcur > 0.0 ? dropped : 0.0;
        u = 0.0;
        std::fill(du.begin(), du.end(), 0.0);
        ws.act.clear();
        if (k1 == k0) {
          for (arma::uword j : ws.ract) {
            const double* r = &ws.res[j * S];
            if (any_nz(r, S)) {
              double* x = &ws.vc[j * S];
              for (int c = 0; c < S; ++c) x[c] = r[c];
              ws.act.push_back(j); u += r[0];
              for (int k = 0; k < np; ++k) du[k] += r[1 + k];
            }
          }
        } else if (!full) {
          for (arma::uword k = k0; k < k1; ++k) {
            arma::uword j = Lri[k];
            if (ws.mark_r[j] != ws.stamp_r) continue;
            const double* r = &ws.res[j * S];
            if (Lv[k] > 0.0 && any_nz(r, S)) {
              double* x = &ws.vc[j * S];
              x[0] = r[0] * Lv[k]; u += x[0];
              for (int c = 1; c < S; ++c) { x[c] = r[c] * Lv[k]; du[c - 1] += x[c]; }
              ws.act.push_back(j);
            }
          }
        } else {
          // outlier mixture: v = res % ((1-p) L + p/ns); d/dlogit_p adds res % (-L + 1/ns) * dp
          for (arma::uword j : ws.ract) ws.vn[j * S] = 0.0;           // vn[.,0] holds L_j temporarily
          for (arma::uword k = k0; k < k1; ++k) {
            arma::uword j = Lri[k];
            if (ws.mark_r[j] == ws.stamp_r) ws.vn[j * S] = Lv[k];
          }
          for (arma::uword j : ws.ract) {
            const double Lj = ws.vn[j * S]; ws.vn[j * S] = 0.0;
            const double* r = &ws.res[j * S];
            const double w = (1.0 - p) * Lj + p / ns;
            if (w > 0.0 && any_nz(r, S)) {
              double* x = &ws.vc[j * S];
              x[0] = r[0] * w; u += x[0];
              for (int c = 1; c < S; ++c) x[c] = r[c] * w;
              if (ip >= 0) x[1 + ip] += r[0] * (-Lj + 1.0 / ns) * dp_dlogit;
              for (int k = 0; k < np; ++k) du[k] += x[1 + k];
              ws.act.push_back(j);
            }
          }
        }
        if (k1 > k0) { double lm = 0.0; for (arma::uword k = k0; k < k1; ++k) lm = std::max(lm, Lv[k]);
        wmax = full ? (1.0 - p) * lm + p / ns : lm; }
        last_wmax = wmax;
        // same decisions as segment_n2ll() (u is the same number there), so value and gradient match
        const bool redo_drop = dcur > 0.0 && dropped > 0.0 && !(dropped * wmax <= kDropRelTolG * u);
        const double cur_tol = (ltol < 0.0) ? std::exp(ltol) : tol;
        const bool redo_trunc = dt(i) > 0.0 && alpha_0 > 0.0 && ltol > kTruncMinLogTol &&
          !(2.0 * cur_tol * wmax <= kTruncRelTol * u);
          if ((redo_drop || redo_trunc) && !*ctl.abort) {
            for (arma::uword j : ws.act) zero_cell(ws.vc, j, S);
            for (arma::uword j : ws.ract) zero_cell(ws.res, j, S);
            ws.act.assign(ws.sav_i.begin(), ws.sav_i.end());
            for (std::size_t q = 0; q < ws.sav_i.size(); ++q)
              for (int c = 0; c < S; ++c) ws.vc[ws.sav_i[q] * S + c] = ws.sav_v[q * S + c];
            if (redo_drop) dcur = 0.0;
            if (redo_trunc) {
              double target = (u > 0.0) ? std::log(kTruncRelTol * u / (20.0 * wmax)) : std::log(cur_tol) - 46.0;
              ltol = std::max(kTruncMinLogTol, std::min(target, std::log(cur_tol) - 2.3));
            }
            continue;
          }
          break;
      }
      if (drop > 0.0 && !hist.empty() && !*ctl.abort) {
        double carry = 0.0, f = 1.0; int jstar = -1;
        for (int h = (int)hist.size() - 1; h >= 0; --h) {
          if (hist[h].e > 0.0) { carry += hist[h].e * f; jstar = h; }
          f *= hist[h].amp;
        }
        if (jstar >= 0 && !(carry * last_wmax <= kCarryRelTolG * u)) {
          for (arma::uword j : ws.act) zero_cell(ws.vc, j, S);
          for (arma::uword j : ws.ract) zero_cell(ws.res, j, S);
          const CkptG& c = hist[jstar];
          ws.act.assign(c.idx.begin(), c.idx.end());
          for (std::size_t q = 0; q < c.idx.size(); ++q)
            for (int cc = 0; cc < S; ++cc) ws.vc[c.idx[q] * S + cc] = c.val[q * S + cc];
          ll = c.ll; for (int k = 0; k < np; ++k) g[k] = c.g[k];
          exact_until = i;
          i = c.step - 1;
          hist.erase(hist.begin() + jstar, hist.end());
          continue;
        }
      }
      if (drop > 0.0) {
        CkptG c; c.step = i; c.idx.swap(ws.sav_i); c.val.swap(ws.sav_v); c.ll = ll_before;
        c.g = g_before;
        c.e = (u > 0.0) ? last_dropped * last_wmax / u : 0.0;
        c.amp = (u > 0.0) ? 2.0 * last_wmax / u : 0.0;
        hist.push_back(std::move(c));
        if ((int)hist.size() > kLookbackG) hist.pop_front();
      }
      if (!(u > 0.0)) {
        for (arma::uword j : ws.act) zero_cell(ws.vc, j, S);
        ws.act.clear();
        zobs.push_back((int)i + 1);
        for (arma::uword j : ws.ract) zero_cell(ws.res, j, S);
        return R_PosInf;
      }
      for (arma::uword j : ws.ract) zero_cell(ws.res, j, S);
      normalise_g(ws, np, u, du, ll, g);
      if (prune > 0.0) {
        double mx = 0.0;
        for (arma::uword j : ws.act) if (ws.vc[j * S] > mx) mx = ws.vc[j * S];
        double thr = prune * mx, kept = 0.0;
        std::vector<double> dk(np, 0.0);
        ws.nact.clear();
        for (arma::uword j : ws.act) {
          if (ws.vc[j * S] >= thr) {
            ws.nact.push_back(j); kept += ws.vc[j * S];
            for (int k = 0; k < np; ++k) dk[k] += ws.vc[j * S + 1 + k];
          } else zero_cell(ws.vc, j, S);
        }
        std::swap(ws.act, ws.nact);
        // renormalising by 'kept' is the same operation (it does not add to log L)
        double dummy = 0.0; std::vector<double> gz(np, 0.0);
        normalise_g(ws, np, kept, dk, dummy, gz.data());
      }
    }
    for (arma::uword j : ws.act) zero_cell(ws.vc, j, S);
    ws.act.clear();
    for (int k = 0; k < np; ++k) g[k] *= -2.0;
    return -2.0 * ll;
  }
  
  // link value and derivative
  inline void link_fn(int link, double x, double a, double L, double U, double& f, double& df) {
    if (link == 1) {                       // soft_plus
      f = std::max(0.0, x) + std::log1p(std::exp(-std::abs(a * x))) / a;
      df = 1.0 / (1.0 + std::exp(-a * x));
    } else if (link == 2) {                // log (trunc_exp as in load_Q)
      f = arma::trunc_exp(x); df = f;
    } else {                               // logit
      double s = 1.0 / (1.0 + arma::trunc_exp(-x));
      f = L + (U - L) * s; df = (U - L) * s * (1.0 - s);
    }
  }
  
}  // namespace

//' Joint CTMC -2 log-likelihood and its gradient (internal)
 //' @return list(n2ll = per-segment -2 log L, grad = n_par x n_seg matrix of derivatives,
 //'   zero_obs, alpha_0)
 //' @keywords internal
 // [[Rcpp::export]]
 Rcpp::List ctmc_n2ll_grad_arma(
     const arma::sp_mat& L, const arma::vec& dt, const arma::uvec& seg_start, int ns,
     const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
     const arma::mat& X_q_r, const arma::mat& X_q_m, double p, bool has_p,
     const arma::rowvec& delta, double eq_prec,
     int link_r, double a_r, double l_r, double u_r, int link_m, double a_m, bool norm,
     int nthreads = 1, double prune = 0.0, double drop = 0.0)
 {
   arma::uword N = dt.n_elem;
   if (L.n_cols != N) Rcpp::stop("ncol(L) must equal length(dt)");
   if ((int)L.n_rows != ns) Rcpp::stop("nrow(L) must equal ns");
   check_seg_start(seg_start, N);
   const arma::uword E = from_to.n_cols;
   const int kr = X_q_r.n_cols, km = X_q_m.n_cols;
   if ((int)X_q_r.n_rows != ns || Xb_q_r.n_elem != (arma::uword)ns) Rcpp::stop("X_q_r / Xb_q_r do not match ns");
   if (Xb_q_m.n_elem != E || (km > 0 && X_q_m.n_rows != E)) Rcpp::stop("X_q_m / Xb_q_m do not match from_to");
   const int npq = kr + km, np = npq + (has_p ? 1 : 0), S = 1 + np;
   const int ip = has_p ? npq : -1;
   
   // ---- rates and their derivatives --------------------------------------
   std::vector<double> r(ns), dr(ns), mv(E), dm(E);
   for (int i = 0; i < ns; ++i) link_fn(link_r, Xb_q_r[i], a_r, l_r, u_r, r[i], dr[i]);
   for (arma::uword e = 0; e < E; ++e) link_fn(link_m, Xb_q_m[e], a_m, 0.0, 1.0, mv[e], dm[e]); // Updated to pass link_m directly
   // rows of P: count out-edges (from = row 0 of from_to)
   std::vector<arma::uword> cnt(ns, 0);
   for (arma::uword e = 0; e < E; ++e) {
     if ((int)from_to(0, e) >= ns || (int)from_to(1, e) >= ns) Rcpp::stop("from_to out of range");
     if (from_to(0, e) == from_to(1, e)) Rcpp::stop("from_to contains self transitions");
     cnt[from_to(0, e)]++;
   }
   PMat P; P.npq = npq;
   P.rp.assign(ns + 1, 0);
   for (int i = 0; i < ns; ++i) P.rp[i + 1] = P.rp[i] + cnt[i] + 1;   // + diagonal
   const arma::uword nnz = P.rp[ns];
   P.ci.assign(nnz, 0); P.v.assign(nnz, 0.0); P.dv.assign(nnz * (std::size_t)npq, 0.0);
   std::vector<arma::uword> pos(ns);
   for (int i = 0; i < ns; ++i) { P.ci[P.rp[i]] = i; pos[i] = P.rp[i] + 1; }
   // row sums of m (norm = TRUE) and of dm * X_q_m
   std::vector<double> Sm(ns, 0.0), Sdm(ns * (std::size_t)std::max(km, 1), 0.0);
   if (norm) for (arma::uword e = 0; e < E; ++e) {
     arma::uword i = from_to(0, e); Sm[i] += std::abs(mv[e]);
     for (int k = 0; k < km; ++k) Sdm[i * km + k] += dm[e] * X_q_m(e, k);
   }
   // off-diagonal entries: Q_ij and dQ_ij (rate per unit time)
   for (arma::uword e = 0; e < E; ++e) {
     arma::uword i = from_to(0, e), j = from_to(1, e), q = pos[i]++;
     P.ci[q] = j;
     double mt = norm ? (Sm[i] > 0 ? mv[e] / Sm[i] : 0.0) : mv[e];
     double Qij = r[i] * mt;
     P.v[q] = Qij;
     double* d = npq ? &P.dv[q * npq] : nullptr;
     for (int k = 0; k < kr; ++k) d[k] = dr[i] * X_q_r(i, k) * mt;
     for (int k = 0; k < km; ++k) {
       double dmt = norm ? (Sm[i] > 0 ? (dm[e] * X_q_m(e, k) * Sm[i] - mv[e] * Sdm[i * km + k]) / (Sm[i] * Sm[i]) : 0.0)
         : dm[e] * X_q_m(e, k);
       d[kr + k] = r[i] * dmt;
     }
   }
   // diagonal = - row sums; uniformization rate; scale to P
   double alpha_0 = 0.0;
   for (int i = 0; i < ns; ++i) {
     double s = 0.0; arma::uword dq = P.rp[i];
     for (arma::uword q = dq + 1; q < P.rp[i + 1]; ++q) {
       s += P.v[q];
       for (int k = 0; k < npq; ++k) P.dv[dq * npq + k] -= P.dv[q * npq + k];
     }
     if (norm) {
       // load_Q(norm = TRUE): diag(Qm) = -1 for every row, so Q_ii = -r_i (also for rows without
       // edges, where the off-diagonal sum is 0)
       s = r[i];
       for (int k = 0; k < kr; ++k) P.dv[dq * npq + k] = -dr[i] * X_q_r(i, k);
       for (int k = kr; k < npq; ++k) P.dv[dq * npq + k] = 0.0;
     }
     P.v[dq] = -s;
     if (std::abs(s) > alpha_0) alpha_0 = std::abs(s);
   }
   if (alpha_0 > 0.0) {
     for (arma::uword q = 0; q < nnz; ++q) P.v[q] /= alpha_0;
     for (int i = 0; i < ns; ++i) P.v[P.rp[i]] += 1.0;
     for (std::size_t q = 0; q < P.dv.size(); ++q) P.dv[q] /= alpha_0;
   }
   const double dp_dlogit = p * (1.0 - p);
   
   arma::sp_mat Lc(L); Lc.sync();
   const arma::uword* Lcp = Lc.col_ptrs; const arma::uword* Lri = Lc.row_indices; const double* Lv = Lc.values;
   
   arma::uword Sg = seg_start.n_elem;
   arma::vec out(Sg, arma::fill::zeros);
   arma::mat G(np, Sg, arma::fill::zeros);
   std::vector< std::vector<int> > zo(Sg);
   int nt = 1;
#ifdef _OPENMP
   nt = std::max(1, nthreads);
   nt = std::min(nt, omp_get_thread_limit());
   nt = std::min(nt, (int)std::max<arma::uword>(1, Sg));
#else
   (void)nthreads;
#endif
   std::vector<WorkG> wss; wss.reserve(nt);
   for (int t = 0; t < nt; ++t) wss.emplace_back(ns, S);
   for (arma::uword s = 0; s < Sg; ++s) zo[s].reserve(4);
   int failed = 0;
   volatile int abort = 0;
   long long next_seg = 0;
   const long long S_ll = (long long)Sg;
   double* Gp = G.memptr();
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
{
  int tid = 0;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#endif
  CtlG ctl{&abort, tid == 0};
  for (;;) {
    long long sl;
#ifdef _OPENMP
#pragma omp atomic capture
#endif
    sl = next_seg++;
    if (sl >= S_ll) break;
    arma::uword s = (arma::uword)sl;
    arma::uword i0 = seg_start(s), i1 = (s + 1 < Sg) ? seg_start(s + 1) : N;
    if (abort) break;
    try {
      out(s) = segment_g(i0, i1, Lcp, Lri, Lv, dt, delta, p, dp_dlogit, ip, ns, P, np, alpha_0,
          eq_prec, prune, drop, ctl, zo[s], wss[tid], Gp + s * np);
      if (!std::isfinite(out(s))) for (int k = 0; k < np; ++k) Gp[s * np + k] = 0.0;
    } catch (...) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
      failed = 1;
      out(s) = R_NaN;
    }
  }
}
if (abort) throw Rcpp::internal::InterruptedException();
if (failed) Rcpp::stop("ctmc_n2ll_grad_arma: out of memory (or other C++ error) in a worker thread.");
std::vector<int> zall;
for (arma::uword s = 0; s < Sg; ++s) zall.insert(zall.end(), zo[s].begin(), zo[s].end());
return Rcpp::List::create(Rcpp::Named("n2ll") = out, Rcpp::Named("grad") = G,
                          Rcpp::Named("zero_obs") = Rcpp::wrap(zall),
                          Rcpp::Named("alpha_0") = alpha_0);
 }