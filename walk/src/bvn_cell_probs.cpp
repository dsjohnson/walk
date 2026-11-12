// Fast bivariate-normal raster-cell likelihoods for walk::proc_Langevin (and proc_telem)
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <climits>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Standard normal CDF (thread-safe, no R API)
// ---------------------------------------------------------------------------
static inline double Phi(double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); }

// ---------------------------------------------------------------------------
// Upper bivariate normal probability P(X > h, Y > k), corr = r.
// Port of Alan Genz's BVNU (TVPACK), the algorithm mvtnorm uses for 2-D
// problems. Absolute accuracy ~1e-15.
// ---------------------------------------------------------------------------
static const double GL_X3[3]  = {-0.9324695142031522, -0.6612093864662647, -0.2386191860831970};
static const double GL_W3[3]  = { 0.1713244923791705,  0.3607615730481384,  0.4679139345726904};
static const double GL_X6[6]  = {-0.9815606342467191, -0.9041172563704750, -0.7699026741943050,
                                 -0.5873179542866171, -0.3678314989981802, -0.1252334085114692};
static const double GL_W6[6]  = { 0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                                  0.2031674267230659,  0.2334925365383547, 0.2491470458134029};
static const double GL_X10[10] = {-0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
                                  -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
                                  -0.5108670019508271, -0.3737060887154196, -0.2277858511416451,
                                  -0.07652652113349733};
static const double GL_W10[10] = { 0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
                                   0.08327674157670475, 0.1019301198172404,  0.1181945319615184,
                                   0.1316886384491766,  0.1420961093183821,  0.1491729864726037,
                                   0.1527533871307259};

static double bvnu(double h, double k, double r) {
  const double inf = R_PosInf;
  if (h == inf || k == inf) return 0.0;
  if (h == -inf) return (k == -inf) ? 1.0 : Phi(-k);
  if (k == -inf) return Phi(-h);
  if (r == 0.0) return Phi(-h) * Phi(-k);

  const double tp = 2.0 * M_PI;
  const double *x, *w; int lg;
  double ar = std::fabs(r);
  if (ar < 0.3)       { x = GL_X3;  w = GL_W3;  lg = 3;  }
  else if (ar < 0.75) { x = GL_X6;  w = GL_W6;  lg = 6;  }
  else                { x = GL_X10; w = GL_W10; lg = 10; }

  double hk = h * k, p = 0.0;
  if (ar < 0.925) {
    double hs = (h * h + k * k) / 2.0;
    double asr = std::asin(r) / 2.0;
    for (int i = 0; i < lg; ++i) {
      double sn = std::sin(asr * (1.0 - x[i]));
      p += w[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
      sn = std::sin(asr * (1.0 + x[i]));
      p += w[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
    }
    p = p * asr / tp + Phi(-h) * Phi(-k);
  } else {
    if (r < 0) { k = -k; hk = -hk; }
    if (ar < 1.0) {
      double as = (1.0 - r) * (1.0 + r);
      double a = std::sqrt(as);
      double bs = (h - k) * (h - k);
      double asr = -(bs / as + hk) / 2.0;
      double c = (4.0 - hk) / 8.0;
      double d = (12.0 - hk) / 80.0;
      if (asr > -100.0)
        p = a * std::exp(asr) * (1.0 - c * (bs - as) * (1.0 - d * bs) / 3.0 + c * d * as * as);
      if (hk > -100.0) {
        double b = std::sqrt(bs);
        double sp = std::sqrt(tp) * Phi(-b / a);
        p -= std::exp(-hk / 2.0) * sp * b * (1.0 - c * bs * (1.0 - d * bs) / 3.0);
      }
      a = a / 2.0;
      double s = 0.0;
      for (int i = 0; i < 2 * lg; ++i) {
        double xi = (i < lg) ? (1.0 + x[i]) : (1.0 - x[i - lg]);
        double wi = (i < lg) ? w[i] : w[i - lg];
        double xs = (a * xi) * (a * xi);
        double asr2 = -(bs / xs + hk) / 2.0;
        if (asr2 > -100.0) {
          double sp = 1.0 + c * xs * (1.0 + 5.0 * d * xs);
          double rs = std::sqrt(1.0 - xs);
          double ep = std::exp(-(hk / 2.0) * xs / ((1.0 + rs) * (1.0 + rs))) / rs;
          s += std::exp(asr2) * (sp - ep) * wi;
        }
      }
      p = (a * s - p) / tp;
    }
    if (r > 0) {
      p += Phi(-std::max(h, k));
    } else if (h >= k) {
      p = -p;
    } else {
      double L = (h < 0) ? (Phi(k) - Phi(h)) : (Phi(-h) - Phi(-k));
      p = L - p;
    }
  }
  return std::max(0.0, std::min(1.0, p));
}

// Probability a standardized BVN(rho) falls in [l1,u1] x [l2,u2]. The
// rectangle is reflected toward the upper quadrant so the 4-term difference is
// taken on small upper-tail values (limits cancellation).
static inline double rect_prob(double l1, double u1, double l2, double u2, double r) {
  if (l1 + u1 < 0) { double t = l1; l1 = -u1; u1 = -t; r = -r; }
  if (l2 + u2 < 0) { double t = l2; l2 = -u2; u2 = -t; r = -r; }
  if (r == 0.0) return (Phi(-l1) - Phi(-u1)) * (Phi(-l2) - Phi(-u2));
  double p = bvnu(l1, l2, r) - bvnu(u1, l2, r) - bvnu(l1, u2, r) + bvnu(u1, u2, r);
  return p > 0 ? p : 0.0;
}

// 1-D probability that N(0,1) (or a point mass at 0 when degenerate) lies in [l,u).
// For a degenerate dimension the caller passes unstandardised offsets and deg=true.
static inline double interval_prob(double l, double u, bool deg) {
  if (deg) return (l <= 0.0 && 0.0 < u) ? 1.0 : 0.0;
  return (l + u < 0) ? Phi(u) - Phi(l) : Phi(-l) - Phi(-u);
}

// ---------------------------------------------------------------------------
// Raster geometry helper
// ---------------------------------------------------------------------------
struct Geom {
  double xmin, ymax, rx, ry; int nr, nc;
  explicit Geom(const NumericVector& g)
    : xmin(g[0]), ymax(g[1]), rx(g[2]), ry(g[3]), nr((int) g[4]), nc((int) g[5]) {}
  // terra::cellFromXY semantics (outer right/bottom edge -> last col/row); 0 = outside
  long long cell_from_xy(double x, double y) const {
    double xmax = xmin + nc * rx, ymin = ymax - nr * ry;
    if (ISNAN(x) || ISNAN(y) || x < xmin || x > xmax || y < ymin || y > ymax) return 0;
    int col = (x == xmax) ? nc : (int) std::floor((x - xmin) / rx) + 1;
    int row = (y == ymin) ? nr : (int) std::floor((ymax - y) / ry) + 1;
    return (long long)(row - 1) * nc + col;
  }
};

// Per-observation result holder
struct ObsOut {
  std::vector<double> cell, lik;   // kept entries
  bool inad = false;               // all admissible mass ~ 0
  bool outside = false;            // neighbourhood completely outside raster
  bool bad_sigma = false;
  bool failed = false;             // exception (e.g. std::bad_alloc) inside the parallel loop
};

// Mask, normalise, truncate in place
static void finish_obs(ObsOut& o, const int* mask, double trunc) {
  double tot = 0.0;
  for (size_t j = 0; j < o.cell.size(); ++j) {
    if (mask && mask[(long long) o.cell[j] - 1]) o.lik[j] = 0.0;
    tot += o.lik[j];
  }
  if (tot <= std::sqrt(2.220446049250313e-16)) { o.inad = true; return; }  // keep raw (masked) values for return_inadmissible
  size_t k = 0;
  for (size_t j = 0; j < o.cell.size(); ++j) {
    double v = o.lik[j] / tot;
    if (v < trunc || !(v > 0)) continue;
    o.cell[k] = o.cell[j]; o.lik[k] = v; ++k;
  }
  o.cell.resize(k); o.lik.resize(k);
}

static List assemble(std::vector<ObsOut>& res) {
  const R_xlen_t n = (R_xlen_t) res.size();
  long long m = 0, mi = 0;
  for (R_xlen_t i = 0; i < n; ++i)
    if (res[i].failed) stop("Memory allocation failed while computing cell likelihoods.");
  for (R_xlen_t i = 0; i < n; ++i) {
    if (res[i].bad_sigma) stop("Invalid (non-positive-definite or NA) error covariance for observation %d.", (int)(i + 1));
    if (res[i].outside) stop("There are error buffered locations completely outside of raster area.");
    if (res[i].inad) mi += res[i].cell.size(); else m += res[i].cell.size();
  }
  NumericVector obs(m), cell(m), lik(m), iobs(mi), icell(mi), ilik(mi);
  long long a = 0, b = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    ObsOut& o = res[i];
    for (size_t j = 0; j < o.cell.size(); ++j) {
      if (o.inad) { iobs[b] = (double)(i + 1); icell[b] = o.cell[j]; ilik[b] = o.lik[j]; ++b; }
      else        { obs[a]  = (double)(i + 1); cell[a]  = o.cell[j]; lik[a]  = o.lik[j]; ++a; }
    }
    std::vector<double>().swap(o.cell); std::vector<double>().swap(o.lik);
  }
  return List::create(_["obs"] = obs, _["cell"] = cell, _["lik"] = lik,
                      _["inad_obs"] = iobs, _["inad_cell"] = icell, _["inad_lik"] = ilik);
}

// Bivariate normal cell likelihoods on a regular raster (disc neighbourhood)
//
// For each observation, every raster cell touching a disc of radius
// `max_err` receives the exact BVN probability mass in that cell. Orthant
// probabilities are evaluated once per grid vertex and shared between the
// neighbouring cells (~1 BVN evaluation per cell instead of 4). Cells with
// `mask == TRUE` get 0, rows are normalised and values below `trunc` are
// dropped, all inside one parallel pass. Observations with `max_err <= 0`
// or `NA` get likelihood 1 in the cell containing the location.
// @keywords internal
// [[Rcpp::export]]
List bvn_lik_disc(NumericVector x, NumericVector y,
                  NumericVector sxx, NumericVector syy, NumericVector sxy,
                  NumericVector max_err, NumericVector geom, LogicalVector mask,
                  double trunc, int nthreads = 1) {
  if (x.size() > INT_MAX) stop("Too many observations.");
  const int n = (int) x.size();
  if (nthreads < 1) nthreads = 1;
  const Geom G(geom);
  const int* pm = mask.size() ? LOGICAL(mask) : nullptr;
  std::vector<ObsOut> res(n);
  const double *px = x.begin(), *py = y.begin(), *pxx = sxx.begin(), *pyy = syy.begin(),
               *pxy = sxy.begin(), *pme = max_err.begin();

#ifdef _OPENMP
  // static (round-robin chunks) avoids __kmpc_dispatch_* symbols, so the package
  // loads even when the libomp runtime is older than the compiler (e.g. Apple clang 17).
  #pragma omp parallel for schedule(static, 16) num_threads(nthreads)
#endif
  for (int i = 0; i < n; ++i) {
    ObsOut& o = res[i];
    try {   // no C++ exception may escape an OpenMP region
    double R = pme[i];
    if (!(R > 0)) {                                  // exact-location observation
      long long c = G.cell_from_xy(px[i], py[i]);
      if (c == 0) { o.outside = true; continue; }
      o.cell.push_back((double) c); o.lik.push_back(1.0);
      finish_obs(o, pm, trunc);
      continue;
    }
    double sx = std::sqrt(pxx[i]), sy = std::sqrt(pyy[i]);
    if (!(sx >= 0) || !(sy >= 0)) { o.bad_sigma = true; continue; }
    const bool degx = (sx == 0.0), degy = (sy == 0.0);
    double rho = (degx || degy) ? 0.0 : pxy[i] / (sx * sy);
    if (!(std::fabs(rho) <= 1.0)) { o.bad_sigma = true; continue; }
    const double dsx = degx ? 1.0 : sx, dsy = degy ? 1.0 : sy;   // degenerate dims kept unstandardised
    if (!std::isfinite(px[i]) || !std::isfinite(py[i]) || !std::isfinite(R)) { o.outside = true; continue; }

    // bounding window of rows / cols (1-based), clipped to raster
    // (clamp in double before converting so out-of-range values never overflow int)
    auto idx = [](double v, int lo, int hi) {
      double f = std::floor(v) + 1.0;
      return (int) std::min(std::max(f, (double) lo - 1.0), (double) hi + 1.0);
    };
    int c0 = std::max(idx((px[i] - R - G.xmin) / G.rx, 1, G.nc), 1);
    int c1 = std::min(idx((px[i] + R - G.xmin) / G.rx, 1, G.nc), G.nc);
    int r0 = std::max(idx((G.ymax - (py[i] + R)) / G.ry, 1, G.nr), 1);
    int r1 = std::min(idx((G.ymax - (py[i] - R)) / G.ry, 1, G.nr), G.nr);
    if (c1 < c0 || r1 < r0) { o.outside = true; continue; }
    const int W = c1 - c0 + 1, H = r1 - r0 + 1;

    // standardized vertex coordinates
    std::vector<double> zx(W + 1), zy(H + 1);
    for (int a = 0; a <= W; ++a) zx[a] = (G.xmin + (c0 - 1 + a) * G.rx - px[i]) / dsx;   // left -> right
    for (int b = 0; b <= H; ++b) zy[b] = (G.ymax - (r0 - 1 + b) * G.ry - py[i]) / dsy;   // top -> bottom

    // cells touching the disc
    std::vector<unsigned char> in((size_t) W * H, 0);
    const double R2 = R * R, hx = G.rx / 2.0, hy = G.ry / 2.0;
    size_t nin = 0;
    for (int b = 0; b < H; ++b) {
      double cy = G.ymax - (r0 - 1 + b + 0.5) * G.ry;
      double dy = std::max(0.0, std::fabs(py[i] - cy) - hy);
      for (int a = 0; a < W; ++a) {
        double cx = G.xmin + (c0 - 1 + a + 0.5) * G.rx;
        double dx = std::max(0.0, std::fabs(px[i] - cx) - hx);
        if (dx * dx + dy * dy <= R2) { in[(size_t) b * W + a] = 1; ++nin; }
      }
    }
    o.cell.reserve(nin); o.lik.reserve(nin);

    if (rho == 0.0) {
      // separable: product of 1-D interval probabilities (tail-accurate form)
      std::vector<double> qx(W), qy(H);
      for (int a = 0; a < W; ++a) qx[a] = interval_prob(zx[a], zx[a + 1], degx);
      for (int b = 0; b < H; ++b) qy[b] = interval_prob(zy[b + 1], zy[b], degy);
      for (int b = 0; b < H; ++b) for (int a = 0; a < W; ++a) if (in[(size_t) b * W + a]) {
        o.cell.push_back((double)((long long)(r0 - 1 + b) * G.nc + (c0 + a)));
        o.lik.push_back(std::max(0.0, qx[a] * qy[b]));
      }
    } else {
      // Orthant probabilities at grid vertices, cached per quadrant orientation
      // (the orientation is chosen per cell so differences are taken on the
      // small tail side, as in rect_prob()).
      const size_t NV = (size_t)(W + 1) * (H + 1);
      std::vector<double> Gpp(NV, -1.0), Gnn(NV, -1.0), Gpn(NV, -1.0), Gnp(NV, -1.0);
      auto get = [](std::vector<double>& cache, size_t v, double h, double k, double r) {
        double& c = cache[v]; if (c < 0) c = bvnu(h, k, r); return c;
      };
      for (int b = 0; b < H; ++b) for (int a = 0; a < W; ++a) if (in[(size_t) b * W + a]) {
        double l1 = zx[a], u1 = zx[a + 1], l2 = zy[b + 1], u2 = zy[b];
        size_t vLL = (size_t)(b + 1) * (W + 1) + a, vUL = (size_t)(b + 1) * (W + 1) + a + 1; // (l1,l2),(u1,l2)
        size_t vLU = (size_t) b * (W + 1) + a,      vUU = (size_t) b * (W + 1) + a + 1;      // (l1,u2),(u1,u2)
        bool fx = (l1 + u1 < 0), fy = (l2 + u2 < 0);
        double p;
        if (!fx && !fy) {        // P(X>.,Y>.)
          p = get(Gpp, vLL, l1, l2, rho) - get(Gpp, vUL, u1, l2, rho)
            - get(Gpp, vLU, l1, u2, rho) + get(Gpp, vUU, u1, u2, rho);
        } else if (fx && fy) {   // P(X<.,Y<.) = P(-X>-.,-Y>-.)
          p = get(Gnn, vUU, -u1, -u2, rho) - get(Gnn, vLU, -l1, -u2, rho)
            - get(Gnn, vUL, -u1, -l2, rho) + get(Gnn, vLL, -l1, -l2, rho);
        } else if (fx && !fy) {  // P(X<.,Y>.) = P(-X>-.,Y>.), corr -rho
          p = get(Gnp, vUL, -u1, l2, -rho) - get(Gnp, vLL, -l1, l2, -rho)
            - get(Gnp, vUU, -u1, u2, -rho) + get(Gnp, vLU, -l1, u2, -rho);
        } else {                 // P(X>.,Y<.) = P(X>.,-Y>-.), corr -rho
          p = get(Gpn, vLU, l1, -u2, -rho) - get(Gpn, vUU, u1, -u2, -rho)
            - get(Gpn, vLL, l1, -l2, -rho) + get(Gpn, vUL, u1, -l2, -rho);
        }
        o.cell.push_back((double)((long long)(r0 - 1 + b) * G.nc + (c0 + a)));
        o.lik.push_back(p > 0 ? p : 0.0);
      }
    }
    finish_obs(o, pm, trunc);
    } catch (...) { o.failed = true; }
  }
  return assemble(res);
}

// BVN cell likelihoods for an explicit list of (obs, cell) pairs
// (used for `neighborhood = "buffer"`, i.e. cells from terra::buffer + terra::cells).
// Pairs must be grouped by `obs` in ascending order. Masking, normalising and
// truncation as in `bvn_lik_disc()`.
// @keywords internal
// [[Rcpp::export]]
List bvn_lik_list(IntegerVector obs, NumericVector cell,
                  NumericVector x, NumericVector y,
                  NumericVector sxx, NumericVector syy, NumericVector sxy,
                  NumericVector max_err, NumericVector geom, LogicalVector mask,
                  double trunc, int nthreads = 1) {
  if (x.size() > INT_MAX) stop("Too many observations.");
  const int n = (int) x.size();
  if (nthreads < 1) nthreads = 1;
  const Geom G(geom);
  const int* pm = mask.size() ? LOGICAL(mask) : nullptr;
  std::vector<ObsOut> res(n);
  std::vector<long long> st(n + 1, 0);
  for (R_xlen_t j = 0; j < obs.size(); ++j) {
    if (obs[j] < 1 || obs[j] > n || (j > 0 && obs[j] < obs[j - 1])) stop("'obs' must be sorted and in 1..n.");
    st[obs[j]]++;
  }
  for (int i = 0; i < n; ++i) st[i + 1] += st[i];
  const double* pc = cell.begin();
  const double *px = x.begin(), *py = y.begin(), *pxx = sxx.begin(), *pyy = syy.begin(),
               *pxy = sxy.begin(), *pme = max_err.begin();

#ifdef _OPENMP
  // static (round-robin chunks) avoids __kmpc_dispatch_* symbols, so the package
  // loads even when the libomp runtime is older than the compiler (e.g. Apple clang 17).
  #pragma omp parallel for schedule(static, 16) num_threads(nthreads)
#endif
  for (int i = 0; i < n; ++i) {
    ObsOut& o = res[i];
    try {   // no C++ exception may escape an OpenMP region
    if (!(pme[i] > 0)) {
      long long c = G.cell_from_xy(px[i], py[i]);
      if (c == 0) { o.outside = true; continue; }
      o.cell.push_back((double) c); o.lik.push_back(1.0);
      finish_obs(o, pm, trunc); continue;
    }
    if (st[i + 1] == st[i]) { o.outside = true; continue; }
    double sx = std::sqrt(pxx[i]), sy = std::sqrt(pyy[i]);
    if (!(sx >= 0) || !(sy >= 0)) { o.bad_sigma = true; continue; }
    const bool degx = (sx == 0.0), degy = (sy == 0.0);
    double rho = (degx || degy) ? 0.0 : pxy[i] / (sx * sy);
    if (!(std::fabs(rho) <= 1.0)) { o.bad_sigma = true; continue; }
    const double dsx = degx ? 1.0 : sx, dsy = degy ? 1.0 : sy;
    o.cell.reserve(st[i + 1] - st[i]); o.lik.reserve(st[i + 1] - st[i]);
    for (long long j = st[i]; j < st[i + 1]; ++j) {
      long long c = (long long) pc[j] - 1, row = c / G.nc, col = c % G.nc;
      double cx = G.xmin + ((double) col + 0.5) * G.rx, cy = G.ymax - ((double) row + 0.5) * G.ry;
      o.cell.push_back(pc[j]);
      double l1 = (cx - G.rx / 2 - px[i]) / dsx, u1 = (cx + G.rx / 2 - px[i]) / dsx;
      double l2 = (cy - G.ry / 2 - py[i]) / dsy, u2 = (cy + G.ry / 2 - py[i]) / dsy;
      o.lik.push_back((degx || degy) ? interval_prob(l1, u1, degx) * interval_prob(l2, u2, degy)
                                     : rect_prob(l1, u1, l2, u2, rho));
    }
    finish_obs(o, pm, trunc);
    } catch (...) { o.failed = true; }
  }
  return assemble(res);
}

// OpenMP status for this build of walk: whether it was compiled with OpenMP,
// the runtime thread limit, and the number of processors OpenMP can see.
// [[Rcpp::export]]
List omp_info() {
#ifdef _OPENMP
  return List::create(_["enabled"] = true, _["max_threads"] = omp_get_max_threads(),
                      _["thread_limit"] = omp_get_thread_limit(), _["num_procs"] = omp_get_num_procs());
#else
  return List::create(_["enabled"] = false, _["max_threads"] = 1, _["thread_limit"] = 1, _["num_procs"] = NA_INTEGER);
#endif
}
