// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// [[Rcpp::export]]
double ctmc_n2ll_arma_precomputed(
    const arma::sp_mat& L, 
    const arma::vec& dt, 
    const int& ns,
    const arma::umat& from_to, 
    const arma::vec& Xb_q_r, 
    const arma::vec& Xb_q_m,
    const double& p,
    const arma::rowvec& delta, 
    const arma::vec& hij,
    const arma::uvec& active_indices, // Pre-computed 1D vector of all indices
    const arma::uvec& active_offsets, // Offsets array of length N+1
    const double& eq_prec = 1.0e-8,
    const int& link_r = 1,
    const double& a_r = 1.0, 
    const double& l_r = 0.0,
    const double& u_r = 0.0,
    const int& link_m = 1,
    const double& a_m = 1.0, 
    const int& form = 1,
    const double& k = 2.0,
    const bool& norm = true,
    const double& clip = 0.0)
{
  int N = dt.size();
  double u = 0.0;
  arma::vec log_lik_v(N);
  
  // 1. Build Q once
  arma::sp_mat Q;
  if(form == 1){
    Q = load_Q(from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, l_r, u_r, link_m, a_m, norm, clip);
  } else {
    Q = load_Q_sde(from_to, Xb_q_r, Xb_q_m, hij, ns, k, clip);
  }
  
  // 2. Initialize step 0
  arma::rowvec v(ns, fill::zeros);
  arma::rowvec phi = delta % ((1.0 - p) * L.row(0)) + (p / ns);
  
  u = accu(phi);
  log_lik_v(0) = log(u);
  phi = phi / u;
  
  // 3. Forward Loop using Pre-computed Active Subgrids
  for(int i = 1; i < N; i++) {
    
    // Get pre-computed active state slice for step i
    uword start_idx = active_offsets(i);
    uword end_idx   = active_offsets(i + 1) - 1;
    arma::uvec sub_set = active_indices.subvec(start_idx, end_idx);
    
    // Slice phi and Q locally
    arma::rowvec phi_sub = phi.cols(sub_set);
    arma::sp_mat Q_sub   = Q.submat(sub_set, sub_set);
    
    // Propagate only on local subgrid
    arma::rowvec v_sub = phi_exp_lnG(phi_sub, Q_sub * dt(i), eq_prec);
    
    // Reconstruct full v
    v.zeros();
    v.cols(sub_set) = v_sub;
    
    // Apply emission likelihood and normalize
    v = v % ((1.0 - p) * L.row(i)) + (p / ns) * v;
    u = accu(v);
    
    if(u < 1e-300) {
      log_lik_v(i) = -1e10;
      phi.fill(1.0 / ns);
    } else {
      log_lik_v(i) = log(u);
      phi = v / u;
    }
  }
  
  return -2.0 * accu(log_lik_v);
}