#ifndef WALK_TYPES_H
#define WALK_TYPES_H

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// Function prototypes
arma::mat v_exp_Q_t(const arma::rowvec& v, const arma::sp_mat& Q, double t=1, double tolerance=1.0e-12); 
  
arma::sp_mat sp_mat_div(const arma::sp_mat& X, const arma::sp_mat& Y);

arma::vec stat_dist(const arma::sp_mat& Q);

arma::vec logit(const arma::vec& x, double L = 0.0, double U = 0.0);

arma::vec soft_plus(const arma::vec& x, double a = 1.0);

arma::vec hard_plus(const arma::vec& x);

arma::sp_mat clip_Q(const arma::sp_mat& Q, double clip);

// arma::mat phi_exp_lnG(const arma::rowvec& v, const arma::sp_mat& Q, double t, double prec);

arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, 
                    int ns, int link_r = 1, double a_r = 1.0, double l_r = 0.0, double u_r = 0.0, 
                    int link_m = 1, double a_m = 1.0, bool norm = true, double clip = 0.0);

arma::sp_mat load_Q_sde(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, 
                        const arma::vec& hij, int ns, double k, double clip = 0.0);

#endif // WALK_TYPES_H