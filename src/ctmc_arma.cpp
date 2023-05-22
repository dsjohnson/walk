// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;



// [[Rcpp::export]]
arma::mat phi_exp_G(const arma::mat& v, const arma::sp_mat&  G, 
                    const double& prec=1.0e-8) {
  arma::mat out = expQ2::sv_exp_Q(v, G, prec, false, true);
  return out;
}


// [[Rcpp::export]]
arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& idx_q, 
                    const arma::vec& Xb_q, const arma::vec& off_q, 
                    const int& ns){
  int n = from_to.n_cols;
  arma::vec qvals(n); 
  for(int i=0; i<n; i++){qvals(i) = exp(off_q(i) + Xb_q(idx_q(i)));}
  arma::sp_mat Q(from_to, qvals, ns, ns);
  //arma::colvec ones(ns,fill::ones);
  arma::colvec row_sums = Q * ones(ns);
  Q.diag() = -1.0*row_sums;
  return Q;
}



// Calculate likelihood ///////////////
// [[Rcpp::export]]
Rcpp::List ctmc_arma(const arma::vec& id, const  arma::vec& period, 
                     const arma::vec& dt, const arma::vec& cell, 
                     const int& ns, const int& np, 
                     const arma::umat& from_to_q, const arma::mat& X_q, 
                     const arma::vec& off_q,
                     const arma::vec& idx_q, const arma::vec& beta_q)
{
  int N = cell.size();
  
  arma::vec Xb_q = X_q * beta_q;
  
  // arma::vec q_vals = exp();
  // arma::vec l_vals = exp();

  double u = 0.0;
  arma::vec log_lik_v(N, fill::zeros);
  
  arma::sp_mat Q = load_Q(from_to_q, idx_q, Xb_q, off_q, ns);
  arma::sp_mat G(ns,ns);

  // Start forward loop
  arma::rowvec v(ns);
  arma::rowvec phi(ns);
  phi(cell(0)) = 1.0;
  arma::sp_mat P(ns,ns);
  arma::sp_mat L(ns,ns);
  
  // Start Foward alg loop (index = i)
  for(int i=1; i<N; i++){
    if(id(i)!=id(i-1)){
      phi = 0.0*phi; phi(cell(i)) = 1.0;
    } else{
      // Q = load_Q(from, to,Q_mat.col(period(i)), ns); #mod here for dynamic movement
      G = Q*dt(i);
      v = phi_exp_G(phi, G);
      P = 0.0*P;
      if(R_finite(cell(i))){
        P(cell(i),cell(i)) = L(cell(i), cell(i));
        v = v * P; 
      }
      u = accu(v);
      log_lik_v(i) = log(u);
      phi = v/u;
    }
  } // end i
  
  double n2ll = -2*accu(log_lik_v);
  
  return Rcpp::List::create(
    // Rcpp::Named("log_lik_v") = log_lik_v,
    // Rcpp::Named("Q") = Q,
    // Rcpp::Named("G") = G,
    // Rcpp::Named("L_mat") = L_mat,
    // Rcpp::Named("L") = L,
    // Rcpp::Named("P") = P,
    // Rcpp::Named("v") = v,
    // Rcpp::Named("phi") = phi,
    // Rcpp::Named("N") = N,
    Rcpp::Named("n2ll") = n2ll
  );
  
}