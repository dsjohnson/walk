// [[Rcpp::depends(RcppArmadillo)]]

#include "walk_types.h"

using namespace Rcpp;
using namespace arma;

// Calculate likelihood ///////////////

// [[Rcpp::export]]
Rcpp::List ctmc_predict_arma(
    const arma::sp_mat& L, 
    const arma::vec& obs, 
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
  if(form==1){
    Q = load_Q(from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, l_r, u_r, link_m, a_m, norm, clip);
   } else {
    Q = load_Q_sde(from_to, Xb_q_r, Xb_q_m, hij, ns, k, clip);
  }
  
  
  // Forward probs
  arma::mat A(N, ns);
  A.row(0) = delta;
  // Backward probs
  arma::mat B(ns, N);
  B.col(N-1).ones();
  
  // State posterior matrix
  arma::sp_mat G(N, ns);
  
  arma::rowvec v(ns);
  arma::rowvec ab(ns);
  
  // Start Forward alg loop (index = i)
  for(int i=1; i<N; i++){
    v = v_exp_Q_t(A.row(i-1), Q, dt(i), eq_prec);
    if(obs(i)==1) v = v % ((1-p)*L.row(i)) + (p/ns)*v;
    A.row(i) = v/accu(v);
  } // end i
  
  G.row(N-1) = A.row(N-1);
  
  // Start backward loop (index i)
  for(int i=N-1; i>0; i--){
    //v = phi_exp_lnG(B.col(i).t(), (Q*dt(i)).t(), eq_prec);
    if(obs(i)==1) v = v % ((1-p)*L.row(i)) + (p/ns)*v;
    B.col(i-1) = (v/accu(v)).t();
    ab =  A.row(i-1) % B.col(i-1).t();
    ab = ab/accu(ab);
    G.row(i-1) = ab.clean(trunc_tol);
  } // end i
  // G = normalise(G, 1, 1);
  
  return Rcpp::List::create(
    Rcpp::Named("local_state_prob") = G,
    Rcpp::Named("alpha") = A,
    Rcpp::Named("beta") = B
  );
  
}