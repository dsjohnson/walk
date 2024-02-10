// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// function prototypes
arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8);
arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const bool& row_sweep=true);
arma::sp_mat load_Q_sp(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const double& a = 1.0, const bool& row_sweep=true);


// Calculate likelihood ///////////////
// [[Rcpp::export]]
Rcpp::List ctmc_predict_arma(
    const arma::sp_mat& L, 
    const arma::vec& obs, 
    const arma::vec& dt, 
    const int& ns, 
    const arma::umat& from_to, 
    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
    const double& p,
    const arma::rowvec& delta, 
    const double& eq_prec = 1.0e-8,
    const double& trunc_tol = 1.0e-8,
    const int& link = 1,
    const double& a = 1.0,
    const bool& row_sweep=true)
{
  int N = dt.size();
  arma::sp_mat Q;
  if(link==1){
    Q = load_Q_sp(from_to, Xb_q_r, Xb_q_m, ns, a, row_sweep);
  } else{
    Q = load_Q(from_to, Xb_q_r, Xb_q_m, ns, row_sweep);
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
    v = phi_exp_lnG(A.row(i-1), Q*dt(i), eq_prec);
    if(obs(i)==1) v = v % ((1-p)*L.row(i)) + (p/ns)*v;
    A.row(i) = v/accu(v);
  } // end i
  
  G.row(N-1) = A.row(N-1);
  
  // Start backward loop (index i)
  for(int i=N-1; i>0; i--){
    v = phi_exp_lnG(B.col(i).t(), (Q*dt(i)).t(), eq_prec);
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