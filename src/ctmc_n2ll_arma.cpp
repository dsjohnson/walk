// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// function prototypes
arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8);
arma::sp_mat load_Q_mult(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const int& link_r=1, const double& a_r=1.0, const double& l_r=0.0, const double& u_r=0.0, const int& link_m=1, const double& a_m=1.0, const bool& norm=true, const double& clip=0.0);
arma::sp_mat load_Q_add(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const int& link_r=1, const double& a_r=1.0, const int& link_m=1, const double& a_m=1.0, const double& clip=0.0);
arma::sp_mat load_Q_sde(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const double& k,const double& a_r=1.0);


// Calculate likelihood ///////////////
// [[Rcpp::export]]
Rcpp::List ctmc_n2ll_arma(
    const arma::sp_mat& L, 
    const arma::vec& dt, 
    const int& ns,
    const arma::umat& from_to, 
    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
    const double& p,
    const arma::rowvec& delta, 
    const double& eq_prec = 1.0e-8,
    const int& link_r = 1,
    const double& a_r = 1.0, 
    const double& l_r = 0.0,
    const double& u_r = 0.0,
    const int& link_m = 1,
    const double& a_m = 1.0, 
    const int& form = 1,
    const double& k = 2.0,
    const bool& norm=true,
    const double& clip=0.0)
{
  int N = dt.size();
  double u = 0.0;
  arma::vec log_lik_v(N);
  arma::sp_mat Q;
  if(form==1){
    Q = load_Q_mult(from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, l_r, u_r, link_m, a_m, norm, clip);
  } else if(form==2){
    Q = load_Q_add(from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, link_m, a_m, clip);
  } else if(form==3){
    Q = load_Q_sde(from_to, Xb_q_r, Xb_q_m, ns, k, a_r);
  }
  double qii = abs(Q.diag()).max();
 arma::vec rho(N);
  
  // Start forward loop
  arma::rowvec v(ns);
  arma::sp_mat P(ns,ns);
  // arma::vec Lt(ns);
  arma::rowvec phi = delta % ((1-p)*L.row(0)) + (p/ns)*v;
  u = accu(phi);
  log_lik_v(0) = log(u);
  phi = phi/u;
  
  // Start Forward alg loop (index = i)
  for(int i=1; i<N; i++){
    rho(i) = R::qpois(eq_prec, qii*dt(i), true, false);
    v = phi_exp_lnG(phi, Q*dt(i), eq_prec);
    v = v % ((1-p)*L.row(i)) + (p/ns)*v;
    u = accu(v);
    // if(u==0){
    // log_lik_v(i) = log(u);
    // v = phi_exp_lnG(phi, Qt);
    // arma::rowvec Lrow;
    // Lrow = L.row(i);
    // return Rcpp::List::create(
    //   Rcpp::Named("i") = i+1,
    //   Rcpp::Named("log_lik_v") = log_lik_v,
    //   Rcpp::Named("v") = v,
    //   Rcpp::Named("Lrow") = Lrow
    // );
    // }
    log_lik_v(i) = log(u);
    phi = v/u;
  } // end i
  
  //Rcout << "max rho: " << rho.max() << "\n";
  
  double n2ll = -2*accu(log_lik_v);
  
  return Rcpp::List::create(
    Rcpp::Named("log_lik_v") = log_lik_v,
    Rcpp::Named("rho") = rho,
    Rcpp::Named("n2ll") = n2ll
  );
  
}