// // [[Rcpp::depends(RcppArmadillo)]]
// 
// #include <RcppArmadillo.h>
// #include <expQ2.h>
// 
// using namespace Rcpp;
// using namespace expQ2;
// using namespace arma;
// 
// // function prototypes
// arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8);
// arma::sp_mat load_Q_mult(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const int& link_r=1, const double& a_r=1.0, const double& l_r=0.0, const double& u_r=0.0, const int& link_m=1, const double& a_m=1.0, const bool& norm=true, const double& clip=0.0);
// arma::sp_mat load_Q_add(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const int& link_r=1, const double& a_r=1.0, const int& link_m=1, const double& a_m=1.0, const double& clip=0.0);
// arma::sp_mat load_Q_sde(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const double& k,const double& a_r=1.0);
// 
// double ctmc_n2ll_arma(const arma::sp_mat& L, const arma::vec& dt, const int& ns,const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const double& p, const arma::rowvec& delta, const double& eq_prec = 1.0e-8, const int& link_r = 1, const double& a_r = 1.0,  const double& l_r = 0.0, const double& u_r = 0.0, const int& link_m = 1, const double& a_m = 1.0, const int& form = 1, const double& k = 2.0, const bool& norm=true, const double& clip=0.0);
// 
