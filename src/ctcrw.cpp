
#include <TMB.hpp>
#include <cmath>

using namespace density;
using std::sqrt;


// Function for creating cov mat for ctcrw SSM
template<class Type>
matrix<Type> make_Q(Type beta, Type sigma2, Type dt){
  matrix<Type> Q(4,4);
  Q.setZero();
  Q(0,0) = sigma2 * (dt -(2/beta)*(1-exp(-beta*dt)) + (1/(2*beta))*(1-exp(-2*beta*dt)));
  Q(1,1) = sigma2*beta*(1-exp(-2*beta*dt))/2;
  Q(0,1) = sigma2*(1 - 2*exp(-beta*dt) + exp(-2*beta*dt))/2;
  Q(1,0) = sigma2*(1 - 2*exp(-beta*dt) + exp(-2*beta*dt))/2;
  Q(2,2) = sigma2 * (dt -(2/beta)*(1-exp(-beta*dt)) + (1/(2*beta))*(1-exp(-2*beta*dt)));
  Q(3,3) = sigma2*beta*(1-exp(-2*beta*dt))/2;
  Q(2,3) = sigma2*(1 - 2*exp(-beta*dt) + exp(-2*beta*dt))/2;
  Q(3,2) = sigma2*(1 - 2*exp(-beta*dt) + exp(-2*beta*dt))/2;
  return Q;
}

// Function for creating projection matrix for ctcrw SSM
template<class Type>
matrix<Type> make_T(Type beta, Type dt){
  matrix<Type> T(4,4);
  T.setZero();
  T(0,0) = 1; 
  T(0,1) = (1-exp(-beta*dt))/beta;
  T(1,1) = exp(-beta*dt);
  T(2,2) = 1;
  T(2,3) = (1-exp(-beta*dt))/beta;
  T(3,3) = exp(-beta*dt);
  return T;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_MATRIX(Y);	            //  (x, y) observations
  DATA_VECTOR(dt);         	//  time diff in some appropriate unit. this should contain dt for both interp and obs positions.
  DATA_VECTOR(alpha0);        //  initial state mean
  DATA_SCALAR(Vmu0);          // Prior variance of first location
  DATA_IVECTOR(isd);          //  indexes observations vs. interpolation points
  DATA_IVECTOR(obs_mod);      //  indicates which obs error model to be used
  
  // for KF observation model
  DATA_VECTOR(m);             //  m is the semi-minor axis length
  DATA_VECTOR(M);             //  M is the semi-major axis length
  DATA_VECTOR(c);             //  c is the orientation of the error ellipse
  // for LS observation model
  DATA_MATRIX(K);                 // error weighting factors for LS obs model
  
  // PROCESS PARAMETERS
  // for CRW
  PARAMETER(log_beta);				
  PARAMETER(log_sigma);
  // random variables
  PARAMETER_MATRIX(alpha); /* State matrix */
  
  // OBSERVATION PARAMETERS
  // for KF OBS MODEL
  PARAMETER(l_psi); 				// error SD scaling parameter to account for possible uncertainty in Argos error ellipse variables
  // for LS OBS MODEL
  PARAMETER_VECTOR(l_tau);     	// error dispersion for LS obs model (log scale)
  PARAMETER(l_rho_o);             // error correlation
  
  // Transform parameters
  Type beta = exp(log_beta);
  Type sigma2 = exp(2*log_sigma);
  Type psi = exp(l_psi);
  vector<Type> tau = exp(l_tau);
  Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);
  
  int timeSteps = dt.size();
  
  /* Define likelihood */
  //parallel_accumulator<Type> jnll(this);
  Type nll_proc=0.0;
  Type nll_obs=0.0;
  Type nll=0.0;
  vector<Type> x;
  
  // CRW
  // Setup object for evaluating multivariate normal likelihood
  matrix<Type> Q(4,4);
  matrix<Type> T(4,4);
  matrix<Type> mu(2,timeSteps);
  matrix<Type> cov_obs(2, 2);
  cov_obs.setZero();
  
  // Set up initial Q
  Q(0,0) = Vmu0;
  Q(1,1) = sigma2*beta/2;
  Q(2,2) = Vmu0;
  Q(3,3) = sigma2*beta/2;
  x = alpha.col(0)-alpha0.matrix();
  nll_proc += MVNORM(Q)(x);

  for(int i = 1; i < timeSteps; i++) {
    Q = make_Q(beta, sigma2, dt(i));
    T = make_T(beta, dt(i));
    x = alpha.col(i) - T * alpha.col(i-1);
    nll_proc += MVNORM(Q)(x);
  }
  
  // OBSERVATION MODEL
  // 2 x 2 covariance matrix for observations
  
  for(int i = 0; i < timeSteps; ++i) {
  mu(0,i) = alpha(0,i);
  mu(1,i) = alpha(2,i);
  if(isd(i) == 1) {
    if(obs_mod(i) == 0) {
      // Argos Least Squares observations
      Type s = tau(0) * K(i,0);
      Type q = tau(1) * K(i,1);
      cov_obs(0,0) = s * s;
      cov_obs(1,1) = q * q;
      //cov_obs(0,1) = s * q * rho_o;
      //cov_obs(1,0) = cov_obs(0,1);
    } else if(obs_mod(i) == 1) {
      // Argos Kalman Filter (or Kalman Smoothed) observations
      Type z = sqrt(Type(2.0));
      Type s2c = sin(c(i)) * sin(c(i));
      Type c2c = cos(c(i)) * cos(c(i));
      Type M2  = (M(i) / z) * (M(i) / z);
      Type m2 = (m(i) * psi / z) * (m(i) * psi / z);
      cov_obs(0,0) = (M2 * s2c + m2 * c2c);
      cov_obs(1,1) = (M2 * c2c + m2 * s2c);
      cov_obs(0,1) = (0.5 * (M(i) * M(i) - (m(i) * psi * m(i) * psi))) * cos(c(i)) * sin(c(i));
      cov_obs(1,0) = cov_obs(0,1);
    }
    nll_obs +=  MVNORM(cov_obs)(Y.col(i) - mu.col(i));   // CRW innovations


  }
  }
  
  nll = nll_proc + nll_obs;
  
  ADREPORT(beta);
  ADREPORT(sigma2);
  ADREPORT(rho_o);
  ADREPORT(tau);
  ADREPORT(psi);
  
  return nll;
}

