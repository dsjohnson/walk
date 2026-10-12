
// [[Rcpp::depends(RcppArmadillo)]]

#include "walk_types.h"

 arma::mat v_exp_Q_t(
     const arma::rowvec& v, const arma::sp_mat& Q, double t, double tolerance) {
   int N = Q.n_rows;
   
   // Handle the trivial edge case instantly
   if (t <= 0.0) {
     return arma::conv_to<arma::mat>::from(v);
   }
   
   arma::spdiagview<double> diag_Q = Q.diag();
   double alpha_0 = 0.0;
   for (arma::uword i = 0; i < diag_Q.n_elem; ++i) {
     double abs_val = std::abs(diag_Q(i));
     if (abs_val > alpha_0) {
       alpha_0 = abs_val;
     }
   }
   
   // 1. Calculate the time-dependent uniformization rate
   double alpha = alpha_0 * t;
   
   // 2. Dynamically calculate the safe upper truncation bound 'm'
   int m = 0;
   double running_right_sum = 0.0;
   while (running_right_sum < (1.0 - tolerance)) {
     double log_p_m = m * std::log(alpha) - alpha - std::lgamma(m + 1);
     running_right_sum += std::exp(log_p_m);
     m++;
     
     // Safety guard for extreme edge cases / infinite loops
     if (m > 1000000) break; 
   }
   
   // 3. Dynamically calculate the lower truncation bound 'm_lower'
   int m_lower = 0;
   double running_left_sum = 0.0;
   for (int n = 0; n <= m; ++n) {
     double log_p_n = n * std::log(alpha) - alpha - std::lgamma(n + 1);
     running_left_sum += std::exp(log_p_n);
     if (running_left_sum < tolerance) {
       m_lower = n; 
     } else {
       break; // Active numerical window reached
     }
   }
   
   // 4. Build the Uniformized Transition Matrix for this specific t scale
   arma::sp_mat I;
   I.eye(N, N);
   arma::sp_mat P = I + (Q / alpha_0); // Matrix remains scale-invariant!
   
   // 5. Setup dense tracking containers (zero reallocations)
   arma::rowvec v_curr = v;
   arma::rowvec v_next(N);
   arma::rowvec result(N, arma::fill::zeros); 
   
   // 6. Bare-Metal Execution Loop
   for (int n = 0; n <= m; ++n) {
     
     // Only accumulate if we are within the active Poisson window
     if (n >= m_lower) {
       double log_p_n = n * std::log(alpha) - alpha - std::lgamma(n + 1);
       double p_n = std::exp(log_p_n);
       
       result += p_n * v_curr;
     }
     
     // Super high-speed Dense Vector * Sparse Matrix multiplication
     if (n < m) {
       v_next = v_curr * P;
       v_curr = v_next; // Blazing fast element assignment
     }
   }  
   
   return result;
 }