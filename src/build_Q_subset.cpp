// Construct a local rate matrix Q_sub with an explicit absorbing sink at index n_sub
arma::sp_mat build_Q_sub_absorbing(const arma::sp_mat& Q, const arma::uvec& sub_set) {
  uword n_sub = sub_set.n_elem;
  arma::sp_mat Q_raw = Q.submat(sub_set, sub_set);
  
  // Calculate total row sums of the raw submatrix
  // For a valid Markov rate matrix, row sums should equal 0.
  // Any positive discrepancy (-row_sum) represents rate lost to out-of-bounds states.
  arma::vec row_sums = arma::vec(Q_raw * arma::ones<arma::vec>(n_sub));
  
  // Create augmented sparse matrix of size (n_sub + 1) x (n_sub + 1)
  arma::sp_mat Q_sub(n_sub + 1, n_sub + 1);
  
  // Copy interior transition rates
  Q_sub.submat(0, 0, n_sub - 1, n_sub - 1) = Q_raw;
  
  // Direct out-of-bounds transition rates into the sink column (index n_sub)
  for(uword i = 0; i < n_sub; ++i) {
    double rate_to_sink = -row_sums(i);
    if(rate_to_sink > 1e-12) {
      Q_sub(i, n_sub) = rate_to_sink;
    }
  }
  // Row n_sub remains 0 (absorbing state)
  
  return Q_sub;
}



/////////
////////
////////
// add to likelihood function 
// 1. Build local subgrid with sink
arma::sp_mat Q_sub = build_Q_sub_absorbing(Q, sub_set);

// 2. Pad phi_sub with 0.0 for the sink state
arma::rowvec phi_sub(sub_set.n_elem + 1, fill::zeros);
phi_sub.head(sub_set.n_elem) = phi.cols(sub_set);

// 3. Propagate uniformization on (n_sub + 1) dimension
arma::rowvec v_sub = phi_exp_lnG(phi_sub, Q_sub * dt(i), eq_prec);

// 4. Map only the interior states back to full v (discard the sink value at v_sub(n_sub))
v.zeros();
v.cols(sub_set) = v_sub.head(sub_set.n_elem);
