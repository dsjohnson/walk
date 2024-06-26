// expm_generator:  v^T %*% exp(M) ->  M( from, to ) AND rowSums(M) = 1
template<class Type>
matrix<Type> make_M( int CTMC_version,
                     int n_g,
                     Type DeltaD,
                     matrix<int> At_zz,
                     Type ln_D,
                     vector<Type> h_g,
                     vector<Type> colsumA_g ){
  
  int n_z = At_zz.rows();
  Type D = exp( ln_D );
  Eigen::SparseMatrix<Type> M_gg( n_g, n_g );
  if( CTMC_version==0 ){
    // Diffusion .. equal rate by cell, spread equally among neighbors
    for(int z=0; z<n_z; z++){
      M_gg.coeffRef( At_zz(z,0), At_zz(z,1) ) += D / colsumA_g( At_zz(z,0) ) / pow(DeltaD,2);
      M_gg.coeffRef( At_zz(z,0), At_zz(z,0) ) -= D / colsumA_g( At_zz(z,0) ) / pow(DeltaD,2);
    }
    // Taxis
    for(int z=0; z<n_z; z++){
      M_gg.coeffRef( At_zz(z,0), At_zz(z,1) ) += (h_g(At_zz(z,1)) - h_g(At_zz(z,0))) / DeltaD;
      M_gg.coeffRef( At_zz(z,0), At_zz(z,0) ) -= (h_g(At_zz(z,1)) - h_g(At_zz(z,0))) / DeltaD;
    }
  }else{
    // Combined taxis and diffusion
    for(int z=0; z<n_z; z++){
      M_gg.coeffRef( At_zz(z,0), At_zz(z,1) ) += D / pow(DeltaD,2) * exp( (h_g(At_zz(z,1)) - h_g(At_zz(z,0))) / DeltaD );
      M_gg.coeffRef( At_zz(z,0), At_zz(z,0) ) -= D / pow(DeltaD,2) * exp( (h_g(At_zz(z,1)) - h_g(At_zz(z,0))) / DeltaD );
    }
  }
  return M_gg;
}