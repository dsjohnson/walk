// 
// // 
// #include <RcppEigen.h>
// 
// // [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::export]]
// Eigen::VectorXd ctmc_limiting_distribution(const Eigen::SparseMatrix<double> &Q) {
//   int n = Q.rows();
//   Eigen::SparseMatrix<double> QT = Q.transpose();  // work with left nullspace
//   Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
// 
//   // Replace last row with ones
//   for (int j = 0; j < n; ++j) {
//     QT.coeffRef(n - 1, j) = 1.0;
//   }
//   b[n - 1] = 1.0;
// 
//   // Solve QT * ud = b
//   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//   solver.compute(QT);
//   if (solver.info() != Eigen::Success) {
//     Rcpp::stop("Decomposition failed");
//   }
// 
//   Eigen::VectorXd pi = solver.solve(b);
//   if (solver.info() != Eigen::Success) {
//     Rcpp::stop("Solving failed");
//   }
//   
//   // Normalize
//   pi /= pi.sum();
// 
//   return pi;
// }
// 
// 
// // [[Rcpp::export]]
// Eigen::VectorXd bicgstab_stationary(const Eigen::SparseMatrix<double> &Q, 
//                                     int max_iter = 1000, 
//                                     double tol = 1e-10) {
//   using namespace Eigen;
//   int n = Q.rows();
//   
//   // Compute transpose of Q
//   SparseMatrix<double> Qt = Q.transpose();
//   
//   // Replace last row of Qt with ones to enforce sum(pi) = 1 constraint
//   std::vector<Triplet<double>> triplets;
//   for (int k = 0; k < Qt.outerSize(); ++k) {
//     for (SparseMatrix<double>::InnerIterator it(Qt, k); it; ++it) {
//       if (it.row() != n - 1) {
//         triplets.emplace_back(it.row(), it.col(), it.value());
//       }
//     }
//   }
//   for (int j = 0; j < n; ++j) {
//     triplets.emplace_back(n - 1, j, 1.0);
//   }
//   SparseMatrix<double> Qt_mod(n, n);
//   Qt_mod.setFromTriplets(triplets.begin(), triplets.end());
//   
//   // RHS vector enforcing sum(pi) = 1
//   VectorXd b = VectorXd::Zero(n);
//   b[n - 1] = 1.0;
//   
//   // Use BiCGSTAB solver from Eigen
//   BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> solver;
//   solver.setMaxIterations(max_iter);
//   solver.setTolerance(tol);
//   solver.compute(Qt_mod);
//   if (solver.info() != Success) {
//     Rcpp::stop("Failed to decompose matrix with BiCGSTAB");
//   }
//   
//   VectorXd pi = solver.solve(b);
//   if (solver.info() != Success) {
//     Rcpp::stop("BiCGSTAB solver failed to converge");
//   }
//   
//   // Normalize in case of numerical drift
//   pi /= pi.sum();
//   
//   return pi;
// }
// 
