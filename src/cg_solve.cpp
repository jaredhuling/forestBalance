#include "forestBalance_types.h"

// Conjugate gradient solver for Z Z^T x = rhs (matrix-free).
// Uses Eigen sparse mat-vecs with the C++ CG loop to avoid R overhead.
// [[Rcpp::export]]
Rcpp::NumericVector cg_solve_cpp(Rcpp::S4 Z_s4, Rcpp::NumericVector rhs_r,
                                  double tol, int maxiter) {
  const Eigen::MappedSparseMatrix<double> Z(Rcpp::as<Eigen::MappedSparseMatrix<double> >(Z_s4));
  const Eigen::Map<Eigen::VectorXd> rhs(rhs_r.begin(), rhs_r.size());

  int n = Z.rows();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd r = rhs;
  Eigen::VectorXd p = r;
  double rs = r.squaredNorm();

  for (int i = 0; i < maxiter; i++) {
    // Matrix-free mat-vec: A*p = Z * (Z^T * p)
    Eigen::VectorXd Ztp = Z.transpose() * p;
    Eigen::VectorXd Ap = Z * Ztp;

    double pAp = p.dot(Ap);
    if (pAp <= 0) break;  // numerical safeguard

    double alpha = rs / pAp;
    x += alpha * p;
    r -= alpha * Ap;

    double rsnew = r.squaredNorm();
    if (std::sqrt(rsnew) < tol) break;

    double beta = rsnew / rs;
    p = r + beta * p;
    rs = rsnew;
  }

  return Rcpp::wrap(x);
}
