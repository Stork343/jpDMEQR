#include <RcppEigen.h>
#include <cmath>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd sb_cluster_curvature(const MatrixXd& Z, const VectorXd& phi,
                             const MatrixXd& Lambda) {
  const int m = Z.rows();
  VectorXd weight = phi / m;
  MatrixXd WZ = weight.asDiagonal() * Z;
  Eigen::LLT<MatrixXd> factor(Z.transpose() * WZ + Lambda);
  if (factor.info() != Eigen::Success) Rcpp::stop("Invalid nuisance curvature.");
  MatrixXd L = weight.asDiagonal().toDenseMatrix();
  L.noalias() -= WZ * factor.solve(WZ.transpose());
  return (L + L.transpose()) / 2;
}

// [[Rcpp::export]]
Eigen::MatrixXd sb_geometry_cpp(
    const Eigen::Map<MatrixXd> X, const Eigen::Map<MatrixXd> Z,
    const Eigen::Map<VectorXd> phi, Rcpp::IntegerVector order,
    Rcpp::IntegerVector pointers, const Eigen::Map<MatrixXd> Lambda,
    bool return_factor = false, Rcpp::Nullable<Rcpp::NumericVector> cluster_weights = R_NilValue) {
  const int N = X.rows(), p = X.cols(), q = Z.cols(), n = pointers.size() - 1;
  if (N < 1 || p < 1 || q < 1 || n < 1 || Z.rows() != N || phi.size() != N ||
      order.size() != N || pointers[0] != 0 || pointers[n] != N ||
      Lambda.rows() != q || Lambda.cols() != q || !X.allFinite() ||
      !Z.allFinite() || !phi.allFinite() || phi.minCoeff() < 0) {
    Rcpp::stop("Invalid profile geometry input.");
  }
  MatrixXd output = MatrixXd::Zero(return_factor ? N : p, p);
  Rcpp::NumericVector weights = cluster_weights.isNotNull() ?
    Rcpp::NumericVector(cluster_weights) : Rcpp::NumericVector(n, 1.0);
  if (weights.size() != n) Rcpp::stop("Invalid cluster weights.");
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(weights[i]) || weights[i] < 0) Rcpp::stop("Invalid cluster weight.");
    const int m = pointers[i + 1] - pointers[i];
    if (m < 1) Rcpp::stop("Empty geometry cluster.");
    MatrixXd Xi(m, p), Zi(m, q);
    VectorXd pi(m);
    for (int j = 0; j < m; ++j) {
      const int r = order[pointers[i] + j];
      if (r < 0 || r >= N) Rcpp::stop("Invalid geometry row.");
      Xi.row(j) = X.row(r);
      Zi.row(j) = Z.row(r);
      pi[j] = phi[r];
    }
    const MatrixXd L = sb_cluster_curvature(Zi, pi, Lambda);
    if (return_factor) {
      Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(L);
      if (eigen.info() != Eigen::Success ||
          eigen.eigenvalues().minCoeff() < -1e-10 * std::max(1.0, L.norm())) {
        Rcpp::stop("Nonpositive profile factor.");
      }
      MatrixXd F = eigen.eigenvectors() *
        eigen.eigenvalues().cwiseMax(0.0).cwiseSqrt().asDiagonal() *
        eigen.eigenvectors().transpose() * Xi;
      for (int j = 0; j < m; ++j) output.row(order[pointers[i] + j]) = std::sqrt(weights[i]) * F.row(j);
    } else {
      output.noalias() += weights[i] * Xi.transpose() * L * Xi / n;
    }
  }
  return output;
}
