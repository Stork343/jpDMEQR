#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>

extern "C" SEXP jpDMEQR_population_crossprod(SEXP xSEXP, SEXP ySEXP) {
  BEGIN_RCPP
  const Rcpp::NumericMatrix x(xSEXP);
  const Eigen::Map<const Eigen::MatrixXd> left(x.begin(), x.nrow(), x.ncol());
  Eigen::MatrixXd product;
  if (Rf_isNull(ySEXP)) {
    product.noalias() = left.transpose() * left;
  } else {
    const Rcpp::NumericMatrix rhs(ySEXP);
    if (rhs.nrow() != x.nrow()) Rcpp::stop("Crossproduct row dimensions differ.");
    const Eigen::Map<const Eigen::MatrixXd> right(rhs.begin(), rhs.nrow(), rhs.ncol());
    product.noalias() = left.transpose() * right;
  }
  return Rcpp::wrap(product);
  END_RCPP
}
