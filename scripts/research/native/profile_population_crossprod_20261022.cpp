#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::NumericMatrix profile_population_crossprod_cpp(Rcpp::NumericMatrix x,
                                                    Rcpp::Nullable<Rcpp::NumericMatrix> y = R_NilValue) {
  const Eigen::Map<const Eigen::MatrixXd> left(x.begin(), x.nrow(), x.ncol());
  Eigen::MatrixXd product;
  if (y.isNull()) {
    product.noalias() = left.transpose() * left;
  } else {
    const Rcpp::NumericMatrix rhs(y.get());
    if (rhs.nrow() != x.nrow()) Rcpp::stop("Crossproduct row dimensions differ.");
    const Eigen::Map<const Eigen::MatrixXd> right(rhs.begin(), rhs.nrow(), rhs.ncol());
    product.noalias() = left.transpose() * right;
  }
  return Rcpp::wrap(product);
}
