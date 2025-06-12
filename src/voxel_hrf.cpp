#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat estimate_hrf_cpp(const arma::mat& X, const arma::mat& Y) {
  // Solve the multi-response linear model X * B = Y for B.
  // This uses Armadillo's solve which will automatically
  // choose the appropriate least squares solver when X is
  // not square.
  return arma::solve(X, Y);
}
