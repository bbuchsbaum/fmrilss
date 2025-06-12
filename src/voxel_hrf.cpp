#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Placeholder C++ function for voxel-wise HRF estimation
// This will be expanded in later sprints.

// [[Rcpp::export]]
arma::mat estimate_hrf_cpp(const arma::mat& X, const arma::mat& Y) {
  arma::mat out(X.n_cols, Y.n_cols, arma::fill::zeros);
  return out;
}
