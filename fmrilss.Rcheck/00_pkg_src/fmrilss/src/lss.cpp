#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// PATCH A: Projection without Q matrix
// [[Rcpp::export]]
List compute_residuals_cpp(const arma::mat& X,          // (n×k)
                           const arma::mat& Y,          // (n×V)
                           const arma::mat& C) {        // (n×T)

    arma::mat XtX   = X.t() * X;                        // k×k
    arma::mat XtXinv= inv_sympd(XtX);                   // k×k  (Cholesky, no SVD)
    arma::mat XtY   = X.t() * Y;                        // k×V
    arma::mat XtC   = X.t() * C;                        // k×T

    arma::mat Y_res = Y - X * (XtXinv * XtY);           // n×V
    arma::mat C_res = C - X * (XtXinv * XtC);           // n×T

    return List::create(Named("Q_dmat_ran") = C_res,
                        Named("residual_data") = Y_res);
}

// PATCH B: Single-pass LS-S solver (vectorised)
// [[Rcpp::export]]
arma::mat lss_compute_cpp(const arma::mat& C,   // projected (n×T)
                          const arma::mat& Y) { // projected (n×V)

    // 1. row-sums of C  (length n, held as rowvec 1×n)
    arma::rowvec sumC = sum(C, 1);                        // 1×n (rowvec)
    double ss         = dot(sumC, sumC);                  // scalar

    // 2. Cᵀ Y and Cᵀ sumC
    arma::mat CtY = C.t() * Y;                            // T×V
    arma::vec Cs  = C.t() * sumC.t();                     // T×1   <-- FIX ❶
    arma::rowvec C2 = sum(square(C), 0);                  // 1×T (column sums)
    
    // 3. sumC Y (1×V)
    arma::rowvec sumCY = sumC * Y;                        // 1×V   <-- FIX ❷
    arma::mat num = CtY * ss - Cs * sumCY;                // T×V

    arma::mat den = repmat(C2.t() * ss, 1, Y.n_cols) - repmat(square(Cs), 1, Y.n_cols);

    den.elem(find(abs(den) < 1e-12)).fill(datum::nan);   // avoid 0-division
    return num / den;
}

// User-friendly exported names
// [[Rcpp::export]]
List project_confounds_cpp(const arma::mat& X_confounds,
                           const arma::mat& Y_data,
                           const arma::mat& C_trials) {
    return compute_residuals_cpp(X_confounds, Y_data, C_trials);
}

// [[Rcpp::export]]
arma::mat lss_beta_cpp(const arma::mat& C_projected,
                       const arma::mat& Y_projected) {
    return lss_compute_cpp(C_projected, Y_projected);
}
