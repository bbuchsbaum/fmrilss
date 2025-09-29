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

    // Handle empty or zero-column nuisance robustly
    if (X.n_cols == 0) {
        return List::create(Named("Q_dmat_ran") = C,
                            Named("residual_data") = Y);
    }

    arma::mat XtX = X.t() * X;                          // k×k
    arma::mat XtXinv;
    bool spd = false;
    try {
        spd = arma::inv_sympd(XtXinv, XtX);             // Cholesky if SPD
    } catch (...) {
        spd = false;
    }
    if (!spd) {
        // Fallback to pseudo-inverse for singular/non-SPD cases
        XtXinv = arma::pinv(XtX);
    }

    arma::mat XtY = X.t() * Y;                          // k×V
    arma::mat XtC = X.t() * C;                          // k×T

    arma::mat Y_res = Y - X * (XtXinv * XtY);           // n×V
    arma::mat C_res = C - X * (XtXinv * XtC);           // n×T

    return List::create(Named("Q_dmat_ran") = C_res,
                        Named("residual_data") = Y_res);
}

// PATCH B: Single-pass LS-S solver (vectorised)
// [[Rcpp::export]]
arma::mat lss_compute_cpp(const arma::mat& C,   // projected (n×T)
                          const arma::mat& Y) { // projected (n×V)

    if (C.n_cols == 1) {                // single-trial guard
        double cc = dot(C, C);
        return (C.t() * Y) / cc;
    }

    arma::vec  total   = sum(C, 1);               // n
    double     ss_tot  = dot(total, total);

    arma::mat  CtY   = C.t() * Y;                 // T×V
    arma::vec  CtC   = sum(square(C), 0).t();     // T
    arma::vec  CtT   = C.t() * total;             // T
    arma::rowvec totalY = total.t() * Y;          // 1×V

    arma::mat  BtY   = repmat(totalY, C.n_cols, 1) - CtY; // T×V
    arma::vec  bt2   = ss_tot - 2*CtT + CtC;
    arma::vec  ctb   = CtT   - CtC;

    arma::vec ctb_bt2 = ctb / bt2;                // T vector
    
    // Memory-efficient numerator calculation (Fix 2, final)
    arma::mat num = CtY;
    for(uword i = 0; i < num.n_cols; ++i) {
        num.col(i) -= BtY.col(i) % ctb_bt2;
    }

    arma::mat den = repmat(CtC - square(ctb)/bt2, 1, Y.n_cols);

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
