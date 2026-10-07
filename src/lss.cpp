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

arma::mat lss_weight_matrix_cpp(const arma::mat& C, const arma::ivec& groups,
                                double eps, double ridge_x, double ridge_b);

// Single-pass LS-S solver: beta = W' Y with W the LSS weight matrix of the
// projected trial regressors. Y may be raw or projected data; the result is
// identical because W lies in the confound residual space.
// [[Rcpp::export]]
arma::mat lss_compute_cpp(const arma::mat& C,   // projected (n×T)
                          const arma::mat& Y) { // data (n×V)
    arma::ivec no_groups;
    arma::mat W = lss_weight_matrix_cpp(C, no_groups, 1e-12, 0.0, 0.0);
    return W.t() * Y;
}

// C++ aliases retained for internal native callers; the public R wrappers call
// compute_residuals_cpp() and lss_compute_cpp() directly.
List project_confounds_cpp(const arma::mat& X_confounds,
                           const arma::mat& Y_data,
                           const arma::mat& C_trials) {
    return compute_residuals_cpp(X_confounds, Y_data, C_trials);
}

arma::mat lss_beta_cpp(const arma::mat& C_projected,
                       const arma::mat& Y_projected) {
    return lss_compute_cpp(C_projected, Y_projected);
}

// Non-allocating finiteness check for numeric (double) or integer arrays.
// [[Rcpp::export]]
bool all_finite_cpp(SEXP x) {
    if (TYPEOF(x) == REALSXP) {
        const double* p = REAL(x);
        const R_xlen_t n = XLENGTH(x);
        // Accumulating x * 0 stays 0 for finite inputs and becomes NaN on
        // any NA, NaN or +/-Inf, which keeps the loop branch-free.
        double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        R_xlen_t i = 0;
        for (; i + 4 <= n; i += 4) {
            a0 += p[i] * 0.0;
            a1 += p[i + 1] * 0.0;
            a2 += p[i + 2] * 0.0;
            a3 += p[i + 3] * 0.0;
        }
        for (; i < n; ++i) a0 += p[i] * 0.0;
        return (a0 + a1 + a2 + a3) == 0.0;
    }
    if (TYPEOF(x) == INTSXP || TYPEOF(x) == LGLSXP) {
        const int* p = INTEGER(x);
        const R_xlen_t n = XLENGTH(x);
        for (R_xlen_t i = 0; i < n; ++i) if (p[i] == NA_INTEGER) return false;
        return true;
    }
    Rcpp::stop("all_finite_cpp expects a numeric array");
}

// Per-voxel autocorrelations at lags 1..max_lag of a residual matrix, using
// the biased (PSD) autocovariance estimator. The mean is removed per run and
// lag products never span a run boundary. `run_starts` are 0-based.
// Returns a max_lag x V matrix.
// [[Rcpp::export]]
arma::mat voxel_acf_cpp(const arma::mat& E, const arma::uvec& run_starts,
                        int max_lag) {
    const uword n = E.n_rows;
    const uword V = E.n_cols;
    const uword L = static_cast<uword>(std::max(max_lag, 1));
    arma::uvec bounds = arma::join_cols(run_starts, arma::uvec({n}));
    arma::mat out(L, V, arma::fill::zeros);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (uword v = 0; v < V; ++v) {
        arma::vec gam(L + 1, arma::fill::zeros);
        for (uword r = 0; r + 1 < bounds.n_elem; ++r) {
            const uword a = bounds[r], b = bounds[r + 1];
            if (b <= a) continue;
            arma::vec e = E.col(v).subvec(a, b - 1);
            e -= arma::mean(e);
            const uword m = e.n_elem;
            for (uword k = 0; k <= L && k < m; ++k) {
                gam[k] += arma::dot(e.head(m - k), e.tail(m - k));
            }
        }
        if (gam[0] > 0) out.col(v) = gam.tail(L) / gam[0];
    }
    return out;
}
