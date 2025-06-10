#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;
using namespace arma;

// --- OPTIMIZATION 1 & 2: Use solver and combine matrices ---
// Replaces the original 'compute_residuals_cpp'
// [[Rcpp::export]]
List compute_residuals_optim(const arma::mat& X,  // (n x k) Confound matrix
                             const arma::mat& Y,  // (n x V) Data matrix
                             const arma::mat& C) {  // (n x T) Trial matrix
    
    // Combine Y and C to perform calculations only once
    arma::mat Z = arma::join_horiz(Y, C);

    // Use arma::solve for superior speed and numerical stability over inv_sympd
    // This solves the normal equation: (X'X)B = X'Z for B
    arma::mat B = arma::solve(X, Z, arma::solve_opts::fast);
    
    // Calculate residuals in one operation
    arma::mat Z_res = Z - X * B;

    // Split the residualized matrix back into Y and C components
    // .head_cols() and .tail_cols() are efficient subview operations (no data copy)
    arma::mat Y_res = Z_res.head_cols(Y.n_cols);
    arma::mat C_res = Z_res.tail_cols(C.n_cols);

    return List::create(Named("Q_dmat_ran") = C_res,
                        Named("residual_data") = Y_res);
}


// --- OPTIMIZATION 3, 4, & 5: No repmat, vectorized, and parallelized ---
// Replaces the original 'lss_compute_cpp'
// [[Rcpp::export]]
arma::mat lss_compute_optim(const arma::mat& C,   // projected (n x T)
                            const arma::mat& Y) { // projected (n x V)

    if (C.n_cols == 1) {
        double cc = dot(C, C);
        if (cc == 0) return arma::mat(1, Y.n_cols, arma::fill::zeros);
        return (C.t() * Y) / cc;
    }

    const uword T = C.n_cols;
    const uword V = Y.n_cols;

    // These calculations are efficient and remain the same
    arma::vec  total   = sum(C, 1);
    double     ss_tot  = dot(total, total);
    arma::mat  CtY     = C.t() * Y;
    arma::vec  CtC     = sum(square(C), 0).t();
    arma::vec  CtT     = C.t() * total;
    arma::rowvec totalY = total.t() * Y;

    // --- Optimization Area ---
    // Calculate intermediate vectors once
    arma::vec bt2 = ss_tot - 2 * CtT + CtC;
    arma::vec ctb = CtT - CtC;
    
    // Add a safeguard for division by zero if a trial's "other" regressor has no variance
    bt2.elem(find(bt2 < 1e-12)).fill(arma::datum::inf);
    
    arma::vec ctb_bt2 = ctb / bt2;
    arma::vec den_vec = CtC - square(ctb) / bt2;
    den_vec.elem(find(den_vec < 1e-12)).fill(arma::datum::inf);

    // OPTIMIZATION 3: Avoid repmat for BtY
    // Create BtY by broadcasting totalY over the rows of -CtY
    arma::mat BtY = -CtY;
    BtY.each_row() += totalY;

    // OPTIMIZATION 4: Vectorize the numerator calculation
    // Multiply each column of BtY by the ctb_bt2 vector and subtract from CtY
    arma::mat num = CtY - BtY.each_col() % ctb_bt2;
    
    // OPTIMIZATION 5: Use broadcasting for division and parallelize with OpenMP
    // The final result matrix
    arma::mat beta_estimates(T, V);

    // This loop is "embarrassingly parallel" as each voxel is independent.
    // OpenMP will distribute the columns (voxels) across available CPU cores.
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (uword i = 0; i < V; ++i) {
        beta_estimates.col(i) = num.col(i) / den_vec;
    }

    return beta_estimates;
}
