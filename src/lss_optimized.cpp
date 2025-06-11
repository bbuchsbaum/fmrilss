#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;
using namespace arma;

//' Fused Single-Pass LSS Solver (C++)
//'
//' This function computes Least Squares-Separate (LSS) beta estimates using
//' a memory-efficient, single-pass algorithm. It fuses the projection and
//' estimation steps, processing voxels in parallel blocks to maximize cache
//' efficiency.
//'
//' @param X The nuisance regressor matrix (confounds).
//' @param Y The data matrix (e.g., fMRI data).
//' @param C The trial-wise design matrix.
//' @param block_size The number of voxels to process in each parallel block.
//' @return A matrix of LSS beta estimates.
//' @keywords internal
// [[Rcpp::export]]
arma::mat lss_fused_optim_cpp(const arma::mat& X,          // (n x k) Nuisance regressors
                              const arma::mat& Y,          // (n x V) Data matrix
                              const arma::mat& C,          // (n x T) Trial matrix
                              int block_size = 96) {       // Voxel block size

    const uword V = Y.n_cols;
    const uword T = C.n_cols;

    // --- 1. Pre-compute confound projection components ---
    // These are small matrices, computed once.
    arma::mat XtX = X.t() * X;
    arma::mat XtXinv;
    bool success = arma::inv_sympd(XtXinv, XtX);
    if (!success) {
        // Fallback to slower but more stable pseudo-inverse if sympd fails
        XtXinv = arma::pinv(XtX);
    }
    arma::mat P_confound_mult = XtXinv * X.t(); // (k x n)

    // --- 2. Project the trial matrix (C) ---
    // C_res is (n x T), also computed once.
    arma::mat C_res = C - X * (P_confound_mult * C);

    // --- 3. Pre-compute LSS components for projected C ---
    // These are small vectors/scalars, computed once.
    if (T == 1) { // Guard for single trial case
        arma::mat Y_res = Y - X * (P_confound_mult * Y);
        double cc = dot(C_res, C_res);
        if (cc == 0) return arma::mat(1, V, arma::fill::zeros);
        return (C_res.t() * Y_res) / cc;
    }

    arma::vec  total_c = sum(C_res, 1);
    double     ss_tot_c = dot(total_c, total_c);
    arma::vec  C_res_t_C_res = sum(square(C_res), 0).t();
    arma::vec  C_res_t_total_c = C_res.t() * total_c;

    arma::vec bt2 = ss_tot_c - 2 * C_res_t_total_c + C_res_t_C_res;
    arma::vec ctb = C_res_t_total_c - C_res_t_C_res;
    
    bt2.elem(find(bt2 < 1e-12)).fill(arma::datum::inf);
    
    arma::vec ctb_bt2 = ctb / bt2;
    arma::vec den_vec = C_res_t_C_res - square(ctb) / bt2;
    den_vec.elem(find(den_vec < 1e-12)).fill(arma::datum::inf);

    // --- 4. Process Y in parallel blocks ---
    arma::mat beta_estimates(T, V);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (uword j = 0; j < V; j += block_size) {
        // Define the current block of voxels
        uword j_end = std::min(j + (uword)block_size, V) - 1;
        arma::mat Y_block = Y.cols(j, j_end);

        // Project the Y block (n x block_size)
        arma::mat Y_res_block = Y_block - X * (P_confound_mult * Y_block);
        
        // --- LSS calculation for the block ---
        arma::mat C_res_t_Y_res_block = C_res.t() * Y_res_block; // (T x block_size)
        arma::rowvec total_c_t_Y_res_block = total_c.t() * Y_res_block; // (1 x block_size)

        arma::mat BtY_block = -C_res_t_Y_res_block;
        BtY_block.each_row() += total_c_t_Y_res_block;

        arma::mat num_block = C_res_t_Y_res_block - BtY_block.each_col() % ctb_bt2;
        num_block.each_col() /= den_vec;

        beta_estimates.cols(j, j_end) = std::move(num_block);
    }

    return beta_estimates;
}
