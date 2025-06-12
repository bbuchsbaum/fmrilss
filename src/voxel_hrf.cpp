#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif
// [[Rcpp::depends(RcppArmadillo)]]

// forward declarations from lss.cpp
Rcpp::List compute_residuals_cpp(const arma::mat& X, const arma::mat& Y,
                                 const arma::mat& C);
arma::mat lss_beta_cpp(const arma::mat& C_projected, const arma::mat& Y_projected);

// [[Rcpp::export]]
arma::mat estimate_hrf_cpp(const arma::mat& X, const arma::mat& Y) {
  // Solve the multi-response linear model X * B = Y for B.
  // This uses Armadillo's solve which will automatically
  // choose the appropriate least squares solver when X is
  // not square.
  return arma::solve(X, Y);
}

// helper to build trial design matrix for a voxel
static arma::mat build_trial_matrix(const arma::vec& hrf_kernel,
                                    const arma::uvec& onset_idx,
                                    const arma::vec& durations,
                                    unsigned int n_time,
                                    unsigned int step) {
  const unsigned int n_trials = onset_idx.n_elem;
  arma::mat C(n_time, n_trials, arma::fill::zeros);

  for (unsigned int t = 0; t < n_trials; ++t) {
    unsigned int start_tr = onset_idx[t] / step; // convert fine index to TR index
    if (start_tr >= n_time) continue;

    unsigned int len = hrf_kernel.n_elem;
    unsigned int end_tr = std::min(start_tr + len - 1, n_time - 1);
    C.submat(start_tr, t, end_tr, t) += hrf_kernel.subvec(0, end_tr - start_tr);
  }
  return C;
}

// [[Rcpp::export]]
arma::mat lss_engine_vox_hrf(const arma::mat& Y,
                             const arma::mat& coeffs,
                             const arma::mat& basis_kernels,
                             const arma::uvec& onset_idx,
                             const arma::vec& durations,
                             const arma::mat& nuisance,
                             const int chunk_size,
                             bool verbose) {
  const unsigned int n_time = Y.n_rows;
  const unsigned int V = Y.n_cols;
  const unsigned int step = 10;              // 1 / 0.1 (fine_dt)
  const arma::uvec tr_idx = arma::regspace<arma::uvec>(0, step, basis_kernels.n_rows - 1);
  const unsigned int n_trials = onset_idx.n_elem;

  arma::mat betas(n_trials, V, arma::fill::zeros);

  for (unsigned int start = 0; start < V; start += chunk_size) {
    unsigned int end = std::min(start + (unsigned int)chunk_size, V) - 1;
    unsigned int chunkV = end - start + 1;
    arma::mat chunk_out(n_trials, chunkV);

    #pragma omp parallel
    {
      arma::vec beta_local(n_trials);
      #pragma omp for schedule(static)
      for (unsigned int idx = 0; idx < chunkV; ++idx) {
        unsigned int v = start + idx;
        arma::vec hrf_fine = basis_kernels * coeffs.col(v);
        arma::vec hrf_tr = hrf_fine.elem(tr_idx);
        arma::mat C = build_trial_matrix(hrf_tr, onset_idx, durations, n_time, step);
        Rcpp::List res = compute_residuals_cpp(nuisance, Y.col(v), C);
        arma::mat beta_v = lss_beta_cpp(res["Q_dmat_ran"], res["residual_data"]);
        beta_local = beta_v.col(0);
        #pragma omp critical
        chunk_out.col(idx) = beta_local;
      }
    }
    betas.cols(start, end) = chunk_out;
  }

  return betas;
}

