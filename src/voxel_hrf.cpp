#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory)]]
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

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
    // Convert fine-grid onset index to TR index (TR assumed to align with grid sampling)
    double start_tr_d = static_cast<double>(onset_idx[t]) / static_cast<double>(step);
    int start_tr = static_cast<int>(std::round(start_tr_d)) - 1;
    if (start_tr < 0) start_tr = 0;
    if (start_tr >= static_cast<int>(n_time)) continue;

    // Duration in TRs (inclusive semantics: length d -> indices [0..d])
    unsigned int d_tr = 0;
    if (t < durations.n_elem) {
      double dv = durations[t];
      if (dv < 0.0) dv = 0.0;
      d_tr = static_cast<unsigned int>(std::round(dv));
    }

    // Add HRF contributions for each TR within the duration window
    const unsigned int klen = static_cast<unsigned int>(hrf_kernel.n_elem);
    for (unsigned int k = 0; k <= d_tr; ++k) {
      unsigned int s = static_cast<unsigned int>(start_tr + static_cast<int>(k));
      if (s >= n_time) break;
      unsigned int max_len = std::min<unsigned int>(klen, n_time - s);
      if (max_len == 0U) break;
      C.submat(s, t, s + max_len - 1U, t) += hrf_kernel.subvec(0, max_len - 1U);
    }
  }
  return C;
}

// [[Rcpp::export]]
void lss_engine_vox_hrf(const arma::mat& Y,
                        const arma::mat& coeffs,
                        const arma::mat& basis_kernels,
                        const arma::uvec& onset_idx,
                        const arma::vec& durations,
                        const arma::mat& nuisance,
                        SEXP betas_ptr,
                        Rcpp::Function progress,
                        const int chunk_size,
                        bool verbose) {
  const unsigned int n_time = Y.n_rows;
  const unsigned int V = Y.n_cols;
  const unsigned int step = 10;              // 1 / 0.1 (fine_dt)
  const arma::uvec tr_idx = arma::regspace<arma::uvec>(0, step, basis_kernels.n_rows - 1);
  const unsigned int n_trials = onset_idx.n_elem;

  Rcpp::XPtr<BigMatrix> bm_ptr(betas_ptr);
  MatrixAccessor<double> bm_acc(*bm_ptr);

  for (unsigned int start = 0; start < V; start += chunk_size) {
    unsigned int end = std::min(start + (unsigned int)chunk_size, V) - 1;
    unsigned int chunkV = end - start + 1;
    arma::mat chunk_out(n_trials, chunkV);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (unsigned int idx = 0; idx < chunkV; ++idx) {
      unsigned int v = start + idx;

      // Build voxel-specific HRF at fine grid and downsample to TRs
      arma::vec coeff_v = coeffs.col(v);
      for (arma::uword k = 0; k < coeff_v.n_elem; ++k) {
        if (!std::isfinite(coeff_v[k])) {
          coeff_v[k] = 0.0;
        }
      }
      arma::vec hrf_fine = basis_kernels * coeff_v;      // L_fine
      for (arma::uword k = 0; k < hrf_fine.n_elem; ++k) {
        if (!std::isfinite(hrf_fine[k])) {
          hrf_fine[k] = 0.0;
        }
      }
      // Downsample by selecting the fine-sample aligned to each TR (matches fmrihrf conv on TR grid)
      arma::vec hrf_tr = hrf_fine.elem(tr_idx);
      for (arma::uword k = 0; k < hrf_tr.n_elem; ++k) {
        if (!std::isfinite(hrf_tr[k])) hrf_tr[k] = 0.0;
      }
      arma::mat C        = build_trial_matrix(hrf_tr, onset_idx, durations, n_time, step); // n_time x n_trials

      // FWL residualization by nuisance: project out nuisance from Y and C
      arma::vec y = Y.col(v);
      arma::mat C_res;
      arma::vec y_res;
      if (nuisance.n_cols == 0) {
        C_res = C;
        y_res = y;
      } else {
        arma::mat Q, R;
        arma::qr_econ(Q, R, nuisance);
        y_res = y - Q * (Q.t() * y);
        C_res = C - Q * (Q.t() * C);
      }

      // Closed-form LSS solution with fallback on degenerate trials
      const unsigned int T = C_res.n_cols;
      arma::mat y_res_mat(C_res.n_rows, 1);
      y_res_mat.col(0) = y_res;
      arma::mat beta_mat = lss_beta_cpp(C_res, y_res_mat);
      arma::vec beta = beta_mat.col(0);

      if (!beta.is_finite()) {
        arma::vec total = arma::sum(C_res, 1);
        double ss_tot = arma::dot(total, total);
        arma::vec CtY = C_res.t() * y_res;
        arma::vec CtC = arma::sum(arma::square(C_res), 0).t();
        arma::vec CtT = C_res.t() * total;
        double totalY = arma::dot(total, y_res);
        arma::vec BtY = totalY - CtY;
        arma::vec bt2 = ss_tot - 2.0 * CtT + CtC;
        const double eps = 1e-10;

        for (arma::uword j = 0; j < T; ++j) {
          if (std::isfinite(beta[j])) continue;
          double a = CtC[j];
          double n1 = CtY[j];
          double b = bt2[j];
          double c = CtT[j] - CtC[j];
          double n2 = BtY[j];
          double det = a * b - c * c;
          double val = 0.0;

          if (std::isfinite(det) && std::fabs(det) >= eps &&
              std::isfinite(a) && std::isfinite(b) && std::isfinite(c) &&
              std::isfinite(n1) && std::isfinite(n2)) {
            double numer = n1 * b - c * n2;
            val = numer / det;
            if (!std::isfinite(val)) {
              val = 0.0;
            }
          } else if (std::isfinite(a) && std::fabs(a) >= eps && std::isfinite(n1)) {
            val = n1 / a;
            if (!std::isfinite(val)) {
              val = 0.0;
            }
          }

          beta[j] = val;
        }
      }

      chunk_out.col(idx) = beta;
    }

    for (unsigned int idx = 0; idx < chunkV; ++idx) {
      unsigned int v = start + idx;
      for (unsigned int t = 0; t < n_trials; ++t) {
        bm_acc[v][t] = chunk_out(t, idx);
      }
    }

    if (verbose) progress(chunkV);
  }

  return;
}
