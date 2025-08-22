#include <Rcpp.h>
using namespace Rcpp;

// Helper: convert List of matrices to std::vector<NumericMatrix>
static std::vector<NumericMatrix> as_matrix_list(const List& L, int expected_rows, int expected_cols) {
  const int K = L.size();
  std::vector<NumericMatrix> out;
  out.reserve(K);
  for (int k = 0; k < K; ++k) {
    NumericMatrix Mk = as<NumericMatrix>(L[k]);
    if (Mk.nrow() != expected_rows || Mk.ncol() != expected_cols) {
      stop("basis_convolved[[%d]] has inconsistent dimensions: expected %d x %d, got %d x %d",
           k + 1, expected_rows, expected_cols, Mk.nrow(), Mk.ncol());
    }
    out.push_back(Mk);
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix lss_engine_vox_hrf_cpp(
    const NumericMatrix& Y,              // n_time x n_vox
    const NumericMatrix& coeffs,         // K x n_vox
    const List& basis_convolved,         // length K; each n_time x n_trials
    const NumericMatrix& Z               // n_time x pz (>=1; intercept allowed)
) {
  const int n_time   = Y.nrow();
  const int n_vox    = Y.ncol();
  const int K        = coeffs.nrow();
  const int n_trials = as<NumericMatrix>(basis_convolved[0]).ncol();
  const int pz       = Z.ncol();

  // Validate list and sizes early
  std::vector<NumericMatrix> Dk = as_matrix_list(basis_convolved, n_time, n_trials);
  if (coeffs.ncol() != n_vox) stop("coeffs must be [K x n_vox] with n_vox matching Y");
  if (pz < 1) stop("Z must have at least one column (e.g., intercept).");

  NumericMatrix betas(n_trials, n_vox);

  // Prepare access to stats::qr and stats::qr.coef once
  Function qr_fn("qr");
  Function qr_coef("qr.coef");

  // Pre-allocate workspaces reused inside loops to reduce allocations
  NumericMatrix Xv(n_time, n_trials);
  std::vector<double> xall(n_time);

  for (int v = 0; v < n_vox; ++v) {
    // Combine per-basis convolved designs using voxel weights
    std::fill(Xv.begin(), Xv.end(), 0.0);
    for (int k = 0; k < K; ++k) {
      const double w = coeffs(k, v);
      if (w == 0.0) continue;
      NumericMatrix& D = Dk[k];
      for (int t = 0; t < n_time; ++t) {
        for (int i = 0; i < n_trials; ++i) {
          Xv(t, i) += w * D(t, i);
        }
      }
    }

    // Precompute "all others" sum per timepoint
    for (int t = 0; t < n_time; ++t) {
      double s = 0.0;
      for (int i = 0; i < n_trials; ++i) s += Xv(t, i);
      xall[t] = s;
    }

    NumericVector yv = Y(_, v);

    // Per-trial GLMs
    for (int i = 0; i < n_trials; ++i) {
      // Build Xdesign = [Z, Xi, (Xother if n_trials>1)]
      const bool have_other = (n_trials > 1);
      const int p = pz + 1 + (have_other ? 1 : 0);
      NumericMatrix Xdesign(n_time, p);

      // Copy Z
      for (int c = 0; c < pz; ++c) {
        for (int t = 0; t < n_time; ++t) {
          Xdesign(t, c) = Z(t, c);
        }
      }
      // Xi
      for (int t = 0; t < n_time; ++t) {
        Xdesign(t, pz) = Xv(t, i);
      }
      // Xother
      if (have_other) {
        for (int t = 0; t < n_time; ++t) {
          Xdesign(t, pz + 1) = xall[t] - Xv(t, i);
        }
      }

      // Solve via QR using R's LAPACK through stats::qr + stats::qr.coef
      SEXP qrobj = qr_fn(Xdesign);
      NumericVector coef = qr_coef(qrobj, yv);

      // Take coefficient of Xi, which is at index pz (0-based) in coef
      if (coef.size() > pz && R_finite(coef[pz])) {
        betas(i, v) = coef[pz];
      } else {
        betas(i, v) = NA_REAL;
      }
    }
  }

  return betas;
}