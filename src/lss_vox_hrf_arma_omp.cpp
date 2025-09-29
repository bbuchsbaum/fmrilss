#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Safe coercion to base double matrix; falls back to as.matrix for S4 Matrix, etc.
static Rcpp::NumericMatrix coerce_dense_matrix(SEXP x, int nr, int nc, const char* label) {
  if (Rf_isMatrix(x) && TYPEOF(x) == REALSXP) {
    Rcpp::NumericMatrix M(x);
    if (M.nrow() != nr || M.ncol() != nc) {
      Rcpp::stop("%s has shape %d x %d, expected %d x %d", label, M.nrow(), M.ncol(), nr, nc);
    }
    return M;
  }
  Rcpp::Function as_matrix("as.matrix");
  Rcpp::NumericMatrix M = Rcpp::as<Rcpp::NumericMatrix>(as_matrix(x));
  if (M.nrow() != nr || M.ncol() != nc) {
    Rcpp::stop("%s has shape %d x %d, expected %d x %d", label, M.nrow(), M.ncol(), nr, nc);
  }
  return M;
}

// Convert an R list of matrices into a vector of Armadillo views (no copies).
static std::vector<arma::Mat<double>> as_arma_list(const List& L, int expected_rows, int expected_cols) {
  const int K = L.size();
  std::vector<arma::Mat<double>> out;
  out.reserve(K);
  for (int k = 0; k < K; ++k) {
    Rcpp::NumericMatrix Mk = coerce_dense_matrix(L[k], expected_rows, expected_cols, "basis_convolved[[k]]");
    arma::Mat<double> Vk(Mk.begin(), Mk.nrow(), Mk.ncol(), /*copy_aux_mem=*/false, /*strict=*/true);
    out.push_back(Vk);
  }
  return out;
}

inline arma::vec ridge_normal_eq(const arma::mat& X, const arma::vec& y) {
  const arma::mat XtX = X.t() * X;
  double lam = 1e-8 * XtX.diag().max();
  arma::mat A = XtX;
  A.diag() += lam;
  arma::vec b = X.t() * y;
  // prefer symmetric positive definite solver if possible
  arma::vec beta;
  bool ok = arma::solve(beta, A, b, arma::solve_opts::likely_sympd + arma::solve_opts::equilibrate);
  if (!ok) {
    // fall back to generic solver
    beta = arma::solve(A, b);
  }
  return beta;
}

inline arma::vec qr_coef(const arma::mat& X, const arma::vec& y) {
  arma::mat Q, R;
  // economical QR: X = Q (n x p) * R (p x p), with p <= n
  arma::qr_econ(Q, R, X);
  // Solve R * beta = Q' * y
  arma::vec qty = Q.t() * y;
  arma::vec beta;
  bool ok = arma::solve(beta, arma::trimatu(R), qty, arma::solve_opts::fast);
  if (!ok) {
    // Rare: fallback to ridge-stabilized normal equations
    beta = ridge_normal_eq(X, y);
  }
  return beta;
}

// Core engine using Armadillo (sequential)
// [[Rcpp::export]]
Rcpp::NumericMatrix lss_engine_vox_hrf_arma(
    const Rcpp::NumericMatrix& Y,       // n_time x n_vox
    const Rcpp::NumericMatrix& coeffs,  // K x n_vox
    const Rcpp::List& basis_convolved,  // length K; each n_time x n_trials
    const Rcpp::NumericMatrix& Z        // n_time x pz (>=1; intercept allowed)
) {
  const int n_time   = Y.nrow();
  const int n_vox    = Y.ncol();
  const int K        = coeffs.nrow();
  const int n_trials = as<Rcpp::NumericMatrix>(basis_convolved[0]).ncol();
  const int pz       = Z.ncol();

  if (coeffs.ncol() != n_vox) stop("coeffs must be [K x n_vox] with n_vox matching Y");
  if (pz < 1) stop("Z must have at least one column (e.g., intercept).");

  // Armadillo views (no copies)
  arma::Mat<double> Ya(const_cast<double*>(Y.begin()), n_time, n_vox, false, true);
  arma::Mat<double> Ca(const_cast<double*>(coeffs.begin()), K, n_vox, false, true);
  arma::Mat<double> Za(const_cast<double*>(Z.begin()), n_time, pz, false, true);
  std::vector<arma::Mat<double>> Dk = as_arma_list(basis_convolved, n_time, n_trials);

  Rcpp::NumericMatrix betas(n_trials, n_vox);
  arma::Mat<double> Betas(betas.begin(), n_trials, n_vox, false, true);

  for (int v = 0; v < n_vox; ++v) {
    // Combine per-basis convolved designs using voxel weights
    arma::mat Xv(n_time, n_trials, arma::fill::zeros);
    for (int k = 0; k < K; ++k) {
      const double w = Ca(k, v);
      if (w == 0.0) continue;
      Xv += w * Dk[k];
    }
    arma::vec xall = arma::sum(Xv, 1);     // row-wise sums (n_time x 1)
    arma::vec yv   = Ya.col(v);            // voxel time series

    for (int i = 0; i < n_trials; ++i) {
      const bool have_other = (n_trials > 1);
      const int p = pz + 1 + (have_other ? 1 : 0);
      arma::mat Xd(n_time, p);
      // Z
      Xd.cols(0, pz - 1) = Za;
      // Xi
      Xd.col(pz) = Xv.col(i);
      // Xother
      if (have_other) Xd.col(pz + 1) = xall - Xv.col(i);

      arma::vec beta = qr_coef(Xd, yv);
      Betas(i, v) = beta(pz);  // coefficient of Xi
    }
  }

  return betas;
}

// OpenMP-parallel engine using Armadillo (parallel over voxels)
// [[Rcpp::export]]
Rcpp::NumericMatrix lss_engine_vox_hrf_omp(
    const Rcpp::NumericMatrix& Y,       // n_time x n_vox
    const Rcpp::NumericMatrix& coeffs,  // K x n_vox
    const Rcpp::List& basis_convolved,  // length K; each n_time x n_trials
    const Rcpp::NumericMatrix& Z        // n_time x pz (>=1; intercept allowed)
) {
  const int n_time   = Y.nrow();
  const int n_vox    = Y.ncol();
  const int K        = coeffs.nrow();
  const int n_trials = as<Rcpp::NumericMatrix>(basis_convolved[0]).ncol();
  const int pz       = Z.ncol();

  if (coeffs.ncol() != n_vox) stop("coeffs must be [K x n_vox] with n_vox matching Y");
  if (pz < 1) stop("Z must have at least one column (e.g., intercept).");

  arma::Mat<double> Ya(const_cast<double*>(Y.begin()), n_time, n_vox, false, true);
  arma::Mat<double> Ca(const_cast<double*>(coeffs.begin()), K, n_vox, false, true);
  arma::Mat<double> Za(const_cast<double*>(Z.begin()), n_time, pz, false, true);
  std::vector<arma::Mat<double>> Dk = as_arma_list(basis_convolved, n_time, n_trials);

  Rcpp::NumericMatrix betas(n_trials, n_vox);
  arma::Mat<double> Betas(betas.begin(), n_trials, n_vox, false, true);

  // Parallelize across voxels; each thread writes to distinct column of Betas
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int v = 0; v < n_vox; ++v) {
    arma::mat Xv(n_time, n_trials, arma::fill::zeros);
    for (int k = 0; k < K; ++k) {
      const double w = Ca(k, v);
      if (w == 0.0) continue;
      Xv += w * Dk[k];
    }
    arma::vec xall = arma::sum(Xv, 1);
    arma::vec yv   = Ya.col(v);

    for (int i = 0; i < n_trials; ++i) {
      const bool have_other = (n_trials > 1);
      const int p = pz + 1 + (have_other ? 1 : 0);
      arma::mat Xd(n_time, p);
      Xd.cols(0, pz - 1) = Za;
      Xd.col(pz) = Xv.col(i);
      if (have_other) Xd.col(pz + 1) = xall - Xv.col(i);

      arma::vec beta = qr_coef(Xd, yv);
      Betas(i, v) = beta(pz);
    }
  }

  return betas;
}
