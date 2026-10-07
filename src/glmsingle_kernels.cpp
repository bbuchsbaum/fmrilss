// Kernels for glmsingle(): per-voxel conversion of ridge fractions to
// penalties (alpha), reproducing fracridge's grid interpolation or solving
// the fraction equation exactly.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Thread count for per-voxel loops: n <= 0 means the OpenMP default.
static int glms_threads(int n) {
#ifdef _OPENMP
  return n > 0 ? n : omp_get_max_threads();
#else
  (void)n;
  return 1;
#endif
}

// numpy.interp(x, xp, fp) for increasing xp (ties resolved to the last match).
static double np_interp(double x, const double* xp, const double* fp, int n) {
  if (n == 0 || !std::isfinite(x) || !std::isfinite(xp[0]) ||
      !std::isfinite(xp[n - 1])) return NA_REAL;
  if (x < xp[0] || n == 1) return fp[0];
  if (x >= xp[n - 1]) return fp[n - 1];
  int j = int(std::upper_bound(xp, xp + n, x) - xp) - 1;  // xp[j] <= x < xp[j+1]
  if (j < 0 || j >= n - 1) return NA_REAL;
  if (xp[j] == x) return fp[j];
  const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
  return slope * (x - xp[j]) + fp[j];
}

// fracridge-compatible alphas.
// newlen: grid x voxels matrix of coefficient norms (not yet normalised),
// grid: alpha grid (first entry 0), fracs: requested fractions.
// Returns fracs x voxels alphas (NaN where the OLS norm is zero).
// [[Rcpp::export]]
NumericMatrix glms_frac_alpha_grid(const NumericMatrix& newlen,
                                   const NumericVector& grid,
                                   const NumericVector& fracs,
                                   int n_threads = 1) {
  const int G = newlen.nrow(), V = newlen.ncol(), F = fracs.size();
  if (G == 0 || grid.size() != G) Rcpp::stop("grid must match the nonempty rows of newlen");
  arma::mat res(F, V);
  std::vector<double> fp(G);
  for (int g = 0; g < G; ++g) fp[g] = std::log(1.0 + grid[G - 1 - g]);
  const arma::mat nl(const_cast<double*>(newlen.begin()), G, V, false, true);
  const std::vector<double> fr(fracs.begin(), fracs.end());
  #pragma omp parallel for num_threads(glms_threads(n_threads)) schedule(static)
  for (int v = 0; v < V; ++v) {
    std::vector<double> xp(G);
    const double top = nl(0, v);
    bool valid = std::isfinite(top) && top > 0.0;
    if (valid) {
      for (int g = 0; g < G; ++g) {
        xp[g] = nl(G - 1 - g, v) / top;
        valid = valid && std::isfinite(xp[g]);
      }
    }
    for (int f = 0; f < F; ++f) {
      double t = valid ? np_interp(fr[f], xp.data(), fp.data(), G) : NA_REAL;
      res(f, v) = std::exp(t) - 1.0;
    }
  }
  return Rcpp::wrap(res);
}

// Exact alphas: solve sum(a2 * (s2/(s2+alpha))^2) / sum(a2) = frac^2 for
// alpha >= 0, per voxel and fraction, by bisection on log(alpha).
// a2: components x voxels squared OLS coefficients (bad components zeroed),
// s2: squared singular values per component.
// [[Rcpp::export]]
NumericMatrix glms_frac_alpha_exact(const NumericMatrix& a2,
                                    const NumericVector& s2,
                                    const NumericVector& fracs,
                                    int n_threads = 1) {
  const int N = a2.nrow(), V = a2.ncol(), F = fracs.size();
  arma::mat out(F, V);
  const arma::mat A2(const_cast<double*>(a2.begin()), N, V, false, true);
  const std::vector<double> S2(s2.begin(), s2.end()), FR(fracs.begin(), fracs.end());
  double smax = 0.0;
  for (int i = 0; i < N; ++i) smax = std::max(smax, S2[i]);
  const double lo0 = std::log(1e-12 * std::max(smax, 1e-300));
  const double hi0 = std::log(1e12 * std::max(smax, 1e-300));
  #pragma omp parallel for num_threads(glms_threads(n_threads)) schedule(static)
  for (int v = 0; v < V; ++v) {
    double tot = 0.0;
    for (int i = 0; i < N; ++i) tot += A2(i, v);
    for (int f = 0; f < F; ++f) {
      const double target = FR[f] * FR[f];
      if (tot <= 0.0) { out(f, v) = arma::datum::nan; continue; }
      if (target >= 1.0) { out(f, v) = 0.0; continue; }
      double lo = lo0, hi = hi0;
      for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double alpha = std::exp(mid);
        double acc = 0.0;
        for (int i = 0; i < N; ++i) {
          if (A2(i, v) == 0.0) continue;
          const double sh = S2[i] / (S2[i] + alpha);
          acc += A2(i, v) * sh * sh;
        }
        if (acc / tot > target) lo = mid; else hi = mid;
        if (hi - lo < 1e-13) break;
      }
      out(f, v) = std::exp(0.5 * (lo + hi));
    }
  }
  return Rcpp::wrap(out);
}

// Cross-validation loss of fractional ridge candidates (fractions x voxels).
//
// For each fraction f and run r, the shrunk eigen coefficients are
// W = s2 / (s2 + alpha_f) * a_r (components x voxels) and the candidate
// betas on the cross-validated trial rows are B = V_r[used_r, ] W. With
// z = (B - mu[session]) * isd[session], the loss accumulates sum d z^2 - 2 sum z M (rows
// split additively across runs) plus a constant. Runs with no used rows are
// skipped. Inputs per run are lists; `row_offset` gives each run's first row
// in the stacked used-row arrays mu, isd and M.
// [[Rcpp::export]]
arma::mat glms_frac_cv_loss(const Rcpp::List& Vu,
                            const Rcpp::List& s2,
                            const Rcpp::List& a,
                            const arma::ivec& row_offset,
                            const arma::mat& alphas,
                            const arma::mat& mu,
                            const arma::mat& isd,
                            const arma::ivec& sess,
                            const arma::vec& d,
                            const arma::mat& M,
                            const arma::rowvec& cst,
                            int n_threads = 1) {
  const arma::uword F = alphas.n_rows, V = alphas.n_cols;
  const int R = Vu.size();
  const int nt = glms_threads(n_threads);
  std::vector<arma::mat> Vs(R), As(R);
  std::vector<arma::vec> Ss(R);
  for (int r = 0; r < R; ++r) {
    Vs[r] = Rcpp::as<arma::mat>(Vu[r]);
    Ss[r] = Rcpp::as<arma::vec>(s2[r]);
    As[r] = Rcpp::as<arma::mat>(a[r]);
  }
  arma::mat out(F, V);
  for (arma::uword f = 0; f < F; ++f) {
    arma::rowvec loss = cst;
    arma::rowvec al = alphas.row(f);
    al.replace(arma::datum::nan, 0.0);
    for (int r = 0; r < R; ++r) {
      const arma::mat& Vr = Vs[r];
      if (Vr.n_rows == 0) continue;
      const arma::vec& sr = Ss[r];
      const arma::mat& ar = As[r];
      arma::mat W(ar.n_rows, V);
      #pragma omp parallel for num_threads(nt) schedule(static)
      for (arma::uword v = 0; v < V; ++v) {
        for (arma::uword i = 0; i < ar.n_rows; ++i) {
          const double den = sr[i] + al[v];
          W(i, v) = den > 0.0 ? sr[i] / den * ar(i, v) : 0.0;
        }
      }
      const arma::mat B = Vr * W;
      const arma::uword off = row_offset[r];
      #pragma omp parallel for num_threads(nt) schedule(static)
      for (arma::uword v = 0; v < V; ++v) {
        double acc = 0.0;
        for (arma::uword i = 0; i < B.n_rows; ++i) {
          const int sv = sess[off + i];
          const double z = (B(i, v) - mu(sv, v)) * isd(sv, v);
          acc += d[off + i] * z * z - 2.0 * z * M(off + i, v);
        }
        loss[v] += acc;
      }
    }
    out.row(f) = loss;
  }
  return out;
}

// Column means and sums of squared deviations (two-pass, no copy of Y).
// [[Rcpp::export]]
Rcpp::List glms_col_mean_ss(const arma::mat& Y, int n_threads = 1) {
  const arma::uword T = Y.n_rows, V = Y.n_cols;
  arma::rowvec mean(V), ss(V);
  #pragma omp parallel for num_threads(glms_threads(n_threads)) schedule(static)
  for (arma::uword v = 0; v < V; ++v) {
    const double* y = Y.colptr(v);
    double m = 0.0;
    for (arma::uword t = 0; t < T; ++t) m += y[t];
    m /= double(T);
    double s = 0.0;
    for (arma::uword t = 0; t < T; ++t) { const double e = y[t] - m; s += e * e; }
    mean[v] = m; ss[v] = s;
  }
  return Rcpp::List::create(Rcpp::_["mean"] = mean, Rcpp::_["ss"] = ss);
}

// B with column v multiplied by w[v].
// [[Rcpp::export]]
arma::mat glms_scale_cols(const arma::mat& B, const arma::vec& w, int n_threads = 1) {
  arma::mat out(B.n_rows, B.n_cols);
  #pragma omp parallel for num_threads(glms_threads(n_threads)) schedule(static)
  for (arma::uword v = 0; v < B.n_cols; ++v) out.col(v) = B.col(v) * w[v];
  return out;
}

// Compile GLMsingle's repeated-trial cross-validation (calcbadness) for a
// reference fit `ref` (trials x voxels).
//
// Betas are z-scored per session with the reference's mean and SD (ddof 1).
// For every fold (a set of held-out runs) and condition, each training trial
// i of that condition is compared with each held-out trial j of the same
// condition. With W_ij counting these pairs, the summed squared error of a
// candidate z is sum_i d_i z_i^2 - 2 sum_i z_i M_i + const, where
// d = rowSums(W), M = W z_ref and const = sum_j colSums(W)_j z_ref_j^2.
// Only trials with d > 0 ("used") are returned. `cond_trials` lists the
// 0-based trials of each condition; `test_runs` is runs x folds (logical).
// [[Rcpp::export]]
Rcpp::List glms_cv_compile(const arma::mat& ref,
                           const arma::ivec& session,
                           const arma::ivec& run,
                           const Rcpp::List& cond_trials,
                           const arma::umat& test_runs) {
  const arma::uword N = ref.n_rows, V = ref.n_cols;
  const int S = session.max() + 1;
  arma::mat mu(S, V, arma::fill::zeros), sd(S, V, arma::fill::zeros);
  arma::vec cnt(S, arma::fill::zeros);
  for (arma::uword i = 0; i < N; ++i) { mu.row(session[i]) += ref.row(i); cnt[session[i]] += 1; }
  for (int s = 0; s < S; ++s) if (cnt[s] > 0) mu.row(s) /= cnt[s];
  for (arma::uword i = 0; i < N; ++i) {
    const arma::rowvec e = ref.row(i) - mu.row(session[i]);
    sd.row(session[i]) += e % e;
  }
  for (int s = 0; s < S; ++s) sd.row(s) = arma::sqrt(sd.row(s) / (cnt[s] - 1.0));
  arma::mat isd_ref(S, V);
  for (int s = 0; s < S; ++s)
    for (arma::uword v = 0; v < V; ++v)
      isd_ref(s, v) = sd(s, v) == 0.0 ? 0.0 : 1.0 / sd(s, v);
  arma::mat zref(N, V);
  for (arma::uword i = 0; i < N; ++i)
    zref.row(i) = (ref.row(i) - mu.row(session[i])) % isd_ref.row(session[i]);

  arma::vec d(N, arma::fill::zeros);
  arma::mat M(N, V, arma::fill::zeros);
  arma::rowvec cst(V, arma::fill::zeros);
  const arma::uword F = test_runs.n_cols;
  arma::rowvec tsum(V), tsq(V);
  for (int c = 0; c < cond_trials.size(); ++c) {
    const Rcpp::IntegerVector tr = cond_trials[c];
    if (tr.size() < 2) continue;
    for (arma::uword f = 0; f < F; ++f) {
      tsum.zeros(); tsq.zeros();
      int n_test = 0, n_train = 0;
      for (int k = 0; k < tr.size(); ++k) {
        const int j = tr[k];
        if (test_runs(run[j], f)) {
          tsum += zref.row(j); tsq += zref.row(j) % zref.row(j); ++n_test;
        } else {
          ++n_train;
        }
      }
      if (!n_test || !n_train) continue;
      for (int k = 0; k < tr.size(); ++k) {
        const int i = tr[k];
        if (test_runs(run[i], f)) continue;
        d[i] += n_test;
        M.row(i) += tsum;
      }
      cst += double(n_train) * tsq;
    }
  }
  const arma::uvec used = arma::find(d > 0);
  const arma::ivec used1 = arma::conv_to<arma::ivec>::from(used) + 1;
  const arma::vec d_used = d.elem(used);
  const arma::mat M_used = M.rows(used);
  const arma::ivec sess_used = session.elem(used);
  return Rcpp::List::create(
    Rcpp::_["used"] = used1, Rcpp::_["d"] = d_used, Rcpp::_["M"] = M_used,
    Rcpp::_["const"] = cst, Rcpp::_["mu"] = mu, Rcpp::_["isd_ref"] = isd_ref,
    Rcpp::_["sd"] = sd, Rcpp::_["session_used"] = sess_used);
}

// Loss (per voxel) of candidate betas on the used rows. `isd` is the
// per-session inverse SD applied to candidates (zero-SD handling differs
// between the reference and candidates in GLMsingle's Python code).
// [[Rcpp::export]]
arma::rowvec glms_cv_loss_cpp(const arma::mat& cand, const arma::mat& mu,
                              const arma::mat& isd, const arma::ivec& sess,
                              const arma::vec& d, const arma::mat& M,
                              const arma::rowvec& cst, int n_threads = 1) {
  const arma::uword N = cand.n_rows, V = cand.n_cols;
  arma::rowvec out = cst;
  #pragma omp parallel for num_threads(glms_threads(n_threads)) schedule(static)
  for (arma::uword v = 0; v < V; ++v) {
    double acc = 0.0;
    for (arma::uword i = 0; i < N; ++i) {
      const int s = sess[i];
      const double z = (cand(i, v) - mu(s, v)) * isd(s, v);
      acc += d[i] * z * z - 2.0 * z * M(i, v);
    }
    out[v] += acc;
  }
  return out;
}
