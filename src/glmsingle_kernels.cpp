// Kernels for glmsingle(): per-voxel conversion of ridge fractions to
// penalties (alpha), reproducing fracridge's grid interpolation or solving
// the fraction equation exactly.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// numpy.interp(x, xp, fp) for increasing xp (ties resolved to the last match).
static double np_interp(double x, const double* xp, const double* fp, int n) {
  if (std::isnan(x)) return NA_REAL;
  if (x <= xp[0]) return (x == xp[0] || n == 1) ? fp[0] : fp[0];
  if (x >= xp[n - 1]) return fp[n - 1];
  int j = int(std::upper_bound(xp, xp + n, x) - xp) - 1;  // xp[j] <= x < xp[j+1]
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
                                   const NumericVector& fracs) {
  const int G = newlen.nrow(), V = newlen.ncol(), F = fracs.size();
  NumericMatrix out(F, V);
  std::vector<double> xp(G), fp(G);
  for (int g = 0; g < G; ++g) fp[g] = std::log(1.0 + grid[G - 1 - g]);
  for (int v = 0; v < V; ++v) {
    const double top = newlen(0, v);
    for (int g = 0; g < G; ++g) xp[g] = newlen(G - 1 - g, v) / top;
    for (int f = 0; f < F; ++f) {
      double t = np_interp(fracs[f], xp.data(), fp.data(), G);
      out(f, v) = std::exp(t) - 1.0;
    }
  }
  return out;
}

// Exact alphas: solve sum(a2 * (s2/(s2+alpha))^2) / sum(a2) = frac^2 for
// alpha >= 0, per voxel and fraction, by bisection on log(alpha).
// a2: components x voxels squared OLS coefficients (bad components zeroed),
// s2: squared singular values per component.
// [[Rcpp::export]]
NumericMatrix glms_frac_alpha_exact(const NumericMatrix& a2,
                                    const NumericVector& s2,
                                    const NumericVector& fracs) {
  const int N = a2.nrow(), V = a2.ncol(), F = fracs.size();
  NumericMatrix out(F, V);
  double smax = 0.0;
  for (int i = 0; i < N; ++i) smax = std::max(smax, s2[i]);
  for (int v = 0; v < V; ++v) {
    double tot = 0.0;
    for (int i = 0; i < N; ++i) tot += a2(i, v);
    for (int f = 0; f < F; ++f) {
      const double target = fracs[f] * fracs[f];
      if (tot <= 0.0) { out(f, v) = NA_REAL; continue; }
      if (target >= 1.0) { out(f, v) = 0.0; continue; }
      auto ratio = [&](double alpha) {
        double acc = 0.0;
        for (int i = 0; i < N; ++i) {
          if (a2(i, v) == 0.0) continue;
          const double sh = s2[i] / (s2[i] + alpha);
          acc += a2(i, v) * sh * sh;
        }
        return acc / tot;
      };
      double lo = std::log(1e-12 * std::max(smax, 1e-300));
      double hi = std::log(1e12 * std::max(smax, 1e-300));
      for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (ratio(std::exp(mid)) > target) lo = mid; else hi = mid;
        if (hi - lo < 1e-13) break;
      }
      out(f, v) = std::exp(0.5 * (lo + hi));
    }
  }
  return out;
}

// Cross-validation loss of fractional ridge candidates (fractions x voxels).
//
// For each fraction f and run r, the shrunk eigen coefficients are
// W = s2 / (s2 + alpha_f) * a_r (components x voxels) and the candidate
// betas on the cross-validated trial rows are B = V_r[used_r, ] W. With
// z = (B - mu) * isd, the loss accumulates sum d z^2 - 2 sum z M (rows
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
                            const arma::vec& d,
                            const arma::mat& M,
                            const arma::rowvec& cst) {
  const arma::uword F = alphas.n_rows, V = alphas.n_cols;
  const int R = Vu.size();
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
      for (arma::uword v = 0; v < V; ++v) {
        for (arma::uword i = 0; i < ar.n_rows; ++i) {
          const double den = sr[i] + al[v];
          W(i, v) = den > 0.0 ? sr[i] / den * ar(i, v) : 0.0;
        }
      }
      const arma::mat B = Vr * W;
      const arma::uword off = row_offset[r];
      for (arma::uword v = 0; v < V; ++v) {
        double acc = 0.0;
        for (arma::uword i = 0; i < B.n_rows; ++i) {
          const double z = (B(i, v) - mu(off + i, v)) * isd(off + i, v);
          acc += d[off + i] * z * z - 2.0 * z * M(off + i, v);
        }
        loss[v] += acc;
      }
    }
    out.row(f) = loss;
  }
  return out;
}
