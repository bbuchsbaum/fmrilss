// File: src/oasis_core.cpp
// OASIS-HRF implementation for fmrilss package
// Supports both single-basis (K=1) and multi-basis (K>1) HRFs
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ============================================================================
// Single-basis (K=1) functions
// ============================================================================

/*
 * Precompute residualized design terms for one condition (single-basis case).
 * Inputs:
 *   X_trials (T x N): trial-specific columns from fmrihrf::evaluate(regressor_set(...))
 *   N_nuis   (T x P): nuisances: confounds + other-condition aggregates (no within-condition aggregate)
 * Outputs: Q, A=R X_trials, s_all=R sum_j x_j, and design-only scalars {d, alpha, s}
 */
// [[Rcpp::export]]
Rcpp::List oasis_precompute_design(const arma::mat& X_trials,
                                   const arma::mat& N_nuis) {
  const uword T = X_trials.n_rows, N = X_trials.n_cols;
  arma::mat Q, Rq;
  if (N_nuis.n_cols > 0) qr_econ(Q, Rq, N_nuis); else Q.set_size(T, 0);

  arma::mat A = X_trials;
  if (Q.n_cols > 0) A -= Q * (Q.t() * A);

  arma::vec s_all = sum(A, 1);
  const double s_all_norm2 = dot(s_all, s_all);
  arma::vec d = sum(square(A), 0).t();
  arma::vec aTs = A.t() * s_all;
  arma::vec alpha = aTs - d;
  arma::vec s = s_all_norm2 + d - 2.0 * aTs;

  return Rcpp::List::create(_["Q"]=Q,_["A"]=A,_["s_all"]=s_all,
                            _["d"]=d,_["alpha"]=alpha,_["s"]=s,
                            _["s_all_norm2"]=s_all_norm2);
}

/*
 * Blocked products: N_Y = A' (R Y), S_Y = s_all' (R Y), where R = I - QQ'
 */
// [[Rcpp::export]]
Rcpp::List oasis_AtY_SY_blocked(const arma::mat& A,
                                const arma::vec& s_all,
                                const arma::mat& Q,
                                const arma::mat& Y,
                                const int block_cols = 4096) {
  const uword V = Y.n_cols, N = A.n_cols;
  arma::mat N_Y(N, V, fill::zeros);
  arma::rowvec S_Y(V, fill::zeros);
  arma::rowvec RY_norm2(V, fill::zeros);

  for (uword start = 0; start < V; start += block_cols) {
    uword end = std::min<uword>(V - 1, start + block_cols - 1);
    arma::mat Yb = Y.cols(start, end);
    if (Q.n_cols > 0) Yb -= Q * (Q.t() * Yb);           // residualize once
    N_Y.cols(start, end) = A.t() * Yb;                  // n1 = A'Ry
    S_Y.cols(start, end) = s_all.t() * Yb;              // s_all'Ry
    RY_norm2.cols(start, end) = sum(square(Yb), 0);    // for SEs
  }
  return Rcpp::List::create(_["N_Y"]=N_Y, _["S_Y"]=S_Y, _["RY_norm2"]=RY_norm2);
}

/*
 * Closed-form LSS (with optional ridge on the [a_j, b_j] 2x2 Gram)
 *   d_j = ||a_j||^2, e_j = ||b_j||^2 = s_j, c_j = a_j' b_j = alpha_j
 *   n1  = a_j'Ry  (row j of N_Y),  n2 = b_j'Ry = S_Y - n1
 *   beta_j = ( n1*(e+λ_b) - c*n2 ) / ( (d+λ_x)*(e+λ_b) - c^2 )
 */
// [[Rcpp::export]]
arma::mat oasis_betas_closed_form(const arma::mat& N_Y,
                                  const arma::rowvec& S_Y,
                                  const arma::vec& d,
                                  const arma::vec& alpha,
                                  const arma::vec& s,
                                  const double ridge_x = 0.0,
                                  const double ridge_b = 0.0,
                                  const double denom_eps = 1e-12) {
  const uword N = N_Y.n_rows, V = N_Y.n_cols;
  arma::mat B(N, V, fill::zeros);
  for (uword j = 0; j < N; ++j) {
    const double dj = d[j] + ridge_x;
    const double ej = s[j] + ridge_b;
    const double cj = alpha[j];
    const double denom = std::max(dj * ej - cj * cj, denom_eps);
    B.row(j) = (N_Y.row(j) * ej - cj * (S_Y - N_Y.row(j))) / denom;
  }
  return B;
}

/*
 * Optional: companion function to also return gamma (aggregator coefficient),
 * useful for computing SSE/SE without re-running products.
 */
// [[Rcpp::export]]
Rcpp::List oasis_betas_gammas(const arma::mat& N_Y,
                              const arma::rowvec& S_Y,
                              const arma::vec& d,
                              const arma::vec& alpha,
                              const arma::vec& s,
                              const double ridge_x = 0.0,
                              const double ridge_b = 0.0,
                              const double denom_eps = 1e-12) {
  const uword N = N_Y.n_rows, V = N_Y.n_cols;
  arma::mat B(N, V, fill::zeros);
  arma::mat G(N, V, fill::zeros);

  for (uword j = 0; j < N; ++j) {
    const double dj = d[j] + ridge_x;
    const double ej = s[j] + ridge_b;
    const double cj = alpha[j];
    const double denom = std::max(dj * ej - cj * cj, denom_eps);

    const arma::rowvec n1 = N_Y.row(j);
    const arma::rowvec n2 = S_Y - n1;

    B.row(j) = ( n1 * ej - cj * n2 ) / denom; // beta
    G.row(j) = ( n2 * dj - cj * n1 ) / denom; // gamma
  }
  return Rcpp::List::create(_["beta"]=B, _["gamma"]=G);
}

// ============================================================================
// Multi-basis (K>1) functions
// ============================================================================

/*
 * Precompute for multi-basis (K>1).
 * X_trials: T x (N*K) with columns grouped by trial: [A1 | A2 | ... | AN], each Aj is T x K.
 * N_nuis:   T x P nuisance design (confounds + other-condition aggregates).
 * Returns: Q, A (residualized X), S (sum_j Aj), and per-trial Gram blocks D_j, C_j, E_j.
 */
// [[Rcpp::export]]
Rcpp::List oasisk_precompute_design(const arma::mat& X_trials,
                                    const arma::mat& N_nuis,
                                    const int K) {
  const uword T = X_trials.n_rows;
  const uword NK = X_trials.n_cols;
  if (NK % K != 0) stop("X_trials ncol not divisible by K.");
  const uword N = NK / K;

  arma::mat Q, Rq;
  if (N_nuis.n_cols > 0) qr_econ(Q, Rq, N_nuis); else Q.set_size(T, 0);

  arma::mat A = X_trials;
  if (Q.n_cols > 0) A -= Q * (Q.t() * A);

  // S = sum_j Aj  (T x K)
  arma::mat S(T, K, fill::zeros);
  for (uword j = 0; j < N; ++j) {
    S += A.cols(j*K, j*K + (K-1));
  }

  // Precompute cross-products
  arma::mat SS = S.t() * S;                 // K x K
  arma::cube D(K, K, N);                    // Aj'Aj
  arma::cube C(K, K, N);                    // Aj'Bj
  arma::cube E(K, K, N);                    // Bj'Bj

  for (uword j = 0; j < N; ++j) {
    arma::mat Aj = A.cols(j*K, j*K + (K-1));   // T x K
    arma::mat Dj = Aj.t() * Aj;                // K x K
    arma::mat Mj = Aj.t() * S;                 // K x K  (Aj'S)
    arma::mat Cj = Mj - Dj;                    // Aj'Bj
    arma::mat Ej = SS - Mj - Mj.t() + Dj;      // Bj'Bj

    D.slice(j) = Dj;
    C.slice(j) = Cj;
    E.slice(j) = Ej;
  }

  return Rcpp::List::create(
    _["Q"] = Q, _["A"] = A, _["S"] = S,
    _["D"] = D, _["C"] = C, _["E"] = E
  );
}

/*
 * Blocked products for multi-basis:
 *   N1 = A' (R Y)  -> (N*K x V)
 *   SY = S' (R Y)  -> (K x V)
 *   RY_norm2       -> (1 x V)
 */
// [[Rcpp::export]]
Rcpp::List oasisk_products(const arma::mat& A,
                           const arma::mat& S,
                           const arma::mat& Q,
                           const arma::mat& Y,
                           const int block_cols = 4096) {
  const uword V  = Y.n_cols;
  const uword NK = A.n_cols;
  const uword K  = S.n_cols;

  arma::mat    N1(NK, V, fill::zeros);
  arma::rowvec RY_norm2(V, fill::zeros);
  arma::mat    SY(K, V, fill::zeros);

  for (uword start = 0; start < V; start += block_cols) {
    const uword end = std::min<uword>(V - 1, start + block_cols - 1);
    arma::mat Yb = Y.cols(start, end);
    if (Q.n_cols > 0) Yb -= Q * (Q.t() * Yb);

    RY_norm2.cols(start, end) = sum(square(Yb), 0);
    N1.cols(start, end) = A.t() * Yb;         // (N*K) x nb
    SY.cols(start, end) = S.t() * Yb;         // K x nb
  }

  return Rcpp::List::create(_["N1"] = N1, _["SY"] = SY, _["RY_norm2"] = RY_norm2);
}

/*
 * Solve for betas with ridge on block Gram (K>1).
 * Inputs:
 *   D, C, E : cubes (K x K x N)
 *   N1      : (N*K x V), stacked by trials
 *   SY      : (K x V)
 *   ridge_x, ridge_b : absolute ridge strengths
 * Output:
 *   B : (N*K x V) betas, stacked by trials (K rows per trial)
 */
// [[Rcpp::export]]
arma::mat oasisk_betas(const arma::cube& D,
                       const arma::cube& C,
                       const arma::cube& E,
                       const arma::mat&  N1,
                       const arma::mat&  SY,
                       const double ridge_x = 0.0,
                       const double ridge_b = 0.0,
                       const double diag_eps = 1e-10) {
  const uword K = D.n_rows;
  const uword N = D.n_slices;
  const uword V = N1.n_cols;

  arma::mat B(N*K, V, fill::zeros);

  arma::mat I_K = eye(K, K);
  arma::mat G(2*K, 2*K), L, RHS, X;

  for (uword j = 0; j < N; ++j) {
    // Assemble block Gram with ridge
    arma::mat Dj = D.slice(j) + ridge_x * I_K;
    arma::mat Cj = C.slice(j);
    arma::mat Ej = E.slice(j) + ridge_b * I_K;

    G.zeros();
    G.submat(0,   0,   K-1,   K-1  ) = Dj;
    G.submat(0,   K,   K-1, 2*K-1 ) = Cj;
    G.submat(K,   0, 2*K-1,   K-1 ) = Cj.t();
    G.submat(K,   K, 2*K-1, 2*K-1 ) = Ej;

    // KxV RHS: N1_j; and KxV RHS: N2_j = SY - N1_j
    const uword r0 = j*K, r1 = r0 + (K-1);
    arma::mat N1j = N1.rows(r0, r1);
    arma::mat N2j = SY - N1j;

    RHS.set_size(2*K, V);
    RHS.rows(0,   K-1)   = N1j;
    RHS.rows(K, 2*K-1)   = N2j;

    // Cholesky solve for many RHS
    bool ok = chol(L, G, "lower");
    if (!ok) {
      G.diag() += diag_eps;         // tiny jitter if needed
      chol(L, G, "lower");
    }
    arma::mat Z = solve(trimatl(L), RHS);
    X           = solve(trimatu(L.t()), Z);

    // Betas are the top K rows
    B.rows(r0, r1) = X.rows(0, K-1);
  }

  return B;
}

// [[Rcpp::export]]
arma::vec oasisk_compute_RY_norm2(const arma::mat& Q, const arma::mat& Y) {
  // Compute ||R*Y||^2 for each column, where R = I - Q*Q'
  arma::mat RY = Y - Q * (Q.t() * Y);
  return sum(RY % RY, 0).t();  // Column-wise squared norms
}

// [[Rcpp::export]]
Rcpp::List oasisk_betas_se(const arma::cube& D,
                           const arma::cube& C,
                           const arma::cube& E,
                           const arma::mat& N1,
                           const arma::mat& SY,
                           const arma::vec& RY_norm2,
                           double ridge_x = 0.0,
                           double ridge_b = 0.0) {
  
  int N = D.n_slices;
  int K = D.n_rows;
  int V = N1.n_cols;
  
  arma::mat B(N * K, V, arma::fill::zeros);
  arma::mat SE(N * K, V, arma::fill::zeros);
  
  double diag_eps = 1e-10;
  
  // Precompute dof (simplified - could be refined)
  int dof = std::max(1, (int)RY_norm2.n_elem - N*K - 2*K);
  
  for (int j = 0; j < N; ++j) {
    int r0 = j * K;
    int r1 = r0 + K - 1;
    
    // Build 2K×2K Gram
    arma::mat G(2*K, 2*K, arma::fill::zeros);
    G.submat(0,   0,   K-1, K-1) = D.slice(j) + ridge_x * arma::eye(K, K);
    G.submat(0,   K,   K-1, 2*K-1) = C.slice(j);
    G.submat(K,   0,   2*K-1, K-1) = C.slice(j).t();
    G.submat(K,   K,   2*K-1, 2*K-1) = E.slice(j) + ridge_b * arma::eye(K, K);
    
    // Prepare many RHS - extract K rows for trial j
    arma::mat N1j = N1.rows(r0, r1);      // K x V
    arma::mat sum_N1j = sum(N1j, 0);      // 1 x V (sum across K basis)
    arma::mat N2j = SY - repmat(sum_N1j, K, 1);  // K x V
    arma::mat RHS(2*K, V);
    RHS.rows(0,   K-1) = N1j;
    RHS.rows(K, 2*K-1) = N2j;
    
    // Cholesky solve
    arma::mat L;
    bool ok = chol(L, G, "lower");
    if (!ok) {
      G.diag() += diag_eps;
      chol(L, G, "lower");
    }
    arma::mat Z = solve(trimatl(L), RHS);
    arma::mat X = solve(trimatu(L.t()), Z);
    
    // Extract betas and gammas
    arma::mat beta_j = X.rows(0, K-1);
    arma::mat gamma_j = X.rows(K, 2*K-1);
    B.rows(r0, r1) = beta_j;
    
    // Compute SSE for each voxel
    for (int v = 0; v < V; ++v) {
      arma::vec bv = beta_j.col(v);
      arma::vec gv = gamma_j.col(v);
      arma::vec n1v = N1j.col(v);
      arma::vec n2v = N2j.col(v);
      
      double sse = RY_norm2(v) - 
                   2*dot(bv, n1v) - 2*dot(gv, n2v) +
                   as_scalar(bv.t() * D.slice(j) * bv) +
                   as_scalar(gv.t() * E.slice(j) * gv) +
                   2*as_scalar(bv.t() * C.slice(j) * gv);
      
      double sigma2 = std::max(sse / dof, 0.0);
      
      // G^{-1} for variance (top-left K×K block)
      arma::mat Ginv = inv_sympd(G);
      arma::mat var_beta = sigma2 * Ginv.submat(0, 0, K-1, K-1);
      
      // Extract diagonal for SEs
      for (int k = 0; k < K; ++k) {
        SE(r0 + k, v) = std::sqrt(var_beta(k, k));
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("beta") = B,
    Rcpp::Named("se") = SE
  );
}