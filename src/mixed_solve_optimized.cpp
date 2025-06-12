#include <RcppArmadillo.h>
#include <Rmath.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

//' Mixed Model Workspace for Optimized Computation
//' 
//' Stores precomputed matrices and decompositions that can be reused
//' across multiple voxels to avoid repeated expensive computations.
//' 
//' @name MixedWorkspace
struct MixedWorkspace {
  mat Q;           // n × (n-p) matrix from QR of nuisance regressors
  mat U;           // Eigenvectors from spectral decomposition  
  vec phi;         // Eigenvalues from spectral decomposition
  vec theta;       // Transformed eigenvalues for REML
  mat X_orig;      // Original fixed effects matrix
  mat Z_orig;      // Original random effects matrix  
  mat XtXinv;      // (X'X)^-1 for fixed effects
  uword n;         // Number of timepoints
  uword p;         // Number of fixed effects
  uword q;         // Number of random effects
  bool use_identity_K; // Whether kinship matrix is identity
  
  // Default constructor
  MixedWorkspace() : n(0), p(0), q(0), use_identity_K(true) {}
};

//' Fast analytical REML estimation for single variance component
//' 
//' For a single variance component model, the REML estimate of λ = σe²/σu²
//' has a closed-form solution that can be computed efficiently.
//' 
//' @param omega Transformed response vector Q'y
//' @param theta Transformed eigenvalues 
//' @param tol Convergence tolerance for Newton iterations
//' @param max_iter Maximum Newton iterations
//' @return Estimated variance ratio λ
//' @name fast_reml_lambda
double fast_reml_lambda(const vec& omega, const vec& theta, 
                       double tol = 1e-8, int max_iter = 10) {
  
  uword n_minus_p = omega.n_elem;
  
  // Initial guess: method of moments estimator
  double ss_total = dot(omega, omega);
  double trace_theta = sum(theta);
  double lambda = std::max(1e-6, (ss_total / n_minus_p - 1.0) / (trace_theta / n_minus_p));
  
  // Newton-Raphson iterations for REML equation
  for (int iter = 0; iter < max_iter; ++iter) {
    
    // Compute current likelihood derivatives
    vec theta_plus_lambda = theta + lambda;
    vec weights = 1.0 / theta_plus_lambda;
    
    // First derivative of REML log-likelihood
    double trace_term = sum(weights);
    double quad_term = dot(omega % weights, omega);
    double first_deriv = -0.5 * trace_term + 0.5 * quad_term / lambda;
    
    // Second derivative  
    vec weights_sq = weights % weights;
    double trace_term2 = sum(weights_sq);
    double second_deriv = 0.5 * trace_term2 - 0.5 * quad_term / (lambda * lambda);
    
    // Newton update
    double delta = first_deriv / second_deriv;
    double lambda_new = lambda - delta;
    
    // Ensure lambda stays positive
    lambda_new = std::max(1e-9, lambda_new);
    
    // Check convergence
    if (std::abs(delta) < tol * lambda) {
      return lambda_new;
    }
    
    lambda = lambda_new;
  }
  
  // If no convergence, return current estimate
  return std::max(1e-9, lambda);
}

//' Precompute workspace for mixed model optimization
//' 
//' Performs all the expensive computations that don't depend on the response
//' vector y, so they can be reused across multiple voxels.
//' 
//' @param X Fixed effects design matrix (n × p)
//' @param Z Random effects design matrix (n × q) 
//' @param K Kinship/covariance matrix for random effects (q × q)
//' @return MixedWorkspace object containing precomputed matrices
// [[Rcpp::export]]
List mixed_precompute_workspace(const arma::mat& X, 
                               const arma::mat& Z,
                               const arma::mat& K) {
  
  MixedWorkspace ws;
  ws.n = X.n_rows;
  ws.p = X.n_cols;
  ws.q = Z.n_cols;
  ws.X_orig = X;
  ws.Z_orig = Z;
  
  // Check if K is identity matrix (common case)
  ws.use_identity_K = (K.n_rows == K.n_cols) && 
                      approx_equal(K, eye(K.n_rows, K.n_cols), "absdiff", 1e-10);
  
  try {
    
    // Step 1: QR decomposition of fixed effects X
    mat Q_full, R;
    qr_econ(Q_full, R, X);
    
    // Check if X is full rank
    vec r_diag = R.diag();
    if (any(abs(r_diag) < 1e-12)) {
      stop("Fixed effects matrix X is not full rank");
    }
    
    ws.XtXinv = inv_sympd(R.t() * R);
    
    // Step 2: Orthogonal complement to X (projection away from fixed effects)
    mat I_n = eye(ws.n, ws.n);
    mat P_X = Q_full * Q_full.t();
    mat P_perp = I_n - P_X;
    
    // Step 3: Project Z onto orthogonal complement  
    mat Z_proj = P_perp * Z;
    
    // Step 4: Form covariance matrix and perform spectral decomposition
    mat Sigma;
    if (ws.use_identity_K) {
      // K = I case: much simpler
      Sigma = Z_proj * Z_proj.t();
    } else {
      // General case: K not identity
      Sigma = Z_proj * K * Z_proj.t();
    }
    
    // Eigendecomposition of projected covariance matrix
    vec eigenvals;
    mat eigenvecs;
    
    if (!eig_sym(eigenvals, eigenvecs, Sigma)) {
      stop("Failed to compute eigendecomposition of covariance matrix");
    }
    
    // Remove near-zero eigenvalues and corresponding eigenvectors
    uvec pos_idx = find(eigenvals > 1e-10);
    if (pos_idx.n_elem == 0) {
      stop("All eigenvalues are effectively zero");
    }
    
    ws.phi = eigenvals(pos_idx);
    ws.U = eigenvecs.cols(pos_idx);
    
    // Step 5: Form Q matrix for transformation (reduced rank)
    // Use only the columns corresponding to positive eigenvalues
    uword effective_rank = pos_idx.n_elem;
    ws.Q = ws.U;  // For simplicity, use eigenvectors directly
    
    // For REML, we need theta = phi (eigenvalues)
    ws.theta = ws.phi;
    
    Rcout << "Workspace precomputed successfully:" << std::endl;
    Rcout << "  - n=" << ws.n << ", p=" << ws.p << ", q=" << ws.q << std::endl;
    Rcout << "  - Effective rank: " << ws.phi.n_elem << std::endl;
    Rcout << "  - Using identity K: " << (ws.use_identity_K ? "yes" : "no") << std::endl;
    
  } catch (const std::exception& e) {
    stop("Error in workspace precomputation: " + std::string(e.what()));
  }
  
  // Return as R list (we'll extract this in the R wrapper)
  return List::create(
    Named("Q") = ws.Q,
    Named("U") = ws.U, 
    Named("phi") = ws.phi,
    Named("theta") = ws.theta,
    Named("X_orig") = ws.X_orig,
    Named("Z_orig") = ws.Z_orig,
    Named("XtXinv") = ws.XtXinv,
    Named("n") = ws.n,
    Named("p") = ws.p,
    Named("q") = ws.q,
    Named("use_identity_K") = ws.use_identity_K
  );
}

//' Convert R list back to MixedWorkspace
//' @name list_to_workspace
MixedWorkspace list_to_workspace(const List& ws_list) {
  MixedWorkspace ws;
  ws.Q = as<mat>(ws_list["Q"]);
  ws.U = as<mat>(ws_list["U"]);
  ws.phi = as<vec>(ws_list["phi"]);
  ws.theta = as<vec>(ws_list["theta"]);
  ws.X_orig = as<mat>(ws_list["X_orig"]);
  ws.Z_orig = as<mat>(ws_list["Z_orig"]);
  ws.XtXinv = as<mat>(ws_list["XtXinv"]);
  ws.n = as<uword>(ws_list["n"]);
  ws.p = as<uword>(ws_list["p"]);
  ws.q = as<uword>(ws_list["q"]);
  ws.use_identity_K = as<bool>(ws_list["use_identity_K"]);
  return ws;
}

//' Fast single-voxel mixed model estimation using precomputed workspace
//' 
//' @param y Response vector for single voxel
//' @param ws_list Precomputed workspace (as R list)
//' @param compute_se Whether to compute standard errors
//' @return List with beta, u, Vu, Ve, and optionally standard errors
// [[Rcpp::export]]
List mixed_single_voxel_cpp(const arma::vec& y,
                           const List& ws_list,
                           bool compute_se = false) {
  
  MixedWorkspace ws = list_to_workspace(ws_list);
  
  if (y.n_elem != ws.n) {
    stop("Response vector length doesn't match workspace dimensions");
  }
  
  try {
    
    // Step 1: Transform response vector
    vec omega = ws.Q.t() * y;
    
    // Ensure omega and theta have compatible dimensions
    if (omega.n_elem != ws.theta.n_elem) {
      // Truncate omega to match theta dimensions
      omega = omega(span(0, ws.theta.n_elem - 1));
    }
    
    // Step 2: Fast REML estimation of λ
    double lambda = fast_reml_lambda(omega, ws.theta);
    
    // Step 3: Compute variance components
    vec theta_plus_lambda = ws.theta + lambda;
    vec weights = omega / theta_plus_lambda;
    double Vu = dot(omega, weights) / (ws.n - ws.p);
    double Ve = lambda * Vu;
    
    // Step 4: Compute BLUP estimates
    // Fixed effects β - use more stable computation
    vec weights_inv = 1.0 / (ws.phi + lambda);
    mat H_inv = ws.U * diagmat(weights_inv) * ws.U.t();
    
    // Add small regularization to ensure numerical stability
    mat XtHinvX = ws.X_orig.t() * H_inv * ws.X_orig;
    XtHinvX.diag() += 1e-12; // Small regularization
    
    mat XtHinvX_inv;
    if (!inv_sympd(XtHinvX_inv, XtHinvX)) {
      // Fallback to regular inverse if symmetric positive definite fails
      XtHinvX_inv = inv(XtHinvX);
    }
    vec beta = XtHinvX_inv * ws.X_orig.t() * H_inv * y;
    
    // Random effects u (BLUP)
    vec resid = y - ws.X_orig * beta;
    mat K_ZtHinv = ws.Z_orig.t() * H_inv; // Assuming K = I
    vec u = K_ZtHinv * resid;
    
    // Standard errors (if requested)
    vec beta_se, u_se;
    if (compute_se) {
      mat beta_cov = Vu * XtHinvX_inv;
      beta_se = sqrt(beta_cov.diag());
      
      // Approximate standard errors for random effects
      mat u_cov_diag = Vu * diagmat(ones(ws.q)); // Simplified
      u_se = sqrt(u_cov_diag.diag());
    }
    
    List result = List::create(
      Named("beta") = beta,
      Named("u") = u,
      Named("Vu") = Vu,
      Named("Ve") = Ve,
      Named("lambda") = lambda
    );
    
    if (compute_se) {
      result["beta_se"] = beta_se;
      result["u_se"] = u_se;
    }
    
    return result;
    
  } catch (const std::exception& e) {
    stop("Error in single voxel computation: " + std::string(e.what()));
  }
}

//' Optimized multi-voxel mixed model estimation
//' 
//' Uses precomputed workspace and parallel processing to efficiently
//' estimate mixed models across many voxels.
//' 
//' @param Y Response matrix (n × V) where V is number of voxels
//' @param ws_list Precomputed workspace (as R list)
//' @param compute_se Whether to compute standard errors
//' @param n_threads Number of OpenMP threads (0 = auto)
//' @return List with matrices of estimates across voxels
// [[Rcpp::export]]
List mixed_multi_voxel_cpp(const arma::mat& Y,
                          const List& ws_list,
                          bool compute_se = false,
                          int n_threads = 0) {
  
  MixedWorkspace ws = list_to_workspace(ws_list);
  uword n_voxels = Y.n_cols;
  
  if (Y.n_rows != ws.n) {
    stop("Response matrix rows don't match workspace dimensions");
  }
  
  // Set up OpenMP
#ifdef _OPENMP
  if (n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
  int actual_threads = omp_get_max_threads();
  Rcout << "Using " << actual_threads << " OpenMP threads" << std::endl;
#endif
  
  // Initialize result matrices
  mat beta_results(ws.p, n_voxels);
  mat u_results(ws.q, n_voxels);
  vec Vu_results(n_voxels);
  vec Ve_results(n_voxels);
  vec lambda_results(n_voxels);
  
  mat beta_se_results, u_se_results;
  if (compute_se) {
    beta_se_results.set_size(ws.p, n_voxels);
    u_se_results.set_size(ws.q, n_voxels);
  }
  
  // Process voxels in parallel
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (uword v = 0; v < n_voxels; ++v) {
    
    try {
      vec y = Y.col(v);
      
      // Check for missing values
      uvec finite_idx = find_finite(y);
      if (finite_idx.n_elem < ws.n) {
        // Handle missing data by setting results to NA
        beta_results.col(v).fill(datum::nan);
        u_results.col(v).fill(datum::nan);
        Vu_results(v) = datum::nan;
        Ve_results(v) = datum::nan;
        lambda_results(v) = datum::nan;
        continue;
      }
      
      // Transform response
      vec omega = ws.Q.t() * y;
      
      // Ensure compatible dimensions
      if (omega.n_elem != ws.theta.n_elem) {
        omega = omega(span(0, ws.theta.n_elem - 1));
      }
      
      // Fast REML estimation
      double lambda = fast_reml_lambda(omega, ws.theta);
      lambda_results(v) = lambda;
      
      // Variance components
      vec theta_plus_lambda = ws.theta + lambda;
      vec weights = omega / theta_plus_lambda;
      double Vu = dot(omega, weights) / (ws.n - ws.p);
      double Ve = lambda * Vu;
      Vu_results(v) = Vu;
      Ve_results(v) = Ve;
      
             // BLUP estimates - use more stable computation
       vec weights_inv = 1.0 / (ws.phi + lambda);
       mat H_inv = ws.U * diagmat(weights_inv) * ws.U.t();
       
       mat XtHinvX = ws.X_orig.t() * H_inv * ws.X_orig;
       XtHinvX.diag() += 1e-12; // Small regularization
       
       mat XtHinvX_inv;
       if (!inv_sympd(XtHinvX_inv, XtHinvX)) {
         XtHinvX_inv = inv(XtHinvX);
       }
       vec beta = XtHinvX_inv * ws.X_orig.t() * H_inv * y;
      beta_results.col(v) = beta;
      
      vec resid = y - ws.X_orig * beta;
      mat K_ZtHinv = ws.Z_orig.t() * H_inv; // Assuming K = I
      vec u = K_ZtHinv * resid;
      u_results.col(v) = u;
      
      // Standard errors (if requested)
      if (compute_se) {
        mat beta_cov = Vu * XtHinvX_inv;
        beta_se_results.col(v) = sqrt(beta_cov.diag());
        
        mat u_cov_diag = Vu * diagmat(ones(ws.q));
        u_se_results.col(v) = sqrt(u_cov_diag.diag());
      }
      
    } catch (const std::exception& e) {
      // Set results to NA for failed voxels
      beta_results.col(v).fill(datum::nan);
      u_results.col(v).fill(datum::nan);
      Vu_results(v) = datum::nan;
      Ve_results(v) = datum::nan;
      lambda_results(v) = datum::nan;
    }
  }
  
  List result = List::create(
    Named("beta") = beta_results.t(), // Transpose to match R convention (trials × voxels)
    Named("u") = u_results.t(),       // Transpose to match R convention  
    Named("Vu") = Vu_results,
    Named("Ve") = Ve_results,
    Named("lambda") = lambda_results
  );
  
  if (compute_se) {
    result["beta_se"] = beta_se_results.t();
    result["u_se"] = u_se_results.t();
  }
  
  return result;
} 