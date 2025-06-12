#' Optimized Mixed Model Solver
#'
#' An optimized implementation of mixed model estimation that precomputes
#' expensive matrix operations and can be reused across multiple voxels
#' for significant performance improvements.
#'
#' @param X Fixed effects design matrix (n × p)
#' @param Z Random effects design matrix (n × q) 
#' @param Y Response data - can be a vector (single voxel) or matrix (n × V for multiple voxels)
#' @param K Kinship/covariance matrix for random effects (q × q). Defaults to identity.
#' @param workspace Precomputed workspace (optional, will compute if NULL)
#' @param compute_se Whether to compute standard errors (default: FALSE)
#' @param n_threads Number of OpenMP threads for multi-voxel (0 = auto)
#' @return List with estimated parameters and variance components
#' @export
mixed_solve_optimized <- function(X, Z, Y, K = NULL, workspace = NULL, 
                                 compute_se = FALSE, n_threads = 0) {
  
  # Input validation
  if (missing(Y)) {
    stop("Y (response data) must be provided")
  }
  
  # Convert Y to matrix if it's a vector
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  # Ensure Y is a matrix
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  
  # Default kinship matrix to identity
  if (is.null(K)) {
    K <- diag(ncol(Z))
  }
  
  # Precompute workspace if not provided
  if (is.null(workspace)) {
    workspace <- mixed_precompute_workspace(X, Z, K)
  }
  
  # Single voxel case (1 column)
  if (ncol(Y) == 1) {
    return(mixed_single_voxel_cpp(Y[, 1], workspace, compute_se))
  }
  
  # Multi-voxel case (multiple columns)
  return(mixed_multi_voxel_cpp(Y, workspace, compute_se, n_threads))
}

#' Precompute Workspace for Optimized Mixed Model
#'
#' Performs expensive matrix computations that don't depend on the response
#' vector, allowing for efficient reuse across multiple voxels.
#'
#' @param X Fixed effects design matrix (n × p)
#' @param Z Random effects design matrix (n × q) 
#' @param K Kinship/covariance matrix for random effects (q × q)
#' @return Workspace object for use with mixed_solve_optimized
#' @export
mixed_precompute <- function(X, Z, K = NULL) {
  if (is.null(K)) {
    K <- diag(ncol(Z))
  }
  mixed_precompute_workspace(X, Z, K)
}
