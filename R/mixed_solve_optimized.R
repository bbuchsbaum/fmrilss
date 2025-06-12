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

#' Benchmark Mixed Model Implementations
#'
#' Compare performance between standard mixed_solve and optimized version.
#'
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix
#' @param K Kinship matrix (optional, defaults to identity)
#' @param Y Response matrix (n × V)
#' @param n_reps Number of repetitions for benchmarking
#' @return Data frame with timing results
#' @export
benchmark_mixed_solve <- function(X, Z, K = NULL, Y, n_reps = 5) {
  
  if (is.null(K)) {
    K <- diag(ncol(Z))
  }
  
  n_voxels <- ncol(Y)
  
  cat("Benchmarking mixed model implementations...\n")
  cat("Data: n =", nrow(X), ", p =", ncol(X), ", q =", ncol(Z), 
      ", voxels =", n_voxels, "\n")
  
  # Benchmark standard mixed_solve (single voxel at a time)
  cat("Testing standard mixed_solve...\n")
  times_standard <- replicate(n_reps, {
    start_time <- Sys.time()
    for (v in 1:min(10, n_voxels)) {  # Limit to 10 voxels for standard version
      result <- mixed_solve(Y[, v], X = X, Z = Z, K = K)
    }
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  })
  
  # Benchmark optimized version
  cat("Testing optimized mixed_solve...\n")
  times_optimized <- replicate(n_reps, {
    start_time <- Sys.time()
    result <- mixed_solve_optimized(X, Z, Y, K)
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  })
  
  # Results
  results <- data.frame(
    method = c("standard", "optimized"),
    mean_time = c(mean(times_standard), mean(times_optimized)),
    median_time = c(stats::median(times_standard), stats::median(times_optimized)),
    min_time = c(min(times_standard), min(times_optimized)),
    max_time = c(max(times_standard), max(times_optimized)),
    sd_time = c(stats::sd(times_standard), stats::sd(times_optimized))
  )
  
  # Speedup calculation (per voxel basis)
  voxels_tested_standard <- min(10, n_voxels)
  time_per_voxel_standard <- results$mean_time[1] / voxels_tested_standard
  time_per_voxel_optimized <- results$mean_time[2] / n_voxels
  speedup <- time_per_voxel_standard / time_per_voxel_optimized
  
  cat("\nResults:\n")
  print(results)
  cat("\nPer-voxel timing:\n")
  cat("Standard:", sprintf("%.4f", time_per_voxel_standard), "sec/voxel\n")
  cat("Optimized:", sprintf("%.4f", time_per_voxel_optimized), "sec/voxel\n") 
  cat("Speedup:", sprintf("%.2fx", speedup), "\n")
  
  return(results)
} 