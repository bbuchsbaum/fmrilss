#' Benchmark Mixed Model Implementations
#'
#' Compare performance between the standard `mixed_solve` implementation and the
#' optimized `mixed_solve_optimized` version.
#'
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix
#' @param K Kinship matrix (optional, defaults to identity)
#' @param Y Response matrix (n \u00d7 V)
#' @param n_reps Number of repetitions for benchmarking
#' @return Data frame with timing results
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 2), 100, 2)
#' Z <- matrix(rnorm(100 * 3), 100, 3)
#' Y <- matrix(rnorm(100 * 5), 100, 5)
#' benchmark_mixed_solve(X, Z, Y = Y, n_reps = 2)
#' }
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
