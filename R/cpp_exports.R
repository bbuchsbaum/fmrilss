#' Project Out Confounds Using C++
#'
#' Fast C++ implementation for projecting out confound variables from data and
#' trial design matrices. This uses Cholesky decomposition for numerical stability
#' and avoids creating large projection matrices.
#'
#' @param X_confounds Confound design matrix (n x k)
#' @param Y_data Data matrix (n x V) where V is number of voxels
#' @param C_trials Trial design matrix (n x T) where T is number of trials
#' @return List with projected data (residual_data) and projected trials (Q_dmat_ran)
#'
#' @details
#' This function computes residuals Y - X(X'X)^(-1)X'Y and C - X(X'X)^(-1)X'C
#' without explicitly forming the projection matrix Q = I - X(X'X)^(-1)X'.
#' This approach uses ~100x less memory for large n and is numerically more stable.
#'
#' @examples
#' \donttest{
#' n <- 100; V <- 50; T <- 10
#' X_confounds <- cbind(1, seq_len(n), matrix(rnorm(n * 3), n, 3))
#' Y_data <- matrix(rnorm(n*V), n, V)
#' C_trials <- matrix(rnorm(n*T), n, T)
#'
#' result <- project_confounds_cpp(X_confounds, Y_data, C_trials)
#' }
#'
#' @export
project_confounds_cpp <- function(X_confounds, Y_data, C_trials) {
  # Use the actual exported C++ function name
  compute_residuals_cpp(X_confounds, Y_data, C_trials)
}

#' Vectorized LSS Beta Computation Using C++
#'
#' Fast C++ implementation of least squares separate (LSS) beta estimation using
#' vectorized matrix operations. Computes all trial betas in a single pass without loops.
#'
#' @param C_projected Projected trial regressors (n x T) from project_confounds_cpp
#' @param Y_projected Projected data (n x V) from project_confounds_cpp
#' @return Beta matrix (T x V) with LSS estimates for each trial and voxel
#'
#' @details
#' This vectorized implementation computes all LSS betas simultaneously using
#' matrix algebra. It's significantly faster than per-trial loops and automatically
#' benefits from BLAS multithreading. The algorithm handles numerical edge cases
#' by applying finite lower bounds to degenerate denominators.
#'
#' For best performance on large datasets, ensure your R installation uses
#' optimized BLAS (like OpenBLAS or Intel MKL).
#'
#' @examples
#' \donttest{
#' n <- 40; V <- 3; T <- 4
#' X_confounds <- cbind(1, seq_len(n))
#' Y_data <- matrix(rnorm(n * V), n, V)
#' C_trials <- matrix(rnorm(n * T), n, T)
#' result <- project_confounds_cpp(X_confounds, Y_data, C_trials)
#' betas <- lss_beta_cpp(result$Q_dmat_ran, result$residual_data)
#' }
#'
#' @export
lss_beta_cpp <- function(C_projected, Y_projected) {
  # Use the actual exported C++ function name  
  lss_compute_cpp(C_projected, Y_projected)
}
