#' Optimized LSS Analysis (Pure R)
#'
#' An optimized version of the LSS analysis that avoids creating large intermediate
#' matrices, providing a significant speedup and lower memory usage for the pure R
#' implementation.
#'
#' @param Y A numeric matrix where rows are timepoints and columns are voxels/features.
#' @param bdes A list containing the design matrices.
#' @param dset Optional dataset object.
#' @param use_cpp Logical. If TRUE (default), uses the C++ implementation. If FALSE,
#'   uses the new optimized R implementation.
#' @return A numeric matrix of LSS beta estimates.
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(16), 8, 2)
#' X_trials <- matrix(0, 8, 2)
#' X_trials[2:3, 1] <- 1
#' X_trials[5:6, 2] <- 1
#' bdes <- list(
#'   dmat_base = matrix(1, 8, 1),
#'   dmat_ran = X_trials,
#'   dmat_fixed = NULL,
#'   fixed_ind = NULL
#' )
#' lss_optimized(Y, bdes, use_cpp = FALSE)
#' @export
lss_optimized <- function(Y = NULL, bdes, dset = NULL, use_cpp = TRUE) {
  # This function acts as a wrapper, calling the optimized engine.
  # The original 'lss_fast' is kept for comparison.
  .lss_engine_optimized(dset = dset, bdes = bdes, Y = Y, use_cpp = use_cpp)
}

#' LSS Engine (Optimized)
#'
#' @keywords internal
#' @noRd
.lss_engine_optimized <- function(dset, bdes, Y = NULL, use_cpp = TRUE) {
  # Data preparation and validation (uses helpers from lss.R)
  Y <- if (is.null(Y)) get_data_matrix(dset) else Y
  Y <- .validate_dims(Y, bdes$dmat_ran)

  # Ensure design matrices are matrices
  dmat_base <- as.matrix(bdes$dmat_base)
  dmat_ran <- as.matrix(bdes$dmat_ran)
  dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
  
  # Check for intercept
  if (!any(apply(dmat_base, 2, function(x) all(abs(x - mean(x)) < 1e-10)))) {
    warning("No intercept detected in dmat_base. Consider adding one for proper baseline modeling.")
  }
  
  # Build combined base design matrix
  X_base_fixed <- if (!is.null(dmat_fixed)) {
    cbind(dmat_base, dmat_fixed)
  } else {
    dmat_base
  }
  
  n_events <- ncol(dmat_ran)
  
  # Both paths residualize only the (small) trial design and apply the LSS
  # weight matrix to the raw data in a single matrix product.
  if (use_cpp) {
    return(lss_compute_cpp(.lss_residualize_trials(dmat_ran, X_base_fixed), Y))
  }

  .lss_kernel_r(Y, dmat_ran, X_base_fixed)
}
