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
#' @export
lss_optimized <- function(Y = NULL, bdes, dset = NULL, use_cpp = TRUE) {
  # This function acts as a wrapper, calling the optimized engine.
  # The original 'lss_fast' is kept for comparison.
  .lss_engine_optimized(dset = dset, bdes = bdes, Y = Y, use_cpp = use_cpp)
}

#' Optimized R Projection via Residualization
#'
#' Computes projected data without forming the n x n projection matrix Q.
#' This is mathematically equivalent to `Y - X %*% solve(t(X) %*% X) %*% t(X) %*% Y`
#' but more numerically stable and efficient.
#'
#' @param X Confound design matrix (n x p).
#' @param Y Data matrix (n x V).
#' @param C Trial design matrix (n x T).
#' @return A list containing the residualized (projected) Y and C matrices.
#' @keywords internal
#' @noRd
.project_R_optimized <- function(X, Y, C) {
  XtX  <- crossprod(X)                # k x k, efficient
  Xi   <- chol2inv(chol(XtX))         # Fast and stable inverse
  XtY  <- crossprod(X, Y)             # k x V
  XtC  <- crossprod(X, C)             # k x T
  
  # Return a list of the residualized matrices
  list(
    Y_res = Y - (X %*% (Xi %*% XtY)),
    C_res = C - (X %*% (Xi %*% XtC))
  )
}

#' Memory-Efficient and Algebraically Optimized LSS Beta Computation
#'
#' This version avoids creating a large T x V temporary matrix for BtY.
#' The numerator is rewritten as: num_i = (1 + alpha_i)*C_i'Y - alpha_i*total'Y
#'
#' @param C Projected trial regressors (n x T).
#' @param Y Projected data (n x V).
#' @param eps Numerical tolerance.
#' @return Beta matrix (T x V).
#' @keywords internal
#' @noRd
.lss_beta_vec_optimized <- function(C, Y, eps = 1e-12) {
  T_trials <- ncol(C)
  V_voxels <- ncol(Y)
  
  # Shared building blocks
  total   <- rowSums(C)
  ss_tot  <- sum(total^2)
  CtY     <- crossprod(C, Y)
  CtC     <- colSums(C^2)
  CtT     <- crossprod(C, total)
  total_Y <- drop(crossprod(total, Y))

  # Per-trial "other" pieces
  bt2 <- ss_tot - 2*CtT + CtC
  ctb <- CtT - CtC
  
  # Guard against near-zero regressors
  bt2[bt2 < eps] <- Inf
  
  # Numerator & denominator calculation using algebraic rewrite
  alpha <- ctb / bt2
  
  # Initialize the numerator matrix
  num <- matrix(0, nrow = T_trials, ncol = V_voxels)
  
  # Loop over voxels to avoid broadcasting issues
  for (j in 1:V_voxels) {
    num[, j] <- (1 + alpha) * CtY[, j] - alpha * total_Y[j]
  }

  den <- CtC - (ctb^2) / bt2
  
  # Final division with guard
  sweep(num, 1, pmax(den, eps), `/`)
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
  
  # --- C++ path is unchanged ---
  if (use_cpp) {
    # The C++ path is already optimized and doesn't form the Q matrix explicitly
    res <- compute_residuals_cpp(X_base_fixed, Y, dmat_ran)
    return(lss_compute_cpp(res$Q_dmat_ran, res$residual_data))
  }
  
  # --- OPTIMIZED PURE R PATH ---
  
  # Hot-path early exit for a single event
  if (n_events == 1) {
    # Use the efficient lm.fit approach for single-event case
    residuals_all <- lm.fit(X_base_fixed, cbind(Y, dmat_ran))$residuals
    Ry <- residuals_all[, 1:ncol(Y), drop = FALSE]
    QC <- residuals_all[, ncol(Y) + 1, drop = FALSE]
    
    beta_matrix <- matrix(crossprod(QC, Ry) / drop(crossprod(QC)), nrow = 1)
    return(beta_matrix)
  }

  # For multiple events, use the new optimized projection and beta functions
  proj <- .project_R_optimized(X_base_fixed, Y, dmat_ran)
  
  .lss_beta_vec_optimized(proj$C_res, proj$Y_res)
} 