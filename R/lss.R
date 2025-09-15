#' Least Squares Separate (LSS) Analysis
#'
#' Computes trial-wise beta estimates using the Least Squares Separate approach
#' of Mumford et al. (2012). This method fits a separate GLM for each trial,
#' with the trial of interest and all other trials as separate regressors.
#'
#' @param Y A numeric matrix of size n × V where n is the number of timepoints
#'   and V is the number of voxels/variables
#' @param X A numeric matrix of size n × T where T is the number of trials.
#'   Each column represents the design for one trial
#' @param Z A numeric matrix of size n × F representing experimental regressors
#'   to include in all trial-wise models. These are regressors we want to model
#'   and get beta estimates for, but not trial-wise (e.g., intercept, condition
#'   effects, block effects). If NULL, an intercept-only design is used. Defaults to NULL
#' @param Nuisance A numeric matrix of size n × N representing nuisance regressors
#'   to be projected out before LSS analysis (e.g., motion parameters, physiological
#'   noise). If NULL, no nuisance projection is performed. Defaults to NULL
#' @param method Character string specifying which implementation to use.
#'   Options are:
#'   \itemize{
#'     \item "r_optimized" - Optimized R implementation (recommended, default)
#'     \item "cpp_optimized" - Optimized C++ implementation with parallel support
#'     \item "r_vectorized" - Standard R vectorized implementation  
#'     \item "cpp" - Standard C++ implementation
#'     \item "naive" - Simple loop-based R implementation (for testing)
#'     \item "oasis" - OASIS method with HRF support and ridge regularization
#'   }
#' @param block_size An integer specifying the voxel block size for parallel
#'   processing, only applicable when `method = "cpp_optimized"`. Defaults to 96.
#' @param oasis A list of options for the OASIS method. See Details for available options.
#'
#' @return A numeric matrix of size T × V containing the trial-wise beta estimates.
#'   Note: Currently only returns estimates for the trial regressors (X). Beta
#'   estimates for the experimental regressors (Z) are computed but not returned.
#'
#' @details
#' The LSS approach fits a separate GLM for each trial, where each model includes:
#' \itemize{
#'   \item The trial of interest (from column i of X)
#'   \item All other trials combined (sum of all other columns of X) 
#'   \item Experimental regressors (Z matrix) - these are modeled to get beta estimates but not trial-wise
#' }
#' 
#' If Nuisance regressors are provided, they are first projected out from both
#' Y and X using standard linear regression residualization.
#' 
#' When using method="oasis", the following options are available in the oasis list:
#' \itemize{
#'   \item \code{design_spec}: A list for building trial-wise designs from event onsets using fmrihrf.
#'     Must contain: \code{sframe} (sampling frame), \code{cond} (list with \code{onsets},
#'     \code{hrf}, and optionally \code{span}), and optionally \code{others} (list of other conditions
#'     to be modeled as nuisances). When provided, X can be NULL and will be constructed automatically.
#'   \item \code{K}: Explicit basis dimension for multi-basis HRF models (e.g., 3 for SPMG3).
#'     If not provided, it's auto-detected from X dimensions or defaults to 1 for single-basis HRFs.
#'   \item \code{ridge_mode}: Either "absolute" (default) or "fractional". In absolute mode,
#'     ridge_x and ridge_b are used directly as regularization parameters. In fractional mode,
#'     they represent fractions of the maximum eigenvalue for adaptive regularization.
#'   \item \code{ridge_x}: Ridge parameter for trial-specific regressors (default 0). Controls
#'     regularization strength for individual trial estimates.
#'   \item \code{ridge_b}: Ridge parameter for the aggregator regressor (default 0). Controls
#'     regularization strength for the sum of all other trials.
#'   \item \code{return_se}: Logical, whether to return standard errors (default FALSE). When TRUE,
#'     returns a list with \code{beta} (trial estimates) and \code{se} (standard errors) components.
#'   \item \code{return_diag}: Logical, whether to return design diagnostics (default FALSE).
#'     When TRUE, includes diagnostic information about the design matrix structure.
#'   \item \code{block_cols}: Integer, voxel block size for memory-efficient processing (default 4096).
#'     Larger values use more memory but may be faster for systems with sufficient RAM.
#'   \item \code{whiten}: Logical, whether to apply AR(1) whitening (default FALSE). When TRUE,
#'     estimates AR(1) coefficients and pre-whitens data to account for temporal autocorrelation.
#'   \item \code{ntrials}: Explicit number of trials (used when K > 1 to determine output dimensions).
#'     If not provided, calculated as ncol(X) / K.
#'   \item \code{hrf_grid}: Vector of HRF indices for grid-based HRF selection (advanced use).
#'     Allows testing multiple HRF shapes simultaneously.
#' }
#'
#' The OASIS method provides a mathematically equivalent but computationally optimized version
#' of standard LSS. It reformulates the per-trial GLM fitting as a single matrix operation,
#' eliminating redundant computations. This is particularly beneficial for designs with many
#' trials or when processing large datasets. When K > 1 (multi-basis HRFs), the output will
#' have K*ntrials rows, with basis functions for each trial arranged sequentially.
#'
#' @references
#' Mumford, J. A., Turner, B. O., Ashby, F. G., & Poldrack, R. A. (2012).
#' Deconvolving BOLD activation in event-related designs for multivoxel pattern
#' classification analyses. NeuroImage, 59(3), 2636-2643.
#'
#' @examples
#' # Generate example data
#' n_timepoints <- 100
#' n_trials <- 10
#' n_voxels <- 50
#' 
#' # Create trial design matrix
#' X <- matrix(0, n_timepoints, n_trials)
#' for(i in 1:n_trials) {
#'   start <- (i-1) * 8 + 1
#'   if(start + 5 <= n_timepoints) {
#'     X[start:(start+5), i] <- 1
#'   }
#' }
#' 
#' # Create data with some signal
#' Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
#' true_betas <- matrix(rnorm(n_trials * n_voxels, 0, 0.5), n_trials, n_voxels)
#' for(i in 1:n_trials) {
#'   Y <- Y + X[, i] %*% matrix(true_betas[i, ], 1, n_voxels)
#' }
#' 
#' # Run LSS analysis
#' beta_estimates <- lss(Y, X)
#' 
#' # With experimental regressors (intercept + condition effects)
#' Z <- cbind(1, scale(1:n_timepoints))
#' beta_estimates_with_regressors <- lss(Y, X, Z = Z)
#' 
#' # With nuisance regression (motion parameters)
#' Nuisance <- matrix(rnorm(n_timepoints * 6), n_timepoints, 6)
#' beta_estimates_clean <- lss(Y, X, Z = Z, Nuisance = Nuisance)
#'
#' \dontrun{
#' # Using OASIS method with ridge regularization
#' beta_oasis <- lss(Y, X, method = "oasis",
#'                   oasis = list(ridge_x = 0.1, ridge_b = 0.1,
#'                               ridge_mode = "fractional"))
#'
#' # OASIS with standard errors
#' result_with_se <- lss(Y, X, method = "oasis",
#'                      oasis = list(return_se = TRUE))
#' beta_estimates <- result_with_se$beta
#' standard_errors <- result_with_se$se
#'
#' # Building design from event onsets using fmrihrf
#'   sframe <- sampling_frame(blocklens = 200, TR = 1.0)
#'
#'   # OASIS with automatic design construction
#'   beta_auto <- lss(Y, X = NULL, method = "oasis",
#'                    oasis = list(
#'                      design_spec = list(
#'                        sframe = sframe,
#'                        cond = list(
#'                          onsets = c(10, 30, 50, 70, 90, 110, 130, 150),
#'                          hrf = HRF_SPMG1,
#'                          span = 25
#'                        ),
#'                        others = list(
#'                          list(onsets = c(20, 40, 60, 80, 100, 120, 140))
#'                        )
#'                      )
#'                    ))
#'
#'   # Multi-basis HRF example (3 basis functions per trial)
#'   beta_multibasis <- lss(Y, X = NULL, method = "oasis",
#'                         oasis = list(
#'                           design_spec = list(
#'                             sframe = sframe,
#'                             cond = list(
#'                               onsets = c(10, 30, 50, 70, 90),
#'                               hrf = HRF_SPMG3,  # 3-basis HRF
#'                               span = 30
#'                             )
#'                           ),
#'                           K = 3  # Explicit basis dimension
#'                         ))
#'   # Returns 15 rows (5 trials * 3 basis functions)
#' }
#'
#' @export
lss <- function(Y, X, Z = NULL, Nuisance = NULL, 
                method = c("r_optimized", "cpp_optimized", "r_vectorized", "cpp", "naive", "oasis"),
                block_size = 96, oasis = list()) {
  
  method <- match.arg(method)
  
  # Fast-path to OASIS: it has its own coercion/validation and supports Matrix inputs
  if (method == "oasis") {
    if (is.null(X) && is.null(oasis$design_spec)) {
      stop("For method='oasis', either X or oasis$design_spec must be provided")
    }
    return(.lss_oasis(Y, X, Z, Nuisance, oasis))
  }
  
  # Coerce S4 Matrix/data.frame to base matrices for non-OASIS methods
  to_mat <- function(M) {
    if (is.null(M)) return(NULL)
    if (inherits(M, "Matrix")) return(as.matrix(M))
    if (is.data.frame(M))     return(as.matrix(M))
    M
  }
  Y        <- to_mat(Y)
  X        <- to_mat(X)
  Z        <- to_mat(Z)
  Nuisance <- to_mat(Nuisance)
  
  # Input validation (non-OASIS)
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a numeric matrix")
  }
  if (!is.null(Z) && (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != nrow(Y))) {
    stop("Z must be a numeric matrix with the same number of rows as Y")
  }
  if (!is.null(Nuisance) && (!is.matrix(Nuisance) || !is.numeric(Nuisance) || nrow(Nuisance) != nrow(Y))) {
    stop("Nuisance must be a numeric matrix with the same number of rows as Y")
  }
  
  # For non-OASIS methods, X is required
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix")
  }
  if (nrow(Y) != nrow(X)) {
    stop("Y and X must have the same number of rows (timepoints)")
  }
  
  # Set up default experimental regressors (intercept) if not provided
  if (is.null(Z)) {
    Z <- matrix(1, nrow(Y), 1)
    colnames(Z) <- "Intercept"
  }
  
  # Check for zero or near-zero regressors
  .check_zero_regressors(X)
  
  # Step 1: Project out nuisance regressors if provided
  if (!is.null(Nuisance)) {
    # Create full nuisance design matrix
    X_nuisance <- cbind(Z, Nuisance)
    
    # Project out nuisance from Y and X
    proj_result <- .project_out_nuisance(Y, X, X_nuisance)
    Y_clean <- proj_result$Y_residual
    X_clean <- proj_result$X_residual
  } else {
    Y_clean <- Y
    X_clean <- X
  }
  
  # Step 2: Run LSS analysis with the chosen method
  result <- switch(method,
    "r_optimized" = .lss_r_optimized(Y_clean, X_clean, Z),
    "cpp_optimized" = .lss_cpp_optimized(Y_clean, X_clean, Z, block_size = block_size),
    "r_vectorized" = .lss_r_vectorized(Y_clean, X_clean, Z),
    "cpp" = .lss_cpp(Y_clean, X_clean, Z),
    "naive" = .lss_naive(Y_clean, X_clean, Z),
    stop("Unknown method: ", method)
  )
  
  # Add row and column names if available
  if (!is.null(colnames(X))) {
    rownames(result) <- colnames(X)
  } else {
    rownames(result) <- paste0("Trial_", 1:ncol(X))
  }
  
  if (!is.null(colnames(Y))) {
    colnames(result) <- colnames(Y)
  } else {
    colnames(result) <- paste0("Voxel_", 1:ncol(Y))
  }
  
  return(result)
}

# Helper function to project out nuisance regressors
.project_out_nuisance <- function(Y, X, X_nuisance) {
  # Handle empty or NULL nuisance matrix early
  if (is.null(X_nuisance) || !is.matrix(X_nuisance) || ncol(X_nuisance) == 0L) {
    return(list(Y_residual = Y, X_residual = X))
  }

  # Numerically stable residualization via QR (avoids explicit inverse)
  qrX <- qr(X_nuisance)
  coefY <- qr.coef(qrX, Y)
  coefX <- qr.coef(qrX, X)
  
  Y_residual <- Y - X_nuisance %*% coefY
  X_residual <- X - X_nuisance %*% coefX
  
  list(Y_residual = Y_residual, X_residual = X_residual)
}

# Helper function to check for zero or near-zero regressors
.check_zero_regressors <- function(X, eps = 1e-12) {
  trial_names <- if (!is.null(colnames(X))) colnames(X) else paste0("Trial_", 1:ncol(X))
  
  # Check each trial regressor
  for (i in 1:ncol(X)) {
    regressor_norm <- sqrt(sum(X[, i]^2))
    regressor_var <- stats::var(X[, i])
    
    # Check for exactly zero regressor first
    if (regressor_norm < eps) {
      warning(sprintf("Trial regressor '%s' appears to be zero (norm = %g). This may cause numerical issues or NaN results.",
                     trial_names[i], regressor_norm))
    } else if (regressor_var < eps && regressor_norm >= eps) {
      # Non-zero but constant (or nearly constant) regressor
      warning(sprintf("Trial regressor '%s' has very low variance (%g) and may cause numerical instability.",
                     trial_names[i], regressor_var))
    }
  }
}

# Implementation functions (these will call the existing optimized functions)
.lss_r_optimized <- function(Y, X, Z) {
  # Create design list for compatibility with existing function
  bdes <- list(
    dmat_base = Z,
    dmat_ran = X,
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  return(lss_optimized(Y, bdes, use_cpp = FALSE))
}

.lss_cpp_optimized <- function(Y, X, Z, block_size = 96) {
  # This now calls the new single-pass C++ function
  
  # Ensure Z (confounds) is a matrix, even if NULL or a vector
  if (is.null(Z)) {
    Z <- matrix(0, nrow = nrow(Y), ncol = 0)
  } else if (!is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  # X = trial regressors, Z = confounds, Y = data
  # The C++ function expects: X=confounds, Y=data, C=trials
  lss_fused_optim_cpp(X = Z, Y = Y, C = X, block_size = block_size)
}

.lss_r_vectorized <- function(Y, X, Z) {
  # Use existing lss_fast function with use_cpp = FALSE
  bdes <- list(
    dmat_base = Z,
    dmat_ran = X,
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  return(lss_fast(dset = NULL, bdes = bdes, Y = Y, use_cpp = FALSE))
}

.lss_cpp <- function(Y, X, Z) {
  # Use existing lss_fast function with use_cpp = TRUE  
  bdes <- list(
    dmat_base = Z,
    dmat_ran = X,
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  return(lss_fast(dset = NULL, bdes = bdes, Y = Y, use_cpp = TRUE))
}

.lss_naive <- function(Y, X, Z) {
  bdes <- list(
    dmat_base = Z,
    dmat_ran = X,
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  return(lss_naive(Y, bdes))
}

#' Orthogonal Projection Matrix
#'
#' Computes Q = I - X(X'X)^(-1)X' using QR decomposition for numerical stability.
#'
#' @param X Design matrix
#' @return Projection matrix Q
#' @keywords internal
#' @noRd
.Q_project <- function(X) {
  ## returns Q = I − X (XᵀX)⁻¹ Xᵀ   without allocating I
  qrX <- qr(X)
  Q   <- diag(nrow(X))              # allocate once
  Q   <- Q - tcrossprod(qr.Q(qrX))   # Q = I – QQᵀ   (orthonormal Q from QR)
  Q
}

#' Vectorized LSS Beta Computation
#'
#' Computes LSS beta estimates without explicit loops using vectorized operations.
#'
#' @param QC Projected trial regressors (n x T)
#' @param Ry Projected data (n x V)
#' @param eps Numerical tolerance
#' @return Beta matrix (T x V)
#' @keywords internal
#' @noRd
.lss_beta_vec <- function(C, Y, eps = 1e-12) {

  T <- ncol(C); V <- ncol(Y)

  # Shared building blocks ----------------------------------------------------
  total  <- rowSums(C)                     # n
  ss_tot <- sum(total^2)                   # scalar

  CtY <- crossprod(C, Y)                   # T × V
  CtC <- colSums(C^2)                      # T
  CtT <- crossprod(C, total)               # T
  total_Y <- drop(crossprod(total, Y))     # 1 × V  (row vector)

  # Per-trial "other" pieces ---------------------------------------------------
  #   b_i^T y  (T × V)
  BtY <- matrix(rep(total_Y, each = T), T, V) - CtY

  #   ||b_i||^2 ,  c_i^T b_i  (length T)
  bt2 <- ss_tot - 2*CtT + CtC
  ctb <- CtT   - CtC

  # Add guard for near-zero other-trial regressors (Fix 3)
  bt2[bt2 < eps] <- Inf

  # Numerator & denominator (memory-efficient version, Fix 1)
  ctb_bt2 <- ctb / bt2                     # T vector
  num <- CtY - sweep(BtY, 1, ctb_bt2, `*`) # Sweep for column-wise multiplication
  den <- CtC - (ctb^2) / bt2               # T vector
  
  # Broadcast den vector for division
  sweep(num, 1, pmax(den, eps), `/`)
}

#' Validate Input Dimensions and Types
#'
#' @param Y Data matrix
#' @param dmat_ran Trial design matrix
#' @keywords internal
#' @noRd
.validate_dims <- function(Y, dmat_ran) {
  if (!is.matrix(Y)) {
    if (is.data.frame(Y)) {
      Y <- as.matrix(Y)
      warning("Converting Y from data.frame to matrix")
    } else {
      stop("Y must be a numeric matrix or data.frame")
    }
  }
  if (nrow(Y) != nrow(dmat_ran)) {
    stop(sprintf("Y has %d timepoints, design has %d.",
                 nrow(Y), nrow(dmat_ran)))
  }
  Y
}

#' @keywords internal
#' @noRd
lss_fast <- function(dset, bdes, Y = NULL, use_cpp = TRUE) {
  # Data preparation and validation
  Y <- if (is.null(Y)) get_data_matrix(dset) else Y
  Y <- .validate_dims(Y, bdes$dmat_ran)

  # Ensure design matrices are matrices
  dmat_base <- as.matrix(bdes$dmat_base)
  dmat_ran <- as.matrix(bdes$dmat_ran)
  dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
  
  # Check for intercept in base design
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
  
  # Hot-path early exit for single event
  if (n_events == 1) {
    if (use_cpp) {
      res <- compute_residuals_cpp(X_base_fixed, Y, dmat_ran)
      return(lss_compute_cpp(res$Q_dmat_ran, res$residual_data))
    } else {
      # Simple regression for single event
      Q <- .Q_project(X_base_fixed)
      Ry <- Q %*% Y
      QC <- Q %*% dmat_ran
      beta_matrix <- matrix(crossprod(QC, Ry) / drop(crossprod(QC)), nrow = 1)
      return(beta_matrix)
    }
  }

  if (use_cpp) {
    res <- compute_residuals_cpp(X_base_fixed, Y, dmat_ran)
    return(lss_compute_cpp(res$Q_dmat_ran, res$residual_data))
  }

  # --- pure R fall-back ---
  Q  <- .Q_project(X_base_fixed)
  Ry <- Q %*% Y
  QC <- Q %*% dmat_ran
  .lss_beta_vec(QC, Ry)
}

#' Extract Data Matrix from Dataset
#'
#' Helper function to extract data matrix from various dataset formats.
#' This is a placeholder that should be customized based on your data format.
#'
#' @param dset Dataset object (format depends on your specific use case)
#' @return A numeric matrix where rows are timepoints and columns are voxels
#' @keywords internal
#' @export
get_data_matrix <- function(dset) {
  # This is a placeholder function - customize based on your data format
  if (is.matrix(dset)) {
    return(dset)
  } else if (is.data.frame(dset)) {
    return(as.matrix(dset))
  } else {
    stop("Unsupported dataset format. Please provide Y matrix directly or customize get_data_matrix function.")
  }
}

#' Project Out Confound Variables
#'
#' Computes the orthogonal projection matrix Q = I - X(X'X)^(-1)X' that projects
#' out the space spanned by confound regressors X. This is useful for advanced
#' users who want to cache and reuse projection matrices.
#'
#' @param X Confound design matrix (n x p) where n is number of timepoints
#'   and p is number of confound regressors
#' @return Projection matrix Q (n x n) that projects out the column space of X
#'
#' @details
#' This function uses QR decomposition for numerical stability instead of
#' computing the Moore-Penrose pseudoinverse directly. The resulting matrix
#' Q can be applied to data to remove the influence of confound regressors.
#'
#' @examples
#' \dontrun{
#' # Create confound matrix (intercept + linear trend)
#' n <- 100
#' X_confounds <- cbind(1, 1:n)
#' 
#' # Get projection matrix
#' Q <- project_confounds(X_confounds)
#' 
#' # Apply to data to remove confounds
#' Y_clean <- Q %*% Y_raw
#' }
#'
#' @export
project_confounds <- function(X) {
  .Q_project(as.matrix(X))
}
