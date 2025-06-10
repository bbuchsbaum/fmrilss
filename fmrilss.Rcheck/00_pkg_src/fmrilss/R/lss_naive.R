#' Naive Least Squares Separate (LSS) Analysis
#'
#' Performs LSS analysis using the naive approach where each trial model is fit
#' separately. This is the conceptually simplest implementation but less efficient
#' than the optimized \code{\link{lss}} function.
#'
#' @param Y A numeric matrix where rows are timepoints and columns are voxels/features.
#'   If NULL, the function will attempt to extract data from \code{dset}.
#' @param bdes A list containing design matrices with components:
#'   \itemize{
#'     \item \code{dmat_base}: Base design matrix (e.g., intercept, drift terms)
#'     \item \code{dmat_fixed}: Fixed effects design matrix (optional)
#'     \item \code{dmat_ran}: Random/trial design matrix for LSS analysis
#'     \item \code{fixed_ind}: Indices for fixed effects (optional)
#'   }
#' @param dset Optional dataset object. If provided and Y is NULL, data will be
#'   extracted using \code{get_data_matrix}.
#'
#' @return A numeric matrix with dimensions (n_events x n_voxels) containing
#'   the LSS beta estimates for each trial and voxel.
#'
#' @details
#' This function implements the naive LSS approach where for each trial, a separate
#' GLM is fitted that includes:
#' \itemize{
#'   \item All base regressors (intercept, drift, etc.)
#'   \item All fixed effects regressors (if any)
#'   \item Only the current trial's regressor from the trial design matrix
#' }
#'
#' While less efficient than the optimized \code{\link{lss}} function, this
#' implementation is conceptually simpler and can serve as a reference or for
#' validation purposes.
#'
#' @examples
#' \dontrun{
#' # Using same setup as lss() examples
#' beta_estimates_naive <- lss_naive(Y = Y, bdes = bdes)
#' 
#' # Compare with optimized version
#' beta_estimates_fast <- lss(Y = Y, bdes = bdes)
#' max(abs(beta_estimates_naive - beta_estimates_fast))
#' }
#'
#' @seealso \code{\link{lss}} for the optimized implementation
#' @export
lss_naive <- function(Y = NULL, bdes, dset = NULL) {
  # Data preparation
  if (is.null(Y)) {
    data_matrix <- get_data_matrix(dset)
  } else {
    data_matrix <- Y
  }

  # Check and validate Y dimensions
  if (!is.null(Y)) {
    expected_rows <- nrow(bdes$dmat_ran) # Use design matrix rows as reference
    if (nrow(Y) != expected_rows) {
      stop(sprintf("Y must have %d rows to match design matrix dimensions, but has %d rows", 
                  expected_rows, nrow(Y)))
    }
  }

  n_timepoints <- nrow(data_matrix)
  n_voxels <- ncol(data_matrix)
  
  # Design matrices
  dmat_base <- as.matrix(bdes$dmat_base)
  dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
  dmat_ran <- as.matrix(bdes$dmat_ran)
  
  n_events <- ncol(dmat_ran)
  
  # Prepare baseline and fixed design matrix
  if (!is.null(dmat_fixed)) {
    X_base_fixed <- cbind(dmat_base, dmat_fixed)
  } else {
    X_base_fixed <- dmat_base
  }
  
  # Initialize beta matrix
  beta_matrix <- matrix(NA, nrow = n_events, ncol = n_voxels)
  
  # Loop over each trial and fit separate model
  for (i in seq_len(n_events)) {
    # Create design matrix for this trial: base + fixed + current trial
    X_trial <- cbind(X_base_fixed, dmat_ran[, i, drop = FALSE])
    
    # Fit GLM using pseudoinverse (more stable than solve)
    # Beta coefficients for all regressors
    beta_all <- MASS::ginv(X_trial) %*% data_matrix
    
    # Extract beta for the trial regressor (last column)
    beta_matrix[i, ] <- beta_all[nrow(beta_all), ]
  }
  
  return(beta_matrix)
} 