#' Least Squares All (LSA) Analysis
#'
#' Performs a standard multiple regression analysis where all trial regressors
#' are fitted simultaneously. This provides a reference comparison to the 
#' Least Squares Separate (LSS) approach.
#'
#' @param Y A numeric matrix where rows are timepoints and columns are 
#'   voxels/features. This is the dependent variable data.
#' @param X A numeric matrix where rows are timepoints and columns are 
#'   trial-specific regressors. Each column represents a single trial or event.
#' @param Z A numeric matrix of nuisance regressors (e.g., motion parameters,
#'   drift terms). Defaults to NULL.
#' @param Nuisance An alias for Z, provided for consistency with LSS interface.
#'   If both Z and Nuisance are provided, Z takes precedence.
#' @param method Character string specifying the computational method:
#'   \itemize{
#'     \item "r" - Pure R implementation using lm.fit
#'     \item "cpp" - C++ implementation for better performance
#'   }
#'
#' @return A numeric matrix of size T × V containing the beta estimates for
#'   each trial regressor (rows) and each voxel (columns).
#'
#' @details
#' LSA fits the model: Y = X*beta + Z*gamma + error, where all trial regressors
#' in X are estimated simultaneously. This is in contrast to LSS, which fits
#' each trial separately while treating other trials as nuisance regressors.
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
#' # Run LSA analysis
#' beta_estimates <- lsa(Y, X)
#'
#' @export
lsa <- function(Y, X, Z = NULL, Nuisance = NULL, 
                method = c("r", "cpp")) {
  
  method <- match.arg(method)
  
  # Input validation (reuse logic from lss function)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # Handle Nuisance parameter (alias for Z)
  if (is.null(Z) && !is.null(Nuisance)) {
    Z <- Nuisance
  }
  
  # Ensure Z is a matrix if provided
  if (!is.null(Z) && !is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  # Dimension checks
  if (nrow(Y) != nrow(X)) {
    stop("Y and X must have the same number of rows (timepoints)")
  }
  
  if (!is.null(Z) && nrow(Y) != nrow(Z)) {
    stop("Y and Z must have the same number of rows (timepoints)")
  }
  
  # Step 1: Clean data by removing any rows with missing values
  if (is.null(Z)) {
    combined_data <- cbind(Y, X)
  } else {
    combined_data <- cbind(Y, X, Z)
  }
  
  complete_rows <- complete.cases(combined_data)
  if (!all(complete_rows)) {
    Y <- Y[complete_rows, , drop = FALSE]
    X <- X[complete_rows, , drop = FALSE]
    if (!is.null(Z)) {
      Z <- Z[complete_rows, , drop = FALSE]
    }
    warning(paste("Removed", sum(!complete_rows), "rows with missing values"))
  }
  
  # Step 2: Run LSA analysis with the chosen method
  result <- switch(method,
    "r" = .lsa_r(Y, X, Z),
    "cpp" = .lsa_cpp(Y, X, Z),
    stop("Unknown method: ", method)
  )
  
  # Step 3: Add dimension names
  rownames(result) <- colnames(X)
  colnames(result) <- colnames(Y)
  
  return(result)
}

# Implementation functions
.lsa_r <- function(Y, X, Z) {
  # Combine trial regressors and confounds into full design matrix
  if (is.null(Z)) {
    design_matrix <- X
  } else {
    design_matrix <- cbind(X, Z)
  }
  
  # Use lm.fit for efficient computation
  fit <- lm.fit(design_matrix, Y)
  
  # Extract only the coefficients for trial regressors (first ncol(X) rows)
  # Handle case where Y has only one column (returns vector instead of matrix)
  if (is.matrix(fit$coefficients)) {
    beta_estimates <- fit$coefficients[1:ncol(X), , drop = FALSE]
  } else {
    # Convert vector to matrix for single response case
    beta_estimates <- matrix(fit$coefficients[1:ncol(X)], ncol(X), 1)
  }
  
  return(beta_estimates)
}

.lsa_cpp <- function(Y, X, Z) {
  # For now, fall back to R implementation
  # Could implement C++ version later if needed
  .lsa_r(Y, X, Z)
} 