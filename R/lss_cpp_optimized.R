#' A wrapper for the optimized C++ LSS implementation
#'
#' @param Y the voxel by time data matrix
#' @param bdes the block design list created by \code{block_design}
#' @return a matrix of beta estimates
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(16), 8, 2)
#' X_trials <- matrix(0, 8, 2)
#' X_trials[2:3, 1] <- 1
#' X_trials[5:6, 2] <- 1
#' bdes <- list(
#'   dmat_base = matrix(1, 8, 1),
#'   dmat_fixed = NULL,
#'   dmat_ran = X_trials
#' )
#' lss_cpp_optimized(Y, bdes)
#' @export
lss_cpp_optimized <- function(Y, bdes) {
  X_confounds <- if (!is.null(bdes$dmat_fixed)) {
    cbind(bdes$dmat_base, bdes$dmat_fixed)
  } else {
    bdes$dmat_base
  }
  
  # Step 1: Project confounds from data
  pres <- project_confounds_cpp(X_confounds, Y, bdes$dmat_ran)

  # Step 2: Compute beta estimates
  betas <- lss_beta_cpp(pres$Q_dmat_ran, pres$residual_data)
  
  return(betas)
} 
