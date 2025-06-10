#' A wrapper for the optimized C++ LSS implementation
#'
#' @param Y the voxel by time data matrix
#' @param bdes the block design list created by \code{block_design}
#' @return a matrix of beta estimates
#' @export
lss_cpp_optimized <- function(Y, bdes) {
  X_confounds <- if (!is.null(bdes$dmat_fixed)) {
    cbind(bdes$dmat_base, bdes$dmat_fixed)
  } else {
    bdes$dmat_base
  }
  
  # Step 1: Project confounds from data
  pres <- compute_residuals_optim(X_confounds, Y, bdes$dmat_ran)
  
  # Step 2: Compute beta estimates
  betas <- lss_compute_optim(pres$Q_dmat_ran, pres$residual_data)
  
  return(betas)
} 