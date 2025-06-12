#' VoxelHRF object
#'
#' S3 class returned by \code{estimate_voxel_hrf} containing voxel-wise HRF basis coefficients.
#' @name VoxelHRF
NULL

#' LSSBeta object
#'
#' S3 class returned by \code{lss_with_hrf} containing trial-wise beta estimates.
#' @name LSSBeta
NULL

#' Estimate Voxel-wise HRF Basis Coefficients
#'
#' Fits a GLM to estimate HRF basis coefficients for every voxel.
#'
#' @param Y Numeric matrix of BOLD data (time \times voxels).
#' @param events Data frame with \code{onset}, \code{duration} and \code{condition} columns.
#' @param basis HRF object from the \code{fmrihrf} package.
#' @param nuisance_regs Optional numeric matrix of nuisance regressors.
#'
#' @return A \link{VoxelHRF} object containing at least:
#'   \item{coefficients}{Matrix of HRF basis coefficients.}
#'   \item{basis}{The HRF basis object used.}
#'   \item{conditions}{Character vector of modeled conditions.}
#' @export
estimate_voxel_hrf <- function(Y, events, basis, nuisance_regs = NULL) {
  # Input validation
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a numeric matrix")
  }

  required_cols <- c("onset", "duration", "condition")
  if (!is.data.frame(events) || !all(required_cols %in% names(events))) {
    stop("events must be a data.frame with columns onset, duration, condition")
  }

  if (!inherits(basis, "HRF")) {
    stop("basis must be an 'HRF' object from fmrihrf")
  }

  if (!is.null(nuisance_regs)) {
    if (!is.matrix(nuisance_regs) || !is.numeric(nuisance_regs)) {
      stop("nuisance_regs must be a numeric matrix")
    }
    if (nrow(nuisance_regs) != nrow(Y)) {
      stop("nuisance_regs must have same number of rows as Y")
    }
  }

  # Construct HRF regressor matrix using fmrihrf
  X_basis <- fmrihrf::regressor_set(events, basis = basis, n = nrow(Y))$X

  X_full <- X_basis
  if (!is.null(nuisance_regs)) {
    X_full <- cbind(X_full, nuisance_regs)
  }

  coef_mat <- estimate_hrf_cpp(X_full, Y)

  result <- list(
    coefficients = coef_mat,
    basis = basis,
    conditions = unique(as.character(events$condition))
  )
  class(result) <- "VoxelHRF"
  result
}

#' Perform LSS using Voxel-wise HRFs
#'
#' Computes trial-wise beta estimates using voxel-specific HRFs.
#'
#' @param Y Numeric matrix of BOLD data (time \times voxels).
#' @param events Data frame with \code{onset}, \code{duration} and \code{condition} columns.
#' @param hrf_estimates A \link{VoxelHRF} object returned by \code{estimate_voxel_hrf}.
#' @param nuisance_regs Optional numeric matrix of nuisance regressors.
#' @param engine Computational engine. Currently only "C++" placeholder.
#' @param chunk_size Number of voxels to process per batch.
#' @param verbose Logical; display progress bar.
#'
#' @return An object of class \link{LSSBeta}.
#' @export
lss_with_hrf <- function(Y, events, hrf_estimates, nuisance_regs = NULL,
                         engine = "C++", chunk_size = 5000, verbose = TRUE) {
  stop("lss_with_hrf not yet implemented")
}
