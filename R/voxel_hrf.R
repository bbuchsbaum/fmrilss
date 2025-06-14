#' VoxelHRF object
#'
#' Simple list-based S3 class returned by \code{estimate_voxel_hrf} containing
#' voxel-wise HRF basis coefficients and related metadata.
#' @name VoxelHRF
NULL

#' LSSBeta object
#'
#' Simple list-based S3 class returned by \code{lss_with_hrf} containing
#' trial-wise beta estimates.
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
#'
#' @examples
#' if (requireNamespace("fmrihrf", quietly = TRUE)) {
#'   set.seed(1)
#'   Y <- matrix(rnorm(100), 50, 2)
#'   events <- data.frame(onset = c(5, 25), duration = 1,
#'                        condition = "A")
#'   basis <- fmrihrf::HRF_GAMMA()
#'   X <- fmrihrf::regressor_set(events, basis, n = nrow(Y))$X
#'   coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
#'   Y <- X %*% coef + Y * 0.1
#'   est <- estimate_voxel_hrf(Y, events, basis)
#'   str(est)
#' }
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

  coef_basis <- coef_mat[seq_len(ncol(X_basis)), , drop = FALSE]

  result <- list(
    coefficients = coef_basis,
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
#' @param backing_dir Directory for bigmemory backing files. If NULL, a
#'   temporary directory is used.
#'
#' @return An object of class \link{LSSBeta}.
#'
#' @examples
#' if (requireNamespace("fmrihrf", quietly = TRUE)) {
#'   set.seed(1)
#'   Y <- matrix(rnorm(100), 50, 2)
#'   events <- data.frame(onset = c(5, 25), duration = 1,
#'                        condition = "A")
#'   basis <- fmrihrf::HRF_GAMMA()
#'   X <- fmrihrf::regressor_set(events, basis, n = nrow(Y))$X
#'   coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
#'   Y <- X %*% coef + Y * 0.1
#'   est <- estimate_voxel_hrf(Y, events, basis)
#'   betas <- lss_with_hrf(Y, events, est, verbose = FALSE, chunk_size = 1)
#'   betas$betas[1:2, ]
#' }
#' @export
lss_with_hrf <- function(Y, events, hrf_estimates, nuisance_regs = NULL,
                         engine = "C++", chunk_size = 5000, verbose = TRUE,
                         backing_dir = NULL) {
  engine <- match.arg(engine, c("C++", "R"))

  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a numeric matrix")
  }

  required_cols <- c("onset", "duration", "condition")
  if (!is.data.frame(events) || !all(required_cols %in% names(events))) {
    stop("events must be a data.frame with columns onset, duration, condition")
  }

  if (!inherits(hrf_estimates, "VoxelHRF") ||
      !is.matrix(hrf_estimates$coefficients)) {
    stop("hrf_estimates must be a VoxelHRF object with coefficient matrix")
  }

  if (!is.null(nuisance_regs)) {
    if (!is.matrix(nuisance_regs) || !is.numeric(nuisance_regs) ||
        nrow(nuisance_regs) != nrow(Y)) {
      stop("nuisance_regs must be a numeric matrix with same number of rows as Y")
    }
  }

  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size <= 0) {
    stop("chunk_size must be a positive integer")
  }

  onsets <- as.numeric(events$onset)
  durations <- as.numeric(events$duration)
  conditions <- as.character(events$condition)

  fine_dt <- 0.1
  max_time <- max(onsets + durations) + fmrihrf::hrf_length(hrf_estimates$basis)
  fine_grid <- seq(0, max_time, by = fine_dt)
  hrf_basis_kernels <- fmrihrf::evaluate_hrf_basis(hrf_estimates$basis, fine_grid)

  onset_idx <- as.integer(round(onsets / fine_dt)) + 1L

  if (is.null(backing_dir)) {
    backing_dir <- tempdir()
  }
  bfile <- tempfile("betas", tmpdir = backing_dir, fileext = ".bin")
  descfile <- tempfile("betas", tmpdir = backing_dir, fileext = ".desc")
  betas <- bigmemory::filebacked.big.matrix(
    nrow = length(onsets), ncol = ncol(Y), type = "double",
    backingfile = basename(bfile), descriptorfile = basename(descfile),
    backingpath = dirname(bfile)
  )

  pb <- NULL
  update_progress <- function(step) {}
  if (verbose) {
    pb <- progress::progress_bar$new(total = ncol(Y))
    update_progress <- function(step) pb$tick(step)
  }

  lss_engine_vox_hrf(Y, hrf_estimates$coefficients, hrf_basis_kernels,
                     onset_idx, durations, nuisance_regs,
                     betas@address, update_progress,
                     as.integer(chunk_size), verbose)

  result <- list(betas = betas)
  class(result) <- "LSSBeta"
  result
}
